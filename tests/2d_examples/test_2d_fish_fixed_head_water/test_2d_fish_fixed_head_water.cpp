// Fixed-head fish in initially still water. Solid baseline copied unchanged.
#include "fixed_head_water.h"
#include "sphinxsys.h"
#include <memory>
#include <vector>

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the SPH system.
    //----------------------------------------------------------------------
    const Real requested_end = readPositive("FISH_END_TIME", 0.1);
    const Real wall_limit = readPositive("FISH_WALL_SECONDS", 600.0);
    const Real water_spacing = readPositive("FISH_WATER_SPACING", 0.0025);
    const Real drive_scale = readPositive("FISH_DRIVE_SCALE", 1.0, true);
    swim_amplitude *= drive_scale;
    SPHSystem sph_system(BoundingBox(
        Vec2d(outer_x_min - 4 * water_spacing, cy - outer_y_half - 4 * water_spacing),
        Vec2d(outer_x_max + 4 * water_spacing, cy + outer_y_half + 4 * water_spacing)),
        particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);

    //----------------------------------------------------------------------
    //	Creating bodies: one fish-shaped muscle + 10 rigid spine segments.
    //----------------------------------------------------------------------
    SolidBody muscle(sph_system, makeShared<FishMuscleBody>("MuscleBody"));
    muscle.defineAdaptationRatios(1.15, 2.0);
    muscle.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    muscle.defineMaterial<FishTissueComposite>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? muscle.generateParticles<BaseParticles, Reload>(muscle.getName())
        : muscle.generateParticles<BaseParticles, Lattice>();

    std::vector<std::unique_ptr<SolidBody>> bones;
    bones.reserve(bone_segments);
    for (size_t i = 0; i < bone_segments; ++i)
    {
        bones.push_back(std::make_unique<SolidBody>(
            sph_system, makeShared<BoneSegmentBody>("BoneBody_" + std::to_string(i), i)));
        bones[i]->defineAdaptationRatios(1.15, 2.0);
        bones[i]->defineBodyLevelSetShape()->writeLevelSet(sph_system);
        bones[i]->defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, young_bone, poisson);
        bones[i]->generateParticles<BaseParticles, Lattice>();
    }

    //----------------------------------------------------------------------
    //	Topology: only inner relations are needed (attachment-only coupling).
    //----------------------------------------------------------------------
    InnerRelation muscle_inner(muscle);
    std::vector<std::unique_ptr<InnerRelation>> bone_inners;
    for (size_t i = 0; i < bone_segments; ++i)
        bone_inners.push_back(std::make_unique<InnerRelation>(*bones[i]));

    //----------------------------------------------------------------------
    //	Run particle relaxation for the muscle body only, then write reload particles.
    //	The spine segments stay on the lattice because they are rigidized by Simbody.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_muscle_particles(muscle);
        BodyStatesRecordingToVtp write_muscle_states(muscle);
        ReloadParticleIO write_muscle_reload_files(muscle);
        RelaxationStepInner muscle_relaxation_step(muscle_inner);

        random_muscle_particles.exec(0.25);
        muscle_relaxation_step.SurfaceBounding().exec();
        write_muscle_states.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            muscle_relaxation_step.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "Relaxation steps for muscle body N = " << ite_p << "\n";
                write_muscle_states.writeToFile(ite_p);
            }
        }
        std::cout << "The particle relaxation process for muscle body is finished." << std::endl;
        write_muscle_reload_files.writeToFile();
        return 0;
    }

    FluidBody water(sph_system, makeShared<SurroundingWater>("WaterBody"));
    water.defineAdaptationRatios(1.3, particle_spacing_ref / water_spacing);
    water.defineClosure<WeaklyCompressibleFluid, Viscosity>(
        ConstructArgs(water_rho, water_sound_speed), water_mu);
    water.generateParticles<BaseParticles, Lattice>();
    InnerRelation water_inner(water);
    ContactRelation water_contact(water, {&muscle});
    ContactRelation muscle_water_contact(muscle, {&water});
    ComplexRelation water_complex(water_inner, water_contact);

    //----------------------------------------------------------------------
    //	Solid dynamics for muscle and bones.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> muscle_normal(muscle);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> muscle_corrected_configuration(muscle_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> muscle_first_half(muscle_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> muscle_second_half(muscle_inner);

    std::vector<std::unique_ptr<SimpleDynamics<NormalDirectionFromBodyShape>>> bone_normals;
    std::vector<std::unique_ptr<InteractionWithUpdate<LinearGradientCorrectionMatrixInner>>> bone_corrected_configurations;
    std::vector<std::unique_ptr<Dynamics1Level<solid_dynamics::Integration1stHalfPK2>>> bone_first_halves;
    std::vector<std::unique_ptr<Dynamics1Level<solid_dynamics::Integration2ndHalf>>> bone_second_halves;
    std::vector<std::unique_ptr<ReduceDynamics<solid_dynamics::AcousticTimeStep>>> bone_dts;
    for (size_t i = 0; i < bone_segments; ++i)
    {
        bone_normals.push_back(std::make_unique<SimpleDynamics<NormalDirectionFromBodyShape>>(*bones[i]));
        bone_corrected_configurations.push_back(
            std::make_unique<InteractionWithUpdate<LinearGradientCorrectionMatrixInner>>(*bone_inners[i]));
        bone_first_halves.push_back(std::make_unique<Dynamics1Level<solid_dynamics::Integration1stHalfPK2>>(*bone_inners[i]));
        bone_second_halves.push_back(std::make_unique<Dynamics1Level<solid_dynamics::Integration2ndHalf>>(*bone_inners[i]));
        bone_dts.push_back(std::make_unique<ReduceDynamics<solid_dynamics::AcousticTimeStep>>(*bones[i]));
    }

    ReduceDynamics<solid_dynamics::AcousticTimeStep> muscle_dt(muscle);
    SimpleDynamics<FishMaterialInitialization> tag_material_ids(muscle);
    SimpleDynamics<ImposingActiveStrain> impose_active_strain(muscle);

    //----------------------------------------------------------------------
    //	Body parts coupled to Simbody.
    //	Bone parts are kept for rigidization; only segment 0 keeps hard muscle attachments.
    //----------------------------------------------------------------------
    std::vector<SharedPtr<GeometricShapeBox>> bone_shapes;
    std::vector<std::unique_ptr<SolidBodyPartForSimbody>> bone_parts;
    for (size_t i = 0; i < bone_segments; ++i)
    {
        bone_shapes.push_back(createBoneRigidShape(i));
        bone_parts.push_back(std::make_unique<SolidBodyPartForSimbody>(*bones[i], bone_shapes[i]));
    }
    SharedPtr<GeometricShapeBox> head_muscle_shape_top = createMuscleAttachmentShape(0);
    SharedPtr<GeometricShapeBox> head_muscle_shape_bottom = createMuscleAttachmentShapeBottom(0);
    auto head_muscle_attachment_top = std::make_unique<SolidBodyPartForSimbody>(muscle, head_muscle_shape_top);
    auto head_muscle_attachment_bottom = std::make_unique<SolidBodyPartForSimbody>(muscle, head_muscle_shape_bottom);

    //----------------------------------------------------------------------
    //	Simbody multibody system: head welded, segments 1..9 chained by pin joints.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    SimTK::GeneralForceSubsystem forces(MBsystem);
    std::vector<SimTK::Body::Rigid> bone_infos;
    std::vector<Vec2d> bone_coms;
    for (size_t i = 0; i < bone_segments; ++i)
    {
        bone_infos.push_back(SimTK::Body::Rigid(*bone_parts[i]->body_part_mass_properties_));
        bone_coms.push_back(bone_parts[i]->initial_mass_center_);
    }

    SimTK::MobilizedBody::Weld left_mob(
        matter.Ground(),
        SimTK::Transform(SimTKVec3(bone_coms[0][0], bone_coms[0][1], 0.0)),
        bone_infos[0],
        SimTK::Transform(SimTKVec3(0.0, 0.0, 0.0)));

    std::vector<SimTK::MobilizedBody::Pin> pin_mobs;
    pin_mobs.reserve(bone_segments - 1);
    for (size_t i = 1; i < bone_segments; ++i)
    {
        Vec2d hinge_xy = hingePoint(i);
        Vec2d r_hinge_on_prev = hinge_xy - bone_coms[i - 1];
        Vec2d r_hinge_on_this = hinge_xy - bone_coms[i];
        if (i == 1)
        {
            pin_mobs.emplace_back(
                left_mob,
                SimTK::Transform(SimTKVec3(r_hinge_on_prev[0], r_hinge_on_prev[1], 0.0)),
                bone_infos[i],
                SimTK::Transform(SimTKVec3(r_hinge_on_this[0], r_hinge_on_this[1], 0.0)));
        }
        else
        {
            pin_mobs.emplace_back(
                pin_mobs[i - 2],
                SimTK::Transform(SimTKVec3(r_hinge_on_prev[0], r_hinge_on_prev[1], 0.0)),
                bone_infos[i],
                SimTK::Transform(SimTKVec3(r_hinge_on_this[0], r_hinge_on_this[1], 0.0)));
        }
        pin_mobs.back().setDefaultAngle(0.0);
    }

    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    std::vector<std::unique_ptr<SimTK::Force::MobilityLinearDamper>> hinge_dampers;
    for (size_t i = 0; i + 1 < bone_segments; ++i)
        hinge_dampers.push_back(
            std::make_unique<SimTK::Force::MobilityLinearDamper>(forces, pin_mobs[i], SimTK::MobilizerUIndex(0), 0.6));

    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1.0e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    //----------------------------------------------------------------------
    //	Force feedback and constraints between SPH parts and Simbody bodies.
    //----------------------------------------------------------------------
    std::vector<std::unique_ptr<ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>>> force_on_bone_parts;
    std::vector<SolidBody *> bone_body_ptrs;
    for (size_t i = 0; i < bone_segments; ++i)
        bone_body_ptrs.push_back(bones[i].get());
    DistributedAttachmentNetwork distributed_attachment_network(muscle, bone_body_ptrs);
    SimpleDynamics<DistributedAttachmentForce> distributed_attachment_force(
        muscle, distributed_attachment_network, MBsystem, pin_mobs, integ);
    std::vector<std::unique_ptr<ReduceDynamics<DistributedAttachmentReactionForSimBody>>> distributed_attachment_reactions;
    std::vector<std::unique_ptr<SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>>> constraint_bones;

    constraint_bones.push_back(std::make_unique<SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>>(
        *bone_parts[0], MBsystem, left_mob, integ));
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_head_muscle_top(
        *head_muscle_attachment_top, MBsystem, left_mob, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_head_muscle_bottom(
        *head_muscle_attachment_bottom, MBsystem, left_mob, integ);

    for (size_t i = 1; i < bone_segments; ++i)
    {
        force_on_bone_parts.push_back(std::make_unique<ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>>(
            *bone_parts[i], MBsystem, pin_mobs[i - 1], integ));
        distributed_attachment_reactions.push_back(std::make_unique<ReduceDynamics<DistributedAttachmentReactionForSimBody>>(
            muscle, distributed_attachment_network, MBsystem, pin_mobs[i - 1], integ, i));
        constraint_bones.push_back(std::make_unique<SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>>(
            *bone_parts[i], MBsystem, pin_mobs[i - 1], integ));
    }

    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> far_surface(water_inner, water_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_step(water_inner, water_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_step(water_inner, water_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> density_sum(water_inner, water_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> fluid_viscosity(water_inner, water_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport(water_inner, water_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_dt(water, water_reference_speed);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_dt(water);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> fluid_viscous_load(muscle_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_step)>> fluid_pressure_load(muscle_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration solid_average(muscle);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> update_muscle_normal(muscle);
    SimpleDynamics<FarFieldSponge> sponge(water);
    InteractionDynamics<fluid_dynamics::VorticityInner> vorticity(water_inner);

    //----------------------------------------------------------------------
    //	Output.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Matd>(muscle, "ActiveStrain");
    write_states.addToWrite<int>(muscle, "MaterialID");
    write_states.addToWrite<Vecd>(muscle, "AttachmentDisplacement");
    write_states.addToWrite<Real>(muscle, "AttachmentDisplacementNorm");
    write_states.addToWrite<Vecd>(muscle, "AttachmentForce");
    write_states.addToWrite<Real>(muscle, "AttachmentForceNorm");
    write_states.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(muscle);
    write_states.addToWrite<Real>(water, "Pressure");
    write_states.addToWrite<Real>(water, "Density");
    write_states.addToWrite<Vecd>(water, "Velocity");
    write_states.addToWrite<Real>(water, "SpongeWeight");
    write_states.addToWrite<Vecd>(muscle, "PressureForceFromFluid");
    write_states.addToWrite<Vecd>(muscle, "ViscousForceFromFluid");
    WriteSimBodyPinData write_pin_data(sph_system, integ, pin_mobs.back());

    //----------------------------------------------------------------------
    //	Prepare the simulation.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    tag_material_ids.exec();
    muscle_normal.exec();
    muscle_corrected_configuration.exec();
    for (size_t i = 0; i < bone_segments; ++i)
    {
        bone_normals[i]->exec();
        bone_corrected_configurations[i]->exec();
        constraint_bones[i]->exec();
    }
    constraint_head_muscle_top.exec();
    constraint_head_muscle_bottom.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");

    size_t ite = 0;
    TickCount t1 = TickCount::now();
    auto &wp = water.getBaseParticles();
    auto &mp = muscle.getBaseParticles();
    Vecd *wpos = wp.getVariableDataByName<Vecd>("Position");
    Vecd *wvel = wp.getVariableDataByName<Vecd>("Velocity");
    Real *wrho = wp.getVariableDataByName<Real>("Density");
    Vecd *mpos = mp.getVariableDataByName<Vecd>("Position");
    Vecd *mvel = mp.getVariableDataByName<Vecd>("Velocity");
    Vecd *pf = mp.getVariableDataByName<Vecd>("PressureForceFromFluid");
    Vecd *vf = mp.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    std::vector<Vecd> mpos0(mpos, mpos + mp.TotalRealParticles());
    size_t tail_index = 0;
    for (size_t i = 1; i < mp.TotalRealParticles(); ++i)
        if (mpos0[i][0] > mpos0[tail_index][0]) tail_index = i;
    std::ofstream meta("run_parameters.txt");
    meta << std::setprecision(12) << "end_time=" << requested_end
         << "\nwater_spacing=" << water_spacing << "\nsolid_reference_spacing=" << particle_spacing_ref
         << "\nsolid_local_spacing=" << muscle.getSPHAdaptation().ReferenceSpacing()
         << "\nwater_particles=" << wp.TotalRealParticles() << "\nmuscle_particles=" << mp.TotalRealParticles()
         << "\nwater_rho=" << water_rho << "\nwater_mu=" << water_mu << "\nsound_speed=" << water_sound_speed
         << "\ndrive_scale=" << drive_scale << "\nfrequency=" << swim_frequency
         << "\ncore_x=" << core_x_min << "," << core_x_max << "\ncore_y_half=" << core_y_half
         << "\nsponge_width=" << sponge_width << "\nouter_x=" << outer_x_min << "," << outer_x_max
         << "\nouter_y_half=" << outer_y_half << "\nroot=weld\ngravity=0\ninflow=0\n"
         << "boundary=finite reservoir with graded sponge; not validated nonreflecting\n";
    meta.close();
    std::ofstream diag("fsi_diagnostics.csv");
    diag << "time,tail_dx,tail_dy,hydro_fx,hydro_fy,max_water_speed,max_density_error,far_band_speed,max_muscle_speed";
    for (size_t j = 1; j < bone_segments; ++j) diag << ",hinge_" << j;
    diag << '\n' << std::setprecision(12);
    auto diagnostic = [&]() {
        Real speed = 0, rhoerr = 0, far_speed = 0, muscle_speed = 0;
        for (size_t i = 0; i < wp.TotalRealParticles(); ++i) {
            if (!wpos[i].allFinite() || !wvel[i].allFinite() || !std::isfinite(wrho[i]))
                throw std::runtime_error("Non-finite water state");
            if (wpos[i][0] < outer_x_min - 2*water_spacing || wpos[i][0] > outer_x_max + 2*water_spacing ||
                std::abs(wpos[i][1]-cy) > outer_y_half + 2*water_spacing)
                throw std::runtime_error("Water escaped finite reservoir; revise boundary/domain");
            speed = std::max(speed, wvel[i].norm());
            rhoerr = std::max(rhoerr, std::abs(wrho[i]/water_rho-1));
            if (wpos[i][0] < core_x_min || wpos[i][0] > core_x_max || std::abs(wpos[i][1]-cy)>core_y_half)
                far_speed = std::max(far_speed, wvel[i].norm());
        }
        Vecd load = Vecd::Zero();
        for (size_t i = 0; i < mp.TotalRealParticles(); ++i) {
            if (!mpos[i].allFinite() || !mvel[i].allFinite() || !pf[i].allFinite() || !vf[i].allFinite())
                throw std::runtime_error("Non-finite muscle state");
            load += pf[i] + vf[i];
            muscle_speed = std::max(muscle_speed, mvel[i].norm());
        }
        diag << physical_time << ',' << mpos[tail_index][0]-mpos0[tail_index][0] << ','
             << mpos[tail_index][1]-mpos0[tail_index][1] << ',' << load[0] << ',' << load[1]
             << ',' << speed << ',' << rhoerr << ',' << far_speed << ',' << muscle_speed;
        for (auto &pin : pin_mobs) diag << ',' << pin.getAngle(integ.getState());
        diag << '\n'; diag.flush();
        if (std::max(speed, muscle_speed) / water_sound_speed > 0.1 || rhoerr > 0.05)
            throw std::runtime_error("Exploratory WCSPH limit exceeded (Mach>0.1 or density deviation>5%)");
        return speed;
    };
    far_surface.exec();
    sponge.exec(0);
    distributed_attachment_network.build();
    std::cout << "Fixed head, no gravity/inflow. Water particles=" << wp.TotalRealParticles()
              << ", muscle particles=" << mp.TotalRealParticles() << ", target time=" << requested_end << std::endl;
    diagnostic();
    write_states.writeToFile(0);
    Real next_output = std::min(output_interval, requested_end);
    const bool freeze_body = readPositive("FISH_FREEZE_BODY", 0.0, true) == 1.0;
    if (freeze_body && drive_scale != 0.0)
        throw std::runtime_error("Frozen control requires zero drive");
    auto solid_step = [&](Real dt) {
        if (freeze_body) {
            physical_time += dt;
            return;
        }
        impose_active_strain.exec();
        distributed_attachment_force.exec(dt);

        muscle_first_half.exec(dt);
        for (size_t i = 0; i < bone_segments; ++i)
            bone_first_halves[i]->exec(dt);
        for (size_t i = 0; i < bone_segments; ++i)
        {
            constraint_bones[i]->exec();
        }
        constraint_head_muscle_top.exec();
        constraint_head_muscle_bottom.exec();

        muscle_second_half.exec(dt);
        for (size_t i = 0; i < bone_segments; ++i)
            bone_second_halves[i]->exec(dt);
        for (size_t i = 0; i < bone_segments; ++i)
        {
            constraint_bones[i]->exec();
        }
        constraint_head_muscle_top.exec();
        constraint_head_muscle_bottom.exec();

        SimTK::State &state_for_update = integ.updAdvancedState();
        force_on_bodies.clearAllBodyForces(state_for_update);
        for (size_t i = 1; i < bone_segments; ++i)
        {
            auto total_force = force_on_bone_parts[i - 1]->exec() + distributed_attachment_reactions[i - 1]->exec();
            force_on_bodies.setOneBodyForce(state_for_update, pin_mobs[i - 1], total_force);
        }
        integ.stepBy(dt);

        for (size_t i = 0; i < bone_segments; ++i)
        {
            constraint_bones[i]->exec();
        }
        constraint_head_muscle_top.exec();
        constraint_head_muscle_bottom.exec();


        physical_time += dt;
        // Preserve the original baseline's solid neighbor-update policy.
        muscle.updateCellLinkedList();
        muscle_inner.updateConfiguration();
        for (size_t i = 0; i < bone_segments; ++i) {
            bones[i]->updateCellLinkedList();
            bone_inners[i]->updateConfiguration();
        }
    };
    try {
        while (physical_time < requested_end - 1e-12) {
            if ((TickCount::now()-t1).seconds() > wall_limit)
                throw std::runtime_error("Wall-time limit reached BEFORE requested simulation end");
            Real Dt = std::min(fluid_advection_dt.exec(), requested_end-physical_time);
            far_surface.exec();
            density_sum.exec();
            fluid_viscosity.exec();
            transport.exec();
            fluid_viscous_load.exec();
            update_muscle_normal.exec();
            Real advection_elapsed = 0;
            while (advection_elapsed < Dt - 1e-12) {
                Real dt = std::min({fluid_acoustic_dt.exec(), Dt-advection_elapsed,
                                    next_output-physical_time, requested_end-physical_time});
                if (!std::isfinite(dt) || dt <= 0) throw std::runtime_error("Invalid fluid dt");
                pressure_step.exec(dt);
                fluid_pressure_load.exec();
                density_step.exec(dt);
                sponge.exec(dt);
                solid_average.initialize_displacement_.exec();
                Real sub_elapsed = 0;
                while (sub_elapsed < dt - 1e-14) {
                    Real ds = std::min(muscle_dt.exec(), dt-sub_elapsed);
                    for (auto &bdt : bone_dts) ds = std::min(ds, bdt->exec());
                    if (!std::isfinite(ds) || ds <= 0) throw std::runtime_error("Invalid solid dt");
                    solid_step(ds);
                    sub_elapsed += ds;
                }
                solid_average.update_averages_.exec(dt);
                advection_elapsed += dt;
                if (physical_time >= next_output - 1e-12) {
                    Real vmax = diagnostic();
                    vorticity.exec();
                    write_states.writeToFile();
                    write_pin_data.writeToFile(physical_time);
                    std::cout << "Time=" << physical_time << " max_water_speed=" << vmax
                              << " last_hinge=" << pin_mobs.back().getAngle(integ.getState())
                              << " wall_s=" << (TickCount::now()-t1).seconds() << std::endl;
                    next_output = std::min(next_output + output_interval, requested_end);
                }
            }
            water.updateCellLinkedList();
            muscle.updateCellLinkedList();
            water_complex.updateConfiguration();
            muscle_water_contact.updateConfiguration();
            ++ite;
        }
    } catch (const std::exception &e) {
        std::cerr << "SHORT RUN FAILED: " << e.what() << " at t=" << physical_time << std::endl;
        return 2;
    }
    std::cout << "SHORT RUN COMPLETED at t=" << physical_time
              << ", wall_s=" << (TickCount::now()-t1).seconds() << std::endl;
    return 0;
}
