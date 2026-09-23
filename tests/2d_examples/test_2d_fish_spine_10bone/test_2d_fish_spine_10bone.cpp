/**
 * @file test_2d_fish_spine_10bone.cpp
 * @brief A streamlined fish-shaped muscle body wrapping a central spine that is cut into
 *        10 rigid segments connected by Simbody pin (hinge) joints. No water: the head is
 *        welded and the tail undulates in place. Bone-muscle coupling is attachment-only
 *        (each muscle patch is constrained to follow its spine segment; no contact forces).
 */
#include "fish_spine_10bone.h"
#include "sphinxsys.h"
#include <memory>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the SPH system.
    //----------------------------------------------------------------------
    auto parameter = [](const char *name, Real fallback, bool zero) {
        const char *raw = std::getenv(name); if (!raw) return fallback;
        char *end = nullptr; Real x = std::strtod(raw, &end);
        if (end == raw || *end || !std::isfinite(x) || (zero ? x < 0 : x <= 0))
            throw std::runtime_error(std::string("Invalid parameter ") + name);
        return x;
    };
    end_time = parameter("FISH_END_TIME", end_time, false);
    swim_amplitude *= parameter("FISH_DRIVE_SCALE", 1, true);
    const Real dt_scale = parameter("FISH_DT_SCALE", 1, false);
    const Real hinge_damping = parameter("FISH_HINGE_DAMPING", 0.6, true);
    std::cout << "Hinge damping coefficient=" << hinge_damping << std::endl;
    if (dt_scale > 1) throw std::runtime_error("FISH_DT_SCALE must be <= 1");
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
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
            std::make_unique<SimTK::Force::MobilityLinearDamper>(forces, pin_mobs[i], SimTK::MobilizerUIndex(0), hinge_damping));

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
    const Real max_wall_time_seconds = parameter("FISH_WALL_SECONDS", 1800.0, false);

    std::ofstream audit("attachment_audit.csv");
    audit << "time,energy,dissipation,max_power_residual,hinge9\n" << std::setprecision(15);
    write_states.writeToFile(0);
    write_pin_data.writeToFile(0.0);

    //----------------------------------------------------------------------
    //	Main loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        TimeInterval elapsed = TickCount::now() - t1;
        if (elapsed.seconds() >= max_wall_time_seconds)
        {
            std::cout << "Stop simulation due to wall-time limit: " << max_wall_time_seconds << " s." << std::endl;
            break;
        }

        Real dt = muscle_dt.exec();
        for (size_t i = 0; i < bone_segments; ++i)
            dt = SMIN(dt, bone_dts[i]->exec());
        if (!std::isfinite(dt) || dt <= 0.0)
        {
            std::cout << "Invalid dt detected (" << dt << "), stopping." << std::endl;
            break;
        }

        dt = std::min(dt * dt_scale, end_time - physical_time);
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
        ite++;

        muscle.updateCellLinkedList();
        muscle_inner.updateConfiguration();
        for (size_t i = 0; i < bone_segments; ++i)
        {
            bones[i]->updateCellLinkedList();
            bone_inners[i]->updateConfiguration();
        }

        if (ite % 10000 == 0)
        {
            Real angle = pin_mobs.back().getAngle(integ.getState());
            std::cout << std::fixed << std::setprecision(6) << "N=" << ite << "  Time=" << physical_time << "  dt=" << dt
                      << "  last_hinge_angle(rad)=" << angle << "\n";
        }

        static Real next_output_time = output_interval;
        if (physical_time >= next_output_time)
        {
            write_states.writeToFile();
            write_pin_data.writeToFile(physical_time);
            Real energy=0, dissipation=0, residual=0;
            for (size_t j=0;j<muscle.getBaseParticles().TotalRealParticles();++j)
                for (const auto &link : distributed_attachment_network.linksForParticle(j)) {
                    energy += link.cached_energy_; dissipation += link.cached_dissipation_;
                    residual = std::max(residual, std::abs(link.cached_power_residual_));
                }
            audit << physical_time << ',' << energy << ',' << dissipation << ',' << residual
                  << ',' << pin_mobs.back().getAngle(integ.getState()) << '\n';
            audit.flush();
            next_output_time += output_interval;
        }
    }

    TimeInterval tt = TickCount::now() - t1;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    return physical_time >= end_time - 1e-10 ? 0 : 2;
}
