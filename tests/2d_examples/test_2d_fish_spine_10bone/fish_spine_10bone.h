#ifndef FISH_SPINE_10BONE_H
#define FISH_SPINE_10BONE_H

#include "2d_fish_and_bones.h"
#include "active_model.h"
#include "sphinxsys.h"

#include "attachment_law.h"
#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <vector>

using namespace SPH;

//----------------------------------------------------------------------
// Basic geometry and numerical setup.
//  - One streamlined fish-shaped muscle body (outer outline from CreatFishShape).
//  - A central horizontal "spine" cut into 10 rigid segments (Simbody pin chain).
//  - No water: in-place undulation only (architecture stays water-friendly).
//----------------------------------------------------------------------
Real particle_spacing_ref = 0.00125;
Real BW = 4.0 * particle_spacing_ref;

Real cx = 0.0; /**< Fish head tip x (CreatFishShape headtip). */
Real cy = 0.0; /**< Fish mid-axis y. */
Real fish_length = 0.2;
Real fish_shape_resolution = particle_spacing_ref * 0.5;

static constexpr size_t bone_segments = 10;

/** Spine spans the thick central body region only (the fish tapers to a point at both
 *  tips, so the bone cannot reach the extremities while keeping muscle on both sides). */
Real spine_x_start = cx + 0.04;
Real spine_x_end = cx + 0.13;
Real bone_half = 0.002;   /**< Bone half-thickness; bone thickness 0.004 ~ 3 particles. */
Real attach_half = 0.001; /**< Half-thickness of the muscle attachment patch on each side. */

Real output_interval = 0.01;
Real end_time = 2.0;

//----------------------------------------------------------------------
// Material properties.
//----------------------------------------------------------------------
Real rho0_s = 1050.0;
Real poisson = 0.49;
Real young_muscle = 0.5e6; /**< Active outer muscle shell. */
Real young_bone = 1.1e6;   /**< Rigid spine segments (elastic backup material). */
Real young_buffer = 0.3e6; /**< Passive buffer/connective tissue between shell and bone. */

/** Layering of the soft-tissue body (one composite body, fish outline minus spine slot):
 *  - material_id 0 = active muscle, a thin shell within shell_thickness of the skin;
 *  - material_id 1 = passive buffer, everything else (and the whole head region).
 *  A minimum buffer gap is always kept between the bone and the active shell, so the
 *  bone-side attachment patches stay in the buffer and the shell is never hard-constrained. */
Real shell_thickness = 0.004;     /**< Active muscle shell thickness (~3 particle rows). */
Real min_buffer_gap = 0.002;      /**< Minimum passive buffer kept between bone and shell. */
Real head_passive_length = 0.05;  /**< x < cx + this is passive only (thin welded head). */

/** Compliant myosepta-like attachment from passive buffer particles to Simbody bones.
 *  omega/zeta are converted per particle to k_i = m_i * omega^2 and c_i = 2*zeta*m_i*omega. */
Real omega_attach = 2.0 * Pi * 40.0;
Real zeta_attach = 0.8;
Real attachment_cutoff_radius = 2.0 * particle_spacing_ref;
size_t max_attachment_links_per_buffer_particle = 2;
size_t min_attach_links_per_segment_for_warning = 3;

//----------------------------------------------------------------------
// Active-strain traveling-wave drive (sin^2 form, dorsal/ventral phase shift).
//  Faithful to test_2d_flow_stream_around_fish, but the amplitude envelope is
//  tail-biased here because the head (x = cx) is welded and the tail is free.
//----------------------------------------------------------------------
Real swim_amplitude = 0.12;             /**< Am: global active-strain gain. */
Real swim_frequency = 4.0;              /**< Beat frequency (Hz). */
Real swim_wave_length = 3.0 * fish_length;
Real swim_start_time = 0.2;             /**< Soft-start time constant. */

BoundingBox system_domain_bounds(
    Vec2d(cx - 0.03, cy - 0.03),
    Vec2d(cx + fish_length + 0.03, cy + 0.03));

//----------------------------------------------------------------------
// Case geometry helpers.
//----------------------------------------------------------------------
inline std::vector<Vecd> createRectangleShape(Real x_min, Real x_max, Real y_min, Real y_max)
{
    std::vector<Vecd> shape;
    shape.push_back(Vecd(x_min, y_min));
    shape.push_back(Vecd(x_min, y_max));
    shape.push_back(Vecd(x_max, y_max));
    shape.push_back(Vecd(x_max, y_min));
    shape.push_back(Vecd(x_min, y_min));
    return shape;
}

inline Real segmentLength() { return (spine_x_end - spine_x_start) / Real(bone_segments); }
inline Real segmentStart(size_t segment_index) { return spine_x_start + Real(segment_index) * segmentLength(); }
inline Real segmentEnd(size_t segment_index) { return spine_x_start + Real(segment_index + 1) * segmentLength(); }
inline Vec2d hingePoint(size_t hinge_index) { return Vec2d(segmentStart(hinge_index), cy); }
inline Real clampUnit(Real value) { return value < 0.0 ? 0.0 : (value > 1.0 ? 1.0 : value); }

//----------------------------------------------------------------------
// Bodies.
//----------------------------------------------------------------------
/** Fish-shaped muscle = fish outline minus the central spine slot, so the muscle and
 *  the bone bodies never overlap (coupling is purely through attachment constraints). */
class FishMuscleBody : public MultiPolygonShape
{
  public:
    explicit FishMuscleBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
        /** Mirror in x so the thick "head" sits on the left (welded segment 0) and the
         *  slender "tail" on the right (free, large-amplitude end). */
        for (auto &p : fish_shape)
            p[0] = 2.0 * cx + fish_length - p[0];
        std::reverse(fish_shape.begin(), fish_shape.end()); // keep polygon winding consistent
        multi_polygon_.addAPolygon(fish_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(
            createRectangleShape(spine_x_start, spine_x_end, cy - bone_half, cy + bone_half),
            ShapeBooleanOps::sub);
    }
};

/** Composite soft tissue: material_id 0 = active muscle shell, material_id 1 = passive buffer. */
class FishTissueComposite : public CompositeSolid
{
  public:
    FishTissueComposite() : CompositeSolid(rho0_s)
    {
        add<ActiveModelSolid>(rho0_s, young_muscle, poisson);          // material_id 0: active muscle
        add<SaintVenantKirchhoffSolid>(rho0_s, young_buffer, poisson); // material_id 1: passive buffer
    }
};

/** A single rigid spine segment (rectangle fully inside the fish belly). */
class BoneSegmentBody : public MultiPolygonShape
{
  public:
    BoneSegmentBody(const std::string &shape_name, size_t segment_index) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(
            createRectangleShape(segmentStart(segment_index), segmentEnd(segment_index), cy - bone_half, cy + bone_half),
            ShapeBooleanOps::add);
    }
};

inline SharedPtr<GeometricShapeBox> createBoneRigidShape(size_t segment_index)
{
    Real x_start = segmentStart(segment_index);
    Real x_end = segmentEnd(segment_index);
    Vec2d center(0.5 * (x_start + x_end), cy);
    Vec2d half(0.5 * (x_end - x_start) * 0.98, bone_half * 0.98);
    return makeShared<GeometricShapeBox>(Transform(center), half, "BoneRigid_" + std::to_string(segment_index));
}

/** Thin muscle patch just above the spine within a segment, tied to that segment. */
inline SharedPtr<GeometricShapeBox> createMuscleAttachmentShape(size_t segment_index)
{
    Real x_start = segmentStart(segment_index);
    Real x_end = segmentEnd(segment_index);
    Vec2d center(0.5 * (x_start + x_end), cy + bone_half + attach_half);
    Vec2d half(0.5 * (x_end - x_start) * 0.9, attach_half);
    return makeShared<GeometricShapeBox>(Transform(center), half, "MuscleAttachTop_" + std::to_string(segment_index));
}

/** Thin muscle patch just below the spine within a segment, tied to that segment. */
inline SharedPtr<GeometricShapeBox> createMuscleAttachmentShapeBottom(size_t segment_index)
{
    Real x_start = segmentStart(segment_index);
    Real x_end = segmentEnd(segment_index);
    Vec2d center(0.5 * (x_start + x_end), cy - bone_half - attach_half);
    Vec2d half(0.5 * (x_end - x_start) * 0.9, attach_half);
    return makeShared<GeometricShapeBox>(Transform(center), half, "MuscleAttachBot_" + std::to_string(segment_index));
}

//----------------------------------------------------------------------
// Tag material ids by position: active muscle shell (0) vs passive buffer (1).
//  Shell = behind the passive head, within shell_thickness of the skin, AND at least
//  min_buffer_gap away from the bone band; everything else is buffer.
//----------------------------------------------------------------------
class FishMaterialInitialization : public MaterialIdInitialization
{
  public:
    explicit FishMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        Real x = pos_[index_i][0] - cx; /**< 0 at head (left), fish_length at tail (right). */
        Real y = pos_[index_i][1];
        /** Geometry is mirrored in x, so the local skin half-thickness is sampled from the tail tip. */
        Real half_thickness = outline(fish_length - x, 0.03, fish_length);

        bool behind_head = (x >= head_passive_length);
        bool is_shell = false;
        if (y > cy) // dorsal side
        {
            bool near_skin = y > (cy + half_thickness - shell_thickness);
            bool clear_of_bone = y > (cy + bone_half + min_buffer_gap);
            is_shell = behind_head && near_skin && clear_of_bone;
        }
        else // ventral side
        {
            bool near_skin = y < (cy - half_thickness + shell_thickness);
            bool clear_of_bone = y < (cy - bone_half - min_buffer_gap);
            is_shell = behind_head && near_skin && clear_of_bone;
        }

        material_id_[index_i] = is_shell ? 0 : 1;
    };
};

//----------------------------------------------------------------------
// Imposing active strain on the fish muscle shell (sin^2 traveling wave).
//  Only material_id 0 (the active shell) is driven; the buffer (1) is fully passive.
//  Dorsal (y > cy) and ventral (y < cy) contraction waves are half a period apart,
//  so one contracts while the other relaxes -> the spine chain bends side to side.
//----------------------------------------------------------------------
class ImposingActiveStrain : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit ImposingActiveStrain(SolidBody &solid_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(solid_body),
          material_id_(particles_->getVariableDataByName<int>("MaterialID")),
          pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
          active_strain_(particles_->getVariableDataByName<Matd>("ActiveStrain")),
          physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        active_strain_[index_i] = Matd::Zero();
        if (material_id_[index_i] != 0)
            return; // only the active muscle shell is driven; buffer stays passive

        Real x = pos0_[index_i][0] - cx; /**< 0 at head, fish_length at tail. */
        Real y = pos0_[index_i][1];

        Real w = 2.0 * Pi * swim_frequency;
        Real wave_number = 2.0 * Pi / swim_wave_length;

        /** Tail-biased envelope: ~0 at the welded head, max at the free tail. */
        Real x_hat = clampUnit(x / fish_length);
        Real envelope = x_hat * x_hat;

        Real current_time = *physical_time_;
        Real strength = 1.0 - exp(-current_time / swim_start_time);

        Real phase_shift = (y > cy) ? 0.0 : 0.5 * Pi;
        Real wave = sin(0.5 * w * current_time - 0.5 * wave_number * x + phase_shift);

        active_strain_[index_i](0, 0) = -swim_amplitude * envelope * strength * wave * wave;
    }

  protected:
    int *material_id_;
    Vecd *pos0_;
    Matd *active_strain_;
    Real *physical_time_;
};

//------------------------------------------------------------------------------
// Distributed spring-damper links from passive buffer particles to nearby bone-surface particles.
//------------------------------------------------------------------------------
struct DistributedAttachmentLink
{
    size_t buffer_index_;
    size_t segment_index_;
    Vecd bone_station_initial_;
    Vecd buffer_station_initial_;
    Vecd cached_force_ = Vecd::Zero();
    Real cached_torque_ = 0;
    Real cached_energy_ = 0, cached_dissipation_ = 0, cached_power_residual_ = 0;
};

class DistributedAttachmentNetwork
{
  public:
    DistributedAttachmentNetwork(SolidBody &muscle_body, const std::vector<SolidBody *> &bone_bodies)
        : muscle_body_(muscle_body), bone_bodies_(bone_bodies), links_are_built_(false) {}

    void build()
    {
        if (links_are_built_)
            return;

        BaseParticles &muscle_particles = muscle_body_.getBaseParticles();
        Vecd *muscle_pos0 = muscle_particles.getVariableDataByName<Vecd>("Position");
        int *material_id = muscle_particles.getVariableDataByName<int>("MaterialID");
        size_t total_muscle_particles = muscle_particles.TotalRealParticles();
        links_by_buffer_particle_.clear();
        links_by_buffer_particle_.resize(total_muscle_particles);
        links_per_segment_.assign(bone_segments, 0);

        std::vector<std::vector<Vecd>> bone_surface_points(bone_segments);
        for (size_t segment_index = 1; segment_index < bone_segments; ++segment_index)
        {
            BaseParticles &bone_particles = bone_bodies_[segment_index]->getBaseParticles();
            Vecd *bone_pos0 = bone_particles.getVariableDataByName<Vecd>("Position");
            for (size_t j = 0; j < bone_particles.TotalRealParticles(); ++j)
            {
                Real distance_to_top_bottom = SMIN(
                    std::abs(bone_pos0[j][1] - (cy + bone_half)),
                    std::abs(bone_pos0[j][1] - (cy - bone_half)));
                Real distance_to_ends = SMIN(
                    std::abs(bone_pos0[j][0] - segmentStart(segment_index)),
                    std::abs(bone_pos0[j][0] - segmentEnd(segment_index)));
                if (SMIN(distance_to_top_bottom, distance_to_ends) <= 0.75 * particle_spacing_ref)
                    bone_surface_points[segment_index].push_back(bone_pos0[j]);
            }
        }

        Real cutoff_squared = attachment_cutoff_radius * attachment_cutoff_radius;
        size_t linked_buffer_particles = 0;
        size_t total_links = 0;
        for (size_t i = 0; i < total_muscle_particles; ++i)
        {
            if (material_id[i] != 1)
                continue;

            std::array<Real, 2> nearest_distance = {
                std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()};
            std::array<DistributedAttachmentLink, 2> nearest_link{};

            for (size_t segment_index = 1; segment_index < bone_segments; ++segment_index)
            {
                for (const Vecd &bone_station : bone_surface_points[segment_index])
                {
                    Real distance_squared = (muscle_pos0[i] - bone_station).squaredNorm();
                    if (distance_squared > cutoff_squared || distance_squared >= nearest_distance[1])
                        continue;

                    DistributedAttachmentLink candidate{i, segment_index, bone_station, muscle_pos0[i]};
                    if (distance_squared < nearest_distance[0])
                    {
                        nearest_distance[1] = nearest_distance[0];
                        nearest_link[1] = nearest_link[0];
                        nearest_distance[0] = distance_squared;
                        nearest_link[0] = candidate;
                    }
                    else
                    {
                        nearest_distance[1] = distance_squared;
                        nearest_link[1] = candidate;
                    }
                }
            }

            size_t links_for_particle = 0;
            for (size_t n = 0; n < max_attachment_links_per_buffer_particle && n < nearest_distance.size(); ++n)
            {
                if (nearest_distance[n] == std::numeric_limits<Real>::max())
                    continue;
                links_by_buffer_particle_[i].push_back(nearest_link[n]);
                links_per_segment_[nearest_link[n].segment_index_]++;
                total_links++;
                links_for_particle++;
            }
            if (links_for_particle > 0)
                linked_buffer_particles++;
        }

        std::cout << "Distributed attachment diagnostics: linked buffer particles = "
                  << linked_buffer_particles << ", total spring links = " << total_links << std::endl;
        for (size_t segment_index = 1; segment_index < bone_segments; ++segment_index)
        {
            std::cout << "  segment " << segment_index << " links = " << links_per_segment_[segment_index] << std::endl;
            if (links_per_segment_[segment_index] == 0)
            {
                std::cout << "WARNING: segment " << segment_index
                          << " has zero distributed attachment links." << std::endl;
            }
            else if (links_per_segment_[segment_index] < min_attach_links_per_segment_for_warning)
            {
                std::cout << "WARNING: segment " << segment_index << " has only "
                          << links_per_segment_[segment_index]
                          << " distributed attachment links; coupling may be weak." << std::endl;
            }
        }

        links_are_built_ = true;
    }

    std::vector<DistributedAttachmentLink> &linksForParticle(size_t index_i)
    { return links_by_buffer_particle_[index_i]; }

    const std::vector<DistributedAttachmentLink> &linksForParticle(size_t index_i) const
    {
        return links_by_buffer_particle_[index_i];
    }

  protected:
    SolidBody &muscle_body_;
    std::vector<SolidBody *> bone_bodies_;
    std::vector<std::vector<DistributedAttachmentLink>> links_by_buffer_particle_;
    std::vector<size_t> links_per_segment_;
    bool links_are_built_;
};

class DistributedAttachmentForce : public BaseForcePrior<SPHBody>
{
  public:
    DistributedAttachmentForce(SolidBody &muscle_body, DistributedAttachmentNetwork &attachment_network,
                               SimTK::MultibodySystem &MBsystem, std::vector<SimTK::MobilizedBody::Pin> &pin_mobs,
                               SimTK::RungeKuttaMersonIntegrator &integ)
        : BaseForcePrior<SPHBody>(muscle_body, "AttachmentForce"),
          attachment_network_(attachment_network), MBsystem_(MBsystem), pin_mobs_(pin_mobs), integ_(integ),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          material_id_(particles_->getVariableDataByName<int>("MaterialID")),
          attachment_displacement_(particles_->registerStateVariableData<Vecd>("AttachmentDisplacement")),
          attachment_displacement_norm_(particles_->registerStateVariableData<Real>("AttachmentDisplacementNorm")),
          attachment_force_norm_(particles_->registerStateVariableData<Real>("AttachmentForceNorm")),
          simbody_states_(bone_segments),
          initial_origin_locations_(bone_segments, Vec3d::Zero())
    {
        const SimTK::State *state = &integ_.getState();
        MBsystem_.realize(*state, SimTK::Stage::Velocity);
        for (size_t segment_index = 1; segment_index < bone_segments; ++segment_index)
        {
            SimTK::MobilizedBody &mobod = pin_mobs_[segment_index - 1];
            initial_origin_locations_[segment_index] = SimTKToEigen(mobod.getBodyOriginLocation(*state));
        }
    }

    void setupDynamics(Real dt = 0.0)
    {
        attachment_network_.build();
        const SimTK::State *state = &integ_.getState();
        MBsystem_.realize(*state, SimTK::Stage::Velocity);
        updateSimbodyStates(*state);
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd total_force = Vecd::Zero();
        Vecd total_displacement = Vecd::Zero();

        if (material_id_[index_i] == 1)
        {
            for (DistributedAttachmentLink &link : attachment_network_.linksForParticle(index_i))
            {
                Vecd displacement, link_force;
                computeLinkKinematicsAndForce(link, displacement, link_force);
                total_displacement += displacement;
                total_force += link_force;
            }
        }

        current_force_[index_i] = total_force;
        attachment_displacement_[index_i] = total_displacement;
        attachment_displacement_norm_[index_i] = total_displacement.norm();
        attachment_force_norm_[index_i] = total_force.norm();
        BaseForcePrior<SPHBody>::update(index_i, dt);
    }

    void computeLinkKinematicsAndForce(DistributedAttachmentLink &link, Vecd &displacement, Vecd &link_force)
    {
        Vec3d target_pos, target_vel, target_acc, target_normal;
        simbody_states_[link.segment_index_].findStationLocationVelocityAndAccelerationInGround(
            upgradeToVec3d(link.buffer_station_initial_), Vec3d::UnitY(), target_pos, target_vel, target_acc, target_normal);

        displacement = degradeToVecd(target_pos) - pos_[link.buffer_index_];
        Vecd velocity_difference = degradeToVecd(target_vel) - vel_[link.buffer_index_];
        Real particle_stiffness = mass_[link.buffer_index_] * omega_attach * omega_attach;
        Real particle_damping = 2.0 * zeta_attach * mass_[link.buffer_index_] * omega_attach;
        const auto &bs = simbody_states_[link.segment_index_];
        const Vecd &x = pos_[link.buffer_index_], &v = vel_[link.buffer_index_];
        auto response = attachmentLaw(x[0], x[1], v[0], v[1], target_pos[0], target_pos[1],
            bs.origin_location_[0], bs.origin_location_[1], bs.origin_velocity_[0],
            bs.origin_velocity_[1], bs.angular_velocity_[2], particle_stiffness, particle_damping);
        link_force = Vecd(response.fx, response.fy);
        link.cached_force_ = link_force;
        link.cached_torque_ = response.torque;
        link.cached_energy_ = response.energy;
        link.cached_dissipation_ = response.dissipation;
        link.cached_power_residual_ = response.power_residual;
    }

  protected:
    DistributedAttachmentNetwork &attachment_network_;
    SimTK::MultibodySystem &MBsystem_;
    std::vector<SimTK::MobilizedBody::Pin> &pin_mobs_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    Vecd *pos_, *vel_;
    Real *mass_;
    int *material_id_;
    Vecd *attachment_displacement_;
    Real *attachment_displacement_norm_, *attachment_force_norm_;
    std::vector<SimbodyState> simbody_states_;
    std::vector<Vec3d> initial_origin_locations_;

    void updateSimbodyStates(const SimTK::State &state)
    {
        for (size_t segment_index = 1; segment_index < bone_segments; ++segment_index)
        {
            SimTK::MobilizedBody &mobod = pin_mobs_[segment_index - 1];
            simbody_states_[segment_index].origin_location_ = SimTKToEigen(mobod.getBodyOriginLocation(state));
            simbody_states_[segment_index].origin_velocity_ = SimTKToEigen(mobod.getBodyOriginVelocity(state));
            simbody_states_[segment_index].angular_velocity_ = SimTKToEigen(mobod.getBodyAngularVelocity(state));
            simbody_states_[segment_index].rotation_ = SimTKToEigen(mobod.getBodyRotation(state));
            simbody_states_[segment_index].initial_origin_location_ = initial_origin_locations_[segment_index];
        }
    }
};

class DistributedAttachmentReactionForSimBody
    : public BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, SPHBody>
{
  public:
    DistributedAttachmentReactionForSimBody(SolidBody &muscle_body, DistributedAttachmentNetwork &attachment_network,
                                            SimTK::MultibodySystem &MBsystem, SimTK::MobilizedBody &mobod,
                                            SimTK::RungeKuttaMersonIntegrator &integ, size_t segment_index)
        : BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, SPHBody>(muscle_body),
          attachment_network_(attachment_network), MBsystem_(MBsystem), mobod_(mobod), integ_(integ),
          segment_index_(segment_index),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          material_id_(particles_->getVariableDataByName<int>("MaterialID"))
    {
        this->quantity_name_ = "DistributedAttachmentReactionForSimBody";
        const SimTK::State *state = &integ_.getState();
        MBsystem_.realize(*state, SimTK::Stage::Velocity);
        initial_origin_location_ = SimTKToEigen(mobod_.getBodyOriginLocation(*state));
    }

    void setupDynamics(Real dt = 0.0)
    {
        attachment_network_.build();
        const SimTK::State *state = &integ_.getState();
        MBsystem_.realize(*state, SimTK::Stage::Velocity);
        simbody_state_.origin_location_ = SimTKToEigen(mobod_.getBodyOriginLocation(*state));
        simbody_state_.origin_velocity_ = SimTKToEigen(mobod_.getBodyOriginVelocity(*state));
        simbody_state_.angular_velocity_ = SimTKToEigen(mobod_.getBodyAngularVelocity(*state));
        simbody_state_.rotation_ = SimTKToEigen(mobod_.getBodyRotation(*state));
        simbody_state_.initial_origin_location_ = initial_origin_location_;
        current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*state);
    }

    SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0)
    {
        if (material_id_[index_i] != 1)
            return SimTK::SpatialVec(SimTKVec3(0), SimTKVec3(0));

        SimTKVec3 total_reaction_force(0);
        SimTKVec3 total_reaction_torque(0);
        for (const DistributedAttachmentLink &link : attachment_network_.linksForParticle(index_i))
        {
            if (link.segment_index_ != segment_index_)
                continue;

            // Reuse the SAME per-link load sampled before the tissue step.
            total_reaction_force -= EigenToSimTK(upgradeToVec3d(link.cached_force_));
            total_reaction_torque += SimTKVec3(0, 0, link.cached_torque_);
        }
        return SimTK::SpatialVec(total_reaction_torque, total_reaction_force);
    }

  protected:
    DistributedAttachmentNetwork &attachment_network_;
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    size_t segment_index_;
    SimbodyState simbody_state_;
    Vec3d initial_origin_location_;
    SimTKVec3 current_mobod_origin_location_;
    Vecd *pos_, *vel_;
    Real *mass_;
    int *material_id_;
};

#endif // FISH_SPINE_10BONE_H
