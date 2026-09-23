
#ifndef FIXED_HEAD_WATER_H
#define FIXED_HEAD_WATER_H
#include "fish_spine_10bone.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// SI assumptions for this exploratory case: length m, time s, mass kg (per unit
// out-of-plane thickness in 2D). No freestream speed is imposed.
const Real water_rho = 1000.0;
const Real water_mu = 0.001;
const Real water_sound_speed = 10.0;
const Real water_reference_speed = 0.5; // time-step estimate only
const Real sponge_width = 0.5 * fish_length;
const Real core_x_min = cx - fish_length;
const Real core_x_max = cx + 4.0 * fish_length;
const Real core_y_half = 1.5 * fish_length + 0.03;
const Real outer_x_min = core_x_min - sponge_width;
const Real outer_x_max = core_x_max + sponge_width;
const Real outer_y_half = core_y_half + sponge_width;

Real readPositive(const char *name, Real fallback, bool allow_zero = false)
{
    const char *raw = std::getenv(name);
    if (!raw) return fallback;
    char *end = nullptr;
    Real value = std::strtod(raw, &end);
    if (end == raw || *end != '\0' || !std::isfinite(value) ||
        (allow_zero ? value < 0 : value <= 0))
        throw std::runtime_error(std::string("Invalid parameter: ") + name);
    return value;
}

class SurroundingWater : public MultiPolygonShape
{
  public:
    explicit SurroundingWater(const std::string &name) : MultiPolygonShape(name)
    {
        multi_polygon_.addAPolygon(createRectangleShape(outer_x_min, outer_x_max,
            cy - outer_y_half, cy + outer_y_half), ShapeBooleanOps::add);
        auto outline_points = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
        for (auto &p : outline_points) p[0] = 2.0 * cx + fish_length - p[0];
        std::reverse(outline_points.begin(), outline_points.end());
        // Subtract the COMPLETE outer fish, not the soft tissue with its spine slot.
        multi_polygon_.addAPolygon(outline_points, ShapeBooleanOps::sub);
    }
};

// Finite ambient reservoir + graded sponge. This is an APPROXIMATE far-field
// truncation, not an exact radiation boundary and not a physical free surface.
// No damping in the fish observation region. Both velocity components and
// density perturbations are exponentially relaxed only in the perimeter band.
class FarFieldSponge : public LocalDynamics
{
    Vecd *pos_, *vel_;
    Real *rho_, *pressure_, *weight_;
  public:
    explicit FarFieldSponge(SPHBody &body) : LocalDynamics(body),
        pos_(particles_->getVariableDataByName<Vecd>("Position")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
        rho_(particles_->getVariableDataByName<Real>("Density")),
        pressure_(particles_->getVariableDataByName<Real>("Pressure")),
        weight_(particles_->registerStateVariableData<Real>("SpongeWeight")) {}
    void update(size_t i, Real dt = 0.0)
    {
        Real penetration = std::max({core_x_min - pos_[i][0], pos_[i][0] - core_x_max,
                                    std::abs(pos_[i][1] - cy) - core_y_half, Real(0)});
        Real q = clampUnit(penetration / sponge_width);
        weight_[i] = q * q;
        Real factor = std::exp(-10.0 * water_sound_speed / sponge_width * q * q * dt);
        vel_[i] *= factor;
        rho_[i] = water_rho + (rho_[i] - water_rho) * factor;
        pressure_[i] = water_sound_speed * water_sound_speed * (rho_[i] - water_rho);
    }
};
#endif
