#ifndef IKM_ACTIVATION_MEMORY_H
#define IKM_ACTIVATION_MEMORY_H
#include <cmath>
#include <algorithm>
#include <stdexcept>

// Optional project model extension. Does not modify delivered-damage history.
// A positive half-life defines an exposure-state ODE driven by damage increments.
// -1 retains the legacy persistent signal; 0 gates that signal on recent delivery.
namespace PhysiCell {
inline double ikm_recent_activity_step(double recent, double before, double delivered,
                                       double dt, double half_life)
{
    if(!std::isfinite(recent) || !std::isfinite(before) || !std::isfinite(delivered)
       || !std::isfinite(dt) || !std::isfinite(half_life) || recent < 0 || before < 0
       || delivered < 0 || dt <= 0 || half_life < -1 || (half_life < 0 && half_life != -1))
        throw std::runtime_error("Invalid IKM activation-memory state or parameter");
    const double tolerance=1e-10*std::max(1.0,std::max(before,delivered));
    if(delivered < before-tolerance)
        throw std::runtime_error("Unexpected delivered-damage reset in active IKM cell");
    const double increment=std::max(0.0,delivered-before);
    if(half_life == -1) return delivered;
    if(half_life == 0) return increment > tolerance ? delivered : 0.0;
    const double x=std::log(2.0)*dt/half_life;
    // Exact update for constant delivery rate within this phenotype interval.
    return recent*std::exp(-x) + increment*(-std::expm1(-x)/x);
}

inline bool ikm_update_activation_memory(Cell* cell, double dt)
{
    const int enabled=cell->custom_data.find_variable_index("ikm_activation_enabled");
    if(enabled < 0 || cell->custom_data[enabled] < 0.5) return false;
    const int signal=cell->custom_data.find_variable_index("ikm_recent_activity");
    const int previous=cell->custom_data.find_variable_index("ikm_previous_delivered");
    const int half=cell->custom_data.find_variable_index("ikm_activation_half_life_min");
    if(signal < 0 || previous < 0 || half < 0)
        throw std::runtime_error("Incomplete IKM activation-memory configuration");
    const double delivered=cell->phenotype.cell_interactions.total_damage_delivered;
    if(cell->phenotype.death.dead) cell->custom_data[signal]=0.0;
    else cell->custom_data[signal]=ikm_recent_activity_step(cell->custom_data[signal],
        cell->custom_data[previous],delivered,dt,cell->custom_data[half]);
    cell->custom_data[previous]=delivered;
    return true;
}
}
#endif
