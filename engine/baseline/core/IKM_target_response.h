#ifndef IKM_TARGET_RESPONSE_H
#define IKM_TARGET_RESPONSE_H
#include "PhysiCell_signal_behavior.h"
#include "IKM_target_response_math.h"
namespace PhysiCell {
inline double& ikm_response_value(Cell* cell,const char* name) {
    int idx=cell->custom_data.find_variable_index(name);
    if(idx<0) throw std::runtime_error("Incomplete IKM target-response configuration");
    return cell->custom_data[idx];
}
// Must run once per phenotype interval, AFTER the native rules reset behaviors.
inline bool ikm_apply_target_response(Cell* cell,double dt) {
    int enabled=cell->custom_data.find_variable_index("ikm_target_response_enabled");
    if(enabled<0 || cell->custom_data[enabled]<.5) return false;
    if(cell->phenotype.death.dead) {
        ikm_response_value(cell,"ikm_combo_added_hazard")=0;
        return true;
    }
    const double gamma=ikm_hill_occupancy(get_single_signal(cell,"IFN-gamma"),
        ikm_response_value(cell,"ikm_ifng_half_response"),
        ikm_response_value(cell,"ikm_response_hill_power"));
    const double tnf=ikm_hill_occupancy(get_single_signal(cell,"TNF"),
        ikm_response_value(cell,"ikm_tnf_half_response"),
        ikm_response_value(cell,"ikm_response_hill_power"));
    double strength=ikm_response_value(cell,"ikm_cytostasis_strength");
    double primed=ikm_response_value(cell,"ikm_ifng_primed_state");
    double competence=ikm_response_value(cell,"ikm_combo_death_competence");
    ikm_unit_interval(strength); ikm_unit_interval(primed); ikm_unit_interval(competence);
    double maximum=ikm_response_value(cell,"ikm_combo_hazard_max");
    if(!std::isfinite(maximum) || maximum<0) throw std::runtime_error("Invalid combo hazard");
    double& first=ikm_response_value(cell,"ikm_combo_stage1");
    double& second=ikm_response_value(cell,"ikm_combo_stage2");
    const auto next=ikm_two_stage_response(first,second,primed*tnf,dt,
        ikm_response_value(cell,"ikm_combo_stage_tau_min"));
    first=next.first; second=next.second;
    const double added=competence*maximum*next.mean_second;
    const double multiplier=1-strength*gamma;
    const double native_cycle=get_single_behavior(cell,"cycle entry");
    const double native_death=get_single_behavior(cell,"apoptosis");
    if(!std::isfinite(native_cycle) || !std::isfinite(native_death) || native_cycle<0 || native_death<0)
        throw std::runtime_error("Invalid native target behavior");
    set_single_behavior(cell,"cycle entry",native_cycle*multiplier);
    set_single_behavior(cell,"apoptosis",native_death+added);
    ikm_response_value(cell,"ikm_growth_multiplier")=multiplier;
    ikm_response_value(cell,"ikm_combo_added_hazard")=added;
    return true;
}
}
#endif
