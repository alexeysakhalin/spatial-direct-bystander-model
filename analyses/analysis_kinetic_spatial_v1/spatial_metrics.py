import math
TIMES=[6,24,48,72,96]

def metrics(a,g):
 assert a['id']==g['id']and len(a['rows'])==len(g['rows'])==97
 out={}
 for h in TIMES:
  ar=a['rows'][h];gr=g['rows'][h];assert ar['time_h']==gr['time_h']==h
  assert ar['living_targets']==gr['counts']['live_targets']and ar['living_effectors']==gr['counts']['live_effectors']
  for key in ['living_effectors_inside_initial_target_hull','fraction_living_effectors_inside_initial_target_hull','living_effectors_touching_living_target_spheres','living_effectors_with_saved_live_target_attack_link']:
   out[f'access/{key}/at_{h}h']=ar[key]
  q=ar['nearest_living_target_centre_distance_quantiles_um'];out[f'access/nearest_living_target_centre_median_um/at_{h}h']=q[2]if q is not None else None
  for key in ['living_effector_local_concentration_quantiles','living_effector_local_gradient_magnitude_quantiles_per_um','living_effector_direction_norm_quantiles']:
   q=gr[key];out[f'CXCL/median_{key.replace("_quantiles","")}/at_{h}h']=q[2]if q is not None else None
  n=gr['counts']['live_effectors'];out[f'CXCL/fraction_direction_norm_ge_0_99/at_{h}h']=gr['effectors_with_direction_norm_at_least_0_99']/n if n else None
 for key in ['living_effectors_inside_initial_target_hull','living_effectors_touching_living_target_spheres','living_effectors_with_saved_live_target_attack_link','living_effectors_with_positive_recent_activity_memory']:
  out['first_saved_positive_h/'+key]=a['first_saved_positive_hours'][key]
 return out
