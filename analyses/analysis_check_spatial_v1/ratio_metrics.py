import math
import spatial_metrics as spatial

def metric_values(a,g):
 values=spatial.metrics(a,g)
 for h in spatial.TIMES:
  r=a['rows'][h];n=r['living_effectors']
  for key in ['living_effectors_touching_living_target_spheres','living_effectors_with_saved_live_target_attack_link']:
   count=r[key];assert 0<=count<=n;values[f'access/fraction_{key}/at_{h}h']=count/n if n else None
 assert len(values)==59
 assert all(v is None or math.isfinite(v)for v in values.values())
 return values

def units(key):
 if key.startswith('first_saved_positive_h/'):return 'h'
 if 'fraction_'in key or 'direction_norm'in key:return 'dimensionless'
 if 'gradient_magnitude'in key:return 'relative concentration / um'
 if 'concentration'in key:return 'relative concentration'
 if '_um/'in key:return 'um'
 return 'cells'
