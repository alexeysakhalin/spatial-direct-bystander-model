import math
import numpy as np
COUNTS=["live_targets","live_Ag_positive","live_Ag_negative","live_effectors","live_primed_targets","live_unprimed_targets"]
FIELDS=["IFN-gamma","TNF","CXCL9-11","IL-2"]
PRIMARY_FIELDS=["IFN-gamma","TNF"]
LOCS=["domain_volume_weighted_mean","target_nearest_voxel_mean"]

def finite(x):return x is not None and math.isfinite(float(x))

def area(y):
 if any(x is None for x in y):return None
 assert all(finite(x) for x in y)
 return float(np.trapezoid(y,dx=1) if hasattr(np,'trapezoid') else np.trapz(y,dx=1))

def screen_scalar(value,baseline,mode):
 assert len(baseline)==3,'Three accepted reference realizations required'
 assert all(x is None or finite(x) for x in [value]+baseline),'Nonfinite analysis value'
 if value is None or any(x is None for x in baseline):
  return dict(status='undefined',flag=None,value=value,baseline_seed_values=baseline,reason='Missing live-target sampling is not zero exposure.')
 b0=baseline[0];lo=min(baseline);hi=max(baseline);delta=value-b0
 d=dict(value=value,baseline_seed_values=baseline,baseline_seed0=b0,baseline_min=lo,baseline_max=hi,delta_from_seed0=delta,departure_outside_reference_envelope=max(lo-value,0,value-hi))
 if mode=='count':limit=max(60.,hi-lo)
 else:
  assert mode=='exposure'
  if b0<=1e-6:
   return dict(d,status='relative_screen_undefined_near_zero',flag=None,reason='Baseline <=1e-6; absolute difference retained without a relative adequacy claim.')
  limit=.2*b0;d['relative_delta_from_seed0']=delta/b0
 # Strict "exceeding" tolerance; floating-point equality is not a material departure.
 flag=abs(delta)>limit and not math.isclose(abs(delta),limit,rel_tol=1e-12,abs_tol=1e-12)
 return dict(d,status='flagged' if flag else 'within_reporting_tolerance',flag=flag,absolute_tolerance=limit)

def enrich_counts(rows,validation_rows):
 assert len(rows)==len(validation_rows)==97
 result=[]
 for h,(r,v) in enumerate(zip(rows,validation_rows)):
  assert abs(r['time_h']-h)<1e-7 and abs(v['time_h']-h)<1e-7
  assert r['live_targets']==v['live_targets'] and r['live_effectors']==v['live_effectors']
  row=dict(r);c=v['live_targets_by_type']
  row['live_Ag_positive']=sum(n for k,n in c.items() if 'AgHi' in k)
  row['live_Ag_negative']=sum(n for k,n in c.items() if 'Ag0' in k)
  row['live_primed_targets']=sum(n for k,n in c.items() if k.endswith(' primed'))
  row['live_unprimed_targets']=sum(n for k,n in c.items() if k.endswith(' unprimed'))
  assert row['live_Ag_positive']+row['live_Ag_negative']==row['live_targets']
  assert row['live_primed_targets']+row['live_unprimed_targets']==row['live_targets']
  result.append(row)
 return result

def metric_values(rows):
 assert len(rows)==97 and all(abs(r['time_h']-h)<1e-7 for h,r in enumerate(rows))
 values={}
 for key in COUNTS:
  for h in [24,48,72,96]:values[f'count/{key}/at_{h}h']=(rows[h][key],'count',True)
  values[f'count/{key}/time_average_0_96h']=(area([r[key]for r in rows])/96,'count',True)
 # All hourly counts retained separately; flags are limited to the declared sampled times and burden.
 for sp in FIELDS:
  for loc in LOCS:
   y=[r['fields'][sp][loc]for r in rows]
   for h in [24,48,72,96]:values[f'exposure/{sp}/{loc}/at_{h}h']=(y[h],'exposure',sp in PRIMARY_FIELDS)
   values[f'exposure/{sp}/{loc}/AUC_0_96h']=(area(y),'exposure',sp in PRIMARY_FIELDS)
  for qi,q in enumerate([0,25,50,75,100]):
   for h in [24,48,72,96]:
    quant=rows[h]['fields'][sp]['target_nearest_voxel_quantiles_0_25_50_75_100']
    values[f'exposure/{sp}/target_quantile_{q}/at_{h}h']=(None if quant is None else quant[qi],'exposure',sp in PRIMARY_FIELDS)
 return values

def selftest():
 tests=[]
 def check(name,condition):assert condition,name;tests.append(name)
 check('strict_count_boundary',screen_scalar(160,[100,100,100],'count')['flag'] is False)
 check('count_above_boundary',screen_scalar(161,[100,100,100],'count')['flag'] is True)
 check('baseline_range_controls_tolerance',screen_scalar(179,[100,180,120],'count')['flag'] is False)
 check('negative_count_departure',screen_scalar(39,[100,100,100],'count')['flag'] is True)
 check('exposure_boundary',screen_scalar(1.2,[1,1,1],'exposure')['flag'] is False)
 check('exposure_above_boundary',screen_scalar(1.21,[1,1,1],'exposure')['flag'] is True)
 check('near_zero_not_passed',screen_scalar(1,[0,0,0],'exposure')['flag'] is None)
 check('undefined_target_dose_not_zero',screen_scalar(None,[1,1,1],'exposure')['status']=='undefined')
 for name,value in [('missing_replicate',[1,1]),('nonfinite',[1,float('nan'),1])]:
  try:screen_scalar(1,value,'count')
  except AssertionError:tests.append(name+'_rejected')
  else:raise AssertionError(name+' accepted')
 check('constant_burden_integral',area([10]*97)==960)
 check('missing_AUC_undefined',area([1,None]+[1]*95)is None)
 return tests
