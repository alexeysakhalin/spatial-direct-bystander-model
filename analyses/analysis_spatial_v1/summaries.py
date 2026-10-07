import math,statistics

def descriptive(vals):
 a=[v for v in vals if v is not None];assert all(math.isfinite(v) for v in a)
 return {'n_available':len(a),'n_expected':len(vals),'mean':statistics.fmean(a) if a else None,'min':min(a) if a else None,'max':max(a) if a else None,'is_confidence_interval':False}

def difference(a,b):
 # Signed TE5 minus TE2; a zero denominator never produces an invented percentage.
 valid=a is not None and b is not None
 return {'TE2':a,'TE5':b,'TE5_minus_TE2':b-a if valid else None,'relative_to_TE2_percent':100*(b-a)/a if valid and a!=0 else None,'relative_difference_unavailable_reason':'missing_input' if not valid else 'zero_TE2_denominator' if a==0 else None}
