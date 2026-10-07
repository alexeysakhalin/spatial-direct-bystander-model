import numpy as np

def stats(vals):
 a=np.array([v for v in vals if v is not None],dtype=float);assert np.isfinite(a).all()
 return dict(n_available=len(a),mean=float(a.mean())if len(a)else None,median=float(np.median(a))if len(a)else None,min=float(a.min())if len(a)else None,max=float(a.max())if len(a)else None)
