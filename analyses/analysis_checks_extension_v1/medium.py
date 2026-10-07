from pathlib import Path
import csv,math,xml.etree.ElementTree as ET
import numpy as np
from scipy.io import loadmat

def classify(records,expected):
 assert records and set(x['substrate'] for x in records)=={'ifng','tnf'},'species'
 result={}
 for key in ('ifng','tnf'):
  a=[{k:float(v) for k,v in x.items() if k!='substrate'} for x in records if x['substrate']==key]
  assert len(a)==len(expected),('row count',key)
  assert all(np.isfinite(list(x.values())).all() for x in a),'nonfinite ledger'
  t=np.array([x['time_min'] for x in a]);assert (np.diff(t)>0).all(),'unordered/duplicate time'
  assert np.allclose(t,expected,rtol=0,atol=1e-6),'unexpected ledger time'
  result[key]=a
 return result

def audit_rows(rows,volume,initial_domain,initial_reservoir,decay,flow,external,max_concentration=None):
 assert volume>0 and all(math.isfinite(x) and x>=0 for x in [initial_domain,initial_reservoir,decay,flow,external]),'invalid audit parameters'
 initial=initial_domain+volume*initial_reservoir
 previous_loss=0.;previous_export=0.;scales=[];residuals=[]
 for i,r in enumerate(rows):
  assert r['reservoir_concentration']>=-1e-12 and r['domain_amount']>=-1e-8 and r['reservoir_amount']>=-1e-8,'negative inventory'
  if max_concentration is not None:assert r['reservoir_concentration']<=max_concentration+1e-8,'reservoir upper bound'
  assert abs(r['reservoir_amount']-volume*r['reservoir_concentration'])<=1e-8*(1+abs(r['reservoir_amount'])),'reservoir concentration/amount mismatch'
  loss=r['reservoir_decay_loss'];export=r['external_net_export']
  assert loss>=-1e-8 and loss>=previous_loss-1e-8*(1+abs(previous_loss)),'nonmonotone reservoir loss'
  if decay==0:assert abs(loss)<=1e-12,'loss under zero decay'
  if flow==0:assert abs(export)<=1e-12,'export under zero flow'
  if external==0:
   assert export>=-1e-8 and export>=previous_export-1e-8*(1+abs(previous_export)),'nonmonotone export into zero-concentration bath'
  # Positive external concentration allows signed net import; no false positivity requirement.
  direct=r['domain_amount']+r['reservoir_amount']-initial-r['native_net_change']+export+loss
  scale=1+abs(initial)+abs(r['native_net_change'])+abs(export)+abs(loss)
  assert abs(direct)<=1e-8*scale,'independent inventory balance'
  assert abs(r['balance_residual'])<=1e-8*scale,'reported balance residual'
  assert abs(r['balance_residual']-direct)<=1e-8*scale,'reported versus independently calculated balance'
  if i==0:
   assert abs(r['time_min'])<1e-8,'initial time'
   assert abs(r['domain_amount']-initial_domain)<=1e-8*(1+initial_domain),'initial domain inventory'
   assert abs(r['reservoir_concentration']-initial_reservoir)<1e-12,'initial reservoir concentration'
   assert all(abs(r[k])<1e-8 for k in ['native_net_change','external_net_export','reservoir_decay_loss']),'nonzero cumulative initial ledger'
  scales.append(scale);residuals.append(abs(direct)/scale);previous_loss=loss;previous_export=export
 return {'initial_total_amount':initial,'max_independent_relative_residual':max(residuals),'last':rows[-1]}

def read_field(xml):
 t=ET.parse(xml);variables=t.findall('.//microenvironment/domain/variables/variable');indices={x.get('name'):4+int(x.get('ID')) for x in variables}
 arrays=[v for k,v in loadmat(Path(xml).with_name(Path(xml).stem+'_microenvironment0.mat')).items() if not k.startswith('__')];assert len(arrays)==1
 a=arrays[0];assert a.ndim==2 and np.isfinite(a).all() and (a[3]>0).all(),'invalid field matrix'
 fields={key:a[indices[name]] for key,name in [('ifng','IFN-gamma'),('tnf','TNF')]}
 return float(t.findtext('.//current_time')),a[3],fields

def medium_audit(out,rec,records_override=None):
 out=Path(out);settings=ET.parse(Path(rec['folder'])/'settings.xml');root=settings.getroot()
 def val(path):return float(root.findtext(path))
 dt=val('./overall/dt_diffusion');end=val('./overall/max_time');interval=val('./save/full_data/interval')
 assert dt>0 and interval>0 and end>=0
 assert abs(end/dt-round(end/dt))<1e-6 and abs(interval/dt-round(interval/dt))<1e-6,'audit expects integer timestep-aligned saves/end'
 numbered=sorted(out.glob('output[0-9]*.xml'));expected=np.arange(int(math.floor(end/interval+1e-9))+1)*interval
 assert len(numbered)==len(expected),'numbered frame count'
 xmls=numbered+[out/'final.xml'];times=[float(ET.parse(x).findtext('.//current_time')) for x in xmls]
 assert np.allclose(times[:-1],expected,rtol=0,atol=1e-6),'snapshot times'
 assert abs(times[-1]-(end+dt))<1e-6,'post-loop terminal time'
 records=records_override if records_override is not None else list(csv.DictReader((out/'medium_balance.csv').open()))
 by=classify(records,times);up=root.find('user_parameters')
 volume=float(up.findtext('ikm_medium_reservoir_volume_um3'));flow=float(up.findtext('ikm_medium_external_flow_um3_per_min'))
 assert abs(volume-rec['reservoir_volume_um3'])<=1e-10*volume,'record reservoir volume mismatch'
 ti,vi,ci=read_field(out/'initial.xml');assert abs(ti)<1e-8
 domain_volume=rec['domain_volume_um3'];assert abs(float(vi.sum())-domain_volume)<=1e-8*domain_volume
 summary={};limit=rec.get('max_concentration',None)
 for key in ['ifng','tnf']:
  initial=float(np.dot(vi,ci[key]));decay=float(up.findtext('ikm_medium_'+key+'_decay_per_min'));external=float(up.findtext('ikm_medium_'+key+'_external'));rinit=float(up.findtext('ikm_medium_'+key+'_initial'))
  assert (ci[key]>=-1e-12).all()
  summary[key]=audit_rows(by[key],volume,initial,rinit,decay,flow,external,limit)
  summary[key].update(ledger_rows=len(by[key]),numbered_rows=len(numbered),reservoir_decay_per_min=decay,domain_decay_per_min=val('./microenvironment_setup/variable[@name="'+('IFN-gamma' if key=='ifng' else 'TNF')+'"]/physical_parameter_set/decay_rate'))
 field_rows=[]
 for j,xml in enumerate(xmls):
  t,vol,fields=read_field(xml);assert abs(float(vol.sum())-domain_volume)<=1e-8*domain_volume,'domain volume mismatch'
  row={'time_min':t,'snapshot':xml.name,'kind':'numbered' if j<len(numbered) else 'post_loop_terminal'}
  for key,c in fields.items():
   assert (c>=-1e-12).all(),'negative field'
   if limit is not None:assert c.max()<=limit+1e-8,'field upper bound'
   amount=float(np.dot(c,vol));assert abs(amount-by[key][j]['domain_amount'])<=1e-8*(1+abs(amount)),'saved field/ledger mismatch'
   row[key]={'min':float(c.min()),'max':float(c.max()),'volume_weighted_mean':amount/float(vol.sum())}
  field_rows.append(row)
 return {'passed':True,'validator_version':3,'numbered_frames':len(numbered),'end_time_min':float(expected[-1]),'terminal_time_min':times[-1],'species':summary,'field_time_series':field_rows,'interpretation':'Accounting and snapshot consistency only. Native net change includes domain decay; reservoir decay is a separate inventory loss. Not a measurement of protein lifetime or a validation of biological response.'}
