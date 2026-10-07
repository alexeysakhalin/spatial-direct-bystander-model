from pathlib import Path
import argparse,subprocess,sys,os,hashlib,json
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
assert not a.output.exists();a.output.mkdir(parents=True)
here=Path(__file__).resolve().parent;root=here.parent
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
for ratio in [2,5]:
 subprocess.run([sys.executable,'-B',str(here/'render_main.py'),'--ratio',str(ratio),'--out',str(a.output/f'ratio_{ratio}')],env=env,check=True)
 subprocess.run([sys.executable,'-B',str(here/'render_curves.py'),'--ratio',str(ratio),'--source',str(root/'source_data'/f'counts_{ratio}to1.json'),'--output',str(a.output/f'curves_{ratio}')],env=env,check=True)
print(json.dumps({'ratios':[2,5],'render_completed':True,'output':str(a.output)}))
