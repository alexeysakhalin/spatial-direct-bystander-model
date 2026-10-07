"""Verify and replay saved simulation outputs; does not run the engine."""
from pathlib import Path
import argparse, os, subprocess, sys
p=argparse.ArgumentParser()
p.add_argument('--raw-root',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--component')
a=p.parse_args(); root=Path(__file__).resolve().parents[1]
parts=sorted(x for x in (root/'analyses').iterdir() if (x/'replay.py').is_file())
if a.component:
 parts=[x for x in parts if x.name==a.component]
 if not parts:p.error('Unknown component')
if a.output.exists():p.error('Use a fresh output directory')
a.output.mkdir(parents=True)
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
for part in parts:
 subprocess.run([sys.executable,'-B',str(part/'replay.py'),'--raw-root',str(a.raw_root.resolve()),'--inputs-root',str(root/'configurations'),'--output',str(a.output.resolve()/part.name)],check=True,env=env)
