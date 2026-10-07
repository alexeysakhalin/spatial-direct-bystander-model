"""Prepare and verify an isolated engine source directory without running it."""
from pathlib import Path
import argparse,hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--variant',choices=['baseline','diffusion_parallel'],required=True)
 parser.add_argument('--destination',type=Path,required=True)
 args=parser.parse_args();dest=args.destination.resolve()
 if dest.exists():raise FileExistsError('Destination already exists; use a fresh directory.')
 manifest=json.loads((ROOT/'metadata/baseline_sha256.json').read_text())
 files={str(p.relative_to(ROOT/'baseline'))for p in(ROOT/'baseline').rglob('*')if p.is_file()}
 if files!=set(manifest):raise ValueError('Baseline file set differs from the recorded manifest.')
 for name,expected in manifest.items():
  path=ROOT/'baseline'/name
  if path.is_symlink()or sha(path)!=expected:raise ValueError('Source checksum mismatch: '+name)
 if args.variant=='diffusion_parallel':
  manifest=json.loads((ROOT/'metadata/diffusion_parallel_sha256.json').read_text())
  if sha(ROOT/'variants/diffusion_parallel/main.cpp')!=manifest['main.cpp']:raise ValueError('Variant checksum mismatch.')
 shutil.copytree(ROOT/'baseline',dest)
 if args.variant=='diffusion_parallel':shutil.copyfile(ROOT/'variants/diffusion_parallel/main.cpp',dest/'main.cpp')
 for name,expected in manifest.items():
  if sha(dest/name)!=expected:raise ValueError('Prepared source checksum mismatch: '+name)
 print(json.dumps({'variant':args.variant,'source_files_verified':len(manifest),'destination':str(dest),'compiled':False,'simulation_run':False}))
if __name__=='__main__':main()
