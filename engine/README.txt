Spatial Direct–Bystander Model: simulation engine

This component contains the PhysiCell 1.14.1 source used for the conditional
comparison of two-dimensional, three-dimensional and local-entry scenarios.
The baseline source files and the separate diffusion-thread variant have
individual SHA-256 manifests. The variant changes only main.cpp. Cellular
updates use one thread; the variant allows a separate diffusion thread count.
Configuration values and random seeds are stored beside this component in
../CONFIGURATIONS.json and ../configurations/.

Prepare a separate build directory with Python 3:
  python3 engine/tools/prepare_engine.py --variant baseline --destination /tmp/ikm_baseline_build
Then build in that directory with GNU Make and a GNU C++ compiler supporting
OpenMP:
  make -C /tmp/ikm_baseline_build -j1
The original build uses the flags retained in baseline/Makefile, including
-march=native. Compiler, flags and target CPU can affect numerical replay;
record them when rebuilding. Build parallelism does not set simulation
parallelism.

To reproduce a selected case, change into that case's configuration directory
and run the prepared executable with settings.xml as its argument. Keep the
configuration directory as the working directory so its relative input and
output paths resolve correctly. Run only into a fresh output directory.
The original baseline execution uses OMP_NUM_THREADS=1 and the supplied XML
thread count. The diffusion_parallel variant additionally uses
IKM_DIFFUSION_THREADS=5 for the qualified five-thread diffusion mode while
preserving the one-thread cellular updates. Matching a short replay does not
establish full-trajectory equivalence or numerical convergence.

The source files retain their original citations and licence notices.
VERSION.txt and the exact source manifest identify version 1.14.1; the retained
upstream CITATION.txt also contains illustrative text referring to 1.14.0.
No executable, compiler object, credentials or private correspondence is
included in this component.

These simulations compare specified model representations. They do not
establish a pure dimensionality effect when other conditions differ, and the
local-entry scenario is not a validated model of vasculature, extravasation
or an in vivo tissue. Parameter justification and numerical robustness are
reported separately. The engine alone does not establish those scientific properties.
