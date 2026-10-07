Spatial analysis of saved central simulations

From the package root, run:

python3 analyses/analysis_spatial_v1/replay.py --raw-root ../data/raw --inputs-root configurations --output replay_spatial

Use a new output directory. The Python packages and tested versions are listed in requirements.txt. The program checks the hashes of all files it consumes, reads all 97 hourly states of each of the 18 central runs, and recomputes spatial access and CXCL field derivatives. It compares these results with the recorded reference values and produces 59 metrics per case, 354 scenario summaries, 531 paired ratio differences and 177 summaries of those differences.

The initial target-to-effector ratios are 2:1 and 5:1. Fractions use the current living-effector population, including descendants. Minima and maxima across three seeds describe the observed range; they are not confidence intervals. Corresponding seed labels do not imply identical stochastic event histories.

Spatial access includes membership within the convex hull of initial target centres, instantaneous overlap of cell spheres, nearest living-target centre distance, saved links to living or dead targets, and the first recorded positive model activity state. The hull is fixed, uses xy in 2D and xyz in the other scenarios, and excludes target radii. Sphere radii are calculated from recorded cell volumes in all scenarios. Geometric overlap does not establish antigen recognition; saved links are not counts of deaths. The first positive hourly record only bounds onset within the sampling interval.

CXCL gradients are reconstructed from the saved Cartesian fields with centred interior differences and first-order edge differences. Effector-local values use the containing voxel. Numerical direction normalization follows the recorded software convention. These quantities do not recover historical cached motility vectors, estimate migration speed or establish a biological sensing threshold.

The local-entry scenario does not validate vascular transport, extravasation or in vivo tissue behaviour. Geometry, grids, placement and medium conditions vary together across the three representations. These descriptive results do not isolate a pure dimensionality effect. ECM and Tregs are absent, and CD4/CD8 agents share the specified model laws. Numerical convergence, biological parameter support and final claim review remain separate requirements.
