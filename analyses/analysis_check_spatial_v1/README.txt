Spatial analysis of controls and numerical checks

From the package root, run:

python3 analyses/analysis_check_spatial_v1/replay.py --raw-root ../data/raw --inputs-root configurations --output replay_check_spatial

Use a new output directory. This component covers 74 technically accepted cases outside the 18 central runs. It reconstructs 97 hourly states per case, giving 59 spatial measures and their signed differences from the central scenario with the same seed. The included cases and exact input hashes are listed in CASES.json. Tests of different changes are kept separate. Controls without effectors have an undefined target-to-effector ratio and use the 2:1 central configuration as their comparison reference. Undefined spatial values remain undefined.

All consumed raw files and configuration files are hash checked. Seventy cases have prior native spatial results; the replay checks the full hourly arrays, initial hull equations and all 49 native summary measures. Four earlier medium-exchange comparisons have native population and balance validation but no prior spatial report. Their spatial results are additional derived observations from the saved states, with this distinction recorded per case. The population counts at every spatial observation are reconciled with native validation.

Access measures include membership in the fixed convex hull of initial target centres, overlap of recorded cell spheres, distance to living target centres, saved target links and first recorded positive activity memory. Fractions use the current living-effector population, including descendants. A first positive hourly record bounds onset within the sampling interval; it does not identify an exact event time. Hull membership, sphere overlap and saved links are distinct observables. They do not measure fractions of biological death mechanisms.

CXCL gradients are reconstructed from saved Cartesian fields using centred interior differences and first-order edge differences. Local values use the containing voxel. These values do not reconstruct historical cached motility vectors or establish a measured sensing threshold. Field concentrations have relative model units.

Geometry, grids, placement and medium conditions differ across the three representations. Cross-scenario comparisons do not isolate pure dimensionality. The local-entry configuration is a conditional local model without validated vascular transport, extravasation, ECM or Tregs; CD4 and CD8 agents share the assigned model laws. Replay verifies the saved-data calculation. Numerical adequacy and biological parameter support are evaluated separately.
