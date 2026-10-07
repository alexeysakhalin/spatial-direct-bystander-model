Population and field analysis of central simulations

From the package root, run:

python3 analyses/analysis_field_population_v1/replay.py --raw-root ../data/raw --inputs-root configurations --output replay_field_population

Use a new output directory. Tested Python package versions are listed in requirements.txt. The analysis verifies all consumed input hashes and recomputes the hourly populations and field statistics for 18 central cases. Each case contains 97 numbered hourly states from 0 to 96 hours, plus initialization and a terminal state saved one diffusion step after 96 hours. Terminal states are checked separately and are excluded from the 0–96 hour integrals. The saved-data readers and scalar definitions retain their qualified implementations.

The output contains 150 scalar measures per case: 30 population measures and 120 field measures. These cover six living populations, four recorded cytokine fields, four reporting times (24, 48, 72 and 96 hours), population time averages and field integrals. Domain averages are weighted by voxel volume. Target-local values sample the nearest voxel centre to each currently living target. Their quartiles describe this sampled target population. Integrals of the hourly population averages do not track the cumulative dose experienced by individual cells. Missing live-target samples remain undefined. Oxygen and debris are checked in every field state but are not included in the 150 measures.

The target-to-effector ratios are 2:1 and 5:1. Three-seed means, medians and ranges describe the available realizations, not confidence intervals. Signed differences are reported for each same-seed scenario pair within each ratio, and for 5:1 minus 2:1 within each scenario. A shared seed label does not imply identical stochastic histories. Zero and reversed differences are retained; no statistical significance or universal ordering follows from these summaries.

All four cytokine fields use relative, dimensionless concentrations. IFN-gamma and TNF share identical central field dynamics and have distinct assigned roles in the model; their summaries do not identify independent biological contributions. Antigen and priming categories are assigned model states. Population changes do not quantify measured death mechanisms.

The central configurations differ in geometry, medium conditions and spatial placement. Their comparison does not isolate dimensionality. The local-entry configuration is a conditional spatial scenario, with no validated vascular transport, extravasation or in vivo interpretation. ECM and Tregs are absent, and CD4/CD8 agents share the specified model laws. Numerical robustness and quantitative biological parameter support require separate evidence.
