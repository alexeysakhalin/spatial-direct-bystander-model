Central saved-data analysis

From the package root, run:

python3 analyses/analysis_central_v1/replay.py --raw-root ../data/raw --inputs-root configurations --output replay_central

Use a new output directory for each invocation. Python 3 and its standard library are sufficient. The program reads the 18 central cases, verifies the recorded raw-file hashes and portable inputs, and reconstructs all 97 saved hourly population states per case. It computes normalized target and effector trajectories, assigned antigen composition, target integrals, six descriptive scenario/ratio summaries and 18 corresponding-seed contrast trajectories. IFNG/TNF fields are reconciled with finite-medium inventory ledgers at every saved state, including the terminal step.

The two ratios are initial target-to-effector ratios of 2:1 and 5:1. The local-entry case is a specified spatial scenario. AUC uses the trapezoidal rule on hourly samples. Dividing the count AUC by 96 h gives the time-average number of living targets; dividing by the initial target count gives normalized AUC in hours. Minima and maxima across three seeds are observed ranges, not confidence intervals. Positive, zero and negative contrasts are retained.

This component reproduces the central saved-data analysis. Spatial-access analysis, numerical refinement and model-adequacy reviews are separate components. Reproducing a result does not establish numerical convergence, a pure dimensionality effect or experimental validation.
