Analysis of saved controls and numerical checks

From the package root, run:

python3 analyses/analysis_checks_v1/replay.py --raw-root ../data/raw --inputs-root configurations --output replay_checks

Use a new output directory. The Python package versions are listed in requirements.txt. This component covers 71 technically accepted cases outside the 18 central runs: medium-exchange comparisons, controls, sensitivity and diffusion checks, numerical refinements and completed seed repeats. The nine completed growth and background-loss checks are provided in the kinetic components. No simulations are launched.

Each case is reconstructed from 97 hourly states and separate initialization and terminal states. All consumed files are hash checked. Initialization is compared with the exact input cell table; cell counts, assigned antigen composition and priming states are checked against native validation. All six fields are checked for finite values, recorded bounds and target-local sample counts. IFN-gamma and TNF domain and reservoir inventories are independently reconciled with the recorded mass balance, including external exchange and reservoir loss. Spatial decay is already included in native domain changes and is not counted twice.

Sixty-seven cases have prior all-field references for all 99 states. The other four are earlier matched medium-exchange cases: their existing population and IFN-gamma/TNF domain summaries are used as references, and the remaining field readouts are explicitly marked as additional derived observations. This difference in reference coverage is retained in the verification output.

The component reports 150 scalar measures per case and signed differences from the corresponding central scenario and seed. Controls without effectors use the 2:1 central scenario as their comparison reference; they have no finite target-to-effector ratio. These differences describe model dependence and do not apply a new acceptance threshold. Existing numerical reviews, prespecified reporting tolerances and seed-followup decisions remain separate evidence. Cases are not pooled across different changes.

Field concentrations are relative model units. Target-local values sample the nearest voxel centre to the currently living target population; their time integrals do not follow the dose received by individual cells. Missing samples remain undefined. Four reporting times (24, 48, 72 and 96 hours) and 0–96 hour integrals use the numbered outputs only. The post-loop terminal state is checked separately. Same-seed comparisons do not imply identical stochastic histories. Zero and reversed differences are retained.

Antigen and priming states and death-related laws are assigned model rules, not measured biological mechanism fractions. Cross-scenario differences do not isolate pure dimensionality. The local-entry configuration does not validate vascular transport, extravasation or in vivo tissue behaviour. ECM and Tregs are absent, and CD4/CD8 agents share the specified model laws. Successful replay does not establish numerical convergence, biological calibration or final publication readiness.
