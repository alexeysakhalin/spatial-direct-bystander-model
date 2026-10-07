Additional numerical comparisons

From the package root, run:

python3 analyses/analysis_checks_extension_v1/replay.py --raw-root ../data/raw --inputs-root configurations --output replay_scalar_extension

Use a new output directory. This component adds 5 completed numerical cases listed in CASES.json to the previously qualified scalar analysis. The calculation uses the same scientific functions and verifies every consumed raw file, exact input and complete native reference. It reconstructs 97 hourly states per case; initialization and the separate post-loop state are also checked in the scalar component. A second code and metadata location reads the same independently copied raw data. This is a component relocation check, not a second transfer of the complete package.

Differences here use the central configuration with the same ratio, representation and seed. For the coupled 2:1 local-entry checks, this differs from the primary numerical comparison against the refined mesh. Both references remain explicit in the numerical review. Scalar fields use relative units; target-local averages describe the current living population, not tracked-cell dose. Spatial overlap, initial-hull membership and saved links are separate observables. Hourly onset is a saved-time interval, not an exact event time. Zero and reversed differences are retained.

Replay does not establish convergence or biological calibration. Geometry, initial placement and medium conditions differ across representations. The local-entry scenario does not validate vascular transport, extravasation or in vivo tissue; ECM and Tregs are absent and CD4/CD8 agents share the assigned model laws.
