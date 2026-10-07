# Spatial Direct–Bystander Model

Spatial agent-based simulations of contact-dependent and cytokine-mediated
bystander responses in an antigen-heterogeneous target population.

This repository contains source code, configurations, analysis results and
figures for one simulation study. The accompanying raw-data archive contains
18 central simulations.
It includes a 2D monolayer, a 3D spheroid and a separately labelled local-entry
scenario at initial target:CAR-T ratios of 2:1 and 5:1.

## Scientific question

How do the specified spatial arrangement, access and medium conditions change
the response of targets governed by the same assigned contact and soluble-signal
rules? The model combines antigen-gated direct attacks, IFNγ-associated growth
suppression and priming, and a delayed TNF-associated response that can also
affect antigen-negative targets. Growth suppression is distinct from death.

The central 2:1 configurations retained mean target counts of 53.7, 84.3 and
156.3 at 96 h in 2D, 3D and local entry. At 5:1, the 2D and 3D endpoints were
nearly equal (75.3 and 75.0), while their mean target burdens over 96 h differed
(402.9 and 365.7). Endpoint agreement therefore does not imply matching response
histories. Numerical refinements materially changed several predictions,
particularly local entry; no universal ranking or general numerical convergence
is claimed. Full results, signed contrasts and applicability conditions are
retained in the evidence and analysis components.

## Contents

- `engine/`: exact PhysiCell-derived source, diffusion-thread variant, build
  environment and upstream notices.
- `configurations/` and `CONFIGURATIONS.json`: exact portable inputs and seeds
  for 103 cases, including 18 central cases and 85 controls or challenges.
- `analyses/`: nine saved-data replay components, with raw-file manifests,
  reference values and dependency versions.
- `reference_results/`: population, field, balance and spatial results.
- `figures/`: central comparisons for both ratios, kinetic checks, selected
  numerical follow-ups and matched medium-exchange comparisons. Captions state
  which cases and seed coverage each figure represents; the companion figures
  do not depict every numerical challenge.
- `evidence/`: 17 parameter groups, exact values, 35 primary-source records,
  24 comparative claim statements and 57 linked evidence records.
- `documentation/`: methods, results and selected references.
- `data_manifests/`: original-payload verification and the hashes of exported raw
  files. Only explicitly recorded filesystem-folder strings in settings.xml
  sidecars have been made relative; scientific output arrays are unchanged.

## Reproduce the saved-data analyses

Download `Spatial_Direct_Bystander_raw_central_v1.zip` from the accompanying
Zenodo record. Extract it into the parent directory of this repository,
producing `../data/raw/`. Preserve the archive's directory structure. It contains
18 central simulations: three representations, two initial target:CAR-T ratios
and three seeds per combination.

Use the Python versions and dependencies recorded by the corresponding
component; the repository does not install dependencies automatically.
The following commands replay the central population, field and spatial analyses:

```sh
python3 tools/reproduce.py --component analysis_central_v1 --raw-root ../data/raw --output ../replay_central
python3 tools/reproduce.py --component analysis_field_population_v1 --raw-root ../data/raw --output ../replay_fields
python3 tools/reproduce.py --component analysis_spatial_v1 --raw-root ../data/raw --output ../replay_spatial
```

Each component verifies its required raw files and portable inputs, recomputes
its observations and checks recorded references in a fresh output directory.
The scripts do not launch new simulations. The other analysis components
document their case requirements in their own manifests and READMEs.

Central figures can be rebuilt from their included saved figure data:

```sh
python3 figures/central/reproduce/replay.py --output ../central_figures_replay
```

The matched medium-exchange panels and their separate legend can be rebuilt with:

```sh
python3 figures/medium_exchange/render.py --output ../medium_exchange_replay
```

To rebuild the engine, follow `engine/README.txt`. Exact compiler flags and
thread settings matter. A successful saved-data replay does not imply bitwise
identity of a new simulation under another compiler or CPU.

## Interpretation and limitations

The simulations are conditional model comparisons, without quantitative fitting
to the experimental system. Thirteen of the 17 material parameter-choice bases
remain unresolved for unrestricted biological interpretation. Their assumptions
and sensitivity coverage are reported explicitly. The 33 numerical challenge
cases comprise ten three-seed groups and three single-seed screens. All nine
growth/background-loss checks use ratio 2:1 and seed 0; they are assumption
challenges, not independent calibration. Seed ranges are descriptive, not
confidence intervals.

Geometry, access and medium exchange vary jointly across central scenarios;
these comparisons do not isolate dimensionality. Matched exchange controls are
reported separately. Local entry does not model validated vessels,
extravasation, extracellular matrix or regulatory T cells. CD4/CD8 labels share
assigned model laws. Central IFNγ and TNF fields are identical dimensionless
fields with different prescribed response roles, so separate empirical mechanism
fractions cannot be inferred. The model does not resolve molecular necroptosis.
Antigen-negative targets represent tumour cells, not healthy tissue.

`evidence/comparative_claims.json` records each retained, conditional or rejected
claim with checks and scope. `evidence/claim_source_map.json` resolves its E01–E57 evidence
identifiers. `evidence/provenance_identifiers.json` supplies identifiers and
checksums for the underlying provenance records.

## Licence

Project-owned software: MIT (`LICENSE`). Project-owned data, figures and
scientific documentation: CC BY 4.0 (`DATA_LICENSE.txt`). Third-party components
retain their original licences and required attribution. See
`THIRD_PARTY_NOTICES.txt`.
