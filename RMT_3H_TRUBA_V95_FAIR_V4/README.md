# RMT_3H_TRUBA_V95_FAIR_V4

Corrected Robust Multigrid Technique (RMT) pressure-solver campaign package for the opposed-flow reacting-flow thesis study.

## Scientific comparison rule

This package is generated from the **frozen V95 SG-RBGS campaign source** used by the SG/MG2V/MG2W/MG3V comparisons. The generator replaces only the pressure linear solver and pressure-solver reporting. `tools/audit_fairness.py` requires byte-identical bodies for the momentum predictor, pressure RHS, velocity/face-mass-flux correction, scalar transport, chemistry, thermophysical update, and the complete timed physical timestep loop.

The campaign therefore keeps the same V95 physical/numerical settings as the existing in-house campaigns:

- same 108x36, 216x72, 432x144, 864x288 meshes;
- same 1/2/4/16 thread cases;
- same physical end time 2 s;
- dt0=1e-7 s, dtMax=1e-4 s;
- same CFL/diffusion timestep control;
- pressure relTol=1e-4, absTol=1e-6;
- scalar cycles=4, sweeps/cycle=5;
- constant campaign rho=1, mu=2e-5, cp=1000, Pr=0.7, Sc=1;
- `perfectGas=0`, `variableCp=0`, `sutherland=0`, matching the frozen SG/MG campaign launchers;
- same flow and Arrhenius 96-case matrix.

## RMT implementation

- factor-three multiple coarse grids;
- independent automatic x/y hierarchy depth;
- cell-centred finite-volume boundary coefficients consistent with the V95 pressure grid;
- arithmetic control-volume restriction;
- residue-class/index-mapped coarse-to-fine correction (no interpolation);
- coarsest-to-finest sawtooth ordering, no presmoothing;
- production default **16 smoothing sweeps and 32 deepest-level sweeps**;
- initial Uzawa/RMT correction factor = 1.0;
- correction damping only by residual-based halving when needed;
- **no hidden RBGS fallback**: if a safeguarded RMT correction cannot reduce the V95 pressure residual, the run aborts and is not marked complete;
- failure to reach the frozen V95 pressure tolerance within the maximum cycle count is fatal;
- reusable full-grid workspace; no per-child recursive allocation tree.

### Why 16/32?

The exact first V95 pressure solve was tested with the same M1/LOW campaign state. `post=8/coarse=16` stagnated and was rejected. `post=12/coarse=24` only just met the frozen `1e-4` relative tolerance (`9.49e-5`). `post=16/coarse=32` reached `3.22e-5` in 7 RMT cycles with full `omega=1`, zero backtracking and zero rejection, and subsequently passed the exact V95 reacting-flow smoke to `t=0.02`. Higher 24/48 and 32/64 variants also converged but add more smoothing work per cycle. The selected 16/32 setting is therefore the minimum tested power-of-two setting with a clear tolerance margin, rather than a per-case timing optimum.

`backtrack_halvings` in diagnostics is the total number of omega halvings, not the number of pressure solves that backtracked.

## Build and local test

```bash
./build.sh
./smoke_local.sh
```

Acceptance requires:

```text
FAIRNESS_AUDIT PASS
SELFTEST PASS
RMT V95 FAIR V4 BUILD PASS
RMT_V95_FAIR_SMOKE PASS
```

The smoke test runs the exact V95 campaign physics/numerics to t=0.02 s, beyond the ~0.0077 s location where the legacy RMT coupling collapsed.

## TRUBA

Generate the exact matched 96-case manifest:

```bash
./rebuild_manifest.sh
```

Build on TRUBA:

```bash
sbatch build_rmt.slurm
```

First full single case (FLOW LOW, M1, 16 threads; array index 3):

```bash
TASK_ID=3 sbatch run_one_case.slurm
```

Do not submit the 96-case array until the single-case run is inspected:

```bash
sbatch rmt_96_array.slurm
```

The array is `0-95%3`, matching the previous campaign concurrency convention.

## Provenance

`src/opposedflow_v95_frozen_sg_reference.c` is copied byte-for-byte from the repository's frozen SG-RBGS 96-case reference source. It is never edited. `tools/generate_rmt_v95.py` plus `tools/finalize_generated_v95.py` create `build/opposedflow_rmt3h_v95_fair.c` deterministically.

The legacy failed RMT source is intentionally not used as the production base because it differed from the frozen V95 projection/timestep sequence.


## Campaign post-processing

After or during the 96-case campaign, build one campaign table with:

```bash
python3 collect_summary.py
```

This writes:

```text
RMT_V95_campaign_summary.csv
```

It always follows the 96-row manifest and reports each case as `complete`, `failed`, `partial`, or `missing`. To require all 96 cases to be complete:

```bash
python3 collect_summary.py --require-complete
```

For 1/2/4/16-thread scaling based on the campaign compute timer:

```bash
python3 analyze_scaling.py
```

This writes:

```text
RMT_V95_scaling_summary.csv
```
