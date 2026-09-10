# RMT_3H_TRUBA_V95_FAIR_V4

Corrected Robust Multigrid Technique (RMT) pressure-solver campaign package for the opposed-flow reacting-flow thesis study.

## Scientific comparison rule

This package is generated from the **frozen V95 SG-RBGS campaign source** used by the SG/MG2V/MG2W/MG3V comparisons.  The generator replaces only the pressure linear solver and reporting labels.  `tools/audit_fairness.py` requires byte-identical bodies for the momentum predictor, pressure RHS, velocity/face-mass-flux correction, scalar transport, chemistry, thermophysical update, and the complete timed physical timestep loop.

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
- default 8 smoothing sweeps and 16 deepest-level sweeps;
- initial Uzawa/RMT correction factor = 1.0;
- correction damping only by residual-based halving when needed;
- **no hidden RBGS fallback**: if a safeguarded RMT correction cannot reduce the V95 pressure residual, the run aborts and is not marked complete;
- reusable full-grid workspace; no per-child recursive allocation tree.

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

`src/opposedflow_v95_frozen_sg_reference.c` is copied byte-for-byte from the repository's frozen SG-RBGS 96-case reference source.  It is never edited.  `tools/generate_rmt_v95.py` creates `build/opposedflow_rmt3h_v95_fair.c` deterministically.

The legacy failed RMT source is intentionally not used as the production base because it differed from the frozen V95 projection/timestep sequence.
