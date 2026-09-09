# RMT_3H_TRUBA_BOOKFIX_V3

Corrected candidate implementation for the production opposed-flow reacting-flow RMT campaign.

## Why this version exists

LEGACY_V2 was stable but computationally infeasible on the large transient cases.  Audit against
Martynenko's RMT description and Sergei Martynenko's public `Robust_Multigrid_Technique_2020`
examples identified implementation differences in the pressure RMT correction:

1. shifted coarse-grid endpoints were treated as physical boundaries even when the mapped point
   was offset from the real boundary;
2. x and y coarsening depths were forced to stop together and the hierarchy was capped early;
3. the entire RMT correction was globally damped with `omega=0.005`;
4. a still-large terminal grid was treated as coarsest and given 48 GS sweeps;
5. recursive child arrays were allocated/freed repeatedly;
6. the recursion was closer to a conventional residual-correction tree than Sergei's finest-grid
   index-mapped coarsest-to-finest sawtooth implementation.

BOOKFIX_V3 is a correction candidate addressing those items before the production campaign is
rerun. It must be validated by the included self-test and then by a single reacting-flow case
before the 96-case campaign is submitted.

## What is deliberately NOT changed (fair-comparison contract)

The following remain the same as LEGACY_V2 campaign cases:

- structured mesh family: M1 108x36, M2 216x72, M3 432x144, M4 864x288;
- FLOW and CHEM parameter sweeps and their task ordering;
- physical end time = 2 s;
- initial fluid timestep = 1e-7 s;
- `maxCo = 0.05` and the same adaptive-timestep law in the production source;
- momentum, species, chemistry, enthalpy and thermodynamic closures;
- maximum pressure cycles = 4;
- scalar iteration control = 4;
- pressure stopping criterion already present in the production source (`1e-6`);
- compiler optimization `-O2` and OpenMP benchmark thread mapping.

Therefore any timing difference between LEGACY_V2 and BOOKFIX_V3 is attributed to the corrected
RMT pressure implementation, not to a shorter simulation, easier mesh, looser CFL, altered
chemistry, or different physical case.

## RMT changes

- factor-three coarsening is retained;
- independent automatic x/y maximum levels are used, targeting <=5 points per 1-D deepest
  shifted grid (the same coarsest-size philosophy used in Sergei's public examples);
- all shifted families are represented by index mapping on finest-grid storage;
- traversal is coarsest -> finest with no presmoothing;
- non-coarsest point-Seidel sweeps = 8 and deepest-level sweeps = 16 by default, following the
  public Example 12 pattern where `NSIL=8` and the coarsest grid uses `2*NSIL`;
- shifted boundary stencils use the actual distance between the mapped point and the physical
  boundary; x pressure-correction BC is homogeneous Dirichlet and y is homogeneous Neumann,
  matching the production pressure boundary type;
- full correction is attempted first (`omega=1`). If and only if the fine residual would grow,
  a transparent 1/2 backtracking safeguard is used. Full accepts, backtracks and rejects are
  reported in `RMT_BOOK_DIAGNOSTICS`;
- recursive per-child `calloc/free` is replaced by a reusable per-mesh workspace;
- `-fieldOutput 0` now truly prevents VTK/centerline writes while progress remains visible;
- progress is flushed every 100 fluid steps by the campaign script.

## Build

```bash
chmod +x build.sh rebuild_manifest.sh run_case.sh
./rebuild_manifest.sh
./build.sh
```

`build.sh` first compiles and runs `tests/rmt_pressure_selftest.c`. The production binary is not
accepted unless the self-test prints `SELFTEST PASS`.

The frozen legacy source is `src/opposedflow_rmt3h_LEGACY_V2.c`. The corrected production source
is generated deterministically into `build/opposedflow_rmt3h_BOOKFIX_V3.c` by
`tools/patch_bookfix_v3.py` and then compiled as `build/opposedflow_rmt3h_bookfix_v3`.

## First reacting-flow run

Do not submit 96 cases first. Use one matched case:

```bash
mkdir -p logs
sbatch --export=ALL,TASK_ID=<task_id> run_one_case.slurm
```

Watch live progress:

```bash
tail -f logs/slurm-RMTV3-one-<jobid>.out
```

The critical diagnostics are:

- physical `t` reached versus wall time;
- `RMTomega` in progress lines (ideally 1.0 most of the time);
- `RMT_BOOK_DIAGNOSTICS full_accepts/backtracks/rejects`;
- `pressure_wall_s / total_wall_s`;
- final Tmax/YPmax/YOmax/YFmax relative to the matched baseline/reference case.

Only after a matched single case is numerically acceptable and substantially faster should the
full campaign be launched.

## Thesis reporting rule

Do not merge LEGACY_V2 timings and BOOKFIX_V3 timings into one solver series. Treat LEGACY_V2 as
an implementation-debugging version and BOOKFIX_V3 as the corrected RMT candidate. Re-run the
RMT side of MMS A-D after BOOKFIX_V3 is accepted. SG/MG2V/MG2W/MG3V and OpenFOAM reference data
do not need to be rerun unless their own code/case definitions change.
