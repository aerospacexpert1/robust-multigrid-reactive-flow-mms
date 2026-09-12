# FINALVOL2

FINALVOL2 is the calibration/optimization revision derived from **RMT FINAL**.

The numerical RMT method is intentionally unchanged.  FINALVOL2 keeps the same
Martynenko-style factor-three hierarchy, coarsest-to-finest sawtooth,
no-presmoothing structure, control-volume residual restriction, direct
coarsest solves, full correction update, no interpolation, and no hidden
fallback.  The frozen V95 reacting-flow physics, timestep loop, pressure
equation, tolerances, meshes, chemistry, scalar iterations and thread counts are
not changed.

## What is optimized

Only implementation overhead that does not change the mathematical RMT cycle:

1. The fine-defect prefix integral is formed **once per RMT cycle**, then reused
   by every coarse level.
2. Shifted-grid finite-volume coefficients are precomputed once for each level.
3. Red/black shifted-grid point maps are precomputed, removing modulo/division
   work from the smoother hot loop.
4. Coarsest shifted-grid matrices are invariant for a fixed mesh/operator, so
   their LU factorizations are cached and reused.  The RHS is still solved
   exactly on every RMT cycle.
5. A regression test compares the correction from FINALVOL2 with the previous
   RMT FINAL implementation on M1 and M4 and fails if they differ beyond roundoff.

These changes reduce cost per RMT cycle; they do not alter the RMT equations or
the convergence target.

## Important: smoothing is NOT retuned yet

Production scripts remain at the already validated robust baseline:

```text
RMT_SMOOTH_SWEEPS=16
```

FINALVOL2 contains a separate calibration package.  The calibration script does
not edit production code and does not use the 96 reacting-flow benchmark
results to choose the parameter.

Run on TRUBA:

```bash
cd FINALVOL2
chmod +x build.sh calibration/*.py calibration/*.slurm tests/*.sh tools/*.py
./build.sh
sbatch calibration/calibrate.slurm
```

or interactively:

```bash
python3 calibration/run_smoothing_calibration.py \
  --threads 16 --candidates 8,10,12,14,16
```

Outputs:

```text
calibration/calibration_results.csv
calibration/calibration_summary.csv
calibration/CALIBRATION_RECOMMENDATION.txt
```

The default calibration uses three deterministic manufactured pressure fields
(LOW, MIXED and COMBINED) on all four campaign mesh sizes.  A candidate must
satisfy the same mixed pressure-correction boundary conditions and pressure
residual target on all 12 calibration problems.  Among robust candidates the
recommendation minimizes total RMT point-update work.  Wall time is recorded
but is not used to choose the numerical parameter.

After the calibration result is reviewed, create a separate frozen production
revision with one globally fixed smoothing count.  Do **not** case-tune by
velocity, chemistry class, mesh or thread count.

## Validation gates

`./build.sh` performs:

- frozen-V95 fairness audit;
- calibration/production separation audit;
- boundary control-volume restriction test;
- 108x36, 216x72, 432x144, 864x288 pressure mesh ladder at the existing
  16-sweep baseline;
- optimization-equivalence test against RMT FINAL on M1 and M4;
- compilation of the independent calibration executable;
- production solver compilation.

Additional tests:

```bash
OMP_NUM_THREADS=4 RMT_SWEEPS=16 ./tests/first_pressure_matrix.sh
OMP_NUM_THREADS=4 RMT_SWEEPS=16 ./smoke_local.sh
```

No new 96-case campaign should be submitted from FINALVOL2 until the calibration
result has been reviewed and the smoothing parameter is frozen in the next
production revision.
