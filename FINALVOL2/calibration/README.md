# FINALVOL2 smoothing calibration

This directory exists specifically to choose the single uniform RMT
postsmoothing count **before** any new 96-case production campaign.

It does not use BC1/BC2/BC3 or chemistry-stiffness production results to choose
the parameter.  The calibration set is deterministic manufactured pressure
data with the exact frozen-V95 mixed pressure-correction boundary conditions.

Three manufactured spectral contents are used on every campaign mesh:
- LOW: low-frequency mode;
- MIXED: higher-frequency mode;
- COMBINED: deterministic low+high mixture.

Meshes:
`108x36, 216x72, 432x144, 864x288`.

Default candidates:
`8,10,12,14,16`.

The recommendation is based on the minimum total point-update work among
candidates that pass all 12 mode/mesh problems and satisfy the declared cycle
and convergence-factor limits.  Measured wall time is reported but is not used
to select the numerical parameter, avoiding hardware/noise-driven tuning.

Run:

```bash
cd FINALVOL2
./build.sh
python3 calibration/run_smoothing_calibration.py
cat calibration/CALIBRATION_RECOMMENDATION.txt
```

TRUBA batch:

```bash
sbatch calibration/calibrate.slurm
```

Important: the script **never changes** `run_case.sh`, the solver default, the
pressure tolerance, the RMT hierarchy, or any physics.  After reviewing the
calibration result, freeze one smoothing count globally and only then create
the next production campaign revision.
