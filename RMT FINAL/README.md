# RMT FINAL

Production-candidate Robust Multigrid Technique (RMT) pressure solver for the
frozen V95 opposed-flow reacting-flow benchmark.

This directory replaces the earlier experimental RMT packages.  The reacting
flow solver is generated from the **same frozen V95 SG-RBGS source** used by the
SG/MG2V/MG2W/MG3V comparison campaigns.  The generator replaces only the
pressure linear solver and pressure-solver reporting.  A byte-level fairness
audit checks the momentum predictor, pressure RHS, velocity/face-flux
correction, scalar transport, chemistry, thermophysical update, and the complete
physical timestep loop.

## Why this package exists

The V4 package passed the 108x36 case but failed the first pressure solve on
864x288.  The failure exposed implementation errors rather than a need for
mesh-dependent tuning:

1. V4 restricted the defect with a symmetric clipped box at shifted coarse-grid
   boundaries, while its coarse finite-volume operator used the true
   boundary-to-midpoint control volume.  The RHS and operator therefore
   represented different control volumes on shifted boundary grids.
2. V4 added a residual-monotonic line search and aborted whenever one full RMT
   correction did not reduce the infinity norm.  That is not the Martynenko
   sawtooth iteration; RMT residual histories need not be monotone cycle by
   cycle.
3. V4 only iterated on the coarsest shifted grids.  RMT FINAL solves each tiny
   coarsest shifted-grid SLAE directly.
4. V4 finest-level smoothing was effectively serial.  RMT FINAL uses
   red-black Gauss-Seidel over all independent shifted grids, including the
   finest grid.

See `DESIGN_NOTES.md` for the derivation and source comparison.

## RMT algorithm used here

- factor-three multiple coarse grids;
- independent automatic x/y hierarchy depth;
- single finest-grid correction storage;
- no presmoothing;
- coarsest-to-finest sawtooth schedule;
- control-volume restriction of the exact fine-grid defect;
- boundary control-volume geometry shared by restriction and coarse operator;
- direct Gaussian elimination on every coarsest shifted grid;
- parallel RBGS postsmoothing on finer levels;
- no interpolation during transfer between levels;
- full correction update `p <- p + c`;
- no monotonic line search and no hidden RBGS fallback;
- frozen V95 pressure tolerances: relTol=1e-4, absTol=1e-6.

The only normal RMT tuning control is the number of postsmoothing sweeps.
The campaign default is `RMT_SMOOTH_SWEEPS=16` on every mesh.

## Build and validation

```bash
chmod +x build.sh smoke_local.sh rebuild_manifest.sh run_case.sh \
         build_rmt.slurm run_one_case.slurm rmt_96_array.slurm \
         submit_gated_96.sh tests/first_pressure_matrix.sh tools/*.py

./build.sh
./tests/first_pressure_matrix.sh
./smoke_local.sh
```

A successful build must contain:

```text
FAIRNESS_AUDIT PASS
BOUNDARY_CV_RESTRICTION PASS
SELFTEST PASS
RMT FINAL BUILD PASS
```

The mesh-ladder self-test exercises 108x36, 216x72, 432x144 and 864x288.
The first-pressure matrix runs the exact frozen V95 first physical timestep for:

- M1 / LOW flow;
- M4 / HIGH flow (the old TASK 47 failure);
- M4 / S3 stiff chemistry (the old TASK 95 failure).

## 96-case campaign

Generate the matched manifest:

```bash
./rebuild_manifest.sh
wc -l campaign_manifest.csv
```

Expected: 97 lines (header + 96 cases).

Recommended gated submission:

```bash
./submit_gated_96.sh
```

This submits:

```text
build
  -> TASK 47 (FLOW high, M4, C16)
  -> TASK 95 (S3 stiff, M4, C16)
  -> 96-case array 0-95%3 only if both sentinels exit 0
```

Direct array submission after independent sentinel inspection:

```bash
sbatch rmt_96_array.slurm
```

## Campaign post-processing

```bash
python3 collect_summary.py
python3 collect_summary.py --require-complete
python3 analyze_scaling.py
```

Outputs:

```text
RMT_FINAL_campaign_summary.csv
RMT_FINAL_scaling_summary.csv
```

The campaign summary also parses `RMT_FINAL_DIAGNOSTICS`, including hierarchy
depth, direct coarse solves, non-monotone RMT cycles, maximum cycle residual
ratio and point-update count.

## Scientific references used for implementation cross-check

- S. I. Martynenko, *The Robust Multigrid Technique: For Black-Box Software*,
  De Gruyter, 2017, especially Sections 2.3.4 and 2.7.
- S. I. Martynenko, "Robust Multigrid Technique for Black Box Software",
  CMAM 6(4), 2006.
- S. I. Martynenko et al., public OpenMP implementation:
  `simartynenko/Robust_Multigrid_Technique_2021_OpenMP`, especially
  `Multigrid_Iterations_M/Multigrid_Iterations_M.f90` and
  `Multigrid_Iterations_M/Multigrid_Iterations_2.f90`.

## Provenance rule

`src/opposedflow_v95_frozen_sg_reference.c` is copied byte-for-byte from the
frozen comparison source and is never edited.  `tools/generate_rmt_final.py`
creates the production source deterministically.  If the fairness audit fails,
the build fails.
