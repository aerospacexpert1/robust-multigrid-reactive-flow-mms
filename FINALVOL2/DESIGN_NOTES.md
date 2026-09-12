# FINALVOL2 design constraints

## Scientific fairness

FINALVOL2 is not allowed to obtain speed-up by changing:
- physical equations or coupling;
- V95 momentum/scalar/chemistry/thermo kernels;
- physical timestep/CFL/diffusion controls;
- pressure equation or pressure tolerances;
- scalar cycle/sweep counts;
- mesh sizes;
- chemistry parameters;
- thread counts used by the comparison campaign.

The byte-level fairness audit retains the same frozen V95 function and complete
physical-timestep-loop hashes used by RMT FINAL.

## Martynenko constraints retained

The implementation remains a multiple-coarse-grid RMT sawtooth:
- factor-three independent shifted grids;
- independent x/y level depth;
- no presmoothing;
- coarsest to finest ordering;
- coarse RHS from fine-defect control-volume averaging;
- correction retained on the fine index space;
- no interpolation between levels;
- direct solution of the tiny coarsest systems;
- one global postsmoothing iteration count for noncoarsest levels;
- full correction application.

FINALVOL2 does not introduce case-specific damping, residual-monotonic
line-search, hidden RBGS fallback, mesh-specific smoothing counts or
physics-specific branches.

## Optimization equivalence

### Prefix reuse
Every level in one RMT cycle restricts the same exact fine-grid defect.
Therefore the two-dimensional prefix integral is invariant within that cycle.
Building it once and sampling different control volumes is algebraically
identical to rebuilding the same prefix before each level.

### Coefficient precomputation
For the frozen constant-coefficient pressure operator, grid geometry and
boundary types are fixed during the run.  Shifted-grid coefficients depend only
on mesh, level, shift and boundary geometry, so computing them once is exactly
equivalent to recomputing them at every point update.

### Red/black maps
Level membership and red/black colour depend only on integer index mapping.
Precomputing those indices removes integer modulo/division operations without
changing update ordering or stencil values.

### Cached coarsest LU
The coarsest matrices do not change between pressure solves.  Reusing an LU
factorization changes only how the same linear system is evaluated:
`A c = r` remains identical and each cycle still performs a direct coarse
solve with its current RHS.

## Calibration firewall

The calibration suite is deliberately outside the production 96-case campaign.
It uses manufactured pressure fields satisfying the same correction BCs and the
four mesh sizes.  It does not inspect the timing/ranking of SG, MG2V, MG2W,
MG3V or OpenFOAM and it does not branch on production case identifiers.

The calibration script is advisory.  It writes a recommendation file but is
forbidden by audit from modifying `run_case.sh` or the production generator.
