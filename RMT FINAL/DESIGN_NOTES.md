# RMT FINAL design and V4 correction notes

## 1. Comparison target

The pressure correction solved by frozen V95 is the constant-coefficient
cell-centred equation

```text
Laplace(p) = rhs
```

with fixed pressure on the two x patches and zeroGradient pressure on the y
patches.  RMT solves the equivalent defect equation with

```text
A = -Laplace
b = -rhs
defect = b - A p = Laplace(p) - rhs.
```

The fine RMT operator therefore has to reproduce the V95 ghost-cell operator
exactly: diagonal 3/hx^2 at a homogeneous Dirichlet x boundary and 1/hy^2 at a
homogeneous Neumann y boundary.

## 2. Martynenko features checked against the sources

The 2017 monograph describes a sawtooth RMT cycle with no presmoothing.  The
multigrid iteration starts on the coarsest level, the coarsest tiny systems are
solved directly, and transfer to finer levels is an index-mapping operation
without interpolation.  The correction is retained as the starting
approximation on the next finer level and is finally added to the finest-grid
approximation.

The 2006 paper derives the coarse RHS as a control-volume average of the
fine-grid residual and emphasizes accurate coarse-grid boundary treatment.

The public 2021 OpenMP implementation follows the same pattern:
`Mini_Cycle` visits levels coarsest to finest, the correction is stored on the
fine index space, and `Vanka_iteration` uses Gaussian elimination for the small
block systems.

RMT FINAL mirrors those algorithmic points for the present 2-D scalar pressure
equation.

## 3. V4 boundary restriction bug

V4 used, for stride `s`,

```text
[i-(s-1)/2, i+(s-1)/2]
```

clipped to the domain to define the coarse residual average.

That is correct only when a same-residue coarse neighbour exists on both sides
or when the representative happens to be the first shifted point.  It is wrong
for general shifted boundary grids.

Example for `s=3`: the shifted grid beginning at fine cell `i=3` has centres

```text
2.5h, 5.5h, 8.5h, ...
```

The first coarse control volume is bounded by the physical boundary `x=0` and
the midpoint `x=4h`.  It therefore contains fine cells 1..4.  V4 averaged
cells 2..4 instead.

At the same time the V4 coarse operator correctly used the boundary-to-midpoint
width `V=d+H/2=4h`.  The restricted RHS and coefficient matrix therefore
described different control volumes.

RMT FINAL defines the coarse CV bounds by topology:

```text
left  = midpoint if west same-residue neighbour exists, else physical boundary
right = midpoint if east same-residue neighbour exists, else physical boundary
```

and uses exactly those same widths in the coarse finite-volume coefficients.
A conservation test checks every shifted grid for strides 3 and 9.

## 4. Removal of the V4 monotonic line search

V4 required every RMT cycle to satisfy

```text
||r_new||_inf < ||r_old||_inf
```

and repeatedly halved the correction when this was not true.  After twelve
trials it aborted with `RMT_PRESSURE_REJECT`.

That rule is not part of the Martynenko sawtooth iteration and can reject a
legitimate non-monotone multigrid iteration.  It was exactly the mechanism that
killed both M4 sentinels at cycle 14.

RMT FINAL always applies the computed correction in full.  A pressure solve is
successful only when the unchanged V95 pressure tolerance is reached.  A
non-finite residual or failure to meet the tolerance within the unchanged
maximum cycle count is still fatal.  Non-monotone cycles are counted and
reported rather than suppressed.

## 5. Coarsest solve

RMT FINAL automatically coarsens each direction by factors of three until each
shifted coarsest grid contains at most five points per direction.  Every such
tiny 2-D system is assembled from the same shifted finite-volume operator and
solved by partial-pivot Gaussian elimination.

This replaces V4's fixed number of point sweeps on the coarsest level.

## 6. Parallel smoother

Every point belongs to exactly one shifted grid on a given level.  Two points
connected by a coarse-grid stencil differ by one coarse index and therefore
have opposite red/black colour.  RMT FINAL performs RBGS simultaneously across
all shifted grids.  This preserves grid independence and makes the finest level
parallel; V4's finest level had only one grid and therefore executed its
lexicographic smoother serially.

## 7. Validation ladder

The build is intentionally stronger than the old 81x27 self-test.

1. Boundary CV restriction conservation for every shifted offset.
2. Manufactured mixed-BC pressure solve on all campaign meshes:
   108x36, 216x72, 432x144, 864x288.
3. Byte-identical V95 physics/timestep fairness audit.
4. Exact frozen-V95 first-pressure test on M1 LOW.
5. Exact first-pressure test on the two former M4 failures.
6. Reacting-flow smoke through t=0.02 on M1.
7. Fresh-extraction rebuild in CI.

The full t=2 M4 sentinels and the 96-case campaign remain TRUBA acceptance tests;
CI does not pretend to replace those production runs.
