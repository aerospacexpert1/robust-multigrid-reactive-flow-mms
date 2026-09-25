#!/usr/bin/env python3
"""Build a self-contained RMT fidelity-corrected TRUBA package from the V2 reference.

The reacting-flow physics and benchmark matrix are intentionally unchanged.  This builder
changes only the RMT pressure-correction implementation and benchmark plumbing needed for
traceable/fair reruns.

Reference basis:
- Martynenko (2017), RMT coarse-grid BC extrapolation, Eqs. 2.17--2.20.
- Direction-wise triple coarsening down to few-point grids.
- Full coarse-grid correction as the default (omega=1), with residual line-search only as
  a safeguard/diagnostic.
- Arithmetic/control-volume residual restriction.
- No presmoothing; post-smoothing retained.
"""
from __future__ import annotations

import csv
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "_reference_sources/RMT_3H_TRUBA_96CASE_TRUBA_FIXED_V2/RMT_3H_TRUBA_96CASE_TRUBA_FIXED"
OUTROOT = ROOT / "dist"
OUT = OUTROOT / "RMT_3H_TRUBA_96CASE_TRUBA_FIDELITY_V3"
SRC = OUT / "src/opposedflow_rmt3h.c"


def must_sub(text: str, pattern: str, repl: str, desc: str, flags=re.S) -> str:
    new, n = re.subn(pattern, repl, text, count=1, flags=flags)
    if n != 1:
        raise RuntimeError(f"Patch failed ({desc}): expected 1 replacement, got {n}")
    return new


def patch_source() -> None:
    s = SRC.read_text()

    # Track the physical location of the first unknown on each shifted grid.
    s = s.replace(
        "    double x1, x2;\n} Grid;",
        "    double x1, x2;\n"
        "    /* RMT geometry: distance of local (1,1) cell centre from physical left/bottom. */\n"
        "    double rmtX0, rmtY0;\n"
        "} Grid;",
        1,
    )

    s = s.replace("static double g_rmtOmega = 0.005;", "static double g_rmtOmega = 1.0;", 1)
    s = s.replace(
        "static long long g_pressure_rmt_calls = 0;",
        "static long long g_pressure_rmt_calls = 0;\n"
        "static long long g_rmt_cycles_attempted = 0;\n"
        "static long long g_rmt_full_correction_accepts = 0;\n"
        "static long long g_rmt_line_search_backtracks = 0;\n"
        "static long long g_rmt_rejected_corrections = 0;\n"
        "static double g_rmt_min_accepted_omega = 1.0;",
        1,
    )

    # Auto hierarchy and small true coarsest solve defaults.  run_case.sh also sets these explicitly.
    s = s.replace("    c->rmtMaxLevels = 4;", "    c->rmtMaxLevels = 0; /* 0 = automatic direction-wise coarsening */", 1)
    s = s.replace("    c->rmtCoarseSweeps = 48;", "    c->rmtCoarseSweeps = 16;", 1)

    # Finest cell-centre positions for the RMT geometry; repeat after CLI mesh changes.
    init_anchor = "    g->x1 = 5.0*g->Lx/11.0; g->x2 = 6.0*g->Lx/11.0;"
    s = s.replace(init_anchor, init_anchor + "\n    g->rmtX0 = 0.5*g->dx; g->rmtY0 = 0.5*g->dy;", 1)
    parse_anchor = "    g->dx = g->Lx / g->Nx; g->dy = g->Ly / g->Ny; g->x1 = 5.0*g->Lx/11.0; g->x2 = 6.0*g->Lx/11.0;"
    s = s.replace(parse_anchor, parse_anchor + "\n    g->rmtX0 = 0.5*g->dx; g->rmtY0 = 0.5*g->dy;", 1)

    new_rmt_block = r'''typedef enum {
    RMT_BC_PRESSURE = 0,
    RMT_BC_SCALAR   = 1
} RMTBCType;

typedef struct {
    int maxLevels;      /* 0 => automatic, each direction stops independently */
    int postSmooth;
    int coarseSweeps;
    int coarseRepeats;
    RMTBCType bcType;
    double ax, ay, ap;
    double refDx, refDy;
} RMTSettings;

/*
 * V3 fidelity notes
 * -----------------
 * The V2 recursive implementation discarded the physical offset of a shifted child grid
 * and consequently imposed every child boundary as if it coincided with the physical
 * boundary.  RMT requires offset-dependent coarse-grid correction boundary conditions.
 *
 * We retain a recursive C realization, but make it algebraically mirror the RMT index
 * mapping: the 3x3 shifted children partition the parent unknowns, no interpolation is
 * performed on the way back, and x/y are allowed to stop coarsening independently.
 */

#define RMT_WS_MAX_DEPTH 16

typedef struct {
    double *bc[RMT_WS_MAX_DEPTH];
    double *xc[RMT_WS_MAX_DEPTH];
    size_t cap_bc[RMT_WS_MAX_DEPTH];
    size_t cap_xc[RMT_WS_MAX_DEPTH];
} RMTThreadWorkspace;

static RMTThreadWorkspace *g_rmt_ws = NULL;
static int g_rmt_ws_threads = 0;

static void rmt_workspace_prepare(void) {
    if (g_rmt_ws) return;
    int nt = 1;
#ifdef _OPENMP
    nt = MAX(1, omp_get_max_threads());
#endif
    g_rmt_ws = (RMTThreadWorkspace*)calloc((size_t)nt, sizeof(RMTThreadWorkspace));
    if (!g_rmt_ws) die("alloc RMT thread workspaces");
    g_rmt_ws_threads = nt;
}

static RMTThreadWorkspace *rmt_workspace_current(void) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if (!g_rmt_ws || tid < 0 || tid >= g_rmt_ws_threads) die("RMT workspace/thread mismatch");
    return &g_rmt_ws[tid];
}

static double *rmt_workspace_buffer(double **p, size_t *cap, size_t n) {
    if (*cap < n) {
        double *q = (double*)realloc(*p, n*sizeof(double));
        if (!q) die("realloc RMT workspace");
        *p = q;
        *cap = n;
    }
    memset(*p, 0, n*sizeof(double));
    return *p;
}

static int rmt_can_coarsen_dim(int N) {
    /* N>=9 gives every 3-shift child at least three unknowns.  The resulting terminal
       grids are typically 3--8 points per active direction, matching the 'few points'
       coarsest-grid intent without prematurely stopping the other direction. */
    return N >= 9;
}

static double rmt_xcoord(const Grid *g, int i) { return g->rmtX0 + (double)(i-1)*g->dx; }
static double rmt_ycoord(const Grid *g, int j) { return g->rmtY0 + (double)(j-1)*g->dy; }

static double rmt_cv_width_x(const Grid *g, int i) {
    double x = rmt_xcoord(g, i);
    double a = MAX(0.0, x - 0.5*g->dx);
    double b = MIN(g->Lx, x + 0.5*g->dx);
    return MAX(0.0, b-a);
}

static double rmt_cv_width_y(const Grid *g, int j) {
    double y = rmt_ycoord(g, j);
    double a = MAX(0.0, y - 0.5*g->dy);
    double b = MIN(g->Ly, y + 0.5*g->dy);
    return MAX(0.0, b-a);
}

static void rmt_dirichlet_left(const Grid *g, double *phi, const RMTSettings *st, int j) {
    const int ny = g->Ny;
    if (g->Nx < 2 || g->dx <= 1.0000001*st->refDx) {
        /* Exact finest-grid cell-centred homogeneous Dirichlet correction. */
        phi[IDX(0,j,ny)] = -phi[IDX(1,j,ny)];
        return;
    }
    const double xi = MAX(g->rmtX0/g->dx, 1.0e-12);
    /* Martynenko (2017), Eq. 2.17, c|boundary=0. */
    phi[IDX(0,j,ny)] = 2.0*(xi-1.0)/xi * phi[IDX(1,j,ny)]
                     - (xi-1.0)/(xi+1.0) * phi[IDX(2,j,ny)];
}

static void rmt_dirichlet_right(const Grid *g, double *phi, const RMTSettings *st, int j) {
    const int nx = g->Nx, ny = g->Ny;
    if (nx < 2 || g->dx <= 1.0000001*st->refDx) {
        phi[IDX(nx+1,j,ny)] = -phi[IDX(nx,j,ny)];
        return;
    }
    const double xlast = rmt_xcoord(g, nx);
    const double xi = MAX((g->Lx-xlast)/g->dx, 1.0e-12);
    /* Martynenko (2017), Eq. 2.19, c|boundary=0. */
    phi[IDX(nx+1,j,ny)] = 2.0*(xi-1.0)/xi * phi[IDX(nx,j,ny)]
                        - (xi-1.0)/(xi+1.0) * phi[IDX(nx-1,j,ny)];
}

static void rmt_neumann_bottom(const Grid *g, double *phi, int i) {
    const int ny = g->Ny;
    if (ny < 2) { phi[IDX(i,0,ny)] = phi[IDX(i,1,ny)]; return; }
    const double xi = MAX(g->rmtY0/g->dy, 1.0e-12);
    /* Martynenko (2017), Eq. 2.18, dc/dn=0. */
    phi[IDX(i,0,ny)] = (4.0*xi/(2.0*xi+1.0))*phi[IDX(i,1,ny)]
                     - ((2.0*xi-1.0)/(2.0*xi+1.0))*phi[IDX(i,2,ny)];
}

static void rmt_neumann_top(const Grid *g, double *phi, int i) {
    const int nx = g->Nx, ny = g->Ny;
    (void)nx;
    if (ny < 2) { phi[IDX(i,ny+1,ny)] = phi[IDX(i,ny,ny)]; return; }
    const double ylast = rmt_ycoord(g, ny);
    const double xi = MAX((g->Ly-ylast)/g->dy, 1.0e-12);
    /* Martynenko (2017), Eq. 2.20, dc/dn=0. */
    phi[IDX(i,ny+1,ny)] = (4.0*xi/(2.0*xi+1.0))*phi[IDX(i,ny,ny)]
                        - ((2.0*xi-1.0)/(2.0*xi+1.0))*phi[IDX(i,ny-1,ny)];
}

static void apply_bc_rmt_correction(const Grid *g, double *phi, const RMTSettings *st) {
    const int nx = g->Nx, ny = g->Ny;
    if (st->bcType == RMT_BC_PRESSURE) {
        for (int j=1; j<=ny; ++j) {
            rmt_dirichlet_left(g, phi, st, j);
            rmt_dirichlet_right(g, phi, st, j);
        }
        for (int i=1; i<=nx; ++i) {
            rmt_neumann_bottom(g, phi, i);
            rmt_neumann_top(g, phi, i);
        }
    } else {
        /* Scalar RMT is not active in the production V3 path; retain conservative
           homogeneous extension if this helper is used in a diagnostic. */
        for (int j=1; j<=ny; ++j) {
            phi[IDX(0,j,ny)] = phi[IDX(1,j,ny)];
            phi[IDX(nx+1,j,ny)] = phi[IDX(nx,j,ny)];
        }
        for (int i=1; i<=nx; ++i) {
            phi[IDX(i,0,ny)] = phi[IDX(i,1,ny)];
            phi[IDX(i,ny+1,ny)] = phi[IDX(i,ny,ny)];
        }
    }
    phi[IDX(0,0,ny)] = phi[IDX(1,1,ny)];
    phi[IDX(0,ny+1,ny)] = phi[IDX(1,ny,ny)];
    phi[IDX(nx+1,0,ny)] = phi[IDX(nx,1,ny)];
    phi[IDX(nx+1,ny+1,ny)] = phi[IDX(nx,ny,ny)];
}

static void rmt_smooth_gs(const Grid *g, double *x, const double *b,
                          const RMTSettings *st, int sweeps) {
    const int nx = g->Nx, ny = g->Ny;
    for (int it=0; it<sweeps; ++it) {
        apply_bc_rmt_correction(g, x, st);
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            x[IDX(i,j,ny)] = (b[IDX(i,j,ny)]
                            + st->ax*(x[IDX(i+1,j,ny)] + x[IDX(i-1,j,ny)])
                            + st->ay*(x[IDX(i,j+1,ny)] + x[IDX(i,j-1,ny)])) / MAX(st->ap, 1.0e-20);
        }
    }
    apply_bc_rmt_correction(g, x, st);
}

static int rmt_shift_count(int N, int s) {
    int first = s + 1;
    if (first > N) return 0;
    return 1 + (N-first)/3;
}

static Grid rmt_make_shifted_coarse_grid(const Grid *fine, int sx, int sy,
                                         int coarsenX, int coarsenY) {
    Grid c = *fine;
    if (coarsenX) {
        c.Nx = rmt_shift_count(fine->Nx, sx);
        c.rmtX0 = fine->rmtX0 + (double)sx*fine->dx;
        c.dx = 3.0*fine->dx;
    }
    if (coarsenY) {
        c.Ny = rmt_shift_count(fine->Ny, sy);
        c.rmtY0 = fine->rmtY0 + (double)sy*fine->dy;
        c.dy = 3.0*fine->dy;
    }
    return c;
}

static void rmt_build_level_coeffs(const Grid *g, const RMTSettings *fineSt, RMTSettings *st) {
    *st = *fineSt;
    if (fineSt->bcType == RMT_BC_PRESSURE) {
        st->ax = 1.0/(g->dx*g->dx);
        st->ay = 1.0/(g->dy*g->dy);
        st->ap = 2.0*st->ax + 2.0*st->ay;
    } else {
        const double dtRhoDiffX = fineSt->ax*fineSt->refDx*fineSt->refDx;
        const double dtRhoDiffY = fineSt->ay*fineSt->refDy*fineSt->refDy;
        const double rhoDiag = fineSt->ap - 2.0*fineSt->ax - 2.0*fineSt->ay;
        st->ax = dtRhoDiffX/(g->dx*g->dx);
        st->ay = dtRhoDiffY/(g->dy*g->dy);
        st->ap = rhoDiag + 2.0*st->ax + 2.0*st->ay;
    }
}

static void rmt_restrict_shifted_cv(const Grid *fine, int sx, int sy,
                                    int coarsenX, int coarsenY,
                                    const double *rf, const Grid *coarse, double *rc) {
    const int nyc = coarse->Ny;
    memset(rc, 0, (size_t)(coarse->Nx+2)*(size_t)(coarse->Ny+2)*sizeof(double));
    for (int ic=1; ic<=coarse->Nx; ++ic) {
        const int i0 = coarsenX ? sx+1+3*(ic-1) : ic;
        const int ilo = coarsenX ? MAX(1,i0-1) : i0;
        const int ihi = coarsenX ? MIN(fine->Nx,i0+1) : i0;
        for (int jc=1; jc<=coarse->Ny; ++jc) {
            const int j0 = coarsenY ? sy+1+3*(jc-1) : jc;
            const int jlo = coarsenY ? MAX(1,j0-1) : j0;
            const int jhi = coarsenY ? MIN(fine->Ny,j0+1) : j0;
            double sum = 0.0, vol = 0.0;
            for (int i=ilo; i<=ihi; ++i) {
                const double wx = rmt_cv_width_x(fine, i);
                for (int j=jlo; j<=jhi; ++j) {
                    const double w = wx*rmt_cv_width_y(fine, j);
                    sum += w*rf[IDX(i,j,fine->Ny)];
                    vol += w;
                }
            }
            rc[IDX(ic,jc,nyc)] = (vol > 0.0) ? sum/vol : 0.0;
        }
    }
}

static void rmt_inject_shifted_add(const Grid *fine, int sx, int sy,
                                   int coarsenX, int coarsenY,
                                   double *xf, const Grid *coarse, const double *xc) {
    const int nyf = fine->Ny;
    for (int ic=1; ic<=coarse->Nx; ++ic) {
        const int i0 = coarsenX ? sx+1+3*(ic-1) : ic;
        for (int jc=1; jc<=coarse->Ny; ++jc) {
            const int j0 = coarsenY ? sy+1+3*(jc-1) : jc;
            xf[IDX(i0,j0,nyf)] += xc[IDX(ic,jc,coarse->Ny)];
        }
    }
}

static void rmt_multiple_recursive(const Grid *g, double *x, const double *b,
                                   const RMTSettings *fineSt, int level) {
    if (level >= RMT_WS_MAX_DEPTH-1) die("RMT hierarchy exceeded workspace depth");
    RMTSettings st;
    rmt_build_level_coeffs(g, fineSt, &st);

    int canX = rmt_can_coarsen_dim(g->Nx);
    int canY = rmt_can_coarsen_dim(g->Ny);
    if (fineSt->maxLevels > 0 && level+1 >= fineSt->maxLevels) canX = canY = 0;

    if (!canX && !canY) {
        rmt_smooth_gs(g, x, b, &st, fineSt->coarseSweeps);
        return;
    }

    const int nsx = canX ? 3 : 1;
    const int nsy = canY ? 3 : 1;

    for (int rep=0; rep<fineSt->coarseRepeats; ++rep) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if(level == 0)
#endif
        for (int sx=0; sx<nsx; ++sx) for (int sy=0; sy<nsy; ++sy) {
            Grid gc = rmt_make_shifted_coarse_grid(g, sx, sy, canX, canY);
            if (gc.Nx < 1 || gc.Ny < 1) continue;
            const size_t nc = (size_t)(gc.Nx+2)*(size_t)(gc.Ny+2);
            RMTThreadWorkspace *ws = rmt_workspace_current();
            double *bc = rmt_workspace_buffer(&ws->bc[level], &ws->cap_bc[level], nc);
            double *xc = rmt_workspace_buffer(&ws->xc[level], &ws->cap_xc[level], nc);
            rmt_restrict_shifted_cv(g, sx, sy, canX, canY, b, &gc, bc);
            rmt_multiple_recursive(&gc, xc, bc, fineSt, level+1);
            rmt_inject_shifted_add(g, sx, sy, canX, canY, x, &gc, xc);
        }
    }

    rmt_smooth_gs(g, x, b, &st, fineSt->postSmooth);
}

'''

    s = must_sub(
        s,
        r"typedef enum \{\s*RMT_BC_PRESSURE = 0,.*?\nstatic double rmt_residual_norm_pressure",
        new_rmt_block + "static double rmt_residual_norm_pressure",
        "replace RMT hierarchy/boundary block",
    )

    new_pressure = r'''static void solve_pressure_poisson(const Grid *g, const Phys *ph, Fields *f, int iters, double relaxP) {
    const double __rmt_pressure_t0 = rmt_wall_now();
    const int nx=g->Nx, ny=g->Ny;
    const double dx2 = g->dx*g->dx, dy2 = g->dy*g->dy;
    const double ae = 1.0/dx2, aw = 1.0/dx2, an = 1.0/dy2, as = 1.0/dy2;
    const double ap = ae + aw + an + as;
    const int sweepsPerCycle = MAX(4, g_rmtPostSmooth + g_rmtCoarseRepeats + 1);
    const double omegaRMT = g_rmtOmega;
    const size_t n = (size_t)(nx+2)*(size_t)(ny+2);

    /* Reuse existing full-grid work arrays: no per-time-step pressure heap churn. */
    double *corr = f->work;
    double *res  = f->rhsScalar;
    rmt_workspace_prepare();

    RMTSettings st;
    st.maxLevels = g_rmtMaxLevels;
    st.postSmooth = g_rmtPostSmooth;
    st.coarseSweeps = g_rmtCoarseSweeps;
    st.coarseRepeats = g_rmtCoarseRepeats;
    st.bcType = RMT_BC_PRESSURE;
    st.ax = ae; st.ay = an; st.ap = ap; st.refDx = g->dx; st.refDy = g->dy;

    for (int cyc=0; cyc<iters; ++cyc) {
        apply_bc_velocity_pressure(g, ph, f);
        const double oldRes = residual_norm_pressure_variable(g, f, f->rhs);
        if (oldRes < 1.0e-6) break;

        memset(corr, 0, n*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            const double Ap = ap*f->p[IDX(i,j,ny)]
                            - ae*(f->p[IDX(i+1,j,ny)] + f->p[IDX(i-1,j,ny)])
                            - an*(f->p[IDX(i,j+1,ny)] + f->p[IDX(i,j-1,ny)]);
            res[IDX(i,j,ny)] = (-f->rhs[IDX(i,j,ny)]) - Ap;
        }

        rmt_multiple_recursive(g, corr, res, &st, 0);
        ++g_rmt_cycles_attempted;

        double omegaTry = omegaRMT;
        int accepted = 0;
        int acceptedLs = -1;
        for (int ls=0; ls<8; ++ls) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] += omegaTry*corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g, ph, f);
            const double newRes = residual_norm_pressure_variable(g, f, f->rhs);
            if (isfinite(newRes) && newRes < oldRes) {
                accepted = 1;
                acceptedLs = ls;
                break;
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] -= omegaTry*corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g, ph, f);
            omegaTry *= 0.5;
        }

        if (accepted) {
            if (acceptedLs == 0) ++g_rmt_full_correction_accepts;
            else g_rmt_line_search_backtracks += acceptedLs;
            if (omegaTry < g_rmt_min_accepted_omega) g_rmt_min_accepted_omega = omegaTry;
        } else {
            ++g_rmt_rejected_corrections;
            apply_bc_velocity_pressure(g, ph, f);
        }

        /* Same finest-grid pressure smoother and tolerance as V2; this is intentionally
           not retuned so the comparison isolates the RMT correction implementation. */
        for (int sweep=0; sweep<sweepsPerCycle; ++sweep) {
            for (int color=0; color<2; ++color) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
                for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
                    if (((i+j)&1) != color) continue;
                    const double pGs = (ae*f->p[IDX(i+1,j,ny)] + aw*f->p[IDX(i-1,j,ny)]
                                      + an*f->p[IDX(i,j+1,ny)] + as*f->p[IDX(i,j-1,ny)]
                                      - f->rhs[IDX(i,j,ny)]) / MAX(ap,1.0e-20);
                    f->p[IDX(i,j,ny)] = (1.0-relaxP)*f->p[IDX(i,j,ny)] + relaxP*pGs;
                }
                apply_bc_velocity_pressure(g, ph, f);
            }
        }
        if (residual_norm_pressure_variable(g, f, f->rhs) < 1.0e-6) break;
    }

    g_pressure_rmt_wall_seconds += rmt_wall_now() - __rmt_pressure_t0;
    g_pressure_rmt_calls += 1;
}

'''

    s = must_sub(
        s,
        r"static void solve_pressure_poisson\(.*?\n\}\n\nstatic void correct_velocity",
        new_pressure + "static void correct_velocity",
        "replace pressure solve",
    )

    # Truly suppress heavy field I/O when benchmarking with -fieldOutput 0.
    s = s.replace(
        "            write_vti(&g,&ph,&f,time); write_centerline_csv(&g,&f,time); nextWrite += c.writeInterval;",
        "            if (g_fieldOutput) { write_vti(&g,&ph,&f,time); write_centerline_csv(&g,&f,time); }\n"
        "            nextWrite += c.writeInterval;",
        1,
    )

    # Line-buffer progress so TRUBA diagnostics are visible while a case is running.
    s = s.replace(
        "int main(int argc, char **argv) {\n    struct timespec wall_t0, wall_t1;",
        "int main(int argc, char **argv) {\n    setvbuf(stdout, NULL, _IOLBF, 0);\n    struct timespec wall_t0, wall_t1;",
        1,
    )

    # Explicit diagnostics; keep the original RMT_RUN_SUMMARY key/value contract intact.
    diag_anchor = '    printf("Pressure RMT calls: %lld\\n",\n           g_pressure_rmt_calls);'
    diag_repl = diag_anchor + r'''

    printf("RMT_DIAGNOSTICS cycles=%lld full_accepts=%lld backtracks=%lld rejected=%lld min_accepted_omega=%.17g\n",
           g_rmt_cycles_attempted,
           g_rmt_full_correction_accepts,
           g_rmt_line_search_backtracks,
           g_rmt_rejected_corrections,
           g_rmt_min_accepted_omega);'''
    if diag_anchor not in s:
        raise RuntimeError("Patch failed: diagnostic anchor not found")
    s = s.replace(diag_anchor, diag_repl, 1)

    SRC.write_text(s)


def write_manifest() -> None:
    meshes = [("M1",108,36),("M2",216,72),("M3",432,144),("M4",864,288)]
    threads = [1,2,4,16]
    rows = []
    task = 0

    # Keep the exact V2 ordering so task IDs remain comparable.
    for level, vin in [("LOW",0.05),("NOMINAL",0.10),("HIGH",0.20)]:
        for mesh,nx,ny in meshes:
            for th in threads:
                rows.append([task,f"FLOW_{level}_{mesh}_T{th}","FLOW",level,mesh,nx,ny,th,
                             format(vin,'.17g'),"10000000","5","1000","2"])
                task += 1
    for level,A,beta,Ta in [("S1",1e7,5,1000),("S2",1e8,5,500),("S3",1e9,5,100)]:
        for mesh,nx,ny in meshes:
            for th in threads:
                rows.append([task,f"CHEM_{level}_{mesh}_T{th}","CHEM",level,mesh,nx,ny,th,
                             format(0.10,'.17g'),format(A,'.17g'),format(beta,'.17g'),format(Ta,'.17g'),"2"])
                task += 1
    if task != 96:
        raise RuntimeError(f"manifest size is {task}, expected 96")
    with (OUT/"campaign_manifest.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["task_id","case_id","family","level","mesh","Nx","Ny","threads","vIn","arrA","arrBeta","arrTa","endTime"])
        w.writerows(rows)


def write_run_case() -> None:
    text = r'''#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ID="${1:?usage: ./run_case.sh TASK_ID}"
MANIFEST="${ROOT}/campaign_manifest.csv"
BIN="${ROOT}/build/opposedflow_rmt3h"
[[ -f "${MANIFEST}" ]] || { echo "ERROR: missing ${MANIFEST}"; exit 1; }
[[ -x "${BIN}" ]] || { echo "ERROR: missing binary; run ./build.sh"; exit 1; }
[[ "${TASK_ID}" =~ ^[0-9]+$ ]] || { echo "TASK_ID must be integer"; exit 1; }
(( TASK_ID >= 0 && TASK_ID <= 95 )) || { echo "TASK_ID must be 0..95"; exit 1; }
ROW="$(sed -n "$((TASK_ID+2))p" "${MANIFEST}")"
IFS=',' read -r task_id case_id family level mesh Nx Ny threads vIn arrA arrBeta arrTa endTime <<< "${ROW}"
[[ "${task_id}" == "${TASK_ID}" ]] || { echo "manifest task mismatch"; exit 1; }
OUTDIR="${ROOT}/results/${family}/${level}/${mesh}/T${threads}"
mkdir -p "${OUTDIR}"
rm -f "${OUTDIR}/RUN_COMPLETE" "${OUTDIR}/RUN_FAILED"
cat > "${OUTDIR}/metadata.txt" <<META
solver=RMT_3H_FIDELITY_V3
task_id=${task_id}
case_id=${case_id}
family=${family}
level=${level}
mesh=${mesh}
Nx=${Nx}
Ny=${Ny}
threads=${threads}
vIn=${vIn}
arrA=${arrA}
arrBeta=${arrBeta}
arrTa=${arrTa}
endTime=${endTime}
rmt_levels=auto_directionwise
rmt_post=3
rmt_coarse_sweeps=16
rmt_coarse_repeats=1
rmt_omega_initial=1.0
pressure_tolerance=1e-6
host=$(hostname)
date_start=$(date --iso-8601=seconds)
slurm_job_id=${SLURM_JOB_ID:-none}
slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-none}
slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-none}
META
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    allocated="${SLURM_CPUS_PER_TASK:-0}"
    [[ "${allocated}" =~ ^[0-9]+$ ]] && (( allocated >= threads )) || { echo "ERROR: allocation < benchmark threads"; exit 20; }
fi
export OMP_NUM_THREADS="${threads}"
export OMP_THREAD_LIMIT="${threads}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores
printf '%s\n' "======================================================================" \
"RMT FIDELITY V3 case" \
"Task      : ${task_id}" \
"Case      : ${case_id}" \
"Grid      : ${Nx} x ${Ny}" \
"Threads   : ${threads}" \
"vIn       : ${vIn}" \
"Arr A     : ${arrA}" \
"End time  : ${endTime}" \
"======================================================================"
cd "${OUTDIR}"
set +e
"${BIN}" \
    -Nx "${Nx}" -Ny "${Ny}" \
    -end "${endTime}" -dt 1.0e-7 -maxCo 0.05 \
    -vIn "${vIn}" -arrA "${arrA}" -arrBeta "${arrBeta}" -arrTa "${arrTa}" \
    -pCycles 4 -sCycles 4 \
    -rmtLevels 0 -rmtPost 3 -rmtCoarseSweeps 16 -rmtCoarseRepeats 1 -rmtOmega 1.0 \
    -fieldOutput 0 \
    2>&1 | tee solver.log
RC=${PIPESTATUS[0]}
set -e
printf 'date_end=%s\nexit_code=%s\n' "$(date --iso-8601=seconds)" "${RC}" >> metadata.txt
if (( RC != 0 )); then touch RUN_FAILED; exit "${RC}"; fi
if ! grep -q 'RMT_RUN_SUMMARY' solver.log; then echo "ERROR: missing RMT_RUN_SUMMARY"; touch RUN_FAILED; exit 50; fi
python3 "${ROOT}/tools/extract_one_summary.py" "${OUTDIR}" "${task_id}" "${case_id}" "${family}" "${level}" "${mesh}" "${Nx}" "${Ny}" "${threads}" "${vIn}" "${arrA}" "${arrBeta}" "${arrTa}"
touch RUN_COMPLETE
printf '\nRUN COMPLETE: %s\n' "${case_id}"
'''
    (OUT/"run_case.sh").write_text(text)


def write_tools() -> None:
    td = OUT/"tools"
    td.mkdir(parents=True, exist_ok=True)
    (td/"extract_one_summary.py").write_text(r'''#!/usr/bin/env python3
import csv,re,sys
from pathlib import Path
out=Path(sys.argv[1])
keys=["task_id","case_id","family","level","mesh","Nx","Ny","threads","vIn","arrA","arrBeta","arrTa"]
vals=sys.argv[2:14]
log=(out/"solver.log").read_text(errors="replace")
lines=[x for x in log.splitlines() if x.startswith("RMT_RUN_SUMMARY ")]
if not lines: raise SystemExit("missing RMT_RUN_SUMMARY")
kv=dict(re.findall(r"([A-Za-z_]+)=([^ ]+)",lines[-1]))
diag=[x for x in log.splitlines() if x.startswith("RMT_DIAGNOSTICS ")]
dkv=dict(re.findall(r"([A-Za-z_]+)=([^ ]+)",diag[-1])) if diag else {}
row=dict(zip(keys,vals)); row.update(kv); row.update({"diag_"+k:v for k,v in dkv.items()}); row["solver_version"]="RMT_3H_FIDELITY_V3"
fields=keys+["solver_version","final_time","total_wall_s","pressure_wall_s","pressure_calls","Tmax","YPmax","YOmax","YFmax","diag_cycles","diag_full_accepts","diag_backtracks","diag_rejected","diag_min_accepted_omega"]
with (out/"summary.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerow({k:row.get(k,"") for k in fields})
''')
    (OUT/"collect_summary.py").write_text(r'''#!/usr/bin/env python3
import csv
from pathlib import Path
root=Path(__file__).resolve().parent
files=list((root/"results").glob("**/summary.csv"))
rows=[]; fields=None
for p in files:
    with p.open() as f:
        r=csv.DictReader(f); one=list(r)
        if one:
            rows.extend(one); fields=r.fieldnames
rows.sort(key=lambda x:int(x["task_id"]))
out=root/"campaign_summary.csv"
if fields:
    with out.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
else:
    out.write_text("")
print(f"Wrote {out} with {len(rows)} rows")
print(f"Missing summaries: {96-len(rows)}")
''')
    (OUT/"check_campaign_status.sh").write_text(r'''#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Completed : $(find "${ROOT}/results" -name RUN_COMPLETE 2>/dev/null | wc -l) / 96"
echo "Failed    : $(find "${ROOT}/results" -name RUN_FAILED 2>/dev/null | wc -l)"
echo "Started   : $(find "${ROOT}/results" -name metadata.txt 2>/dev/null | wc -l) / 96"
echo "Summaries : $(find "${ROOT}/results" -name summary.csv 2>/dev/null | wc -l) / 96"
echo; echo "Slurm queue:"; squeue -u "${USER}" || true
''')
    (OUT/"submit_one.sh").write_text(r'''#!/usr/bin/env bash
set -euo pipefail
TASK_ID="${1:?usage: ./submit_one.sh TASK_ID}"
[[ "${TASK_ID}" =~ ^[0-9]+$ ]] && (( TASK_ID>=0 && TASK_ID<=95 )) || { echo "TASK_ID must be 0..95"; exit 2; }
sbatch --job-name="RMT3H_V3_${TASK_ID}" --array="${TASK_ID}-${TASK_ID}%1" rmt_96_array.slurm
''')


def write_slurm() -> None:
    (OUT/"rmt_96_array.slurm").write_text(r'''#!/bin/bash
#SBATCH --job-name=RMT3H_V3_96
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-95%3
#SBATCH --output=logs/slurm-RMT3H_V3-%A_%a.out
#SBATCH --error=logs/slurm-RMT3H_V3-%A_%a.err
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"
mkdir -p logs
./run_case.sh "${SLURM_ARRAY_TASK_ID}"
''')


def write_docs() -> None:
    (OUT/"VERSION.txt").write_text("RMT_3H_TRUBA_96CASE_TRUBA_FIDELITY_V3\n")
    patch_manifest = {
        "base": "RMT_3H_TRUBA_96CASE_TRUBA_FIXED_V2",
        "base_source_blob_sha": "2b7d7cacfc00a030a007c171d4024f2358cb876b",
        "solver_version": "RMT_3H_FIDELITY_V3",
        "physics_changed": False,
        "benchmark_matrix_changed": False,
        "numerical_changes": [
            "Martynenko Eq. 2.17-2.20 offset-aware homogeneous coarse correction BC",
            "direction-wise automatic triple coarsening to few-point grids",
            "control-volume weighted shifted residual restriction",
            "full RMT correction default omega=1.0; residual-decrease line search retained only as safeguard",
            "coarsest GS sweeps 48 -> 16 after hierarchy is allowed to reach true few-point grids",
            "recursive buffers reused per OpenMP thread; per-child calloc/free removed",
            "fieldOutput=0 now suppresses periodic VTK/centerline I/O",
            "RMT line-search/cycle diagnostics added",
        ],
        "unchanged_for_fair_comparison": [
            "Nx/Ny mesh family and 96 task ordering",
            "endTime=2 s",
            "dt0=1e-7 and maxCo=0.05",
            "vIn levels and Arrhenius parameter levels",
            "pressure max cycles=4 and pressure tolerance=1e-6",
            "scalar solver path and scalar cycles=4",
            "thermo, transport, chemistry, momentum and low-Mach projection equations",
            "OpenMP benchmark thread counts 1,2,4,16",
            "gcc -O2 -std=c11 -fopenmp build flags inherited from V2",
        ],
        "validation_gate": "Do not mix V2 and V3 RMT rows. Revalidate RMT MMS and one matched production case before the 96-case thesis campaign.",
    }
    (OUT/"PATCH_MANIFEST.json").write_text(json.dumps(patch_manifest,indent=2)+"\n")
    (OUT/"README_FAIR_COMPARISON.md").write_text(r'''# RMT Fidelity V3 — fair-comparison rules

This package is a **new solver version**, not a retuning of selected benchmark cases.
All 96 cases use one fixed V3 algorithm and the same physical/campaign inputs as V2.

## What is intentionally unchanged
- mesh and task matrix (`M1..M4`, `T1/T2/T4/T16`)
- transient interval `endTime=2 s`
- initial timestep `1e-7` and `maxCo=0.05`
- flow-rate and chemistry-stiffness definitions
- momentum/species/enthalpy/chemistry/thermo paths
- maximum pressure cycles and pressure tolerance
- compiler optimization flags

## What changed (RMT implementation only)
1. Shifted coarse-grid correction boundary conditions use Martynenko's offset parameter `xi` (Eqs. 2.17–2.20).
2. x and y coarsen independently; a short direction no longer truncates the other direction.
3. The hierarchy automatically reaches few-point terminal grids.
4. Coarse residual restriction is control-volume weighted near partial shifted boundaries.
5. The RMT correction starts at full strength (`omega=1`); halving is only a residual-decrease safeguard and is recorded.
6. Per-child recursive heap allocation was replaced with reusable per-thread workspaces.
7. `-fieldOutput 0` now actually disables periodic VTK/centerline writes.

## Thesis provenance
- Never append V3 RMT timings to a file labelled as V2.
- SG/MG2V/MG2W/MG3V/OpenFOAM baselines are not rerun merely because RMT changed.
- RMT MMS results that used the V2 correction path must be regenerated before final thesis claims.
- Use `PATCH_MANIFEST.json` in the thesis/repository record to state exactly what changed.

## Required gate before 96 cases
1. `./build.sh`
2. Run one matched case with `./submit_one.sh TASK_ID`.
3. Check `RMT_DIAGNOSTICS`: ideally full correction acceptance should dominate; repeated deep backtracking is a red flag.
4. Compare final fields and wall time against the exact same SG/MG2V/MG2W/MG3V case.
5. Re-run the RMT branch of MMS A–D with this solver logic.
6. Only then submit the complete 96-case array.
''')


def chmod_scripts() -> None:
    for p in [OUT/"build.sh",OUT/"run_case.sh",OUT/"rmt_96_array.slurm",OUT/"check_campaign_status.sh",OUT/"submit_one.sh"]:
        p.chmod(0o755)
    for p in (OUT/"tools").glob("*.py"):
        p.chmod(0o755)
    (OUT/"collect_summary.py").chmod(0o755)


def validate_manifest() -> None:
    rows=list(csv.DictReader((OUT/"campaign_manifest.csv").open()))
    assert len(rows)==96
    checks={
        28:("FLOW_NOMINAL_M4_T1","864","288","1"),
        40:("FLOW_HIGH_M3_T1","432","144","1"),
        44:("FLOW_HIGH_M4_T1","864","288","1"),
        9:("FLOW_LOW_M3_T2","432","144","2"),
        13:("FLOW_LOW_M4_T2","864","288","2"),
        25:("FLOW_NOMINAL_M3_T2","432","144","2"),
    }
    for idx, expected in checks.items():
        r=rows[idx]
        got=(r["case_id"],r["Nx"],r["Ny"],r["threads"])
        if got!=expected: raise RuntimeError(f"manifest task {idx}: {got} != {expected}")


def main() -> None:
    if not BASE.exists(): raise SystemExit(f"missing base package {BASE}")
    OUTROOT.mkdir(exist_ok=True)
    if OUT.exists(): shutil.rmtree(OUT)
    shutil.copytree(BASE,OUT)
    (OUT/"src/opposedflow_rmt3h_ORIGINAL_V2.c").write_text(SRC.read_text())
    patch_source()
    write_manifest()
    write_run_case()
    write_tools()
    write_slurm()
    write_docs()
    chmod_scripts()
    validate_manifest()
    print(OUT)

if __name__ == "__main__": main()
