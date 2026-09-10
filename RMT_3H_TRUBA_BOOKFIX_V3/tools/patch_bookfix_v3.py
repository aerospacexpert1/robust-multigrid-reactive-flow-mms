#!/usr/bin/env python3
"""Generate the BOOKFIX_V3 reacting-flow source from the frozen LEGACY_V2 source.

This is intentionally a deterministic source transformation.  The legacy file is retained
unchanged for thesis auditability.  The generated source changes only the RMT pressure
implementation and runtime diagnostics/output behaviour; governing equations, chemistry,
mesh, time stepping and convergence tolerance are not changed here.
"""
from pathlib import Path
import re
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: patch_bookfix_v3.py LEGACY_SOURCE OUTPUT_SOURCE")

src_path = Path(sys.argv[1])
out_path = Path(sys.argv[2])
s = src_path.read_text()

required = [
    "static void solve_pressure_poisson(",
    "static void correct_velocity(",
    "time += dt; double Co = compute_max_co",
    "Pressure RMT calls:",
]
for token in required:
    if token not in s:
        raise SystemExit(f"BOOKFIX patch refused: expected token not found: {token}")

marker = "static void solve_pressure_poisson(const Grid *g, const Phys *ph, Fields *f, int iters, double relaxP) {"
if '#include "rmt_book2d_impl.h"' not in s:
    s = s.replace(marker, '#include "rmt_book2d_impl.h"\n\n' + marker, 1)

new_pressure = r'''static void solve_pressure_poisson(const Grid *g, const Phys *ph, Fields *f, int iters, double relaxP) {
    const double __rmt_pressure_t0 = rmt_wall_now();
    const int nx=g->Nx, ny=g->Ny;
    const double dx2=g->dx*g->dx, dy2=g->dy*g->dy;
    const double ae=1.0/dx2, aw=1.0/dx2, an=1.0/dy2, as=1.0/dy2;
    const double ap=ae+aw+an+as;
    RMTBookWorkspace *ws=rmt_book_workspace(g);
    double *corr=ws->corr;
    double *res=ws->res;
    const double omegaStart=(g_rmtOmega > 0.0) ? g_rmtOmega : 1.0;
    RMTBookConfig cfg;
    cfg.requestedLevels=g_rmtMaxLevels;
    cfg.smoothSweeps=MAX(1,g_rmtPostSmooth);
    cfg.coarseSweeps=MAX(cfg.smoothSweeps,g_rmtCoarseSweeps);
    cfg.nCoarsest=5;
    (void)ph;

    for (int cyc=0; cyc<iters; ++cyc) {
        apply_bc_velocity_pressure(g,ph,f);
        const double oldRes=residual_norm_pressure_variable(g,f,f->rhs);
        if (oldRes < 1.0e-6) break;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            const double Ap=ap*f->p[IDX(i,j,ny)]
                          -ae*(f->p[IDX(i+1,j,ny)] + f->p[IDX(i-1,j,ny)])
                          -an*(f->p[IDX(i,j+1,ny)] + f->p[IDX(i,j-1,ny)]);
            res[IDX(i,j,ny)]=(-f->rhs[IDX(i,j,ny)])-Ap;
        }

        RMTBookDiag d={0,0,0};
        rmt_book_build_correction(g,res,corr,&cfg,&d);

        double omegaTry=omegaStart;
        int accepted=0, acceptedLS=-1;
        for (int ls=0; ls<8; ++ls) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] += omegaTry*corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g,ph,f);
            const double newRes=residual_norm_pressure_variable(g,f,f->rhs);
            if (isfinite(newRes) && newRes < oldRes) {
                accepted=1;
                acceptedLS=ls;
                g_rmtBookLastAcceptedOmega=omegaTry;
                if (ls==0) ++g_rmtBookFullAccepts;
                else g_rmtBookBacktracks += ls;
                break;
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] -= omegaTry*corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g,ph,f);
            omegaTry *= 0.5;
        }

        if (!accepted) {
            ++g_rmtBookRejects;
            g_rmtBookLastAcceptedOmega=0.0;
            for (int sweep=0; sweep<2; ++sweep) {
                for (int color=0; color<2; ++color) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
                    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
                        if (((i+j)&1) != color) continue;
                        const double pGs=(ae*f->p[IDX(i+1,j,ny)] + aw*f->p[IDX(i-1,j,ny)]
                                        +an*f->p[IDX(i,j+1,ny)] + as*f->p[IDX(i,j-1,ny)]
                                        -f->rhs[IDX(i,j,ny)])/MAX(ap,1.0e-20);
                        f->p[IDX(i,j,ny)]=(1.0-relaxP)*f->p[IDX(i,j,ny)] + relaxP*pGs;
                    }
                    apply_bc_velocity_pressure(g,ph,f);
                }
            }
        }
        (void)acceptedLS;
        if (residual_norm_pressure_variable(g,f,f->rhs) < 1.0e-6) break;
    }

    g_pressure_rmt_wall_seconds += rmt_wall_now()-__rmt_pressure_t0;
    g_pressure_rmt_calls += 1;
}
'''

pat = re.compile(r"static void solve_pressure_poisson\(const Grid \*g, const Phys \*ph, Fields \*f, int iters, double relaxP\) \{.*?\n\}\n\nstatic void correct_velocity", re.S)
m = pat.search(s)
if not m:
    raise SystemExit("BOOKFIX patch refused: could not isolate solve_pressure_poisson")
s = s[:m.start()] + new_pressure + "\nstatic void correct_velocity" + s[m.end():]

progress_pat = re.compile(
    r"time \+= dt; double Co = compute_max_co\(&g,&f,dt\); double dtCo=\(Co>1e-12\)\?\(0\.95\*c\.maxCo\*dt/Co\):\(1\.2\*dt\); double dtGrow=1\.2\*dt, dtMax=1e-4; dt=MIN\(dtMax,MIN\(dtCo,dtGrow\)\);\n"
    r"\s*if \(time >= nextWrite - 1e-14\) \{.*?nextWrite \+= c\.writeInterval;\n\s*\}", re.S)
progress_new = r'''time += dt; double Co = compute_max_co(&g,&f,dt); double dtCo=(Co>1e-12)?(0.95*c.maxCo*dt/Co):(1.2*dt); double dtGrow=1.2*dt, dtMax=1e-4; dt=MIN(dtMax,MIN(dtCo,dtGrow));
        if (c.progressEvery > 0 && (step % c.progressEvery) == 0) {
            double Tmax=field_max(&g,f.T), YPmax=field_max(&g,f.YP), YOmax=field_max(&g,f.YO), YFmax=field_max(&g,f.YF);
            printf("step=%6d t=%9.6f dt=%9.3e Co=%8.3e Tmax=%10.3f YPmax=%10.6f YOmax=%10.6f YFmax=%10.6f RMTomega=%9.3e RMTlevels=%d/%d\n",
                   step,time,dt,Co,Tmax,YPmax,YOmax,YFmax,g_rmtBookLastAcceptedOmega,g_rmtBookLastLevelX,g_rmtBookLastLevelY);
            fflush(stdout);
        }
        if (time >= nextWrite - 1e-14) {
            if (g_fieldOutput) {
                write_vti(&g,&ph,&f,time);
                write_centerline_csv(&g,&f,time);
            }
            nextWrite += c.writeInterval;
        }'''
s2, n = progress_pat.subn(lambda _m: progress_new, s, count=1)
if n != 1:
    raise SystemExit(f"BOOKFIX patch refused: progress/output block matches={n}")
s = s2

diag_anchor = '''    printf("Pressure RMT calls: %lld\\n",
           g_pressure_rmt_calls);
'''
diag_repl = diag_anchor + '''    printf("RMT_BOOK_DIAGNOSTICS levelsX=%d levelsY=%d full_accepts=%lld backtracks=%lld rejects=%lld last_omega=%.17g point_updates=%lld\\n",
           g_rmtBookLastLevelX, g_rmtBookLastLevelY,
           g_rmtBookFullAccepts, g_rmtBookBacktracks, g_rmtBookRejects,
           g_rmtBookLastAcceptedOmega, g_rmtBookPointUpdatesTotal);
'''
if diag_anchor not in s:
    raise SystemExit("BOOKFIX patch refused: diagnostics anchor not found")
s = s.replace(diag_anchor, diag_repl, 1)

banner = "// BOOKFIX_V3_GENERATED: Martynenko-style index-mapped pressure RMT; legacy source preserved separately.\n"
s = banner + s
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(s)
print(f"Generated {out_path} ({len(s)} bytes)")
