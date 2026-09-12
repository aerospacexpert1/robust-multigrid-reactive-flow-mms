#!/usr/bin/env python3
"""Generate FINALVOL2 from the frozen V95 SG-RBGS production source.

Only the pressure linear solver, its controls, and reporting are replaced.  The
physical timestep loop and all reacting-flow kernels remain byte-identical to
the frozen V95 comparison source.
"""
from pathlib import Path
import re
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: generate_finalvol2.py FROZEN_V95.c OUTPUT.c")
base = Path(sys.argv[1])
out = Path(sys.argv[2])
s = base.read_text()

required = [
    "static PressureStats solve_pressure_poisson(",
    "static void correct_velocity_and_face_flux(",
    "Non-incremental projection: pressure is intentionally excluded here",
    "memcpy(f.rhoOld, f.rho",
    "correct_velocity_and_face_flux(&g,&ph,&f,dt);",
]
for token in required:
    if token not in s:
        raise SystemExit(f"generation refused: baseline invariant missing: {token}")

# Reuse legacy mg* storage only as an ABI-preserving carrier for RMT controls.
old_defaults = """c->mgPreSmooth = 3;
    c->mgPostSmooth = 3;
    c->mgCoarseSweeps = 40;
    c->mgOmega = 1.0;"""
new_defaults = """c->mgPreSmooth = 0;
    c->mgPostSmooth = 16;
    c->mgCoarseSweeps = 0;
    c->mgOmega = 1.0;"""
if old_defaults not in s:
    raise SystemExit("generation refused: baseline control defaults not found")
s = s.replace(old_defaults, new_defaults, 1)

aliases = [
    (
        'else if (!strcmp(argv[i], "-mgPre") && i+1<argc) c->mgPreSmooth = atoi(argv[++i]);',
        'else if ((!strcmp(argv[i], "-mgPre") || !strcmp(argv[i], "-rmtLevels")) && i+1<argc) c->mgPreSmooth = atoi(argv[++i]);',
    ),
    (
        'else if (!strcmp(argv[i], "-mgPost") && i+1<argc) c->mgPostSmooth = atoi(argv[++i]);',
        'else if ((!strcmp(argv[i], "-mgPost") || !strcmp(argv[i], "-rmtPost") || !strcmp(argv[i], "-rmtSweeps")) && i+1<argc) c->mgPostSmooth = atoi(argv[++i]);',
    ),
    (
        'else if (!strcmp(argv[i], "-mgCoarseSweeps") && i+1<argc) c->mgCoarseSweeps = atoi(argv[++i]);',
        'else if ((!strcmp(argv[i], "-mgCoarseSweeps") || !strcmp(argv[i], "-rmtCoarseSweeps")) && i+1<argc) c->mgCoarseSweeps = atoi(argv[++i]);',
    ),
    (
        'else if ((!strcmp(argv[i], "-mgOmega") || !strcmp(argv[i], "-rbgsOmega")) && i+1<argc) c->mgOmega = atof(argv[++i]);',
        'else if ((!strcmp(argv[i], "-mgOmega") || !strcmp(argv[i], "-rbgsOmega") || !strcmp(argv[i], "-rmtOmega")) && i+1<argc) c->mgOmega = atof(argv[++i]);',
    ),
]
for old,new in aliases:
    if old not in s:
        raise SystemExit(f"generation refused: CLI anchor missing: {old}")
    s=s.replace(old,new,1)

marker = "static PressureStats solve_pressure_poisson("
s = s.replace(marker, '#include "rmt_finalvol2_impl.h"\n\n' + marker, 1)

new_pressure = r'''static PressureStats solve_pressure_poisson(const Grid *g, const Phys *ph, Fields *f,
                                             const Controls *c,
                                             CycleHistory *history) {
    PressureStats stats = {0, 0.0, 0.0, 0.0, 0.0, 1.0};
    const int nx=g->Nx, ny=g->Ny;
    const double ax=1.0/(g->dx*g->dx), ay=1.0/(g->dy*g->dy);
    const double ap=2.0*ax+2.0*ay;
    RMTV2Workspace *ws=rmt_v2_workspace(g);
    double *corr=ws->corr, *defect=ws->res;

    cycle_history_prepare(history, c->poissonIters + 1);
    history->count = 0;
    apply_bc_velocity_pressure(g, ph, f);

    stats.rhsNorm = max_abs_interior(g, f->rhs);
    stats.initialResidual = residual_norm_pressure_variable(g, f, f->rhs);
    stats.absResidual = stats.initialResidual;
    stats.relResidual = stats.absResidual / MAX(stats.rhsNorm, 1.0e-30);
    history->absResidual[0] = stats.absResidual;
    history->relResidual[0] = stats.relResidual;
    history->count = 1;

    if (stats.absResidual <= c->pressureAbsTol || stats.relResidual <= c->pressureRelTol)
        return stats;

    if (fabs(c->mgOmega-1.0) > 1.0e-14)
        die("FINALVOL2 requires the full Martynenko correction (rmtOmega=1)");

    RMTV2Config cfg;
    cfg.requestedLevels = MAX(0, c->mgPreSmooth);
    cfg.smoothSweeps = MAX(1, c->mgPostSmooth);
    cfg.nCoarsest = 5;

    for (int cyc=0; cyc<c->poissonIters; ++cyc) {
        apply_bc_velocity_pressure(g, ph, f);
        const double oldRes = residual_norm_pressure_variable(g, f, f->rhs);

        /* A=-L and b=-rhs, hence defect b-Ap = Lp-rhs. */
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            const double Ap = ap*f->p[IDX(i,j,ny)]
                            - ax*(f->p[IDX(i+1,j,ny)] + f->p[IDX(i-1,j,ny)])
                            - ay*(f->p[IDX(i,j+1,ny)] + f->p[IDX(i,j-1,ny)]);
            defect[IDX(i,j,ny)] = -f->rhs[IDX(i,j,ny)] - Ap;
        }

        RMTV2Diag d = {0,0,0,0,0};
        rmt_v2_build_correction(g, defect, corr, &cfg, &d);

        /* Martynenko sawtooth update: transfer is index remapping; add the
         * computed correction in full.  RMT does not require monotonic
         * residual decrease at every multigrid iteration. */
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
            f->p[IDX(i,j,ny)] += corr[IDX(i,j,ny)];

        apply_bc_velocity_pressure(g, ph, f);
        const double newRes = residual_norm_pressure_variable(g, f, f->rhs);
        if (!isfinite(newRes))
            die("FINALVOL2 produced a non-finite pressure residual");

        const double ratio = newRes/MAX(oldRes,1.0e-300);
        g_rmtV2LastCycleRatio = ratio;
        g_rmtV2MaxCycleRatio = MAX(g_rmtV2MaxCycleRatio,ratio);
        if (newRes > oldRes) ++g_rmtV2NonmonotoneCycles;

        stats.cycles = cyc + 1;
        stats.absResidual = newRes;
        stats.relResidual = stats.absResidual / MAX(stats.rhsNorm, 1.0e-30);
        history->absResidual[history->count] = stats.absResidual;
        history->relResidual[history->count] = stats.relResidual;
        ++history->count;

        if (stats.absResidual <= c->pressureAbsTol || stats.relResidual <= c->pressureRelTol)
            break;
    }

    if (stats.absResidual > c->pressureAbsTol && stats.relResidual > c->pressureRelTol) {
        fprintf(stderr,
                "RMT_PRESSURE_NOT_CONVERGED cycles=%d absResidual=%.17g relResidual=%.17g absTol=%.17g relTol=%.17g\n",
                stats.cycles, stats.absResidual, stats.relResidual,
                c->pressureAbsTol, c->pressureRelTol);
        die("FINALVOL2 pressure solve hit max cycles without meeting the frozen V95 tolerance");
    }

    if (stats.cycles > 0 && stats.initialResidual > 0.0 && stats.absResidual > 0.0)
        stats.convergenceFactor =
            pow(stats.absResidual/stats.initialResidual, 1.0/(double)stats.cycles);
    else
        stats.convergenceFactor = 0.0;

    return stats;
}
'''
pat = re.compile(r"static PressureStats solve_pressure_poisson\(.*?\n\}\n\nstatic void correct_velocity_and_face_flux", re.S)
m = pat.search(s)
if not m:
    raise SystemExit("generation refused: could not isolate frozen SG pressure solve")
s = s[:m.start()] + new_pressure + "\nstatic void correct_velocity_and_face_flux" + s[m.end():]

# Reporting-only changes outside the byte-audited physical timestep loop.
s = s.replace("SG_RBGS", "RMT_FINALVOL2")
s = s.replace("SG-RBGS", "FINALVOL2")
s = s.replace("single-grid RBGS pressure solve",
              "Martynenko-style multiple-coarse-grid RMT pressure solve")
s = s.replace("Pressure grid: finest grid only, %dx%d (no coarse grids)",
              "RMT pressure grid: finest %dx%d; automatic independent factor-3 hierarchy")
s = s.replace(
    'printf("FINALVOL2: omega=%g, maxSweeps=%d, threads=%d\\n", c.mgOmega, c.poissonIters, c.threads);',
    'printf("FINALVOL2: smoothSweeps=%d, maxCycles=%d, threads=%d\\n", c.mgPostSmooth, c.poissonIters, c.threads);'
)
s = s.replace("rbgs_omega,rbgs_max_sweeps", "rmt_smoothing_sweeps,rmt_max_cycles")
s = s.replace(
    "total_rbgs_sweeps,avg_rbgs_sweeps_per_pressure_solve,max_rbgs_sweeps_per_pressure_solve",
    "total_rmt_cycles,avg_rmt_cycles_per_pressure_solve,max_rmt_cycles_per_pressure_solve"
)
s = s.replace("final_pressure_sweeps", "final_pressure_cycles")
s = s.replace("Final pressure: sweeps=", "Final pressure: cycles=")
s = s.replace(
    "Pressure totals: solves=%lld sweeps=%lld avgSweeps=%.6f maxSweeps=%d",
    "Pressure totals: solves=%lld cycles=%lld avgCycles=%.6f maxCycles=%d"
)
s = s.replace("c->mgOmega, c->poissonIters, finalTime, steps",
              "(double)c->mgPostSmooth, c->poissonIters, finalTime, steps", 1)

# Preserve the exact frozen progress label inside the timestep loop.
s = s.replace("rhoRMT", "rhoRBGS")

needle = "    const double pressureFraction = computeElapsed > 0.0\n"
if needle not in s:
    raise SystemExit("generation refused: summary pressureFraction anchor missing")
insert = (
    "    const double rmtEquivalentFineSweeps = (g->Nx > 0 && g->Ny > 0)\n"
    "        ? (double)g_rmtV2PointUpdatesTotal/((double)g->Nx*(double)g->Ny) : 0.0;\n"
    "    const double rmtEquivalentFineSweepsPerSolve = pressureSolves > 0\n"
    "        ? rmtEquivalentFineSweeps/(double)pressureSolves : 0.0;\n"
)
s = s.replace(needle, insert + needle, 1)
s = s.replace(
    'fprintf(fp, ",%.17g,%.17g", (double)totalPressureSweeps, avgSweeps);',
    'fprintf(fp, ",%.17g,%.17g", rmtEquivalentFineSweeps, rmtEquivalentFineSweepsPerSolve);',
    1,
)

anchor = '    printf("Final mass residual:'
diag = (
    '    printf("RMT_FINALVOL2_DIAGNOSTICS levelsX=%d levelsY=%d direct_solves=%lld direct_unknowns=%lld '
    'lu_factorizations=%lld nonmonotone_cycles=%lld last_cycle_ratio=%.17g max_cycle_ratio=%.17g point_updates=%lld\\n",\n'
    '           g_rmtV2LastLevelX, g_rmtV2LastLevelY,\n'
    '           g_rmtV2DirectSolvesTotal, g_rmtV2DirectUnknownsTotal, g_rmtV2LUFactorizationsTotal,\n'
    '           g_rmtV2NonmonotoneCycles, g_rmtV2LastCycleRatio,\n'
    '           g_rmtV2MaxCycleRatio, g_rmtV2PointUpdatesTotal);\n'
)
if anchor not in s:
    raise SystemExit("generation refused: final diagnostics anchor missing")
s = s.replace(anchor, diag + anchor, 1)

banner = "// RMT_FINALVOL2_GENERATED: frozen V95 physics/timestep; pressure linear solver only replaced by optimized Martynenko-style RMT.\n"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(banner + s)
print(f"Generated {out} from frozen V95 SG reference ({len(s)} bytes)")
