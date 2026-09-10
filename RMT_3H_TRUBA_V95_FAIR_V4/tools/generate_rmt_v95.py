#!/usr/bin/env python3
"""Generate the fair V95 reacting-flow RMT solver from the frozen SG-RBGS baseline.

The frozen V95 source is the thesis campaign reference.  This transformation replaces
ONLY the pressure linear solver and pressure-solver reporting labels.  The physical
model, timestep loop, momentum predictor, scalar/enthalpy transport, chemistry,
thermodynamics, pressure RHS, velocity/face-flux correction and adaptive timestep are
left unchanged.
"""
from pathlib import Path
import re, sys

if len(sys.argv) != 3:
    raise SystemExit("usage: generate_rmt_v95.py BASELINE.c OUTPUT.c")
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

# RMT controls are stored in the legacy mg* slots only to avoid changing the frozen
# V95 timestep/control ABI: mgPreSmooth=requested level cap (0=auto),
# mgPostSmooth=RMT smoothing sweeps, mgCoarseSweeps=coarsest sweeps, mgOmega=Uzawa start.
s = s.replace(
    "c->mgPreSmooth = 3;\n    c->mgPostSmooth = 3;\n    c->mgCoarseSweeps = 40;\n    c->mgOmega = 1.0;",
    "c->mgPreSmooth = 0;\n    c->mgPostSmooth = 8;\n    c->mgCoarseSweeps = 16;\n    c->mgOmega = 1.0;",
    1,
)

# Add transparent RMT CLI aliases while retaining legacy aliases for auditability.
s = s.replace(
    'else if (!strcmp(argv[i], "-mgPre") && i+1<argc) c->mgPreSmooth = atoi(argv[++i]);',
    'else if ((!strcmp(argv[i], "-mgPre") || !strcmp(argv[i], "-rmtLevels")) && i+1<argc) c->mgPreSmooth = atoi(argv[++i]);',
    1,
)
s = s.replace(
    'else if (!strcmp(argv[i], "-mgPost") && i+1<argc) c->mgPostSmooth = atoi(argv[++i]);',
    'else if ((!strcmp(argv[i], "-mgPost") || !strcmp(argv[i], "-rmtPost")) && i+1<argc) c->mgPostSmooth = atoi(argv[++i]);',
    1,
)
s = s.replace(
    'else if (!strcmp(argv[i], "-mgCoarseSweeps") && i+1<argc) c->mgCoarseSweeps = atoi(argv[++i]);',
    'else if ((!strcmp(argv[i], "-mgCoarseSweeps") || !strcmp(argv[i], "-rmtCoarseSweeps")) && i+1<argc) c->mgCoarseSweeps = atoi(argv[++i]);',
    1,
)
s = s.replace(
    'else if ((!strcmp(argv[i], "-mgOmega") || !strcmp(argv[i], "-rbgsOmega")) && i+1<argc) c->mgOmega = atof(argv[++i]);',
    'else if ((!strcmp(argv[i], "-mgOmega") || !strcmp(argv[i], "-rbgsOmega") || !strcmp(argv[i], "-rmtOmega")) && i+1<argc) c->mgOmega = atof(argv[++i]);',
    1,
)

marker = "static PressureStats solve_pressure_poisson("
s = s.replace(marker, '#include "rmt_book2d_impl.h"\n\n' + marker, 1)

new_pressure = r'''static PressureStats solve_pressure_poisson(const Grid *g, const Phys *ph, Fields *f,
                                             const Controls *c,
                                             CycleHistory *history) {
    PressureStats stats = {0, 0.0, 0.0, 0.0, 0.0, 1.0};
    const int nx=g->Nx, ny=g->Ny;
    const double ax=1.0/(g->dx*g->dx), ay=1.0/(g->dy*g->dy);
    const double ap=2.0*ax+2.0*ay;
    RMTBookWorkspace *ws=rmt_book_workspace(g);
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

    RMTBookConfig cfg;
    cfg.requestedLevels = MAX(0, c->mgPreSmooth); /* 0 = automatic independent x/y depth */
    cfg.smoothSweeps = MAX(1, c->mgPostSmooth);
    cfg.coarseSweeps = MAX(cfg.smoothSweeps, c->mgCoarseSweeps);
    cfg.nCoarsest = 5;

    for (int cyc=0; cyc<c->poissonIters; ++cyc) {
        apply_bc_velocity_pressure(g, ph, f);
        const double oldRes = residual_norm_pressure_variable(g, f, f->rhs);

        /* A=-L and b=-rhs, therefore defect b-Ap = Lp-rhs. */
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            const double Ap = ap*f->p[IDX(i,j,ny)]
                            - ax*(f->p[IDX(i+1,j,ny)] + f->p[IDX(i-1,j,ny)])
                            - ay*(f->p[IDX(i,j+1,ny)] + f->p[IDX(i,j-1,ny)]);
            defect[IDX(i,j,ny)] = -f->rhs[IDX(i,j,ny)] - Ap;
        }

        RMTBookDiag d = {0,0,0};
        rmt_book_build_correction(g, defect, corr, &cfg, &d);

        const double omegaStart = (c->mgOmega > 0.0) ? c->mgOmega : 1.0;
        double omegaTry = omegaStart;
        int accepted = 0;
        double newRes = oldRes;
        for (int ls=0; ls<12; ++ls) {
            #pragma omp parallel for collapse(2) schedule(static)
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] += omegaTry*corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g, ph, f);
            newRes = residual_norm_pressure_variable(g, f, f->rhs);
            if (isfinite(newRes) && newRes < oldRes) {
                accepted = 1;
                g_rmtBookLastAcceptedOmega = omegaTry;
                if (ls == 0) ++g_rmtBookFullAccepts;
                else g_rmtBookBacktracks += ls; /* total halvings, not cycle count */
                break;
            }
            #pragma omp parallel for collapse(2) schedule(static)
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] -= omegaTry*corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g, ph, f);
            omegaTry *= 0.5;
        }

        if (!accepted) {
            ++g_rmtBookRejects;
            g_rmtBookLastAcceptedOmega = 0.0;
            fprintf(stderr, "RMT_PRESSURE_REJECT cycle=%d oldResidual=%.17g after 12 safeguarded trials\n", cyc, oldRes);
            die("RMT pressure correction could not reduce the V95 pressure residual; refusing hidden fallback");
        }

        stats.cycles = cyc + 1;
        stats.absResidual = newRes;
        stats.relResidual = stats.absResidual / MAX(stats.rhsNorm, 1.0e-30);
        history->absResidual[history->count] = stats.absResidual;
        history->relResidual[history->count] = stats.relResidual;
        ++history->count;

        if (stats.absResidual <= c->pressureAbsTol || stats.relResidual <= c->pressureRelTol)
            break;
    }

    if (stats.cycles > 0 && stats.initialResidual > 0.0 && stats.absResidual > 0.0)
        stats.convergenceFactor = pow(stats.absResidual/stats.initialResidual, 1.0/(double)stats.cycles);
    else
        stats.convergenceFactor = 0.0;

    return stats;
}
'''
pat = re.compile(r"static PressureStats solve_pressure_poisson\(.*?\n\}\n\nstatic void correct_velocity_and_face_flux", re.S)
m = pat.search(s)
if not m:
    raise SystemExit("generation refused: could not isolate SG pressure solve")
s = s[:m.start()] + new_pressure + "\nstatic void correct_velocity_and_face_flux" + s[m.end():]

# Reporting-only transformations.  Timed physical loop is not edited.
s = s.replace("SG_RBGS", "RMT_3H_BOOK_V95")
s = s.replace("SG-RBGS", "RMT-3H-BOOK-V95")
s = s.replace("single-grid RBGS pressure solve", "Martynenko-style multiple-coarse-grid RMT pressure solve")
s = s.replace("Pressure grid: finest grid only, %dx%d (no coarse grids)",
              "RMT pressure grid: finest %dx%d; automatic independent factor-3 hierarchy")
s = s.replace("rbgs_omega,rbgs_max_sweeps", "rmt_initial_omega,rmt_max_cycles")
s = s.replace("total_rbgs_sweeps,avg_rbgs_sweeps_per_pressure_solve,max_rbgs_sweeps_per_pressure_solve",
              "total_rmt_cycles,avg_rmt_cycles_per_pressure_solve,max_rmt_cycles_per_pressure_solve")
s = s.replace("final_pressure_sweeps", "final_pressure_cycles")
s = s.replace("Final pressure: sweeps=", "Final pressure: cycles=")
s = s.replace("Pressure totals: solves=%lld sweeps=%lld avgSweeps=%.6f maxSweeps=%d", 
              "Pressure totals: solves=%lld cycles=%lld avgCycles=%.6f maxCycles=%d")
s = s.replace("rhoRBGS", "rhoRMT")
s = s.replace("meanRhoRBGS", "meanRhoRMT")

# Make the two existing equivalent-work fields meaningful for RMT by using point updates/N.
needle = "    const double pressureFraction = computeElapsed > 0.0\n"
if needle not in s:
    raise SystemExit("generation refused: summary pressureFraction anchor missing")
insert = ("    const double rmtEquivalentFineSweeps = (g->Nx > 0 && g->Ny > 0)\n"
          "                                         ? (double)g_rmtBookPointUpdatesTotal/((double)g->Nx*(double)g->Ny) : 0.0;\n"
          "    const double rmtEquivalentFineSweepsPerSolve = pressureSolves > 0\n"
          "                                                 ? rmtEquivalentFineSweeps/(double)pressureSolves : 0.0;\n")
s = s.replace(needle, insert + needle, 1)
s = s.replace('fprintf(fp, ",%.17g,%.17g", (double)totalPressureSweeps, avgSweeps);',
              'fprintf(fp, ",%.17g,%.17g", rmtEquivalentFineSweeps, rmtEquivalentFineSweepsPerSolve);', 1)

# Add explicit RMT diagnostics to stdout without changing the timestep loop.
anchor = '    printf("Final mass residual:'
diag = ('    printf("RMT_BOOK_DIAGNOSTICS levelsX=%d levelsY=%d full_accepts=%lld backtrack_halvings=%lld rejects=%lld last_omega=%.17g point_updates=%lld\\n",\n'
        '           g_rmtBookLastLevelX, g_rmtBookLastLevelY, g_rmtBookFullAccepts,\n'
        '           g_rmtBookBacktracks, g_rmtBookRejects, g_rmtBookLastAcceptedOmega,\n'
        '           g_rmtBookPointUpdatesTotal);\n')
if anchor not in s:
    raise SystemExit("generation refused: final diagnostics anchor missing")
s = s.replace(anchor, diag + anchor, 1)

banner = "// RMT_V95_FAIR_V4_GENERATED: frozen V95 physics/timestep, pressure solver only replaced by RMT.\n"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(banner + s)
print(f"Generated {out} from frozen V95 SG reference ({len(s)} bytes)")
