//   gcc -O2 -std=c11 opposedflowV91_RMT_variableDensityLowMach.c -lm -o opposedflowRMT

//  ./opposedflowRMT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define IDX(i,j,ny) ((i)*((ny)+2) + (j))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

typedef struct {
    int Nx, Ny;
    double Lx, Ly;
    double dx, dy;
    double x1, x2;
} Grid;

typedef struct {
    const char *name;
    double W;
    double Tlow, Thigh, Tcommon;
    double lowCpCoeffs[7];
    double highCpCoeffs[7];
    double hRef;
    double As, Ts;
} JanafSpecies;

typedef struct {
    double rho, mu, cp, Pr, Sc;
    double nu, alpha, D;
    bool variableRhoOn;
    bool sutherlandOn;
    bool hydroPressureInDensityOn;
    double rhoRelax;

    double W_F, W_O, W_P, W_N2;
    double Hf_F, Hf_O, Hf_P, Hf_N2;
    double A, beta, Ta;

    double T_in;
    double vFuel, vOx;
    double YF_fuel, YO_fuel, YP_fuel, YN2_fuel;
    double YF_ox,   YO_ox,   YP_ox,   YN2_ox;

    bool chemistryOn;
    double initialChemicalTimeStep;
    double maxChemicalTimeStep;
    double dTchemMax;
    double reactedFracMax;
    double p0; /* absolute total pressure at outlet */
    double Tref; /* sensible enthalpy reference temperature */
    bool variableCpOn;
    JanafSpecies spF, spO, spP, spN2;
} Phys;

typedef struct {
    double startTime, endTime, dt0, writeInterval, maxCo;
    int nOuterCorrectors, nCorrectors;
    int poissonIters;   /* now interpreted as number of RMT correction cycles */
    int scalarIters;    /* now interpreted as number of RMT correction cycles */
    int rmtMaxLevels;
    int rmtPostSmooth;
    int rmtCoarseSweeps;
    int rmtCoarseRepeats;
    double relaxP, relaxU, relaxSc, relaxH;
    int progressEvery;
} Controls;

typedef struct {
    double *p, *u, *v, *T, *h, *YF, *YO, *YP, *YN2;
    double *rho, *rhoOld, *muMix, *nuMix, *alphaMix, *DMix, *cpMix;
    double *us, *vs, *rhs, *work, *phiOld, *rhsScalar;
    double *ufx, *vfy, *mfx, *mfy; /* collocated face-normal volume and mass fluxes */
} Fields;

static double u_face_x(const Grid *g, const Fields *f, int iface, int j);
static double v_face_y(const Grid *g, const Phys *ph, const Fields *f, int i, int jface);
static double rho_face_x(const Grid *g, const Fields *f, int iface, int j);
static double rho_face_y(const Grid *g, const Fields *f, int i, int jface);
static double face_avg_x_coeff(const Grid *g, const double *a, int i, int j);
static double face_avg_y_coeff(const Grid *g, const double *a, int i, int j);

static int g_rmtMaxLevels = 6;
static int g_rmtPostSmooth = 3;
static int g_rmtCoarseSweeps = 48;
static int g_rmtCoarseRepeats = 2;
static double g_rmtOmega = 0.005;

static double g_pressure_rmt_wall_seconds = 0.0;
static long long g_pressure_rmt_calls = 0;
static int g_fieldOutput = 1;

static double rmt_wall_now(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec;
#endif
}



static void die(const char *msg) { fprintf(stderr, "FATAL: %s\n", msg); exit(1); }
static double clamp(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }
static int ensure_dir(const char *path) { if (mkdir(path, 0777) == 0 || errno == EEXIST) return 0; return -1; }

static const double RU_J_PER_KMOLK = 8314.46261815324;

static void set_janaf_species(JanafSpecies *sp, const char *name, double W, double Tlow, double Thigh, double Tcommon,
                              const double low[7], const double high[7], double As, double Ts) {
    sp->name = name;
    sp->W = W;
    sp->Tlow = Tlow;
    sp->Thigh = Thigh;
    sp->Tcommon = Tcommon;
    for (int k=0; k<7; ++k) {
        sp->lowCpCoeffs[k] = low[k];
        sp->highCpCoeffs[k] = high[k];
    }
    sp->hRef = 0.0;
    sp->As = As;
    sp->Ts = Ts;
}

static inline double janaf_limit_T(const JanafSpecies *sp, double T) {
    return clamp(T, sp->Tlow, sp->Thigh);
}

static inline const double *janaf_coeffs(const JanafSpecies *sp, double T) {
    double Tc = janaf_limit_T(sp, T);
    return (Tc < sp->Tcommon) ? sp->lowCpCoeffs : sp->highCpCoeffs;
}

static double janaf_Cp_mass(const JanafSpecies *sp, double T) {
    T = janaf_limit_T(sp, T);
    const double *a = janaf_coeffs(sp, T);
    return (RU_J_PER_KMOLK / sp->W) * ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
}

static double janaf_Ha_mass(const JanafSpecies *sp, double T) {
    T = janaf_limit_T(sp, T);
    const double *a = janaf_coeffs(sp, T);
    return (RU_J_PER_KMOLK / sp->W) * (((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T + a[5]);
}

static void init_janaf_defaults(Phys *ph) {
    const double AsCO = 1.663e-06, TsCO = 136.0;
    const double AsO2 = 1.6934e-06, TsO2 = 127.0;
    const double AsCO2 = 1.370e-06, TsCO2 = 222.0;
    const double AsN2 = 1.406e-06, TsN2 = 111.0;

    static const double CO_low[7]  = { 3.57953347, -0.00061035368,  1.01681433e-06,  9.07005884e-10, -9.04424499e-13, -14344.0860,   3.50840928 };
    static const double CO_high[7] = { 2.71518561,  0.00206252743, -9.98825771e-07,  2.30053008e-10, -2.03647716e-14, -14151.8724,   7.81868772 };

    static const double O2_low[7]  = { 3.78245636, -0.00299673416,  9.84730201e-06, -9.68129509e-09,  3.24372837e-12,  -1063.94356,  3.65767573 };
    static const double O2_high[7] = { 3.28253784,  0.00148308754, -7.57966669e-07,  2.09470555e-10, -2.16717794e-14,  -1088.45772,  5.45323129 };

    static const double CO2_low[7]  = { 2.35677352,  0.00898459677, -7.12356269e-06,  2.45919022e-09, -1.43699548e-13, -48371.9697,  9.90105222 };
    static const double CO2_high[7] = { 3.85746029,  0.00441437026, -2.21481404e-06,  5.23490188e-10, -4.72084164e-14, -48759.1660,  2.27163806 };

    static const double N2_low[7]  = { 3.29867700,  0.00140824040, -3.96322200e-06,  5.64151500e-09, -2.44485400e-12,  -1020.89990,  3.95037200 };
    static const double N2_high[7] = { 2.92664000,  0.00148797680, -5.68476000e-07,  1.00970380e-10, -6.75335100e-15,   -922.79770,  5.98052800 };

    set_janaf_species(&ph->spF,  "F/CO",  28.01055, 200.0, 3500.0, 1000.0, CO_low,  CO_high,  AsCO,  TsCO);
    set_janaf_species(&ph->spO,  "O/O2",  31.99880, 200.0, 3500.0, 1000.0, O2_low,  O2_high,  AsO2,  TsO2);
    set_janaf_species(&ph->spP,  "P/CO2", 44.00995, 200.0, 3500.0, 1000.0, CO2_low, CO2_high, AsCO2, TsCO2);
    set_janaf_species(&ph->spN2, "N2",    28.01340, 250.0, 5000.0, 1000.0, N2_low,  N2_high,  AsN2,  TsN2);
}

static void update_janaf_references(Phys *ph) {
    ph->spF.hRef  = janaf_Ha_mass(&ph->spF,  ph->Tref);
    ph->spO.hRef  = janaf_Ha_mass(&ph->spO,  ph->Tref);
    ph->spP.hRef  = janaf_Ha_mass(&ph->spP,  ph->Tref);
    ph->spN2.hRef = janaf_Ha_mass(&ph->spN2, ph->Tref);
}

static inline double species_h_sensible_rel(const JanafSpecies *sp, double T) {
    return janaf_Ha_mass(sp, T) - sp->hRef;
}

static inline double mixture_cp(const Phys *ph, double YF, double YO, double YP, double YN2, double T) {
    if (!ph->variableCpOn) return ph->cp;
    return YF*janaf_Cp_mass(&ph->spF, T) + YO*janaf_Cp_mass(&ph->spO, T) + YP*janaf_Cp_mass(&ph->spP, T) + YN2*janaf_Cp_mass(&ph->spN2, T);
}

static inline double mixture_h_sensible_rel(const Phys *ph, double YF, double YO, double YP, double YN2, double T) {
    if (!ph->variableCpOn) return ph->cp*(T - ph->Tref);
    return YF*species_h_sensible_rel(&ph->spF, T) + YO*species_h_sensible_rel(&ph->spO, T) + YP*species_h_sensible_rel(&ph->spP, T) + YN2*species_h_sensible_rel(&ph->spN2, T);
}

static inline double mixture_mw(const Phys *ph, double YF, double YO, double YP, double YN2) {
    double invW = YF/ph->W_F + YO/ph->W_O + YP/ph->W_P + YN2/ph->W_N2;
    return (invW > 1.0e-16) ? 1.0/invW : ph->W_N2;
}

static inline double mixture_R_mass(const Phys *ph, double YF, double YO, double YP, double YN2) {
    return RU_J_PER_KMOLK / mixture_mw(ph, YF, YO, YP, YN2);
}

static inline double sutherland_mu_species(const JanafSpecies *sp, double T) {
    T = MAX(T, 150.0);
    return sp->As*sqrt(T)/(1.0 + sp->Ts/T);
}

static inline double mixture_mu(const Phys *ph, double YF, double YO, double YP, double YN2, double T) {
    if (!ph->sutherlandOn) return ph->mu;
    return YF*sutherland_mu_species(&ph->spF, T)
         + YO*sutherland_mu_species(&ph->spO, T)
         + YP*sutherland_mu_species(&ph->spP, T)
         + YN2*sutherland_mu_species(&ph->spN2, T);
}

static inline double mixture_density(const Phys *ph, double YF, double YO, double YP, double YN2, double T, double pAbs) {
    if (!ph->variableRhoOn) return ph->rho;
    double Rmix = mixture_R_mass(ph, YF, YO, YP, YN2);
    double pUse = MAX(pAbs, 1.0e3);
    return pUse / MAX(Rmix*MAX(T, 200.0), 1.0);
}

static double T_from_h_mix(const Phys *ph, double YF, double YO, double YP, double YN2, double h, double Tguess) {
    if (!ph->variableCpOn) return ph->Tref + h/ph->cp;
    double Tmin = 250.0, Tmax = 5000.0;
    if (Tguess < Tmin || Tguess > Tmax || !isfinite(Tguess)) Tguess = ph->Tref + h / MAX(ph->cp, 1.0);
    double hLo = mixture_h_sensible_rel(ph, YF, YO, YP, YN2, Tmin);
    double hHi = mixture_h_sensible_rel(ph, YF, YO, YP, YN2, Tmax);
    if (h <= hLo) return Tmin;
    if (h >= hHi) return Tmax;
    double T = clamp(Tguess, Tmin, Tmax);
    for (int it=0; it<12; ++it) {
        double hEval = mixture_h_sensible_rel(ph, YF, YO, YP, YN2, T);
        double cpMix = MAX(mixture_cp(ph, YF, YO, YP, YN2, T), 1.0);
        double dT = (hEval - h) / cpMix;
        T -= dT;
        if (T <= Tmin || T >= Tmax || !isfinite(T)) break;
        if (fabs(dT) < 1.0e-10*MAX(1.0, fabs(T))) return T;
    }
    double lo = Tmin, hi = Tmax;
    for (int it=0; it<60; ++it) {
        double mid = 0.5*(lo + hi);
        double hMid = mixture_h_sensible_rel(ph, YF, YO, YP, YN2, mid);
        if (hMid < h) lo = mid; else hi = mid;
    }
    return 0.5*(lo + hi);
}

static inline double h_from_T_mix(const Phys *ph, double YF, double YO, double YP, double YN2, double T) {
    return mixture_h_sensible_rel(ph, YF, YO, YP, YN2, T);
}


static void update_thermo_transport_fields(const Grid *g, const Phys *ph, Fields *f) {
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) {
        double yF = clamp(f->YF[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yO = clamp(f->YO[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yP = clamp(f->YP[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yN2 = MAX(0.0, 1.0 - yF - yO - yP);
        double T = MAX(f->T[IDX(i,j,g->Ny)], 200.0);
        double cpMix = mixture_cp(ph, yF, yO, yP, yN2, T);
        double pUse = ph->hydroPressureInDensityOn ? f->p[IDX(i,j,g->Ny)] : ph->p0;
        double rhoRaw = mixture_density(ph, yF, yO, yP, yN2, T, pUse);
        double rhoOld = f->rho[IDX(i,j,g->Ny)];
        double rhoMix = rhoRaw;
        if (ph->variableRhoOn && isfinite(rhoOld) && rhoOld > 0.0 && ph->rhoRelax < 0.999999) {
            rhoMix = ph->rhoRelax*rhoRaw + (1.0 - ph->rhoRelax)*rhoOld;
        }
        double muMix = mixture_mu(ph, yF, yO, yP, yN2, T);
        f->YN2[IDX(i,j,g->Ny)] = yN2;
        f->cpMix[IDX(i,j,g->Ny)] = cpMix;
        f->rho[IDX(i,j,g->Ny)] = rhoMix;
        f->muMix[IDX(i,j,g->Ny)] = muMix;
        f->nuMix[IDX(i,j,g->Ny)] = muMix / MAX(rhoMix, 1.0e-12);
        f->alphaMix[IDX(i,j,g->Ny)] = muMix / MAX(rhoMix*ph->Pr, 1.0e-12);
        f->DMix[IDX(i,j,g->Ny)] = muMix / MAX(rhoMix*ph->Sc, 1.0e-12);
    }
}

static double interior_average(const Grid *g, const double *a) {
    double s = 0.0; int n = 0;
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) { s += a[IDX(i,j,g->Ny)]; ++n; }
    return (n > 0) ? s/(double)n : 0.0;
}

static inline double reaction_heat_release_per_extent(const Phys *ph, double T) {
    if (ph->variableCpOn) {
        return 2.0*ph->spF.W*janaf_Ha_mass(&ph->spF, T) + 1.0*ph->spO.W*janaf_Ha_mass(&ph->spO, T) - 2.0*ph->spP.W*janaf_Ha_mass(&ph->spP, T);
    }
    return (2.0*ph->W_F*ph->Hf_F + ph->W_O*ph->Hf_O - 2.0*ph->W_P*ph->Hf_P);
}

static void refresh_temperature_from_enthalpy(const Grid *g, const Phys *ph, Fields *f) {
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) {
        double yF = clamp(f->YF[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yO = clamp(f->YO[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yP = clamp(f->YP[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yN2 = MAX(0.0, 1.0 - yF - yO - yP);
        f->YN2[IDX(i,j,g->Ny)] = yN2;
        f->T[IDX(i,j,g->Ny)] = T_from_h_mix(ph, yF, yO, yP, yN2, f->h[IDX(i,j,g->Ny)], f->T[IDX(i,j,g->Ny)]);
    }
}

static void init_case(Grid *g, Phys *ph, Controls *c) {
    g->Nx = 220; g->Ny = 80;
    g->Lx = 0.02; g->Ly = 0.02;
    g->dx = g->Lx / g->Nx; g->dy = g->Ly / g->Ny;
    g->x1 = 5.0*g->Lx/11.0; g->x2 = 6.0*g->Lx/11.0;

    ph->rho = 1.0; ph->mu = 2.0e-5; ph->cp = 1000.0; ph->Pr = 0.7; ph->Sc = 1.0;
    ph->variableRhoOn = true; ph->sutherlandOn = true;
    ph->hydroPressureInDensityOn = false; ph->rhoRelax = 0.35;
    ph->nu = ph->mu/ph->rho; ph->alpha = ph->nu/ph->Pr; ph->D = ph->nu;
    init_janaf_defaults(ph);
    ph->W_F = ph->spF.W; ph->W_O = ph->spO.W; ph->W_P = ph->spP.W; ph->W_N2 = ph->spN2.W;
    ph->Hf_F = 1.0e8 / ph->W_F; ph->Hf_O = 0.0; ph->Hf_P = 0.0; ph->Hf_N2 = 0.0;
    ph->A = 1.0e9; ph->beta = 5.0; ph->Ta = 100.0;
    ph->T_in = 293.0; ph->Tref = ph->T_in; ph->vFuel = +0.1; ph->vOx = -0.1;
    ph->variableCpOn = true;
    update_janaf_references(ph);
    ph->YF_fuel = 1.0; ph->YO_fuel = 0.0;  ph->YP_fuel = 0.0; ph->YN2_fuel = 0.0;
    ph->YF_ox   = 0.0; ph->YO_ox   = 0.21; ph->YP_ox   = 0.0; ph->YN2_ox   = 0.79;
    ph->chemistryOn = true;
    ph->initialChemicalTimeStep = 1e-7;
    ph->maxChemicalTimeStep     = 1e-6;
    ph->dTchemMax               = 20.0;
    ph->reactedFracMax          = 1.0;
    ph->p0 = 1.0e5;

    c->startTime = 0.0; c->endTime = 1.5; c->dt0 = 1e-7;
    c->writeInterval = 0.1; c->maxCo = 0.05;
    c->nOuterCorrectors = 1; c->nCorrectors = 1;
    c->poissonIters = 4; c->scalarIters = 4;
    c->rmtMaxLevels = 4;
    c->rmtPostSmooth = 3;
    c->rmtCoarseSweeps = 48;
    c->rmtCoarseRepeats = 1;
    c->relaxP = 0.3; c->relaxU = 0.7; c->relaxSc = 0.8; c->relaxH = 0.8;
    c->progressEvery = 100;
}

static void parse_args(int argc, char **argv, Grid *g, Phys *ph, Controls *c) {
    for (int i=1; i<argc; ++i) {
        if (!strcmp(argv[i], "-Nx") && i+1<argc) g->Nx = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-Ny") && i+1<argc) g->Ny = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-end") && i+1<argc) c->endTime = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dt") && i+1<argc) c->dt0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-write") && i+1<argc) c->writeInterval = atof(argv[++i]);
        else if (!strcmp(argv[i], "-maxCo") && i+1<argc) c->maxCo = atof(argv[++i]);
        else if (!strcmp(argv[i], "-outerCorr") && i+1<argc) c->nOuterCorrectors = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-corr") && i+1<argc) c->nCorrectors = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-pCycles") && i+1<argc) c->poissonIters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-sCycles") && i+1<argc) c->scalarIters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-progressEvery") && i+1<argc) c->progressEvery = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmtLevels") && i+1<argc) c->rmtMaxLevels = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmtPost") && i+1<argc) c->rmtPostSmooth = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmtCoarseSweeps") && i+1<argc) c->rmtCoarseSweeps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmtCoarseRepeats") && i+1<argc) c->rmtCoarseRepeats = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmtOmega") && i+1<argc) g_rmtOmega = atof(argv[++i]);
        else if (!strcmp(argv[i], "-vIn") && i+1<argc) {
            double vmag = fabs(atof(argv[++i]));
            ph->vFuel = +vmag;
            ph->vOx   = -vmag;
        }
        else if (!strcmp(argv[i], "-arrA") && i+1<argc)
            ph->A = atof(argv[++i]);
        else if (!strcmp(argv[i], "-arrBeta") && i+1<argc)
            ph->beta = atof(argv[++i]);
        else if (!strcmp(argv[i], "-arrTa") && i+1<argc)
            ph->Ta = atof(argv[++i]);
        else if (!strcmp(argv[i], "-fieldOutput") && i+1<argc)
            g_fieldOutput = atoi(argv[++i]) != 0;
        else if (!strcmp(argv[i], "-perfectGas") && i+1<argc) ph->variableRhoOn = atoi(argv[++i]) != 0;
        else if (!strcmp(argv[i], "-sutherland") && i+1<argc) ph->sutherlandOn = atoi(argv[++i]) != 0;
        else if (!strcmp(argv[i], "-hydroPInRho") && i+1<argc) ph->hydroPressureInDensityOn = atoi(argv[++i]) != 0;
        else if (!strcmp(argv[i], "-rhoRelax") && i+1<argc) ph->rhoRelax = atof(argv[++i]);
    }
    g->dx = g->Lx / g->Nx; g->dy = g->Ly / g->Ny; g->x1 = 5.0*g->Lx/11.0; g->x2 = 6.0*g->Lx/11.0;
}

static void alloc_fields(const Grid *g, Fields *f) {
    size_t n = (size_t)(g->Nx+2) * (size_t)(g->Ny+2);
    #define ALLOC(name) do { f->name=(double*)calloc(n,sizeof(double)); if(!f->name) die("alloc " #name); } while(0)
    ALLOC(p); ALLOC(u); ALLOC(v); ALLOC(T); ALLOC(h); ALLOC(YF); ALLOC(YO); ALLOC(YP); ALLOC(YN2);
    ALLOC(rho); ALLOC(rhoOld); ALLOC(muMix); ALLOC(nuMix); ALLOC(alphaMix); ALLOC(DMix); ALLOC(cpMix);
    ALLOC(us); ALLOC(vs); ALLOC(rhs); ALLOC(work); ALLOC(phiOld); ALLOC(rhsScalar); ALLOC(ufx); ALLOC(vfy); ALLOC(mfx); ALLOC(mfy);
    #undef ALLOC
}
static void free_fields(Fields *f) {
    free(f->p); free(f->u); free(f->v); free(f->T); free(f->h);
    free(f->YF); free(f->YO); free(f->YP); free(f->YN2);
    free(f->rho); free(f->rhoOld); free(f->muMix); free(f->nuMix); free(f->alphaMix); free(f->DMix); free(f->cpMix);
    free(f->us); free(f->vs); free(f->rhs); free(f->work); free(f->phiOld); free(f->rhsScalar); free(f->ufx); free(f->vfy); free(f->mfx); free(f->mfy);
}

static int slot_cell(const Grid *g, int i) {
    double x = (i - 0.5) * g->dx;
    return (x >= g->x1 && x < g->x2);
}

static void init_fields(const Grid *g, const Phys *ph, Fields *f) {
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) {
        f->p [IDX(i,j,g->Ny)] = ph->p0;
        f->u [IDX(i,j,g->Ny)] = 0.0;
        f->v [IDX(i,j,g->Ny)] = 0.0;
        f->T [IDX(i,j,g->Ny)] = ph->T_in;
        f->h [IDX(i,j,g->Ny)] = 0.0;
        f->YF[IDX(i,j,g->Ny)] = 0.0;
        f->YO[IDX(i,j,g->Ny)] = 0.0;
        f->YP[IDX(i,j,g->Ny)] = 0.0;
        f->YN2[IDX(i,j,g->Ny)] = 1.0;
        f->rho[IDX(i,j,g->Ny)] = ph->rho;
        f->rhoOld[IDX(i,j,g->Ny)] = ph->rho;
        f->muMix[IDX(i,j,g->Ny)] = ph->mu;
        f->nuMix[IDX(i,j,g->Ny)] = ph->nu;
        f->alphaMix[IDX(i,j,g->Ny)] = ph->alpha;
        f->DMix[IDX(i,j,g->Ny)] = ph->D;
        f->cpMix[IDX(i,j,g->Ny)] = ph->cp;
    }
}

static inline double T_from_h(const Phys *ph, double h) { return ph->Tref + h/ph->cp; }
static inline double h_from_T(const Phys *ph, double T) { return ph->cp*(T - ph->Tref); }

static void relax_field(const Grid *g, double *phi, const double *old, double alpha) {
    if (alpha >= 0.999999) return;
    const int nx=g->Nx, ny=g->Ny;
    for (int i=1;i<=nx;++i) for (int j=1;j<=ny;++j) {
        double v = alpha*phi[IDX(i,j,ny)] + (1.0-alpha)*old[IDX(i,j,ny)];
        phi[IDX(i,j,ny)] = isfinite(v) ? v : old[IDX(i,j,ny)];
    }
}


/* ---- OpenFOAM-style patch-face approximations on a collocated grid ---- */
static double p_patch_left(const Grid *g, const Phys *ph, const Fields *f, int j) {
    (void)g;
    double uN = f->u[IDX(1,j,g->Ny)];
    if (uN < 0.0) return ph->p0; /* outflow */
    return ph->p0 - 0.5*f->rho[IDX(1,j,g->Ny)]*(uN*uN); /* inflow: use normal component only */
}
static double p_patch_right(const Grid *g, const Phys *ph, const Fields *f, int j) {
    (void)g;
    double uN = f->u[IDX(g->Nx,j,g->Ny)];
    if (uN > 0.0) return ph->p0; /* outflow */
    return ph->p0 - 0.5*f->rho[IDX(g->Nx,j,g->Ny)]*(uN*uN);
}

static void apply_bc_velocity_pressure(const Grid *g, const Phys *ph, Fields *f) {
    const int nx = g->Nx, ny = g->Ny;
    for (int j=1; j<=ny; ++j) {
        /* outlet sides: totalPressure + pressureInletOutletVelocity(value=(0 0)) */
        double pL = p_patch_left(g, ph, f, j), pR = p_patch_right(g, ph, f, j);
        f->p[IDX(0,    j,ny)] = 2.0*pL - f->p[IDX(1, j,ny)];
        f->p[IDX(nx+1, j,ny)] = 2.0*pR - f->p[IDX(nx,j,ny)];

        if (f->u[IDX(1,j,ny)] < 0.0) { /* left outflow */
            f->u[IDX(0,j,ny)] = f->u[IDX(1,j,ny)];
            f->v[IDX(0,j,ny)] = f->v[IDX(1,j,ny)];
        } else { /* left inflow */
            f->u[IDX(0,j,ny)] = -f->u[IDX(1,j,ny)];
            f->v[IDX(0,j,ny)] = -f->v[IDX(1,j,ny)];
        }
        if (f->u[IDX(nx,j,ny)] > 0.0) { /* right outflow */
            f->u[IDX(nx+1,j,ny)] = f->u[IDX(nx,j,ny)];
            f->v[IDX(nx+1,j,ny)] = f->v[IDX(nx,j,ny)];
        } else {
            f->u[IDX(nx+1,j,ny)] = -f->u[IDX(nx,j,ny)];
            f->v[IDX(nx+1,j,ny)] = -f->v[IDX(nx,j,ny)];
        }
    }

    for (int i=1; i<=nx; ++i) {
        int slot = slot_cell(g, i);
        /* p: zeroGradient on fuel/oxidizer/walls */
        f->p[IDX(i,0,    ny)] = f->p[IDX(i,1, ny)];
        f->p[IDX(i,ny+1, ny)] = f->p[IDX(i,ny,ny)];
        if (slot) {
            /* fixedValue U=(0,vIn) -> ghost = 2*value - interior */
            f->u[IDX(i,0,ny)]    = -f->u[IDX(i,1,ny)];
            f->v[IDX(i,0,ny)]    = 2.0*ph->vFuel - f->v[IDX(i,1,ny)];
            f->u[IDX(i,ny+1,ny)] = -f->u[IDX(i,ny,ny)];
            f->v[IDX(i,ny+1,ny)] = 2.0*ph->vOx   - f->v[IDX(i,ny,ny)];
        } else {
            /* no-slip wall */
            f->u[IDX(i,0,ny)]    = -f->u[IDX(i,1,ny)];
            f->v[IDX(i,0,ny)]    = -f->v[IDX(i,1,ny)];
            f->u[IDX(i,ny+1,ny)] = -f->u[IDX(i,ny,ny)];
            f->v[IDX(i,ny+1,ny)] = -f->v[IDX(i,ny,ny)];
        }
    }
}

static void apply_bc_scalar_of(const Grid *g, const Fields *f, double *s,
                               double bottomSlot, double topSlot, double outletInletValue) {
    const int nx = g->Nx, ny = g->Ny;
    for (int j=1; j<=ny; ++j) {
        /* outlet = inletOutlet */
        if (f->u[IDX(1,j,ny)] < 0.0)
            s[IDX(0,j,ny)] = s[IDX(1,j,ny)]; /* outflow zeroGradient */
        else
            s[IDX(0,j,ny)] = 2.0*outletInletValue - s[IDX(1,j,ny)]; /* inflow fixed inletValue */
        if (f->u[IDX(nx,j,ny)] > 0.0)
            s[IDX(nx+1,j,ny)] = s[IDX(nx,j,ny)];
        else
            s[IDX(nx+1,j,ny)] = 2.0*outletInletValue - s[IDX(nx,j,ny)];
    }
    for (int i=1; i<=nx; ++i) {
        int slot = slot_cell(g, i);
        if (slot) {
            s[IDX(i,0,ny)]    = 2.0*bottomSlot - s[IDX(i,1,ny)];
            s[IDX(i,ny+1,ny)] = 2.0*topSlot    - s[IDX(i,ny,ny)];
        } else {
            s[IDX(i,0,ny)]    = s[IDX(i,1,ny)];
            s[IDX(i,ny+1,ny)] = s[IDX(i,ny,ny)];
        }
    }
}

static double limited_linear_psi(double r) { return clamp(MIN(2.0*r, 1.0), 0.0, 1.0); }

static double bounded_between(double val, double a, double b) {
    double lo = a < b ? a : b;
    double hi = a > b ? a : b;
    if (val < lo) return lo;
    if (val > hi) return hi;
    return val;
}

static double scalar_face_x_raw(const Grid *g, const Fields *f, const double *phi,
                            int iface, int j, double inletValue) {
    const int ny = g->Ny, nx = g->Nx;
    if (iface == 0) {
        /* left outlet patch */
        return (f->u[IDX(1,j,ny)] < 0.0) ? phi[IDX(1,j,ny)] : inletValue;
    }
    if (iface == nx) {
        return (f->u[IDX(nx,j,ny)] > 0.0) ? phi[IDX(nx,j,ny)] : inletValue;
    }
    double uf = 0.5*(f->u[IDX(iface,j,ny)] + f->u[IDX(iface+1,j,ny)]);
    const double eps = 1e-30;
    if (uf >= 0.0) {
        double phiUU = phi[IDX(MAX(iface-1,0), j, ny)];
        double phiU  = phi[IDX(iface, j, ny)];
        double phiD  = phi[IDX(iface+1, j, ny)];
        double denom = phiD - phiU;
        double r = (phiU - phiUU) / (fabs(denom) > eps ? denom : copysign(eps, denom + eps));
        double psi = limited_linear_psi(r);
        return phiU + 0.5*psi*(phiD - phiU);
    } else {
        double phiUU = phi[IDX(MIN(iface+2,nx+1), j, ny)];
        double phiU  = phi[IDX(iface+1, j, ny)];
        double phiD  = phi[IDX(iface, j, ny)];
        double denom = phiD - phiU;
        double r = (phiU - phiUU) / (fabs(denom) > eps ? denom : copysign(eps, denom + eps));
        double psi = limited_linear_psi(r);
        return phiU + 0.5*psi*(phiD - phiU);
    }
}

static double scalar_face_y_raw(const Grid *g, const Phys *ph, const Fields *f, const double *phi,
                            int i, int jface, double bottomSlot, double topSlot) {
    const int ny = g->Ny;
    if (jface == 0) {
        return slot_cell(g,i) ? bottomSlot : phi[IDX(i,1,ny)];
    }
    if (jface == ny) {
        return slot_cell(g,i) ? topSlot : phi[IDX(i,ny,ny)];
    }
    double vf = 0.5*(f->v[IDX(i,jface,ny)] + f->v[IDX(i,jface+1,ny)]);
    const double eps = 1e-30;
    if (vf >= 0.0) {
        double phiUU = phi[IDX(i, MAX(jface-1,0), ny)];
        double phiU  = phi[IDX(i, jface, ny)];
        double phiD  = phi[IDX(i, jface+1, ny)];
        double denom = phiD - phiU;
        double r = (phiU - phiUU) / (fabs(denom) > eps ? denom : copysign(eps, denom + eps));
        double psi = limited_linear_psi(r);
        return phiU + 0.5*psi*(phiD - phiU);
    } else {
        double phiUU = phi[IDX(i, MIN(jface+2,ny+1), ny)];
        double phiU  = phi[IDX(i, jface+1, ny)];
        double phiD  = phi[IDX(i, jface, ny)];
        double denom = phiD - phiU;
        double r = (phiU - phiUU) / (fabs(denom) > eps ? denom : copysign(eps, denom + eps));
        double psi = limited_linear_psi(r);
        return phiU + 0.5*psi*(phiD - phiU);
    }
}


static double scalar_face_x_species(const Grid *g, const Fields *f, const double *phi,
                                    int iface, int j, double inletValue) {
    double raw = scalar_face_x_raw(g, f, phi, iface, j, inletValue);
    const int ny=g->Ny, nx=g->Nx;
    if (iface == 0) {
        return bounded_between(raw, phi[IDX(1,j,ny)], inletValue);
    }
    if (iface == nx) {
        return bounded_between(raw, phi[IDX(nx,j,ny)], inletValue);
    }
    double a = phi[IDX(iface,   j, ny)];
    double b = phi[IDX(iface+1, j, ny)];
    return bounded_between(raw, a, b);
}

static double scalar_face_y_species(const Grid *g, const Phys *ph, const Fields *f, const double *phi,
                                    int i, int jface, double bottomSlot, double topSlot) {
    double raw = scalar_face_y_raw(g, ph, f, phi, i, jface, bottomSlot, topSlot);
    const int ny=g->Ny;
    if (jface == 0) {
        double faceVal = slot_cell(g,i) ? bottomSlot : phi[IDX(i,1,ny)];
        return bounded_between(raw, phi[IDX(i,1,ny)], faceVal);
    }
    if (jface == ny) {
        double faceVal = slot_cell(g,i) ? topSlot : phi[IDX(i,ny,ny)];
        return bounded_between(raw, phi[IDX(i,ny,ny)], faceVal);
    }
    double a = phi[IDX(i, jface,   ny)];
    double b = phi[IDX(i, jface+1, ny)];
    return bounded_between(raw, a, b);
}

static double scalar_face_x_h(const Grid *g, const Fields *f, const double *phi,
                              int iface, int j) {
    const int ny=g->Ny, nx=g->Nx;
    if (iface == 0)  return phi[IDX(1, j, ny)];
    if (iface == nx) return phi[IDX(nx,j, ny)];
    return scalar_face_x_raw(g, f, phi, iface, j, 0.0);
}

static double scalar_face_y_h(const Grid *g, const Phys *ph, const Fields *f, const double *phi,
                              int i, int jface) {
    return scalar_face_y_raw(g, ph, f, phi, i, jface, 0.0, 0.0);
}

static void apply_bc_h_exact(const Grid *g, const double *vfield, double *h) {
    const int nx=g->Nx, ny=g->Ny;
    for (int j=1; j<=ny; ++j) {
        h[IDX(0,    j,ny)] = h[IDX(1, j,ny)];
        h[IDX(nx+1, j,ny)] = h[IDX(nx,j,ny)];
    }
    for (int i=1; i<=nx; ++i) {
        int slot = slot_cell(g, i);
        if (slot) {
            h[IDX(i,0,ny)]    = -h[IDX(i,1,ny)];
            h[IDX(i,ny+1,ny)] = -h[IDX(i,ny,ny)];
        } else {
            h[IDX(i,0,ny)]    = h[IDX(i,1,ny)];
            h[IDX(i,ny+1,ny)] = h[IDX(i,ny,ny)];
        }
    }
    (void)vfield;
}

static void build_h_rhs_openfoam_like(const Grid *g, const Phys *ph, const Fields *f,
                                      const double *phi, double *rhs, double dt) {
    const int nx=g->Nx, ny=g->Ny; const double dx=g->dx, dy=g->dy;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double Fe=u_face_x(g,f,i,j), Fw=u_face_x(g,f,i-1,j), Fn=v_face_y(g,ph,f,i,j), Fs=v_face_y(g,ph,f,i,j-1);
        double phiE=scalar_face_x_h(g,f,phi,i,j);
        double phiW=scalar_face_x_h(g,f,phi,i-1,j);
        double phiN=scalar_face_y_h(g,ph,f,phi,i,j);
        double phiS=scalar_face_y_h(g,ph,f,phi,i,j-1);
        rhs[IDX(i,j,ny)] = phi[IDX(i,j,ny)] - dt*((Fe*phiE - Fw*phiW)/dx + (Fn*phiN - Fs*phiS)/dy);
    }
}

static void solve_h_helmholtz(const Grid *g, const Phys *ph, const Fields *f, double *phi, const double *rhs,
                              double dt, double diff, int iters) {
    const int nx=g->Nx, ny=g->Ny; const double ax=dt*diff/(g->dx*g->dx), ay=dt*diff/(g->dy*g->dy); const double ap=1.0+2.0*ax+2.0*ay;
    for (int it=0; it<iters; ++it) {
        apply_bc_h_exact(g, f->v, phi);
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            double val = (rhs[IDX(i,j,ny)] + ax*(phi[IDX(i+1,j,ny)] + phi[IDX(i-1,j,ny)]) + ay*(phi[IDX(i,j+1,ny)] + phi[IDX(i,j-1,ny)]))/ap;
            phi[IDX(i,j,ny)] = MAX(val, 0.0);
        }
    }
    apply_bc_h_exact(g, f->v, phi);
}


typedef enum {
    RMT_BC_PRESSURE = 0,
    RMT_BC_SCALAR   = 1
} RMTBCType;

typedef struct {
    int maxLevels;
    int postSmooth;
    int coarseSweeps;
    int coarseRepeats;
    RMTBCType bcType;
    double ax, ay, ap;
    double refDx, refDy;
} RMTSettings;

static void apply_bc_rmt_correction(const Grid *g, double *phi, RMTBCType bcType) {
    const int nx = g->Nx, ny = g->Ny;
    for (int j=1; j<=ny; ++j) {
        if (bcType == RMT_BC_PRESSURE) {
            phi[IDX(0,    j, ny)] = -phi[IDX(1,  j, ny)];
            phi[IDX(nx+1, j, ny)] = -phi[IDX(nx, j, ny)];
        } else {
            phi[IDX(0,    j, ny)] = phi[IDX(1,  j, ny)];
            phi[IDX(nx+1, j, ny)] = phi[IDX(nx, j, ny)];
        }
    }
    for (int i=1; i<=nx; ++i) {
        if (bcType == RMT_BC_PRESSURE) {
            phi[IDX(i, 0,    ny)] = phi[IDX(i, 1,  ny)];
            phi[IDX(i, ny+1, ny)] = phi[IDX(i, ny, ny)];
        } else {
            int slot = slot_cell(g, i);
            if (slot) {
                phi[IDX(i, 0,    ny)] = -phi[IDX(i, 1,  ny)];
                phi[IDX(i, ny+1, ny)] = -phi[IDX(i, ny, ny)];
            } else {
                phi[IDX(i, 0,    ny)] = phi[IDX(i, 1,  ny)];
                phi[IDX(i, ny+1, ny)] = phi[IDX(i, ny, ny)];
            }
        }
    }
    phi[IDX(0,      0,      ny)] = phi[IDX(1,   1,   ny)];
    phi[IDX(0,      ny+1,   ny)] = phi[IDX(1,   ny,  ny)];
    phi[IDX(nx+1,   0,      ny)] = phi[IDX(nx,  1,   ny)];
    phi[IDX(nx+1,   ny+1,   ny)] = phi[IDX(nx,  ny,  ny)];
}

static void rmt_smooth_gs(const Grid *g, double *x, const double *b,
                          const RMTSettings *st, int sweeps) {
    const int nx = g->Nx, ny = g->Ny;
    for (int it=0; it<sweeps; ++it) {
        apply_bc_rmt_correction(g, x, st->bcType);
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            x[IDX(i,j,ny)] = (b[IDX(i,j,ny)]
                            + st->ax*(x[IDX(i+1,j,ny)] + x[IDX(i-1,j,ny)])
                            + st->ay*(x[IDX(i,j+1,ny)] + x[IDX(i,j-1,ny)])) / MAX(st->ap, 1.0e-20);
        }
    }
    apply_bc_rmt_correction(g, x, st->bcType);
}

static void rmt_compute_residual(const Grid *g, const double *x, const double *b,
                                 const RMTSettings *st, double *res) {
    const int nx = g->Nx, ny = g->Ny;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double Ax = st->ap*x[IDX(i,j,ny)]
                  - st->ax*(x[IDX(i+1,j,ny)] + x[IDX(i-1,j,ny)])
                  - st->ay*(x[IDX(i,j+1,ny)] + x[IDX(i,j-1,ny)]);
        res[IDX(i,j,ny)] = b[IDX(i,j,ny)] - Ax;
    }
}

static int rmt_shift_count(int N, int s) {
    int first = s + 1;
    if (first > N) return 0;
    return 1 + (N - first) / 3;
}

static Grid rmt_make_shifted_coarse_grid(const Grid *fine, int sx, int sy) {
    Grid c = *fine;
    c.Nx = rmt_shift_count(fine->Nx, sx);
    c.Ny = rmt_shift_count(fine->Ny, sy);
    c.dx = 3.0 * fine->dx;
    c.dy = 3.0 * fine->dy;
    return c;
}

static void rmt_build_level_coeffs(const Grid *g, const RMTSettings *fineSt, RMTSettings *st) {
    *st = *fineSt;
    if (fineSt->bcType == RMT_BC_PRESSURE) {
        st->ax = 1.0/(g->dx*g->dx);
        st->ay = 1.0/(g->dy*g->dy);
        st->ap = 2.0*st->ax + 2.0*st->ay;
    } else {
        const double dtRhoDiffX = fineSt->ax * fineSt->refDx * fineSt->refDx;
        const double dtRhoDiffY = fineSt->ay * fineSt->refDy * fineSt->refDy;
        const double rhoDiag = fineSt->ap - 2.0*fineSt->ax - 2.0*fineSt->ay;
        st->ax = dtRhoDiffX/(g->dx*g->dx);
        st->ay = dtRhoDiffY/(g->dy*g->dy);
        st->ap = rhoDiag + 2.0*st->ax + 2.0*st->ay;
    }
}

static void rmt_restrict_shifted_avg(const Grid *fine, int sx, int sy,
                                     const double *rf, const Grid *coarse, double *rc) {
    const int nyc = coarse->Ny;
    memset(rc, 0, (size_t)(coarse->Nx+2)*(size_t)(coarse->Ny+2)*sizeof(double));
    for (int ic=1; ic<=coarse->Nx; ++ic) {
        int i0 = sx + 1 + 3*(ic-1);
        for (int jc=1; jc<=coarse->Ny; ++jc) {
            int j0 = sy + 1 + 3*(jc-1);
            double sum = 0.0, cnt = 0.0;
            for (int i=MAX(1, i0-1); i<=MIN(fine->Nx, i0+1); ++i)
                for (int j=MAX(1, j0-1); j<=MIN(fine->Ny, j0+1); ++j) {
                    sum += rf[IDX(i,j,fine->Ny)];
                    cnt += 1.0;
                }
            rc[IDX(ic,jc,nyc)] = (cnt > 0.0 ? sum/cnt : 0.0);
        }
    }
}

static void rmt_inject_shifted_add(const Grid *fine, int sx, int sy,
                                   double *xf, const Grid *coarse, const double *xc) {
    const int nyf = fine->Ny;
    for (int ic=1; ic<=coarse->Nx; ++ic) {
        int i0 = sx + 1 + 3*(ic-1);
        for (int jc=1; jc<=coarse->Ny; ++jc) {
            int j0 = sy + 1 + 3*(jc-1);
            xf[IDX(i0,j0,nyf)] += xc[IDX(ic,jc,coarse->Ny)];
        }
    }
}

static void rmt_multiple_recursive(const Grid *g, double *x, const double *b,
                                   const RMTSettings *fineSt, int level) {
    RMTSettings st;
    rmt_build_level_coeffs(g, fineSt, &st);
    const int coarsest = (level+1 >= fineSt->maxLevels || g->Nx < 4 || g->Ny < 4);
    if (coarsest) {
        rmt_smooth_gs(g, x, b, &st, fineSt->coarseSweeps);
        return;
    }

    size_t n = (size_t)(g->Nx+2) * (size_t)(g->Ny+2);
    double *res = (double*)calloc(n, sizeof(double));
    if (!res) die("alloc rmt residual");

    apply_bc_rmt_correction(g, x, st.bcType);
    rmt_compute_residual(g, x, b, &st, res);

    for (int rep=0; rep<fineSt->coarseRepeats; ++rep) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if(level == 0)
#endif
        for (int sx=0; sx<3; ++sx) for (int sy=0; sy<3; ++sy) {
            Grid gc = rmt_make_shifted_coarse_grid(g, sx, sy);
            if (gc.Nx < 1 || gc.Ny < 1) continue;
            size_t nc = (size_t)(gc.Nx+2) * (size_t)(gc.Ny+2);
            double *bc = (double*)calloc(nc, sizeof(double));
            double *xc = (double*)calloc(nc, sizeof(double));
            if (!bc || !xc) die("alloc rmt shifted coarse arrays");
            rmt_restrict_shifted_avg(g, sx, sy, res, &gc, bc);
            rmt_multiple_recursive(&gc, xc, bc, fineSt, level+1);
            rmt_inject_shifted_add(g, sx, sy, x, &gc, xc);
            free(bc);
            free(xc);
        }
    }

    rmt_smooth_gs(g, x, b, &st, fineSt->postSmooth);
    free(res);
}

static double rmt_residual_norm_pressure(const Grid *g, const Fields *f) {
    const int nx = g->Nx, ny = g->Ny;
    const double ax = 1.0/(g->dx*g->dx), ay = 1.0/(g->dy*g->dy), ap = 2.0*ax + 2.0*ay;
    double sum = 0.0;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double Ax = ap*f->p[IDX(i,j,ny)] - ax*(f->p[IDX(i+1,j,ny)] + f->p[IDX(i-1,j,ny)])
                                   - ay*(f->p[IDX(i,j+1,ny)] + f->p[IDX(i,j-1,ny)]);
        double r = (-f->rhs[IDX(i,j,ny)]) - Ax;
        sum += r*r;
    }
    return sqrt(sum/((double)nx*(double)ny));
}

static double rmt_residual_norm_scalar(const Grid *g, const double *phi, const double *rhs,
                                       double dt, double diff) {
    const int nx = g->Nx, ny = g->Ny;
    const double ax = dt*diff/(g->dx*g->dx), ay = dt*diff/(g->dy*g->dy), ap = 1.0 + 2.0*ax + 2.0*ay;
    double sum = 0.0;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double Ax = ap*phi[IDX(i,j,ny)] - ax*(phi[IDX(i+1,j,ny)] + phi[IDX(i-1,j,ny)])
                                    - ay*(phi[IDX(i,j+1,ny)] + phi[IDX(i,j-1,ny)]);
        double r = rhs[IDX(i,j,ny)] - Ax;
        sum += r*r;
    }
    return sqrt(sum/((double)nx*(double)ny));
}

static inline double gradp_x_cell(const Grid *g, const Fields *f, int i, int j) {
    const int ny=g->Ny, nx=g->Nx;
    if (i==1)  return (f->p[IDX(2,j,ny)]    - f->p[IDX(1,j,ny)]) / g->dx;
    if (i==nx) return (f->p[IDX(nx,j,ny)]   - f->p[IDX(nx-1,j,ny)]) / g->dx;
    return (f->p[IDX(i+1,j,ny)] - f->p[IDX(i-1,j,ny)]) / (2.0*g->dx);
}
static inline double gradp_y_cell(const Grid *g, const Fields *f, int i, int j) {
    const int ny=g->Ny;
    if (j==1)  return (f->p[IDX(i,2,ny)]    - f->p[IDX(i,1,ny)]) / g->dy;
    if (j==ny) return (f->p[IDX(i,ny,ny)]   - f->p[IDX(i,ny-1,ny)]) / g->dy;
    return (f->p[IDX(i,j+1,ny)] - f->p[IDX(i,j-1,ny)]) / (2.0*g->dy);
}

static void compute_face_flux_rc(const Grid *g, const Phys *ph, Fields *f, double dt, int usePredicted) {
    const int nx=g->Nx, ny=g->Ny;
    double *uc = usePredicted ? f->us : f->u;
    double *vc = usePredicted ? f->vs : f->v;
    const int applyRC = usePredicted ? 0 : 1;
    for (int j=1; j<=ny; ++j) {
        f->ufx[IDX(0,j,ny)] = uc[IDX(1,j,ny)];
        f->mfx[IDX(0,j,ny)] = rho_face_x(g, f, 0, j) * f->ufx[IDX(0,j,ny)];
        f->ufx[IDX(nx,j,ny)] = uc[IDX(nx,j,ny)];
        f->mfx[IDX(nx,j,ny)] = rho_face_x(g, f, nx, j) * f->ufx[IDX(nx,j,ny)];
        for (int i=1; i<nx; ++i) {
            double uf = 0.5*(uc[IDX(i,j,ny)] + uc[IDX(i+1,j,ny)]);
            if (applyRC) {
                double rhoFace = MAX(rho_face_x(g, f, i, j), 1.0e-12);
                double dpface = (f->p[IDX(i+1,j,ny)] - f->p[IDX(i,j,ny)]) / g->dx;
                double gpbar  = 0.5*(gradp_x_cell(g,f,i,j) + gradp_x_cell(g,f,i+1,j));
                uf -= (dt/rhoFace)*(dpface - gpbar);
            }
            f->ufx[IDX(i,j,ny)] = uf;
            f->mfx[IDX(i,j,ny)] = rho_face_x(g, f, i, j) * uf;
        }
    }
    for (int i=1; i<=nx; ++i) {
        f->vfy[IDX(i,0,ny)] = slot_cell(g,i) ? ph->vFuel : 0.0;
        f->mfy[IDX(i,0,ny)] = rho_face_y(g, f, i, 0) * f->vfy[IDX(i,0,ny)];
        f->vfy[IDX(i,ny,ny)] = slot_cell(g,i) ? ph->vOx : 0.0;
        f->mfy[IDX(i,ny,ny)] = rho_face_y(g, f, i, ny) * f->vfy[IDX(i,ny,ny)];
        for (int j=1; j<ny; ++j) {
            double vf = 0.5*(vc[IDX(i,j,ny)] + vc[IDX(i,j+1,ny)]);
            if (applyRC) {
                double rhoFace = MAX(rho_face_y(g, f, i, j), 1.0e-12);
                double dpface = (f->p[IDX(i,j+1,ny)] - f->p[IDX(i,j,ny)]) / g->dy;
                double gpbar  = 0.5*(gradp_y_cell(g,f,i,j) + gradp_y_cell(g,f,i,j+1));
                vf -= (dt/rhoFace)*(dpface - gpbar);
            }
            f->vfy[IDX(i,j,ny)] = vf;
            f->mfy[IDX(i,j,ny)] = rho_face_y(g, f, i, j) * vf;
        }
    }
}

static double u_face_x(const Grid *g, const Fields *f, int iface, int j) {
    (void)g;
    return f->ufx[IDX(iface,j,g->Ny)];
}
static double v_face_y(const Grid *g, const Phys *ph, const Fields *f, int i, int jface) {
    (void)g; (void)ph;
    return f->vfy[IDX(i,jface,g->Ny)];
}
static double rho_face_x(const Grid *g, const Fields *f, int iface, int j) {
    if (iface <= 0) return f->rho[IDX(1,j,g->Ny)];
    if (iface >= g->Nx) return f->rho[IDX(g->Nx,j,g->Ny)];
    return 0.5*(f->rho[IDX(iface,j,g->Ny)] + f->rho[IDX(iface+1,j,g->Ny)]);
}
static double rho_face_y(const Grid *g, const Fields *f, int i, int jface) {
    if (jface <= 0) return f->rho[IDX(i,1,g->Ny)];
    if (jface >= g->Ny) return f->rho[IDX(i,g->Ny,g->Ny)];
    return 0.5*(f->rho[IDX(i,jface,g->Ny)] + f->rho[IDX(i,jface+1,g->Ny)]);
}
static double v_face_x(const Grid *g, const Fields *f, int iface, int j) {
    const int ny = g->Ny, nx = g->Nx;
    /* Same narrow V64C patch for tangential velocity transport at the open side:
       use the internal boundary-cell value at the face instead of clipping to 0. */
    if (iface == 0)  return f->v[IDX(1, j, ny)];
    if (iface == nx) return f->v[IDX(nx,j, ny)];
    return 0.5*(f->v[IDX(iface,j,ny)] + f->v[IDX(iface+1,j,ny)]);
}
static double u_face_y(const Grid *g, const Fields *f, int i, int jface) {
    const int ny = g->Ny;
    if (jface == 0 || jface == ny) return 0.0;
    return 0.5*(f->u[IDX(i,jface,ny)] + f->u[IDX(i,jface+1,ny)]);
}

static void momentum_predictor(const Grid *g, const Phys *ph, Fields *f, double dt, double relaxU) {
    const int nx = g->Nx, ny = g->Ny;
    const double dx = g->dx, dy = g->dy;
    const double dx2 = dx*dx, dy2 = dy*dy;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double uij=f->u[IDX(i,j,ny)], vij=f->v[IDX(i,j,ny)];
        double rhoi = MAX(f->rho[IDX(i,j,ny)], 1.0e-12);

        double ue=scalar_face_x_raw(g,f,f->u,i,j,0.0),   uw=scalar_face_x_raw(g,f,f->u,i-1,j,0.0);
        double un=scalar_face_y_raw(g,ph,f,f->u,i,j,0.0,0.0), us=scalar_face_y_raw(g,ph,f,f->u,i,j-1,0.0,0.0);
        double ve=scalar_face_x_raw(g,f,f->v,i,j,0.0),   vw=scalar_face_x_raw(g,f,f->v,i-1,j,0.0);
        double vn=scalar_face_y_raw(g,ph,f,f->v,i,j,0.0,0.0), vs=scalar_face_y_raw(g,ph,f,f->v,i,j-1,0.0,0.0);

        double mE=f->mfx[IDX(i,  j,ny)], mW=f->mfx[IDX(i-1,j,ny)];
        double mN=f->mfy[IDX(i,j,  ny)], mS=f->mfy[IDX(i,j-1,ny)];

        double du_adv = -((mE*ue - mW*uw)/dx + (mN*un - mS*us)/dy) / rhoi;
        double dv_adv = -((mE*ve - mW*vw)/dx + (mN*vn - mS*vs)/dy) / rhoi;

        double muE = face_avg_x_coeff(g, f->muMix, i,   j);
        double muW = face_avg_x_coeff(g, f->muMix, i-1, j);
        double muN = face_avg_y_coeff(g, f->muMix, i,   j);
        double muS = face_avg_y_coeff(g, f->muMix, i, j-1);
        double du_diff = ((muE*(f->u[IDX(i+1,j,ny)] - uij) - muW*(uij - f->u[IDX(i-1,j,ny)]))/dx2
                        + (muN*(f->u[IDX(i,j+1,ny)] - uij) - muS*(uij - f->u[IDX(i,j-1,ny)]))/dy2) / rhoi;
        double dv_diff = ((muE*(f->v[IDX(i+1,j,ny)] - vij) - muW*(vij - f->v[IDX(i-1,j,ny)]))/dx2
                        + (muN*(f->v[IDX(i,j+1,ny)] - vij) - muS*(vij - f->v[IDX(i,j-1,ny)]))/dy2) / rhoi;

        double dpdx = (f->p[IDX(i+1,j,ny)] - f->p[IDX(i-1,j,ny)])/(2.0*dx);
        double dpdy = (f->p[IDX(i,j+1,ny)] - f->p[IDX(i,j-1,ny)])/(2.0*dy);

        f->us[IDX(i,j,ny)] = relaxU*(uij + dt*(du_adv + du_diff - dpdx/rhoi)) + (1.0-relaxU)*uij;
        f->vs[IDX(i,j,ny)] = relaxU*(vij + dt*(dv_adv + dv_diff - dpdy/rhoi)) + (1.0-relaxU)*vij;
    }
}

static void build_pressure_rhs(const Grid *g, const Phys *ph, Fields *f, double dt) {
    const int nx = g->Nx, ny = g->Ny;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double mE = f->mfx[IDX(i,  j,ny)];
        double mW = f->mfx[IDX(i-1,j,ny)];
        double mN = f->mfy[IDX(i,j,  ny)];
        double mS = f->mfy[IDX(i,j-1,ny)];
        double divm = (mE - mW)/g->dx + (mN - mS)/g->dy;
        double drhoDt = ph->variableRhoOn ? (f->rho[IDX(i,j,ny)] - f->rhoOld[IDX(i,j,ny)]) / MAX(dt,1.0e-20) : 0.0;
        f->rhs[IDX(i,j,ny)] = (divm + drhoDt) / MAX(dt,1.0e-20);
    }
}

static inline double face_avg_x_coeff(const Grid *g, const double *a, int i, int j) {
    if (i <= 0) return a[IDX(1,j,g->Ny)];
    if (i >= g->Nx) return a[IDX(g->Nx,j,g->Ny)];
    return 0.5*(a[IDX(i,j,g->Ny)] + a[IDX(i+1,j,g->Ny)]);
}

static inline double face_avg_y_coeff(const Grid *g, const double *a, int i, int j) {
    if (j <= 0) return a[IDX(i,1,g->Ny)];
    if (j >= g->Ny) return a[IDX(i,g->Ny,g->Ny)];
    return 0.5*(a[IDX(i,j,g->Ny)] + a[IDX(i,j+1,g->Ny)]);
}

static double residual_norm_pressure_variable(const Grid *g, const Fields *f, const double *rhs) {
    const double dx2 = g->dx*g->dx, dy2 = g->dy*g->dy;
    double maxRes = 0.0;
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) {
        double lap = (f->p[IDX(i+1,j,g->Ny)] - 2.0*f->p[IDX(i,j,g->Ny)] + f->p[IDX(i-1,j,g->Ny)])/dx2
                   + (f->p[IDX(i,j+1,g->Ny)] - 2.0*f->p[IDX(i,j,g->Ny)] + f->p[IDX(i,j-1,g->Ny)])/dy2;
        double r = fabs(lap - rhs[IDX(i,j,g->Ny)]);
        if (r > maxRes) maxRes = r;
    }
    return maxRes;
}

static double residual_norm_scalar_variable(const Grid *g, const Fields *f, const double *phi, const double *rhs,
                                            const double *diffField, double dt) {
    const double dx2 = g->dx*g->dx, dy2 = g->dy*g->dy;
    double maxRes = 0.0;
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) {
        double rhoP = MAX(f->rho[IDX(i,j,g->Ny)], 1.0e-12);
        double gE = rho_face_x(g, f, i,   j) * face_avg_x_coeff(g, diffField, i,   j) / dx2;
        double gW = rho_face_x(g, f, i-1, j) * face_avg_x_coeff(g, diffField, i-1, j) / dx2;
        double gN = rho_face_y(g, f, i,   j) * face_avg_y_coeff(g, diffField, i,   j) / dy2;
        double gS = rho_face_y(g, f, i, j-1) * face_avg_y_coeff(g, diffField, i, j-1) / dy2;
        double ap = rhoP + dt*(gE + gW + gN + gS);
        double Aphi = ap*phi[IDX(i,j,g->Ny)]
                    - dt*(gE*phi[IDX(i+1,j,g->Ny)] + gW*phi[IDX(i-1,j,g->Ny)]
                        + gN*phi[IDX(i,j+1,g->Ny)] + gS*phi[IDX(i,j-1,g->Ny)]);
        double r = fabs(Aphi - rhs[IDX(i,j,g->Ny)]);
        if (r > maxRes) maxRes = r;
    }
    return maxRes;
}

static void solve_pressure_poisson(const Grid *g, const Phys *ph, Fields *f, int iters, double relaxP) {
    const double __rmt_pressure_t0 = rmt_wall_now();
    const int nx=g->Nx, ny=g->Ny;
    const double dx2 = g->dx*g->dx, dy2 = g->dy*g->dy;
    const double ae = 1.0/dx2, aw = 1.0/dx2, an = 1.0/dy2, as = 1.0/dy2;
    const double ap = ae + aw + an + as;
    const int sweepsPerCycle = MAX(4, g_rmtPostSmooth + g_rmtCoarseRepeats + 1);
    const double omegaRMT = g_rmtOmega;
    size_t n = (size_t)(nx+2) * (size_t)(ny+2);
    double *corr = (double*)calloc(n, sizeof(double));
    double *res  = (double*)calloc(n, sizeof(double));
    if (!corr || !res) die("alloc pressure RMT workspace");
    (void)ph;

    RMTSettings st;
    st.maxLevels = g_rmtMaxLevels;
    st.postSmooth = g_rmtPostSmooth;
    st.coarseSweeps = g_rmtCoarseSweeps;
    st.coarseRepeats = g_rmtCoarseRepeats;
    st.bcType = RMT_BC_PRESSURE;
    st.ax = ae; st.ay = an; st.ap = ap; st.refDx = g->dx; st.refDy = g->dy;

    for (int cyc=0; cyc<iters; ++cyc) {
        apply_bc_velocity_pressure(g, ph, f);
        double oldRes = residual_norm_pressure_variable(g, f, f->rhs);
        memset(corr, 0, n*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
            double Ap = ap*f->p[IDX(i,j,ny)]
                      - ae*(f->p[IDX(i+1,j,ny)] + f->p[IDX(i-1,j,ny)])
                      - an*(f->p[IDX(i,j+1,ny)] + f->p[IDX(i,j-1,ny)]);
            res[IDX(i,j,ny)] = (-f->rhs[IDX(i,j,ny)]) - Ap;
        }
        rmt_multiple_recursive(g, corr, res, &st, 0);
        double omegaTry = omegaRMT;
        int accepted = 0;
        for (int ls=0; ls<8; ++ls) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] += omegaTry * corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g, ph, f);
            double newRes = residual_norm_pressure_variable(g, f, f->rhs);
            if (isfinite(newRes) && (oldRes <= 0.0 || newRes < oldRes)) { accepted = 1; break; }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j)
                f->p[IDX(i,j,ny)] -= omegaTry * corr[IDX(i,j,ny)];
            apply_bc_velocity_pressure(g, ph, f);
            omegaTry *= 0.5;
        }
        if (!accepted) {
            apply_bc_velocity_pressure(g, ph, f);
        }

        for (int sweep=0; sweep<sweepsPerCycle; ++sweep) {
            for (int color=0; color<2; ++color) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
                for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
                    if (((i + j) & 1) != color) continue;
                    double pGs = (ae*f->p[IDX(i+1,j,ny)] + aw*f->p[IDX(i-1,j,ny)]
                                + an*f->p[IDX(i,j+1,ny)] + as*f->p[IDX(i,j-1,ny)]
                                - f->rhs[IDX(i,j,ny)]) / MAX(ap, 1.0e-20);
                    f->p[IDX(i,j,ny)] = (1.0-relaxP)*f->p[IDX(i,j,ny)] + relaxP*pGs;
                }
                apply_bc_velocity_pressure(g, ph, f);
            }
        }
        if (residual_norm_pressure_variable(g, f, f->rhs) < 1.0e-6) break;
    }
    free(corr);
    free(res);

    g_pressure_rmt_wall_seconds +=
        rmt_wall_now() - __rmt_pressure_t0;

    g_pressure_rmt_calls += 1;
}

static void correct_velocity(const Grid *g, const Phys *ph, Fields *f, double dt) {
    const int nx=g->Nx, ny=g->Ny;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double dpdx;
        if (i == 1) dpdx = (f->p[IDX(i+1,j,ny)] - f->p[IDX(i,j,ny)]) / g->dx;
        else if (i == nx) dpdx = (f->p[IDX(i,j,ny)] - f->p[IDX(i-1,j,ny)]) / g->dx;
        else dpdx = (f->p[IDX(i+1,j,ny)] - f->p[IDX(i-1,j,ny)])/(2.0*g->dx);
        double dpdy;
        if (j == 1) dpdy = (f->p[IDX(i,j+1,ny)] - f->p[IDX(i,j,ny)]) / g->dy;
        else if (j == ny) dpdy = (f->p[IDX(i,j,ny)] - f->p[IDX(i,j-1,ny)]) / g->dy;
        else dpdy = (f->p[IDX(i,j+1,ny)] - f->p[IDX(i,j-1,ny)])/(2.0*g->dy);
        double rhoi = f->rho[IDX(i,j,ny)];
        f->u[IDX(i,j,ny)] = f->us[IDX(i,j,ny)] - dt*dpdx/MAX(rhoi,1.0e-12);
        f->v[IDX(i,j,ny)] = f->vs[IDX(i,j,ny)] - dt*dpdy/MAX(rhoi,1.0e-12);
    }
}

static void build_scalar_rhs_openfoam_like(const Grid *g, const Phys *ph, const Fields *f,
                                           const double *phi, double *rhs, double dt,
                                           double bottomSlot, double topSlot, double outletInletValue) {
    const int nx=g->Nx, ny=g->Ny; const double dx=g->dx, dy=g->dy;
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double mE=f->mfx[IDX(i,  j,ny)], mW=f->mfx[IDX(i-1,j,ny)];
        double mN=f->mfy[IDX(i,j,  ny)], mS=f->mfy[IDX(i,j-1,ny)];
        double phiE=scalar_face_x_species(g,f,phi,i,j,outletInletValue);
        double phiW=scalar_face_x_species(g,f,phi,i-1,j,outletInletValue);
        double phiN=scalar_face_y_species(g,ph,f,phi,i,j,bottomSlot,topSlot);
        double phiS=scalar_face_y_species(g,ph,f,phi,i,j-1,bottomSlot,topSlot);
        rhs[IDX(i,j,ny)] = f->rho[IDX(i,j,ny)]*phi[IDX(i,j,ny)]
                         - dt*((mE*phiE - mW*phiW)/dx + (mN*phiN - mS*phiS)/dy);
    }
}

static void solve_scalar_helmholtz(const Grid *g, const Phys *ph, const Fields *f, double *phi, const double *rhs,
                                   double dt, const double *diffField, int iters,
                                   double lower, double upper,
                                   double bottomSlot, double topSlot, double outletInletValue, double relaxA) {
    const int nx=g->Nx, ny=g->Ny;
    const double dx2 = g->dx*g->dx, dy2 = g->dy*g->dy;
    const int sweepsPerCycle = MAX(4, g_rmtPostSmooth + g_rmtCoarseRepeats + 1);
    (void)ph;

    /*
       Pressure keeps the safeguarded 2D 9-shift RMT correction.  For the reacting
       scalar/enthalpy blocks we return to the stable V91 variable-coefficient RBGS solve.
       In V93 the multi-shift correction used a very approximate coarse scalar operator
       (rhoBar,diffBar), which systematically damped YP/T and drove Tmax below the
       OpenFOAM reference even when the run stayed stable.
    */
    for (int cyc=0; cyc<iters; ++cyc) {
        apply_bc_scalar_of(g, f, phi, bottomSlot, topSlot, outletInletValue);
        for (int sweep=0; sweep<sweepsPerCycle; ++sweep) {
            for (int color=0; color<2; ++color) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
                for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
                    if (((i + j) & 1) != color) continue;
                    double rhoP = MAX(f->rho[IDX(i,j,ny)], 1.0e-12);
                    double gE = rho_face_x(g, f, i,   j) * face_avg_x_coeff(g, diffField, i,   j) / dx2;
                    double gW = rho_face_x(g, f, i-1, j) * face_avg_x_coeff(g, diffField, i-1, j) / dx2;
                    double gN = rho_face_y(g, f, i,   j) * face_avg_y_coeff(g, diffField, i,   j) / dy2;
                    double gS = rho_face_y(g, f, i, j-1) * face_avg_y_coeff(g, diffField, i, j-1) / dy2;
                    double ap = rhoP + dt*(gE + gW + gN + gS);
                    double phiGs = (rhs[IDX(i,j,ny)] + dt*(gE*phi[IDX(i+1,j,ny)] + gW*phi[IDX(i-1,j,ny)]
                                                        + gN*phi[IDX(i,j+1,ny)] + gS*phi[IDX(i,j-1,ny)])) / MAX(ap, 1.0e-20);
                    double val = (1.0-relaxA)*phi[IDX(i,j,ny)] + relaxA*phiGs;
                    phi[IDX(i,j,ny)] = clamp(val, lower, upper);
                }
                apply_bc_scalar_of(g, f, phi, bottomSlot, topSlot, outletInletValue);
            }
        }
        if (residual_norm_scalar_variable(g, f, phi, rhs, diffField, dt) < 1.0e-8) break;
    }
}

static double arrhenius_k(const Phys *ph, double T) {
    T = MAX(T, 200.0);
    double logk = log(ph->A) + ph->beta*log(T) - ph->Ta/T;
    if (logk > 700.0)
        logk = 700.0;

    if (logk < -700.0)
        return 0.0;

    return exp(logk);
}

static double implicit_extent_be(const Phys *ph, double dts, double T, double CF0, double CO0) {
    double k = arrhenius_k(ph, T);
    double xiMax = MIN(0.5*CF0, CO0);
    if (xiMax <= 0.0 || k <= 0.0) return 0.0;
    double xi = MIN(dts*k*CF0*CF0*CO0, 0.5*xiMax);
    xi = clamp(xi, 0.0, xiMax*(1.0-1e-12));
    for (int it=0; it<12; ++it) {
        double a = MAX(CF0 - 2.0*xi, 0.0), b = MAX(CO0 - xi, 0.0);
        double F = xi - dts*k*a*a*b;
        if (fabs(F) < 1e-14*MAX(1.0, xiMax)) break;
        double dF = 1.0 + dts*k*(4.0*a*b + a*a); /* derivative of -a^2 b gives +(4ab+a^2) */
        double xnew = xi - F/dF;
        if (!(xnew >= 0.0 && xnew <= xiMax) || !isfinite(xnew)) xnew = 0.5*(xi + clamp(xnew,0.0,xiMax));
        xnew = clamp(xnew, 0.0, xiMax);
        if (fabs(xnew - xi) < 1e-14*MAX(1.0, xiMax)) { xi = xnew; break; }
        xi = xnew;
    }
    return clamp(xi, 0.0, xiMax);
}

static void chemistry_cell_step_implicit(const Phys *ph, double dtOuter, double pCell, double *YF, double *YO, double *YP, double *h) {
    double yF=clamp(*YF,0.0,1.0), yO=clamp(*YO,0.0,0.21), yP=clamp(*YP,0.0,1.0);
    double hloc = MAX(*h, 0.0);
    double yN2 = MAX(0.0, 1.0 - yF - yO - yP);
    /* Use outer-step lagged temperature to avoid runaway self-acceleration within the same fluid step. */
    const double tempLag = MAX(T_from_h_mix(ph, yF, yO, yP, yN2, hloc, ph->Tref + hloc/MAX(ph->cp,1.0)), 200.0);
    if (yF <= 1e-16 || yO <= 1e-16 || !ph->chemistryOn) { *YF=yF; *YO=yO; *YP=yP; *h=hloc; return; }

    double rhoInit = mixture_density(ph, yF, yO, yP, yN2, tempLag, pCell);
    double CFinit = rhoInit * yF / ph->W_F;
    double COinit = rhoInit * yO / ph->W_O;
    double xiMaxInit = MIN(0.5*CFinit, COinit);
    if (xiMaxInit <= 0.0) { *YF=yF; *YO=yO; *YP=yP; *h=hloc; return; }

    /* Cumulative chemistry budget over the whole outer step: reactedFracMax applies to the entire step,
       not each accepted chemistry substep. */
    double xiBudget = ph->reactedFracMax * xiMaxInit;
    double xiCum = 0.0;
    double tloc=0.0;
    double dts = MIN(ph->initialChemicalTimeStep, dtOuter);
    double qcoeff = reaction_heat_release_per_extent(ph, tempLag);
    int guard=0;
    while (tloc < dtOuter - 1e-15 && ++guard < 400) {
        if (tloc + dts > dtOuter) dts = dtOuter - tloc;
        yN2 = MAX(0.0, 1.0 - yF - yO - yP);
        double rhoLoc = mixture_density(ph, yF, yO, yP, yN2, tempLag, pCell);
        double cpLoc = MAX(mixture_cp(ph, yF, yO, yP, yN2, tempLag), 1.0);
        double CF0 = rhoLoc * yF / ph->W_F;
        double CO0 = rhoLoc * yO / ph->W_O;
        double xiMax = MIN(0.5*CF0, CO0);
        if (xiMax <= 0.0) break;
        double xi = implicit_extent_be(ph, dts, tempLag, CF0, CO0);
        double xiRemain = xiBudget - xiCum;
        if (xiRemain <= 1e-20) break;
        if (xi > xiRemain) xi = xiRemain;
        double dh = qcoeff * xi / MAX(rhoLoc, 1.0e-12);
        double dT = dh / cpLoc;
        if (dT > ph->dTchemMax && dts > 1e-14) { dts *= 0.25; continue; }
        yF = clamp(yF - (2.0*ph->W_F*xi)/MAX(rhoLoc,1.0e-12), 0.0, 1.0);
        yO = clamp(yO - (1.0*ph->W_O*xi)/MAX(rhoLoc,1.0e-12), 0.0, 0.21);
        yP = clamp(yP + (2.0*ph->W_P*xi)/MAX(rhoLoc,1.0e-12), 0.0, 1.0);
        hloc = MAX(0.0, hloc + dh);
        xiCum += xi;
        /* Do not renormalize all active species together; preserve YF/YO shapes and only limit product to available remainder. */
        {
            double ypMax = MAX(0.0, 1.0 - yF - yO);
            if (ypMax < 0.0) ypMax = 0.0;
            yP = clamp(yP, 0.0, ypMax);
        }
        tloc += dts;
        dts = MIN(ph->maxChemicalTimeStep, dts*2.0);
    }
    *YF=yF; *YO=yO; *YP=yP; *h=hloc;
}

static void chemistry_subcycle(const Grid *g, const Phys *ph, Fields *f, double dt) {
    if (!ph->chemistryOn) return;
    for (int i=1;i<=g->Nx;++i) for (int j=1;j<=g->Ny;++j) {
        chemistry_cell_step_implicit(ph, dt, f->p[IDX(i,j,g->Ny)], &f->YF[IDX(i,j,g->Ny)], &f->YO[IDX(i,j,g->Ny)], &f->YP[IDX(i,j,g->Ny)], &f->h[IDX(i,j,g->Ny)]);
        double yF = clamp(f->YF[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yO = clamp(f->YO[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yP = clamp(f->YP[IDX(i,j,g->Ny)], 0.0, 1.0);
        double yN2 = MAX(0.0, 1.0 - yF - yO - yP);
        f->T[IDX(i,j,g->Ny)] = T_from_h_mix(ph, yF, yO, yP, yN2, f->h[IDX(i,j,g->Ny)], f->T[IDX(i,j,g->Ny)]);
    }
}

static void update_inert(const Grid *g, Fields *f) {
    for (int i=1;i<=g->Nx;++i) for (int j=1;j<=g->Ny;++j) {
        f->YF[IDX(i,j,g->Ny)] = clamp(f->YF[IDX(i,j,g->Ny)], 0.0, 1.0);
        f->YO[IDX(i,j,g->Ny)] = clamp(f->YO[IDX(i,j,g->Ny)], 0.0, 0.21);

        /* Preserve fuel/oxidizer fields and only trim product to the available active-species remainder. */
        {
            double ypMax = MAX(0.0, 1.0 - f->YF[IDX(i,j,g->Ny)] - f->YO[IDX(i,j,g->Ny)]);
            f->YP[IDX(i,j,g->Ny)] = clamp(f->YP[IDX(i,j,g->Ny)], 0.0, ypMax);
        }

        double sum = f->YF[IDX(i,j,g->Ny)] + f->YO[IDX(i,j,g->Ny)] + f->YP[IDX(i,j,g->Ny)];
        f->YN2[IDX(i,j,g->Ny)] = MAX(0.0, 1.0-sum);
    }
}

static double compute_max_co(const Grid *g, const Fields *f, double dt) {
    double m=0.0; for (int i=1;i<=g->Nx;++i) for (int j=1;j<=g->Ny;++j) {
        double co = fabs(f->u[IDX(i,j,g->Ny)])*dt/g->dx + fabs(f->v[IDX(i,j,g->Ny)])*dt/g->dy; if (co>m) m=co; }
    return m;
}
static double field_max(const Grid *g, const double *a) {
    double m=-1e300;
    for(int i=1;i<=g->Nx;++i) for(int j=1;j<=g->Ny;++j) {
        double v=a[IDX(i,j,g->Ny)];
        if (isfinite(v) && v>m) m=v;
    }
    return m;
}

static void make_time_tag(double time, char *buf, size_t n) { snprintf(buf,n,"t%.6f",time); for(size_t k=0; buf[k]; ++k) if(buf[k]=='.') buf[k]='p'; }


static double point_u_value(const Grid *g, const Phys *ph, const Fields *f, int ip, int jp) {
    const int nx=g->Nx, ny=g->Ny;
    double x = ip * g->dx;
    int onSlot = (x >= g->x1 - 1e-12 && x <= g->x2 + 1e-12);

    if (jp == 0 || jp == ny) {
        if (onSlot) return 0.0; /* inlet velocity is purely vertical on slot faces */
        return 0.0;             /* walls */
    }

    if (ip == 0)  return f->u[IDX(1,  jp, ny)];
    if (ip == nx) return f->u[IDX(nx, jp, ny)];

    double u00 = f->u[IDX(ip,   jp,   ny)];
    double u10 = f->u[IDX(ip+1, jp,   ny)];
    double u01 = f->u[IDX(ip,   jp+1, ny)];
    double u11 = f->u[IDX(ip+1, jp+1, ny)];
    return 0.25*(u00 + u10 + u01 + u11);
}

static double point_v_value(const Grid *g, const Phys *ph, const Fields *f, int ip, int jp) {
    const int nx=g->Nx, ny=g->Ny;
    double x = ip * g->dx;
    int onSlot = (x >= g->x1 - 1e-12 && x <= g->x2 + 1e-12);

    if (jp == 0)  return onSlot ? ph->vFuel : 0.0;
    if (jp == ny) return onSlot ? ph->vOx   : 0.0;

    if (ip == 0)  return f->v[IDX(1,  jp, ny)];
    if (ip == nx) return f->v[IDX(nx, jp, ny)];

    double v00 = f->v[IDX(ip,   jp,   ny)];
    double v10 = f->v[IDX(ip+1, jp,   ny)];
    double v01 = f->v[IDX(ip,   jp+1, ny)];
    double v11 = f->v[IDX(ip+1, jp+1, ny)];
    return 0.25*(v00 + v10 + v01 + v11);
}

static void write_vti(const Grid *g, const Phys *ph, const Fields *f, double time) {
    if (ensure_dir("vtk") != 0) die("mkdir vtk");
    char tag[64], fname[256]; make_time_tag(time,tag,sizeof(tag)); snprintf(fname,sizeof(fname),"vtk/field_%s.vti",tag);
    FILE *fp=fopen(fname,"w"); if(!fp) die("open vtk");
    fprintf(fp,"<?xml version=\"1.0\"?>\n");
    fprintf(fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp,"  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 -0.01 0\" Spacing=\"%.16e %.16e 1\">\n",g->Nx,g->Ny,g->dx,g->dy);
    fprintf(fp,"    <FieldData>\n");
    fprintf(fp,"      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">%.16e</DataArray>\n",time);
    fprintf(fp,"    </FieldData>\n");
    fprintf(fp,"    <Piece Extent=\"0 %d 0 %d 0 0\">\n",g->Nx,g->Ny);
    fprintf(fp,"      <CellData Scalars=\"T\">\n");
    #define WRITE_SCALAR(name,arr) do{ fprintf(fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",name); \
        for(int j=1;j<=g->Ny;++j){ for(int i=1;i<=g->Nx;++i) fprintf(fp,"%.16e ",(arr)[IDX(i,j,g->Ny)]); fprintf(fp,"\n"); } \
        fprintf(fp,"        </DataArray>\n"); }while(0)
    WRITE_SCALAR("p",f->p); WRITE_SCALAR("T",f->T); WRITE_SCALAR("h",f->h); WRITE_SCALAR("rho",f->rho); WRITE_SCALAR("mu",f->muMix); WRITE_SCALAR("cpMix",f->cpMix); WRITE_SCALAR("YF",f->YF); WRITE_SCALAR("YO",f->YO); WRITE_SCALAR("YP",f->YP); WRITE_SCALAR("YN2",f->YN2); WRITE_SCALAR("u",f->u); WRITE_SCALAR("v",f->v);
    #undef WRITE_SCALAR
    fprintf(fp,"      </CellData>\n");
    fprintf(fp,"      <PointData Vectors=\"U_point\">\n");
    fprintf(fp,"        <DataArray type=\"Float64\" Name=\"u_point\" format=\"ascii\">\n");
    for(int jp=0;jp<=g->Ny;++jp){ for(int ip=0;ip<=g->Nx;++ip) fprintf(fp,"%.16e ", point_u_value(g, ph, f, ip, jp)); fprintf(fp,"\n"); }
    fprintf(fp,"        </DataArray>\n");
    fprintf(fp,"        <DataArray type=\"Float64\" Name=\"v_point\" format=\"ascii\">\n");
    for(int jp=0;jp<=g->Ny;++jp){ for(int ip=0;ip<=g->Nx;++ip) fprintf(fp,"%.16e ", point_v_value(g, ph, f, ip, jp)); fprintf(fp,"\n"); }
    fprintf(fp,"        </DataArray>\n");
    fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Umag_point\" format=\"ascii\">\n");
    for(int jp=0;jp<=g->Ny;++jp){
        for(int ip=0;ip<=g->Nx;++ip){
            double up = point_u_value(g, ph, f, ip, jp);
            double vp = point_v_value(g, ph, f, ip, jp);
            fprintf(fp,"%.16e ", sqrt(up*up + vp*vp));
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"        </DataArray>\n");
    fprintf(fp,"      </PointData>\n");
    fprintf(fp,"    </Piece>\n  </ImageData>\n</VTKFile>\n");
    fclose(fp);
}

static void write_centerline_csv(const Grid *g, const Fields *f, double time) {
    char tag[64], fname[256]; make_time_tag(time,tag,sizeof(tag)); snprintf(fname,sizeof(fname),"centerline_%s.csv",tag);
    FILE *fp=fopen(fname,"w"); if(!fp) return; fprintf(fp,"y,T,YF,YO,YP,u,v,p\n"); int i=g->Nx/2;
    for(int j=1;j<=g->Ny;++j){ double y=-0.01 + (j-0.5)*g->dy; fprintf(fp,"%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n", y,f->T[IDX(i,j,g->Ny)],f->YF[IDX(i,j,g->Ny)],f->YO[IDX(i,j,g->Ny)],f->YP[IDX(i,j,g->Ny)],f->u[IDX(i,j,g->Ny)],f->v[IDX(i,j,g->Ny)],f->p[IDX(i,j,g->Ny)]); }
    fclose(fp);
}

int main(int argc, char **argv) {
    struct timespec wall_t0, wall_t1;
    timespec_get(&wall_t0, TIME_UTC);
    Grid g; Phys ph; Controls c; Fields f;
    init_case(&g,&ph,&c); parse_args(argc,argv,&g,&ph,&c);
    g_rmtMaxLevels = c.rmtMaxLevels;
    g_rmtPostSmooth = c.rmtPostSmooth;
    g_rmtCoarseSweeps = c.rmtCoarseSweeps;
    g_rmtCoarseRepeats = c.rmtCoarseRepeats;
    alloc_fields(&g,&f); init_fields(&g,&ph,&f); refresh_temperature_from_enthalpy(&g, &ph, &f); update_thermo_transport_fields(&g, &ph, &f);
    printf("RMT-3h/9-shift pressure solver: OpenMP 2D multiple-coarse-grid RMT with stable scalar/enthalpy blocks\n");
    printf("Grid Nx=%d Ny=%d dx=%.6e dy=%.6e\n",g.Nx,g.Ny,g.dx,g.dy);
    printf("Time end=%g dt0=%g write=%g maxCo=%g\n",c.endTime,c.dt0,c.writeInterval,c.maxCo);
    printf("Chemistry initialChemicalTimeStep=%g maxChemicalTimeStep=%g dTchemMax=%g reactedFracMax=%g p0=%g\n",ph.initialChemicalTimeStep,ph.maxChemicalTimeStep,ph.dTchemMax,ph.reactedFracMax,ph.p0);
    printf("Thermo variableCp=%s sutherland=%s variableDensity=%s hydroPInRho=%s rhoRelax=%g referenceT=%g\n",
           ph.variableCpOn ? "on" : "off",
           ph.sutherlandOn ? "on" : "off",
           ph.variableRhoOn ? "on" : "off",
           ph.hydroPressureInDensityOn ? "on" : "off",
           ph.rhoRelax, ph.Tref);
    if (ph.variableRhoOn && !ph.hydroPressureInDensityOn)
        printf("Using low-Mach thermodynamic density rho=p0/(Rmix*T) with mass-flux projection.\n");
    if (ph.variableRhoOn && ph.hydroPressureInDensityOn)
        printf("WARNING: hydroPInRho uses local hydrodynamic pressure in density and may be less stable.\n");
    printf("RMT cycles pressure=%d scalar=%d levels=%d post=%d coarseSweeps=%d coarseRepeats=%d omega=%.3g progressEvery=%d\n",c.poissonIters,c.scalarIters,c.rmtMaxLevels,c.rmtPostSmooth,c.rmtCoarseSweeps,c.rmtCoarseRepeats,g_rmtOmega,c.progressEvery);

    double time=c.startTime, dt=c.dt0, nextWrite=c.writeInterval; int step=0; if (g_fieldOutput) {
                write_vti(&g,&ph,&f,time);
                write_centerline_csv(&g,&f,time);
            }
    while (time < c.endTime - 1e-15) {
        ++step; if (time + dt > c.endTime) dt = c.endTime - time;
        apply_bc_velocity_pressure(&g,&ph,&f);
        compute_face_flux_rc(&g,&ph,&f,dt,0);
        apply_bc_scalar_of(&g,&f,f.YF,ph.YF_fuel,ph.YF_ox,0.0);
        apply_bc_scalar_of(&g,&f,f.YO,ph.YO_fuel,ph.YO_ox,0.0);
        apply_bc_scalar_of(&g,&f,f.YP,ph.YP_fuel,ph.YP_ox,0.0);
        apply_bc_h_exact(&g, f.v, f.h);
        refresh_temperature_from_enthalpy(&g, &ph, &f);
        update_thermo_transport_fields(&g, &ph, &f);
        memcpy(f.rhoOld, f.rho, (size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));

        for (int oc=0; oc<c.nOuterCorrectors; ++oc) {
            momentum_predictor(&g,&ph,&f,dt,c.relaxU);
            apply_bc_velocity_pressure(&g,&ph,&f);
            compute_face_flux_rc(&g,&ph,&f,dt,1);
            build_pressure_rhs(&g,&ph,&f,dt);
            for (int pc=0; pc<c.nCorrectors; ++pc) {
                solve_pressure_poisson(&g,&ph,&f,c.poissonIters,c.relaxP);
                apply_bc_velocity_pressure(&g,&ph,&f);
            }
            correct_velocity(&g,&ph,&f,dt);
            apply_bc_velocity_pressure(&g,&ph,&f);
            update_thermo_transport_fields(&g, &ph, &f);
            compute_face_flux_rc(&g,&ph,&f,dt,0);

            memcpy(f.phiOld,f.YF,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            build_scalar_rhs_openfoam_like(&g,&ph,&f,f.phiOld,f.rhsScalar,dt,ph.YF_fuel,ph.YF_ox,0.0);
            memcpy(f.YF,f.rhsScalar,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            solve_scalar_helmholtz(&g,&ph,&f,f.YF,f.rhsScalar,dt,f.DMix,c.scalarIters,0.0,1.0,ph.YF_fuel,ph.YF_ox,0.0,c.relaxSc);

            memcpy(f.phiOld,f.YO,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            build_scalar_rhs_openfoam_like(&g,&ph,&f,f.phiOld,f.rhsScalar,dt,ph.YO_fuel,ph.YO_ox,0.0);
            memcpy(f.YO,f.rhsScalar,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            solve_scalar_helmholtz(&g,&ph,&f,f.YO,f.rhsScalar,dt,f.DMix,c.scalarIters,0.0,0.21,ph.YO_fuel,ph.YO_ox,0.0,c.relaxSc);

            memcpy(f.phiOld,f.YP,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            build_scalar_rhs_openfoam_like(&g,&ph,&f,f.phiOld,f.rhsScalar,dt,ph.YP_fuel,ph.YP_ox,0.0);
            memcpy(f.YP,f.rhsScalar,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            solve_scalar_helmholtz(&g,&ph,&f,f.YP,f.rhsScalar,dt,f.DMix,c.scalarIters,0.0,1.0,ph.YP_fuel,ph.YP_ox,0.0,c.relaxSc);

            chemistry_subcycle(&g,&ph,&f,dt); update_inert(&g,&f); refresh_temperature_from_enthalpy(&g, &ph, &f); update_thermo_transport_fields(&g, &ph, &f);

            memcpy(f.phiOld,f.h,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            build_scalar_rhs_openfoam_like(&g,&ph,&f,f.phiOld,f.rhsScalar,dt,0.0,0.0,0.0);
            memcpy(f.h,f.rhsScalar,(size_t)(g.Nx+2)*(size_t)(g.Ny+2)*sizeof(double));
            solve_scalar_helmholtz(&g,&ph,&f,f.h,f.rhsScalar,dt,f.alphaMix,c.scalarIters,0.0,1.0e12,0.0,0.0,0.0,c.relaxH);
            refresh_temperature_from_enthalpy(&g, &ph, &f);
        }

        time += dt; double Co = compute_max_co(&g,&f,dt); double dtCo=(Co>1e-12)?(0.95*c.maxCo*dt/Co):(1.2*dt); double dtGrow=1.2*dt, dtMax=1e-4; dt=MIN(dtMax,MIN(dtCo,dtGrow));
        if (time >= nextWrite - 1e-14) {
            double Tmax=field_max(&g,f.T), YPmax=field_max(&g,f.YP), YOmax=field_max(&g,f.YO), YFmax=field_max(&g,f.YF);
            printf("step=%6d t=%9.6f dt=%9.3e Co=%8.3e Tmax=%10.3f YPmax=%10.6f YOmax=%10.6f YFmax=%10.6f\n",step,time,dt,Co,Tmax,YPmax,YOmax,YFmax);
            write_vti(&g,&ph,&f,time); write_centerline_csv(&g,&f,time); nextWrite += c.writeInterval;
        }
    }
    timespec_get(&wall_t1, TIME_UTC);
    double elapsed = (double)(wall_t1.tv_sec - wall_t0.tv_sec) + 1.0e-9*(double)(wall_t1.tv_nsec - wall_t0.tv_nsec);
    printf("Done. Final time=%g\n",time);
    printf("Elapsed computational time: %.6f s\n", elapsed);

    double finalTmax  = field_max(&g, f.T);
    double finalYPmax = field_max(&g, f.YP);
    double finalYOmax = field_max(&g, f.YO);
    double finalYFmax = field_max(&g, f.YF);

    printf("Pressure RMT wall time: %.9f s\n",
           g_pressure_rmt_wall_seconds);

    printf("Pressure RMT calls: %lld\n",
           g_pressure_rmt_calls);

#ifdef _OPENMP
    printf("OpenMP max threads: %d\n",
           omp_get_max_threads());
#else
    printf("OpenMP max threads: 1\n");
#endif

    printf(
        "RMT_RUN_SUMMARY "
        "final_time=%.17g "
        "total_wall_s=%.17g "
        "pressure_wall_s=%.17g "
        "pressure_calls=%lld "
        "Tmax=%.17g "
        "YPmax=%.17g "
        "YOmax=%.17g "
        "YFmax=%.17g\n",
        time,
        elapsed,
        g_pressure_rmt_wall_seconds,
        g_pressure_rmt_calls,
        finalTmax,
        finalYPmax,
        finalYOmax,
        finalYFmax
    );

    free_fields(&f); return 0;
}
