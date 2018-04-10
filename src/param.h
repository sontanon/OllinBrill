#include <libconfig.h>

#ifdef MAIN_FILE
/* CONFIG FILE */
config_t cfg;

/* GRID */
double dr = 0.01;
double dz = 0.01;
int NrInterior = 800;
int NzInterior = 800;
int NrTotal = 0;
int NzTotal = 0;
int ghost = 0;
int DIM = 0;

/* TIME STEPPING */
double t = 0.0;
double dtfac_old = 0.2;
double dtfac_new = 0.2;
double FinalTime = 1.0;

/* OUTPUT */
int constraintCheck = 1;
int Noutput2D = 1;
const char *dirname = "Test";

/* BOUNDARIES */
const char *boundtype = "radiative";

/* SLICING */
const char *slicing = "maximal"; /* maximal, 1+log, one */
const char *ilapse = "one"; /* one, isotropic */
double gauge_f = 1.0;
double gl_alpha0 = 0.01;
double gl_rr0 = 4.0;

/* SHIFT */
const char *shift = "none";
const char *ishift = "zero";
const char *bssn_flavor = "lagrangian";
double sigma = 1.0;
double driver_c = 0.0;
double driver_eta = 0.0;
double gs_beta0 = 0.01;
double gs_rr0 = 4.0;

/* INITIAL DATA */
const char *idata = "BrillWave";

/* HORIZON FINDER */
int ahfind_flag = 0;
int ahfind_freq = 1;
int ahfind_bicubic = 0;
int ahfind_mask_filter = 0;
double ahfind_start_time = 0.0;
double ahfind_maxr = 2.5;
double ahfind_minr = 0.25;
int ahfind_bis_iter = 50;

/* EVOLUTION */
const char *spacetime = "dynamic";
const char *integrator = "rk4";
const char *order = "four";
int interior_norm = 1;
int adaptive_boundary_norm = 0;
double eta = 2.0;
double diss = 0.0;
double damp = 0.0;

/* ELLIPTIC SOLVER */
const char *ell_freq = "every";
int use_permutation = 0;
int use_preconditioner = 0;
int use_low_rank = 1;
int perm_flag = 0;
int perm_count = 0;
int perm_count_max = 0;
int precond_flag = 0;
int precond_L = 0;
int low_rank_flag = 0;
int nrobin = 1;

/* BRILL WAVES */
double brill_a0 = 6.0;
double brill_r0 = 0.0;
double brill_z0 = 0.0;
double brill_sr0 = 1.0;
double brill_sz0 = 1.0;

/* SCHWARSCHILD INITIAL DATA */
double schwar_mass = 1.0;

/* CFL ANALYSIS */
int CFL_n_old = 4;
int CFL_n_new = 4;
int CFL_count = 0;
int CFL_stop = 16;
#else
/* CONFIG FILE */
extern config_t cfg;

/* GRID */
extern double dr;
extern double dz;
extern int NrInterior;
extern int NzInterior;
extern int NrTotal;
extern int NzTotal;
extern int ghost;
extern int DIM;

/* TIME STEPPING */
extern double t;
extern double dtfac_old;
extern double dtfac_new;
extern double FinalTime;

/* OUTPUT */
extern int constraintCheck;
extern int Noutput2D;
extern const char *dirname;

/* BOUNDARIES */
extern const char *boundtype;

/* SLICING */
extern const char *slicing;
extern const char *ilapse;
extern double gauge_f;
extern double gl_alpha0;
extern double gl_rr0;

/* SHIFT */
extern const char *shift;
extern const char *ishift;
extern const char *bssn_flavor;
extern double sigma;
extern double driver_c;
extern double driver_eta;
extern double gs_beta0;
extern double gs_rr0;

/* INITIAL DATA */
extern const char *idata;

/* HORIZON FINDER */
extern int ahfind_flag;
extern int ahfind_freq;
extern int ahfind_bicubic;
extern int ahfind_mask_filter;
extern double ahfind_start_time;
extern double ahfind_maxr;
extern double ahfind_minr;
extern int ahfind_bis_iter;

/* EVOLUTION */
extern const char *spacetime;
extern const char *integrator;
extern const char *order;
extern int interior_norm;
extern int adaptive_boundary_norm;
extern double eta;
extern double diss;
extern double damp;

/* ELLIPTIC SOLVER */
extern const char *ell_freq;
extern int use_permutation;
extern int use_preconditioner;
extern int use_low_rank;
extern int perm_flag;
extern int perm_count;
extern int perm_count_max;
extern int precond_flag;
extern int precond_L;
extern int low_rank_flag;
extern int nrobin;

/* BRILL WAVES */
extern double brill_a0;
extern double brill_r0;
extern double brill_z0;
extern double brill_sr0;
extern double brill_sz0;

/* SCHWARSCHILD */
extern double schwar_mass;

/* CFL ANALYSIS */
extern int CFL_n_old;
extern int CFL_n_new;
extern int CFL_count;
extern int CFL_stop;
#endif
