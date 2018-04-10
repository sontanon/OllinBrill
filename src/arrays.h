#ifdef MAIN_FILE
// AUXILIARY ARRAYS.
double *auxarray;
double *auxarray2;
double *auxarray3;
double *auxarray4;

// COORDINATE GRIDS.
double *r;
double *z;
double *rr;
double *r2;

// LASPE.
double *alpha;
double *Dr_alpha;
double *Dz_alpha;
double *Drr_alpha;
double *Dzz_alpha;
double *Drz_alpha;

// BONA-MASSO SOURCE.
double *falpha;

// CONFORMAL METRIC.
double *a;
double *b;
double *h;
double *c;
double *lambda;

// PHYSICAL METRIC
double *phys_a;
double *phys_b;
double *phys_h;
double *phys_c;
double *phys_lambda;

// PHYSICAL SPHERICAL METRIC
double *g_rr;
double *g_rth;
double *g_thth;
double *g_phph;
double *ig_rr;
double *ig_rth;
double *ig_thth;
double *ig_phph;

// SPEED OF LIGHT.
double *c_rr;
double *c_th;
double *c_ph;

// COURANT CONDITION.
double *courant_rr;
double *courant_th;
double *courant_ph;

double *Dr_a;
double *Dr_b;
double *Dr_h;
double *Dr_c;
double *Dr_lambda;

double *Dz_a;
double *Dz_b;
double *Dz_h;
double *Dz_c;
double *Dz_lambda;

double *Drr_a;
double *Drr_b;
double *Drr_h;
double *Drr_c;
double *Drr_lambda;

double *Dzz_a;
double *Dzz_b;
double *Dzz_h;
double *Dzz_c;
double *Dzz_lambda;

double *Drz_a;
double *Drz_b;
double *Drz_h;
double *Drz_c;
double *Drz_lambda;

// INVERSE CONORMAL METRIC.
double *g_a;
double *g_b;
double *g_h;
double *g_c;
double *g_lambda;

double *Dr_g_a;
double *Dr_g_b;
double *Dr_g_h;
double *Dr_g_c;

double *Dz_g_a;
double *Dz_g_b;
double *Dz_g_h;
double *Dz_g_c;

double *hdet;
double *Dr_hdet;
double *Dz_hdet;

// TRACELESS EXTRINSIC CURVATURE.
double *A_a;
double *A_b;
double *A_h;
double *A_c;
double *A_lambda;

double *Dr_A_a;
double *Dr_A_b;
double *Dr_A_h;
double *Dr_A_c;
double *Dr_A_lambda;

double *Dz_A_a;
double *Dz_A_b;
double *Dz_A_h;
double *Dz_A_c;
double *Dz_A_lambda;

// CONFORMAL FACTOR.
double *phi;
double *psi;
double *psi4;
double *Dr_phi;
double *Dz_phi;
double *Drr_phi;
double *Dzz_phi;
double *Drz_phi;

// TRACE OF EXTRINSIC CURVATURE.
double *K;
double *Dr_K;
double *Dz_K;

// DELTAS.
double *Deltar;
double *Dr_Deltar;
double *Dz_Deltar;

double *Deltaz;
double *Dr_Deltaz;
double *Dz_Deltaz;

// Regularization.
double *DDeltar_r;
double *DDalpha_r;
double *DDphi_r;

// RICCI TENSOR.
double *R_a;
double *R_b;
double *R_h;
double *R_c;
double *R_lambda;
double *RSCAL;

// SECOND COVARIANT DERIVATIVE OF THE LAPSE
double *D2alpha_a;
double *D2alpha_b;
double *D2alpha_h;
double *D2alpha_c;
double *D2alpha_lambda;
double *LaplaAlpha;

// QUADRATIC QUANTITIES OF A.
double *A2_a;
double *A2_b;
double *A2_h;
double *A2_c;
double *A2_lambda;
double *A2;

// CONTRAVARIANT A.
double *upA_a;
double *upA_b;
double *upA_h;
double *upA_c;

// DIVERGENCE OF A.
double *divA_r;
double *divA_z;

// TERM IN DELTA SOURCES.
double *ADelta_r;
double *ADelta_z;

// CONSTRAINTS.
double *ham;
double *mom_r;
double *mom_z;
double *CDeltar;
double *CDeltaz;
double *Clambda;
double *CA_lambda;

// OLD ARRAYS.
double *alpha_p;
double *phi_p;
double *a_p;
double *b_p;
double *h_p;
double *c_p;
double *lambda_p;
double *K_p;
double *A_a_p;
double *A_b_p;
double *A_h_p;
double *A_c_p;
double *A_lambda_p;
double *Deltar_p;
double *Deltaz_p;

// ACCUMULATION ARRAYS.
double *alpha_a;
double *phi_a;
double *a_a;
double *b_a;
double *h_a;
double *c_a;
double *lambda_a;
double *K_a;
double *A_a_a;
double *A_b_a;
double *A_h_a;
double *A_c_a;
double *A_lambda_a;
double *Deltar_a;
double *Deltaz_a;

// SOURCES.
double *salpha;
double *sphi;
double *sa;
double *sb;
double *sh;
double *sc;
double *slambda;
double *sK;
double *sA_a;
double *sA_b;
double *sA_h;
double *sA_c;
double *sA_lambda;
double *sDeltar;
double *sDeltaz;

// ELLIPTIC COEFFICIENTS.
double *ell_a;
double *ell_b;
double *ell_c;
double *ell_d;
double *ell_e;
double *ell_s;
double *ell_f;

// TIME CONSTRAINTS.
double *t_ham_norm;
double *t_mom_r_norm;
double *t_mom_z_norm;
double *t_CDeltar_norm;
double *t_CDeltaz_norm;
double *t_Clambda_norm;
double *t_CA_lambda_norm;

// MIXED DERIVATIVES FOR BICUBIC INTERPOLATION.
double *Drz_A_a;
double *Drz_A_b;
double *Drz_A_h;
double *Drz_A_c;
double *Drz_K;
double *Drrz_a;
double *Drrz_b;
double *Drrz_h;
double *Drrz_c;
double *Drrz_phi;
double *Drzz_a;
double *Drzz_b;
double *Drzz_h;
double *Drzz_c;
double *Drzz_phi;

// SHIFT.
double *beta_r;
double *beta_z;
double *dtbeta_r;
double *dtbeta_z;
double *Dr_beta_r;
double *Dz_beta_r;
double *Drr_beta_r;
double *Dzz_beta_r;
double *Drz_beta_r;
double *Dr_beta_z;
double *Dz_beta_z;
double *Drr_beta_z;
double *Dzz_beta_z;
double *Drz_beta_z;
double *DD_beta_rr;
double *Dr_dtbeta_r;
double *Dz_dtbeta_r;
double *Dr_dtbeta_z;
double *Dz_dtbeta_z;
double *div_beta;
double *Dr_div_beta;
double *Dz_div_beta;
double *sbeta_r;
double *sbeta_z;
double *sdtbeta_r;
double *sdtbeta_z;
double *beta_r_p;
double *beta_z_p;
double *dtbeta_r_p;
double *dtbeta_z_p;
double *beta_r_a;
double *beta_z_a;
double *dtbeta_r_a;
double *dtbeta_z_a;

// ADVECTIVE DERIVATIVES.
double *DAr_alpha;
double *DAz_alpha;
double *DAr_phi;
double *DAz_phi;
double *DAr_a;
double *DAz_a;
double *DAr_b;
double *DAz_b;
double *DAr_h;
double *DAz_h;
double *DAr_c;
double *DAz_c;
double *DAr_lambda;
double *DAz_lambda;
double *DAr_A_lambda;
double *DAz_A_lambda;
double *DAr_beta_r;
double *DAz_beta_r;
double *DAr_beta_z;
double *DAz_beta_z;
double *DAr_dtbeta_r;
double *DAz_dtbeta_r;
double *DAr_dtbeta_z;
double *DAz_dtbeta_z;
double *DAr_K;
double *DAz_K;
double *DAr_A_a;
double *DAz_A_a;
double *DAr_A_b;
double *DAz_A_b;
double *DAr_A_h;
double *DAz_A_h;
double *DAr_A_c;
double *DAz_A_c;
double *DAr_Deltar;
double *DAz_Deltar;
double *DAr_Deltaz;
double *DAz_Deltaz;
#else
// AUXILIARY ARRAYS.
extern double *auxarray;
extern double *auxarray2;
extern double *auxarray3;
extern double *auxarray4;

// COORDINATE GRIDS.
extern double *r;
extern double *z;
extern double *rr;
extern double *r2;

// LASPE.
extern double *alpha;
extern double *Dr_alpha;
extern double *Dz_alpha;
extern double *Drr_alpha;
extern double *Dzz_alpha;
extern double *Drz_alpha;

// BONA-MASSO SOURCE.
extern double *falpha;

// CONFORMAL METRIC.
extern double *a;
extern double *b;
extern double *h;
extern double *c;
extern double *lambda;

// PHYSICAL METRIC
extern double *phys_a;
extern double *phys_b;
extern double *phys_h;
extern double *phys_c;
extern double *phys_lambda;

// PHYSICAL SPHERICAL METRIC
extern double *g_rr;
extern double *g_rth;
extern double *g_thth;
extern double *g_phph;
extern double *ig_rr;
extern double *ig_rth;
extern double *ig_thth;
extern double *ig_phph;

// SPEED OF LIGHT.
extern double *c_rr;
extern double *c_th;
extern double *c_ph;

// COURANT CONDITION.
extern double *courant_rr;
extern double *courant_th;
extern double *courant_ph;

extern double *Dr_a;
extern double *Dr_b;
extern double *Dr_h;
extern double *Dr_c;
extern double *Dr_lambda;

extern double *Dz_a;
extern double *Dz_b;
extern double *Dz_h;
extern double *Dz_c;
extern double *Dz_lambda;

extern double *Drr_a;
extern double *Drr_b;
extern double *Drr_h;
extern double *Drr_c;
extern double *Drr_lambda;

extern double *Dzz_a;
extern double *Dzz_b;
extern double *Dzz_h;
extern double *Dzz_c;
extern double *Dzz_lambda;

extern double *Drz_a;
extern double *Drz_b;
extern double *Drz_h;
extern double *Drz_c;
extern double *Drz_lambda;

// INVERSE CONORMAL METRIC.
extern double *g_a;
extern double *g_b;
extern double *g_h;
extern double *g_c;
extern double *g_lambda;

extern double *Dr_g_a;
extern double *Dr_g_b;
extern double *Dr_g_h;
extern double *Dr_g_c;

extern double *Dz_g_a;
extern double *Dz_g_b;
extern double *Dz_g_h;
extern double *Dz_g_c;

extern double *hdet;
extern double *Dr_hdet;
extern double *Dz_hdet;

// TRACELESS EXTRINSIC CURVATURE.
extern double *A_a;
extern double *A_b;
extern double *A_h;
extern double *A_c;
extern double *A_lambda;

extern double *Dr_A_a;
extern double *Dr_A_b;
extern double *Dr_A_h;
extern double *Dr_A_c;
extern double *Dr_A_lambda;

extern double *Dz_A_a;
extern double *Dz_A_b;
extern double *Dz_A_h;
extern double *Dz_A_c;
extern double *Dz_A_lambda;

// CONFORMAL FACTOR.
extern double *phi;
extern double *psi;
extern double *psi4;
extern double *Dr_phi;
extern double *Dz_phi;
extern double *Drr_phi;
extern double *Dzz_phi;
extern double *Drz_phi;

// TRACE OF EXTRINSIC CURVATURE.
extern double *K;
extern double *Dr_K;
extern double *Dz_K;

// DELTAS.
extern double *Deltar;
extern double *Dr_Deltar;
extern double *Dz_Deltar;

extern double *Deltaz;
extern double *Dr_Deltaz;
extern double *Dz_Deltaz;

// Regularization.
extern double *DDeltar_r;
extern double *DDalpha_r;
extern double *DDphi_r;

// RICCI TENSOR.
extern double *R_a;
extern double *R_b;
extern double *R_h;
extern double *R_c;
extern double *R_lambda;
extern double *RSCAL;

// SECOND COVARIANT DERIVATIVE OF THE LAPSE
extern double *D2alpha_a;
extern double *D2alpha_b;
extern double *D2alpha_h;
extern double *D2alpha_c;
extern double *D2alpha_lambda;
extern double *LaplaAlpha;

// QUADRATIC QUANTITIES OF A.
extern double *A2_a;
extern double *A2_b;
extern double *A2_h;
extern double *A2_c;
extern double *A2_lambda;
extern double *A2;

// CONTRAVARIANT A.
extern double *upA_a;
extern double *upA_b;
extern double *upA_h;
extern double *upA_c;

// DIVERGENCE OF A.
extern double *divA_r;
extern double *divA_z;

// TERM IN DELTA SOURCES.
extern double *ADelta_r;
extern double *ADelta_z;

// CONSTRAINTS.
extern double *ham;
extern double *mom_r;
extern double *mom_z;
extern double *CDeltar;
extern double *CDeltaz;
extern double *Clambda;
extern double *CA_lambda;

// OLD ARRAYS.
extern double *alpha_p;
extern double *phi_p;
extern double *a_p;
extern double *b_p;
extern double *h_p;
extern double *c_p;
extern double *lambda_p;
extern double *K_p;
extern double *A_a_p;
extern double *A_b_p;
extern double *A_h_p;
extern double *A_c_p;
extern double *A_lambda_p;
extern double *Deltar_p;
extern double *Deltaz_p;

// ACCUMULATION ARRAYS.
extern double *alpha_a;
extern double *phi_a;
extern double *a_a;
extern double *b_a;
extern double *h_a;
extern double *c_a;
extern double *lambda_a;
extern double *K_a;
extern double *A_a_a;
extern double *A_b_a;
extern double *A_h_a;
extern double *A_c_a;
extern double *A_lambda_a;
extern double *Deltar_a;
extern double *Deltaz_a;

// SOURCES.
extern double *salpha;
extern double *sphi;
extern double *sa;
extern double *sb;
extern double *sh;
extern double *sc;
extern double *slambda;
extern double *sK;
extern double *sA_a;
extern double *sA_b;
extern double *sA_h;
extern double *sA_c;
extern double *sA_lambda;
extern double *sDeltar;
extern double *sDeltaz;

// ELLIPTIC COEFFICIENTS.
extern double *ell_a;
extern double *ell_b;
extern double *ell_c;
extern double *ell_d;
extern double *ell_e;
extern double *ell_s;
extern double *ell_f;

// TIME CONSTRAINTS.
extern double *t_ham_norm;
extern double *t_mom_r_norm;
extern double *t_mom_z_norm;
extern double *t_CDeltar_norm;
extern double *t_CDeltaz_norm;
extern double *t_Clambda_norm;
extern double *t_CA_lambda_norm;

// MIXED DERIVATIVES FOR BICUBIC INTERPOLATION.
extern double *Drz_A_a;
extern double *Drz_A_b;
extern double *Drz_A_h;
extern double *Drz_A_c;
extern double *Drz_K;
extern double *Drrz_a;
extern double *Drrz_b;
extern double *Drrz_h;
extern double *Drrz_c;
extern double *Drrz_phi;
extern double *Drzz_a;
extern double *Drzz_b;
extern double *Drzz_h;
extern double *Drzz_c;
extern double *Drzz_phi;

// SHIFT.
extern double *beta_r;
extern double *beta_z;
extern double *dtbeta_r;
extern double *dtbeta_z;
extern double *Dr_beta_r;
extern double *Dz_beta_r;
extern double *Drr_beta_r;
extern double *Dzz_beta_r;
extern double *Drz_beta_r;
extern double *Dr_beta_z;
extern double *Dz_beta_z;
extern double *Drr_beta_z;
extern double *Dzz_beta_z;
extern double *Drz_beta_z;
extern double *DD_beta_rr;
extern double *Dr_dtbeta_r;
extern double *Dz_dtbeta_r;
extern double *Dr_dtbeta_z;
extern double *Dz_dtbeta_z;
extern double *div_beta;
extern double *Dr_div_beta;
extern double *Dz_div_beta;
extern double *sbeta_r;
extern double *sbeta_z;
extern double *sdtbeta_r;
extern double *sdtbeta_z;
extern double *beta_r_p;
extern double *beta_z_p;
extern double *dtbeta_r_p;
extern double *dtbeta_z_p;
extern double *beta_r_a;
extern double *beta_z_a;
extern double *dtbeta_r_a;
extern double *dtbeta_z_a;

// ADVECTIVE DERIVATIVES.
extern double *DAr_alpha;
extern double *DAz_alpha;
extern double *DAr_phi;
extern double *DAz_phi;
extern double *DAr_a;
extern double *DAz_a;
extern double *DAr_b;
extern double *DAz_b;
extern double *DAr_h;
extern double *DAz_h;
extern double *DAr_c;
extern double *DAz_c;
extern double *DAr_lambda;
extern double *DAz_lambda;
extern double *DAr_A_lambda;
extern double *DAz_A_lambda;
extern double *DAr_beta_r;
extern double *DAz_beta_r;
extern double *DAr_beta_z;
extern double *DAz_beta_z;
extern double *DAr_dtbeta_r;
extern double *DAz_dtbeta_r;
extern double *DAr_dtbeta_z;
extern double *DAz_dtbeta_z;
extern double *DAr_K;
extern double *DAz_K;
extern double *DAr_A_a;
extern double *DAz_A_a;
extern double *DAr_A_b;
extern double *DAz_A_b;
extern double *DAr_A_h;
extern double *DAz_A_h;
extern double *DAr_A_c;
extern double *DAz_A_c;
extern double *DAr_Deltar;
extern double *DAz_Deltar;
extern double *DAr_Deltaz;
extern double *DAz_Deltaz;
#endif
