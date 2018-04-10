#include "tools.h"

#include "param.h"
#include "write_flags.h"

void parser(const char *fname)
{
	// Initialize config.
	config_init(&cfg);

	// Read the file. If there is an error, report and exit.
	if (!config_read_file(&cfg, fname))
	{
		fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
			config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		exit(-1);
	}


	// GRID.
	config_lookup_float(&cfg, "dr", &dr);
	config_lookup_float(&cfg, "dz", &dz);
	config_lookup_int(&cfg, "NrInterior", &NrInterior);
	config_lookup_int(&cfg, "NzInterior", &NzInterior);
	config_lookup_float(&cfg, "FinalTime", &FinalTime);

	// OUTPUT.
	config_lookup_int(&cfg, "constraintCheck", &constraintCheck);
	config_lookup_int(&cfg, "Noutput2D", &Noutput2D);
	config_lookup_string(&cfg, "dirname", &dirname);

	// BOUNDARY.
	config_lookup_string(&cfg, "boundtype", &boundtype);

	// SLICING.
	config_lookup_string(&cfg, "slicing", &slicing);
	config_lookup_string(&cfg, "ilapse", &ilapse);
	config_lookup_float(&cfg, "gauge_f", &gauge_f);
	config_lookup_float(&cfg, "gl_alpha0", &gl_alpha0);
	config_lookup_float(&cfg, "gl_rr0", &gl_rr0);

	// SHIFT.
	config_lookup_string(&cfg, "shift", &shift);
	config_lookup_string(&cfg, "ishift", &ishift);
	config_lookup_string(&cfg, "bssn_flavor", &bssn_flavor);
	config_lookup_float(&cfg, "driver_c", &driver_c);
	config_lookup_float(&cfg, "driver_eta", &driver_eta);
	config_lookup_float(&cfg, "gs_beta0", &gs_beta0);
	config_lookup_float(&cfg, "gs_rr0", &gs_rr0);

	// EVOLUTION.
	config_lookup_string(&cfg, "order", &order);
	config_lookup_string(&cfg, "integrator", &integrator);
	config_lookup_int(&cfg, "interior_norm", &interior_norm);
	config_lookup_int(&cfg, "adaptive_boundary_norm", &adaptive_boundary_norm);
	config_lookup_float(&cfg, "eta", &eta);
	config_lookup_float(&cfg, "diss", &diss);
	config_lookup_float(&cfg, "damp", &damp);

	/* ELLIPTIC SOLVER */
	config_lookup_string(&cfg, "ell_freq", &ell_freq);
	config_lookup_int(&cfg, "use_permutation", &use_permutation);
	config_lookup_int(&cfg, "perm_count_max", &perm_count_max);
	config_lookup_int(&cfg, "use_preconditioner", &use_preconditioner);
	config_lookup_int(&cfg, "precond_L", &precond_L);
	config_lookup_int(&cfg, "use_low_rank", &use_low_rank);

	// INITIAL DATA.
	config_lookup_string(&cfg, "idata", &idata);

	// BRILL PARAMETERS.
	config_lookup_float(&cfg, "brill_a0", &brill_a0);
	config_lookup_float(&cfg, "brill_sr0", &brill_sr0);
	config_lookup_float(&cfg, "brill_sz0", &brill_sz0);

	// SCHWARSCHILD MASS.
	config_lookup_float(&cfg, "schwar_mass", &schwar_mass);

	// AH FINDER.
	config_lookup_int(&cfg, "ahfind_flag", &ahfind_flag);
	config_lookup_int(&cfg, "ahfind_bicubic", &ahfind_bicubic);
	config_lookup_int(&cfg, "ahfind_freq", &ahfind_freq);
	config_lookup_int(&cfg, "ahfind_mask_filter", &ahfind_mask_filter);
	config_lookup_float(&cfg, "ahfind_start_time", &ahfind_start_time);
	config_lookup_float(&cfg, "ahfind_maxr", &ahfind_maxr);
	config_lookup_float(&cfg, "ahfind_minr", &ahfind_minr);
	config_lookup_int(&cfg, "ahfind_bis_iter", &ahfind_bis_iter);

	// WRITE FLAGS.
	config_lookup_int(&cfg, "alpha", &w_alpha);
	config_lookup_int(&cfg, "a", &w_a);
	config_lookup_int(&cfg, "b", &w_b);
	config_lookup_int(&cfg, "h", &w_h);
	config_lookup_int(&cfg, "c", &w_c);
	config_lookup_int(&cfg, "lambda", &w_lambda);
	config_lookup_int(&cfg, "K", &w_K);
	config_lookup_int(&cfg, "A_a", &w_A_a);
	config_lookup_int(&cfg, "A_b", &w_A_b);
	config_lookup_int(&cfg, "A_h", &w_A_h);
	config_lookup_int(&cfg, "A_c", &w_A_c);
	config_lookup_int(&cfg, "A_lambda", &w_A_lambda);
	config_lookup_int(&cfg, "psi", &w_psi);
	config_lookup_int(&cfg, "phi", &w_phi);
	config_lookup_int(&cfg, "Deltar", &w_Deltar);
	config_lookup_int(&cfg, "Deltaz", &w_Deltaz);
	config_lookup_int(&cfg, "g_a", &w_g_a);
	config_lookup_int(&cfg, "g_b", &w_g_b);
	config_lookup_int(&cfg, "g_h", &w_g_h);
	config_lookup_int(&cfg, "g_c", &w_g_c);
	config_lookup_int(&cfg, "g_lambda", &w_g_lambda);
	config_lookup_int(&cfg, "hdet", &w_hdet);
	config_lookup_int(&cfg, "upA_a", &w_upA_a);
	config_lookup_int(&cfg, "upA_b", &w_upA_b);
	config_lookup_int(&cfg, "upA_h", &w_upA_h);
	config_lookup_int(&cfg, "upA_c", &w_upA_c);
	config_lookup_int(&cfg, "A2_a", &w_A2_a);
	config_lookup_int(&cfg, "A2_b", &w_A2_b);
	config_lookup_int(&cfg, "A2_h", &w_A2_h);
	config_lookup_int(&cfg, "A2_c", &w_A2_c);
	config_lookup_int(&cfg, "A2_lambda", &w_A2_lambda);
	config_lookup_int(&cfg, "A2", &w_A2);
	config_lookup_int(&cfg, "D2alpha_a", &w_D2alpha_a);
	config_lookup_int(&cfg, "D2alpha_b", &w_D2alpha_b);
	config_lookup_int(&cfg, "D2alpha_h", &w_D2alpha_h);
	config_lookup_int(&cfg, "D2alpha_c", &w_D2alpha_c);
	config_lookup_int(&cfg, "D2alpha_lambda", &w_D2alpha_lambda);
	config_lookup_int(&cfg, "LaplaAlpha", &w_LaplaAlpha);
	config_lookup_int(&cfg, "R_a", &w_R_a);
	config_lookup_int(&cfg, "R_b", &w_R_b);
	config_lookup_int(&cfg, "R_h", &w_R_h);
	config_lookup_int(&cfg, "R_c", &w_R_c);
	config_lookup_int(&cfg, "R_lambda", &w_R_lambda);
	config_lookup_int(&cfg, "RSCAL", &w_RSCAL);
	config_lookup_int(&cfg, "divA_r", &w_divA_r);
	config_lookup_int(&cfg, "divA_z", &w_divA_z);
	config_lookup_int(&cfg, "ADeltar", &w_ADeltar);
	config_lookup_int(&cfg, "ADeltaz", &w_ADeltaz);
	config_lookup_int(&cfg, "ham", &w_ham);
	config_lookup_int(&cfg, "mom_r", &w_mom_r);
	config_lookup_int(&cfg, "mom_z", &w_mom_z);
	config_lookup_int(&cfg, "CDeltar", &w_CDeltar);
	config_lookup_int(&cfg, "CDeltaz", &w_CDeltaz);
	config_lookup_int(&cfg, "Clambda", &w_Clambda);
	config_lookup_int(&cfg, "CA_lambda", &w_CA_lambda);
	config_lookup_int(&cfg, "courant_rr", &w_courant_rr);
	config_lookup_int(&cfg, "courant_th", &w_courant_th);
	config_lookup_int(&cfg, "courant_ph", &w_courant_ph);
	config_lookup_int(&cfg, "beta_r", &w_beta_r);
	config_lookup_int(&cfg, "beta_z", &w_beta_z);
}
