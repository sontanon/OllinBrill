#include "tools.h"

#include "param.h"
#include "arrays.h"
#include "files.h"
#include "write_flags.h"

void write_time_file(FILE *f_u, const double *u, const double t)
{
	// Auxiliary integers.
	MKL_INT i, j;

	// Print commented time at first line.
	fprintf(f_u, "# t = %9.18E\n", t);

	// Print main body.
	for (i = 0; i < NrTotal; i++)
	{
		for (j = 0; j < NzTotal - 1; j++)
		{
			fprintf(f_u, "%9.18E\t", u[IDX(i, j)]);
		}
		// Print new line.
		j = NzTotal - 1;
		fprintf(f_u, "%9.18E\n", u[IDX(i, j)]);
	}
	// Print blank line.
	fprintf(f_u, "\n");
}

void write_files(const double t)
{
	if (w_alpha)
		write_time_file(f_alpha, alpha, t);

	if (w_a)
		write_time_file(f_a, a, t);
	if (w_b)
		write_time_file(f_b, b, t);
	if (w_c)
		write_time_file(f_c, c, t);
	if (w_h)
		write_time_file(f_h, h, t);
	if (w_lambda)
		write_time_file(f_lambda, lambda, t);

	if (w_K)
		write_time_file(f_K, K, t);
	if (w_A_a)
		write_time_file(f_A_a, A_a, t);
	if (w_A_b)
		write_time_file(f_A_b, A_b, t);
	if (w_A_c)
		write_time_file(f_A_c, A_c, t);
	if (w_A_h)
		write_time_file(f_A_h, A_h, t);
	if (w_A_lambda)
		write_time_file(f_A_lambda, A_lambda, t);

	if (w_psi)
		write_time_file(f_psi, psi, t);
	if (w_phi)
		write_time_file(f_phi, phi, t);

	if (w_Deltar)
		write_time_file(f_Deltar, Deltar, t);
	if (w_Deltaz)
		write_time_file(f_Deltaz, Deltaz, t);

	if (w_g_a)
		write_time_file(f_g_a, g_a, t);
	if (w_g_b)
		write_time_file(f_g_b, g_b, t);
	if (w_g_c)
		write_time_file(f_g_c, g_c, t);
	if (w_g_h)
		write_time_file(f_g_h, g_h, t);
	if (w_g_lambda)
		write_time_file(f_g_lambda, g_lambda, t);
	if (w_hdet)
		write_time_file(f_hdet, hdet, t);

	if (w_upA_a)
		write_time_file(f_upA_a, upA_a, t);
	if (w_upA_b)
		write_time_file(f_upA_b, upA_b, t);
	if (w_upA_c)
		write_time_file(f_upA_c, upA_c, t);
	if (w_upA_h)
		write_time_file(f_upA_h, upA_h, t);

	if (w_A2_a)
		write_time_file(f_A2_a, A2_a, t);
	if (w_A2_b)
		write_time_file(f_A2_b, A2_b, t);
	if (w_A2_c)
		write_time_file(f_A2_c, A2_c, t);
	if (w_A2_h)
		write_time_file(f_A2_h, A2_h, t);
	if (w_A2_lambda)
		write_time_file(f_A2_lambda, A2_lambda, t);
	if (w_A2)
		write_time_file(f_A2, A2, t);

	if (w_ell_a)
		write_time_file(f_ell_a, ell_a, t);
	if (w_ell_b)
		write_time_file(f_ell_b, ell_b, t);
	if (w_ell_c)
		write_time_file(f_ell_c, ell_c, t);
	if (w_ell_d)
		write_time_file(f_ell_d, ell_d, t);
	if (w_ell_e)
		write_time_file(f_ell_e, ell_e, t);
	if (w_ell_s)
		write_time_file(f_ell_s, ell_s, t);
	if (w_ell_f)
		write_time_file(f_ell_f, ell_f, t);

	if (w_D2alpha_a)
		write_time_file(f_D2alpha_a, D2alpha_a, t);
	if (w_D2alpha_b)
		write_time_file(f_D2alpha_b, D2alpha_b, t);
	if (w_D2alpha_h)
		write_time_file(f_D2alpha_h, D2alpha_h, t);
	if (w_D2alpha_c)
		write_time_file(f_D2alpha_c, D2alpha_c, t);
	if (w_D2alpha_lambda)
		write_time_file(f_D2alpha_lambda, D2alpha_lambda, t);

	if (w_LaplaAlpha)
		write_time_file(f_LaplaAlpha, LaplaAlpha, t);

	if (w_R_a)
		write_time_file(f_R_a, R_a, t);
	if (w_R_b)
		write_time_file(f_R_b, R_b, t);
	if (w_R_h)
		write_time_file(f_R_h, R_h, t);
	if (w_R_c)
		write_time_file(f_R_c, R_c, t);
	if (w_R_lambda)
		write_time_file(f_R_lambda, R_lambda, t);
	if (w_RSCAL)
		write_time_file(f_RSCAL, RSCAL, t);

	if (w_divA_r)
		write_time_file(f_divA_r, divA_r, t);
	if (w_divA_z)
		write_time_file(f_divA_z, divA_z, t);

	if (w_ADeltar)
		write_time_file(f_ADeltar, ADelta_r, t);
	if (w_ADeltaz)
		write_time_file(f_ADeltaz, ADelta_z, t);

	if (w_ham)
		write_time_file(f_ham, ham, t);
	if (w_mom_r)
		write_time_file(f_mom_r, mom_r, t);
	if (w_mom_z)
		write_time_file(f_mom_z, mom_z, t);
	if (w_CDeltar)
		write_time_file(f_CDeltar, CDeltar, t);
	if (w_CDeltaz)
		write_time_file(f_CDeltaz, CDeltaz, t);
	if (w_Clambda)
		write_time_file(f_Clambda, Clambda, t);
	if (w_CA_lambda)
		write_time_file(f_CA_lambda, CA_lambda, t);

	if (w_courant_rr)
		write_time_file(f_courant_rr, courant_rr, t);
	if (w_courant_th)
		write_time_file(f_courant_th, courant_th, t);
	if (w_courant_ph)
		write_time_file(f_courant_ph, courant_ph, t);

	if (w_beta_r)
		write_time_file(f_beta_r, beta_r, t);
	if (w_beta_z)
		write_time_file(f_beta_z, beta_z, t);
	if (w_dtbeta_r)
		write_time_file(f_dtbeta_r, dtbeta_r, t);
	if (w_dtbeta_z)
		write_time_file(f_dtbeta_z, dtbeta_z, t);
}

void open_files(void)
{
	if (w_alpha)
		f_alpha = fopen("alpha.asc", "w");

	if (w_a)
		f_a = fopen("a.asc", "w");
	if (w_b)
		f_b = fopen("b.asc", "w");
	if (w_c)
		f_c = fopen("c.asc", "w");
	if (w_h)
		f_h = fopen("h.asc", "w");
	if (w_lambda)
		f_lambda = fopen("lambda.asc", "w");

	if (w_K)
		f_K = fopen("K.asc", "w");
	if (w_A_a)
		f_A_a = fopen("A_a.asc", "w");
	if (w_A_b)
		f_A_b = fopen("A_b.asc", "w");
	if (w_A_c)
		f_A_c = fopen("A_c.asc", "w");
	if (w_A_h)
		f_A_h = fopen("A_h.asc", "w");
	if (w_A_lambda)
		f_A_lambda = fopen("A_lambda.asc", "w");

	if (w_psi)
		f_psi = fopen("psi.asc", "w");
	if (w_phi)
		f_phi = fopen("phi.asc", "w");

	if (w_Deltar)
		f_Deltar = fopen("Deltar.asc", "w");
	if (w_Deltaz)
		f_Deltaz = fopen("Deltaz.asc", "w");

	if (w_g_a)
		f_g_a = fopen("g_a.asc", "w");
	if (w_g_b)
		f_g_b = fopen("g_b.asc", "w");
	if (w_g_c)
		f_g_c = fopen("g_c.asc", "w");
	if (w_g_h)
		f_g_h = fopen("g_h.asc", "w");
	if (w_g_lambda)
		f_g_lambda = fopen("g_lambda.asc", "w");
	if (w_hdet)
		f_hdet = fopen("hdet.asc", "w");

	if (w_upA_a)
		f_upA_a = fopen("upA_a.asc", "w");
	if (w_upA_b)
		f_upA_b = fopen("upA_b.asc", "w");
	if (w_upA_c)
		f_upA_c = fopen("upA_c.asc", "w");
	if (w_upA_h)
		f_upA_h = fopen("upA_h.asc", "w");

	if (w_A2_a)
		f_A2_a = fopen("A2_a.asc", "w");
	if (w_A2_b)
		f_A2_b = fopen("A2_b.asc", "w");
	if (w_A2_c)
		f_A2_c = fopen("A2_c.asc", "w");
	if (w_A2_h)
		f_A2_h = fopen("A2_h.asc", "w");
	if (w_A2_lambda)
		f_A2_lambda = fopen("A2_lambda.asc", "w");
	if (w_A2)
		f_A2 = fopen("A2.asc", "w");

	if (w_ell_a)
		f_ell_a = fopen("ell_a.asc", "w");
	if (w_ell_b)
		f_ell_b = fopen("ell_b.asc", "w");
	if (w_ell_c)
		f_ell_c = fopen("ell_c.asc", "w");
	if (w_ell_d)
		f_ell_d = fopen("ell_d.asc", "w");
	if (w_ell_e)
		f_ell_e = fopen("ell_e.asc", "w");
	if (w_ell_s)
		f_ell_s = fopen("ell_s.asc", "w");
	if (w_ell_f)
		f_ell_f = fopen("ell_f.asc", "w");

	if (w_D2alpha_a)
		f_D2alpha_a = fopen("D2alpha_a.asc", "w");
	if (w_D2alpha_b)
		f_D2alpha_b = fopen("D2alpha_b.asc", "w");
	if (w_D2alpha_h)
		f_D2alpha_h = fopen("D2alpha_h.asc", "w");
	if (w_D2alpha_c)
		f_D2alpha_c = fopen("D2alpha_c.asc", "w");
	if (w_D2alpha_lambda)
		f_D2alpha_lambda = fopen("D2alpha_lambda.asc", "w");

	if (w_LaplaAlpha)
		f_LaplaAlpha = fopen("LaplaAlpha.asc", "w");

	if (w_R_a)
		f_R_a = fopen("R_a.asc", "w");
	if (w_R_b)
		f_R_b = fopen("R_b.asc", "w");
	if (w_R_h)
		f_R_h = fopen("R_h.asc", "w");
	if (w_R_c)
		f_R_c = fopen("R_c.asc", "w");
	if (w_R_lambda)
		f_R_lambda = fopen("R_lambda.asc", "w");
	if (w_RSCAL)
		f_RSCAL = fopen("RSCAL.asc", "w");

	if (w_divA_r)
		f_divA_r = fopen("divA_r.asc", "w");
	if (w_divA_z)
		f_divA_z = fopen("divA_z.asc", "w");

	if (w_ADeltar)
		f_ADeltar = fopen("ADeltar.asc", "w");
	if (w_ADeltaz)
		f_ADeltaz = fopen("ADeltaz.asc", "w");

	if (w_ham)
		f_ham = fopen("ham.asc", "w");
	if (w_mom_r)
		f_mom_r = fopen("mom_r.asc", "w");
	if (w_mom_z)
		f_mom_z = fopen("mom_z.asc", "w");
	if (w_CDeltar)
		f_CDeltar = fopen("CDeltar.asc", "w");
	if (w_CDeltaz)
		f_CDeltaz = fopen("CDeltaz.asc", "w");
	if (w_Clambda)
		f_Clambda = fopen("Clambda.asc", "w");
	if (w_CA_lambda)
		f_CA_lambda = fopen("CA_lambda.asc", "w");

	if (w_t_ham_norm)
		f_t_ham_norm = fopen("t_ham_norm.asc", "w");
	if (w_t_mom_r_norm)
		f_t_mom_r_norm = fopen("t_mom_r_norm.asc", "w");
	if (w_t_mom_z_norm)
		f_t_mom_z_norm = fopen("t_mom_z_norm.asc", "w");
	if (w_t_CDeltar_norm)
		f_t_CDeltar_norm = fopen("t_CDeltar_norm.asc", "w");
	if (w_t_CDeltaz_norm)
		f_t_CDeltaz_norm = fopen("t_CDeltaz_norm.asc", "w");
	if (w_t_Clambda_norm)
		f_t_Clambda_norm = fopen("t_Clambda_norm.asc", "w");
	if (w_t_CA_lambda_norm)
		f_t_CA_lambda_norm = fopen("t_CA_lambda_norm.asc", "w");

	if (w_courant_rr)
		f_courant_rr = fopen("courant_rr.asc", "w");
	if (w_courant_th)
		f_courant_th = fopen("courant_th.asc", "w");
	if (w_courant_ph)
		f_courant_ph = fopen("courant_ph.asc", "w");

	if (w_beta_r)
		f_beta_r = fopen("beta_r.asc", "w");
	if (w_beta_z)
		f_beta_z = fopen("beta_z.asc", "w");
	if (w_dtbeta_r)
		f_dtbeta_r = fopen("dtbeta_r.asc", "w");
	if (w_dtbeta_z)
		f_dtbeta_z = fopen("dtbeta_z.asc", "w");
}

void close_files(void)
{
	if (w_alpha)
		fclose(f_alpha);

	if (w_a)
		fclose(f_a);
	if (w_b)
		fclose(f_b);
	if (w_c)
		fclose(f_c);
	if (w_h)
		fclose(f_h);
	if (w_lambda)
		fclose(f_lambda);

	if (w_K)
		fclose(f_K);
	if (w_A_a)
		fclose(f_A_a);
	if (w_A_b)
		fclose(f_A_b);
	if (w_A_c)
		fclose(f_A_c);
	if (w_A_h)
		fclose(f_A_h);
	if (w_A_lambda)
		fclose(f_A_lambda);

	if (w_psi)
		fclose(f_psi);
	if (w_phi)
		fclose(f_phi);

	if (w_Deltar)
		fclose(f_Deltar);
	if (w_Deltaz)
		fclose(f_Deltaz);

	if (w_g_a)
		fclose(f_g_a);
	if (w_g_b)
		fclose(f_g_b);
	if (w_g_c)
		fclose(f_g_c);
	if (w_g_h)
		fclose(f_g_h);
	if (w_g_lambda)
		fclose(f_g_lambda);
	if (w_hdet)
		fclose(f_hdet);

	if (w_upA_a)
		fclose(f_upA_a);
	if (w_upA_b)
		fclose(f_upA_b);
	if (w_upA_c)
		fclose(f_upA_c);
	if (w_upA_h)
		fclose(f_upA_h);

	if (w_A2_a)
		fclose(f_A2_a);
	if (w_A2_b)
		fclose(f_A2_b);
	if (w_A2_c)
		fclose(f_A2_c);
	if (w_A2_h)
		fclose(f_A2_h);
	if (w_A2_lambda)
		fclose(f_A2_lambda);
	if (w_A2)
		fclose(f_A2);

	if (w_ell_a)
		fclose(f_ell_a);
	if (w_ell_b)
		fclose(f_ell_b);
	if (w_ell_c)
		fclose(f_ell_c);
	if (w_ell_d)
		fclose(f_ell_d);
	if (w_ell_e)
		fclose(f_ell_e);
	if (w_ell_s)
		fclose(f_ell_s);
	if (w_ell_f)
		fclose(f_ell_f);

	if (w_D2alpha_a)
		fclose(f_D2alpha_a);
	if (w_D2alpha_b)
		fclose(f_D2alpha_b);
	if (w_D2alpha_h)
		fclose(f_D2alpha_h);
	if (w_D2alpha_c)
		fclose(f_D2alpha_c);
	if (w_D2alpha_lambda)
		fclose(f_D2alpha_lambda);

	if (w_LaplaAlpha)
		fclose(f_LaplaAlpha);

	if (w_R_a)
		fclose(f_R_a);
	if (w_R_b)
		fclose(f_R_b);
	if (w_R_h)
		fclose(f_R_h);
	if (w_R_c)
		fclose(f_R_c);
	if (w_R_lambda)
		fclose(f_R_lambda);
	if (w_RSCAL)
		fclose(f_RSCAL);

	if (w_divA_r)
		fclose(f_divA_r);
	if (w_divA_z)
		fclose(f_divA_z);

	if (w_ADeltar)
		fclose(f_ADeltar);
	if (w_ADeltaz)
		fclose(f_ADeltaz);

	if (w_ham)
		fclose(f_ham);
	if (w_mom_r)
		fclose(f_mom_r);
	if (w_mom_z)
		fclose(f_mom_z);
	if (w_CDeltar)
		fclose(f_CDeltar);
	if (w_CDeltaz)
		fclose(f_CDeltaz);
	if (w_Clambda)
		fclose(f_Clambda);
	if (w_CA_lambda)
		fclose(f_CA_lambda);

	if (w_t_ham_norm)
		fclose(f_t_ham_norm);
	if (w_t_mom_r_norm)
		fclose(f_t_mom_r_norm);
	if (w_t_mom_z_norm)
		fclose(f_t_mom_z_norm);
	if (w_t_CDeltar_norm)
		fclose(f_t_CDeltar_norm);
	if (w_t_CDeltaz_norm)
		fclose(f_t_CDeltaz_norm);
	if (w_t_Clambda_norm)
		fclose(f_t_Clambda_norm);
	if (w_t_CA_lambda_norm)
		fclose(f_t_CA_lambda_norm);

	if (w_courant_rr)
		fclose(f_courant_rr);
	if (w_courant_th)
		fclose(f_courant_th);
	if (w_courant_ph)
		fclose(f_courant_ph);

	if (w_beta_r)
		fclose(f_beta_r);
	if (w_beta_z)
		fclose(f_beta_z);
	if (w_dtbeta_r)
		fclose(f_dtbeta_r);
	if (w_dtbeta_z)
		fclose(f_dtbeta_z);
}
