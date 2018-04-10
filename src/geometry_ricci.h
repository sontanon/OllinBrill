#pragma once
void calculate_ricci(double *R_a, double *R_b, double *R_h, double *R_c, double *R_lambda, double *R,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double g_a, const double g_b, const double g_h, const double g_c, const double g_lambda,
	const double Dr_a, const double Dr_b, const double Dr_h, const double Dr_c, const double Dr_lambda,
	const double Dz_a, const double Dz_b, const double Dz_h, const double Dz_c, const double Dz_lambda,
	const double Drr_a, const double Drr_b, const double Drr_h, const double Drr_c, const double Drr_lambda,
	const double Dzz_a, const double Dzz_b, const double Dzz_h, const double Dzz_c, const double Dzz_lambda,
	const double Drz_a, const double Drz_b, const double Drz_h, const double Drz_c, const double Drz_lambda,
	const double Deltar, const double Deltaz, const double Dr_Deltar,
	const double Dr_Deltaz, const double Dz_Deltar, const double Dz_Deltaz, const double DDeltar_r);

void calculate_phys_ricci(double *Ra, double *Rb, double *Rh, double *Rc, double *Rlambda, double *R,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double ga, const double gb, const double gh, const double gc, const double glambda,
	const double dRa, const double dRb, const double dRh, const double dRc, const double dRlambda,
	const double dZa, const double dZb, const double dZh, const double dZc, const double dZlambda,
	const double dRphi, const double dZphi, const double dRRphi, const double dZZphi, const double dRZphi, const double ddRphi, const double psi4);
