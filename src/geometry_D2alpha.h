void calculate_D2alpha(double *D2alpha_a, double *D2alpha_b, double *D2alpha_h, double *D2alpha_c,
	double *D2alpha_lambda, double *LaplaAlpha,
	const double r, const double dRalpha, const double dZalpha,
	const double dRRalpha, const double dZZalpha, const double dRZalpha, const double ddRalpha,
	const double ga, const double gb, const double gh, const double gc, const double glambda,
	const double dRa, const double dRb, const double dRh, const double dRc,
	const double dZa, const double dZb, const double dZh, const double dZc,
	const double c, const double h, const double lambda, const double dZlambda);

void calculate_phys_D2alpha(double *D2alpha_a, double *D2alpha_b, double *D2alpha_h,
	double *D2alpha_c, double *D2alpha_lambda, double *LaplaAlpha,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double ga, const double gb, const double gh, const double gc,
	const double dRalpha, const double dZalpha, const double dRphi, const double dZphi,
	const double psi4);