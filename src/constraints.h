void constraints(const int l);

void calculate_constraints(double *ham, double *mom_r, double *mom_z,
	const double R, const double A2, const double divA_r, const double divA_z,
	const double upA_a, const double upA_b, const double upA_c,
	const double r, const double Dr_phi, const double Dz_phi, const double psi4,
	const double K, const double Dr_K, const double Dz_K,
	const double g_a, const double g_b, const double g_c,
	const bool maximal);
