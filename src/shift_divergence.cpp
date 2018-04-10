void beta_div(double *div_beta, const double r, const double beta_r, const double beta_z,
	const double Dr_beta_r, const double Dz_beta_z, const double hdet,
	const double Dr_hdet, const double Dz_hdet)
{
	*div_beta = Dr_beta_r + Dz_beta_z + beta_r / r + 0.5 * (beta_r * Dr_hdet + beta_z * Dz_hdet) / hdet;
}