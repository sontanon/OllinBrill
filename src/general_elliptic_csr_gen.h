int nnz_general_elliptic(const int NrInterior, const int NzInterior, const int order, const int robin);

void csr_gen_general_elliptic(csr_matrix A,
	const int NrInterior,
	const int NzInterior,
	const int order,
	const double dr,
	const double dz,
	const double *ell_a,
	const double *ell_b,
	const double *ell_c,
	const double *ell_d,
	const double *ell_e,
	const double *ell_s,
	double *ell_f,
	const double uInf,
	const int robin,
	const int r_sym,
	const int z_sym);
