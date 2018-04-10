int nnz_flat_laplacian(const int NrInterior, const int NzInterior, const int order, const int robin);

void csr_gen_flat_laplacian(csr_matrix A,
	const int NrInterior,
	const int NzInterior,
	const int order,
	const double dr,
	const double dz,
	const double *s,
	double *f,
	const double uInf,
	const int robin,
	const int r_sym,
	const int z_sym);
