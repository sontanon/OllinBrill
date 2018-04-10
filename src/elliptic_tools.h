// Reduce array to the elliptic solver size.
void ghost_reduce(const double *u, double *g_u, const int NrInterior, const int NzInterior, const int ghost);
// Fill array from elliptic solver sized-array.
void ghost_fill(const double *g_u, double *u, const int r_sym, const int z_sym, const int NrInterior, const int NzInterior, const int ghost);
