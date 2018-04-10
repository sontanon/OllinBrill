#include "tools.h"

#include "param.h"

#undef VERBOSE

// Bilinear interpolator.
double bilinear(const int i1, const int j1, const double di, const double dj, const double *u)
{
	double p0 = u[IDX(i1, j1)] + (u[IDX(i1, j1 + 1)] - u[IDX(i1, j1)]) * dj;
	double p1 = u[IDX(i1 + 1, j1)] + (u[IDX(i1 + 1, j1 + 1)] - u[IDX(i1 + 1, j1)]) * dj;

	return p0 + (p1 - p0) * di;
}

void ah_interpolator(int *flag, const double rr, const double th, 
		double *i_a, double *i_b, double *i_h, 
		double *i_c, double *i_phi,
		double *i_Aa, double *i_Ab, 
		double *i_Ah, double *i_Ac, 
		double *i_K,
		double *i_Dr_a, double *i_Dr_b, 
		double *i_Dr_h, double *i_Dr_c, double *i_Dr_phi,
		double *i_Dz_a, double *i_Dz_b, 
		double *i_Dz_h, double *i_Dz_c, double *i_Dz_phi,
		const double *a, const double *b, const double *h, 
		const double *c, const double *phi,
		const double *Aa, const double *Ab, 
		const double *Ah, const double *Ac, const double *K,
		const double *Dr_a, const double *Dr_b, 
		const double *Dr_h, const double *Dr_c, const double *Dr_phi,
		const double *Dz_a, const double *Dz_b, 
		const double *Dz_h, const double *Dz_c, const double *Dz_phi)
{
	// First calculate the r, z coordinates.
	double t_r  = rr * sin(th);
	double t_z  = rr * cos(th);

	// Now we must find where these coordinates are.
	// First find the floating point coordinate of each.
	double i = t_r / dr - 0.5 + ghost;
	double j = t_z / dz - 0.5 + ghost;

	// Check if interpolation is in-bounds.
	if ((i < 0.) || (i >= NrTotal - 1))
	{
#ifdef VERBOSE
		fprintf(stderr, "INTERPOLATION OUT OF GRID RANGE! rr = %3.3E\n", rr);
#endif
		*flag = 1;
		return;
	}
	else if ((j < 0.) || (j >= NzTotal - 1))
	{
#ifdef VERBOSE
		fprintf(stderr, "INTERPOLATION OUT OF GRID RANGE! rr = %3.3E\n", rr);
#endif
		*flag = 1;
		return;
	}
	else
	{
		*flag = 0;
	}

	// We need the floor and ceiling integers.
	int i1 = (int)floor(i);
	int j1 = (int)floor(j);

	// Get separation from floor integer.
	double di = i - (double)i1;
	double dj = j - (double)j1;

	// With these coorinates, we can fetch the four values involved in
	// the bilinear interpolation.
	*i_a = bilinear(i1, j1, di, dj, a);
	*i_b = bilinear(i1, j1, di, dj, b);
	*i_h = bilinear(i1, j1, di, dj, h);
	*i_c = bilinear(i1, j1, di, dj, c);
	*i_phi = bilinear(i1, j1, di, dj, phi);
	*i_Aa = bilinear(i1, j1, di, dj, Aa);
	*i_Ab = bilinear(i1, j1, di, dj, Ab);
	*i_Ah = bilinear(i1, j1, di, dj, Ah);
	*i_Ac = bilinear(i1, j1, di, dj, Ac);
	*i_K = bilinear(i1, j1, di, dj, K);
	*i_Dr_a = bilinear(i1, j1, di, dj, Dr_a);
	*i_Dr_b = bilinear(i1, j1, di, dj, Dr_b);
	*i_Dr_h = bilinear(i1, j1, di, dj, Dr_h);
	*i_Dr_c = bilinear(i1, j1, di, dj, Dr_c);
	*i_Dr_phi = bilinear(i1, j1, di, dj, Dr_phi);
	*i_Dz_a = bilinear(i1, j1, di, dj, Dz_a);
	*i_Dz_b = bilinear(i1, j1, di, dj, Dz_b);
	*i_Dz_h = bilinear(i1, j1, di, dj, Dz_h);
	*i_Dz_c = bilinear(i1, j1, di, dj, Dz_c);
	*i_Dz_phi = bilinear(i1, j1, di, dj, Dz_phi);
}
