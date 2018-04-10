#include "tools.h"

#include "param.h"

#undef VERBOSE
#define SIMPLE

// Bicubic interpolation using CSR matrix multiplication.
// Base index for CSR multiplication is one.
double bicubic_csr_interpolator(const double f00,
		 const double f01,
		 const double f10,
		 const double f11,
		 const double Dx_f00,
		 const double Dx_f01,
		 const double Dx_f10,
		 const double Dx_f11,
		 const double Dy_f00,
		 const double Dy_f01,
		 const double Dy_f10,
		 const double Dy_f11,
		 const double Dxy_f00,
		 const double Dxy_f01,
		 const double Dxy_f10,
		 const double Dxy_f11,
		 const double dx,
		 const double dy)
{
#ifdef SIMPLE
	double p = 0.0;
	int i, j, k, l;

	const double M[4][4] = {{ 1.0,  0.0,  0.0,  0.0},
				{ 0.0,  0.0,  1.0,  0.0},
				{-3.0,  3.0, -2.0, -1.0},
				{ 2.0, -2.0,  1.0,  1.0}};

	double F[4][4] = {{   f00,    f01,  Dy_f00,  Dy_f01},
			  {   f10,    f11,  Dy_f10,  Dy_f11},
			  {Dx_f00, Dx_f01, Dxy_f00, Dxy_f01},
			  {Dx_f10, Dx_f11, Dxy_f10, Dxy_f11}};

	double x[4] = {1.0, dx, dx * dx, dx * dx * dx};
	double y[4] = {1.0, dy, dy * dy, dy * dy * dy};

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				for (l = 0; l < 4; l++)
				{
					p += x[i] * M[i][j] * F[j][k] * M[l][k] * y[l];
				}
			}
		}
	}

	return p;
#else
	// Matrix system properties.
	const int NNZ = 100;
	const int NROWS = 16;
	const int NCOLS = 16;

	// Matrix multiplication type.
	const char T = 'N';

	// Nonzero elements of the matrix.
	const double  A[NNZ] = { 1.,                                                         //  1 
		                           1.,                                               //  1
			  -3., 3.,        -2.,-1.,                                           //  4
			   2.,-2.,         1., 1.,                                           //  4
	                                                   1.,                               //  1
	                                                                   1.,               //  1
				                          -3., 3.,        -2.,-1.,           //  4
							   2.,-2.,         1., 1.,           //  4
			  -3.,     3.,			  -2.,    -1.,                       //  4
			                  -3.,     3.,                    -2.,    -1.,       //  4
			   9.,-9.,-9., 9., 6., 3.,-6.,-3., 6.,-6., 3.,-3., 4., 2., 2., 1.,   // 16
			  -6., 6., 6.,-6.,-3.,-3., 3., 3.,-4., 4.,-2., 2.,-2.,-2.,-1.,-1.,   // 16
			   2.,    -2.,                     1.,     1.,                       //  4
			                   2.,    -2.,                     1.,     1.,       //  4
			  -6., 6., 6.,-6.,-4.,-2., 4., 2.,-3., 3.,-3., 3.,-2.,-1.,-2.,-1.,   // 16
			   4.,-4.,-4., 4., 2., 2.,-2.,-2., 2.,-2., 2.,-2., 1., 1., 1., 1. }; // 16

	// Column indices.
	const int jA[NNZ] = {  1,                                                            //  1
		                            5,                                               //  2
			    1,  2,          5,  6,                                           //  3
			    1,  2,          5,  6,                                           //  7
	                                                    9,                               // 11
	                                                                   13,               // 12
				                            9, 10,         13, 14,           // 13
							    9, 10,         13, 14,           // 17
			    1,      3,			    9,     11,                       // 21
			                    5,      7,                     13,     15,       // 25
	                    1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,   // 29
	                    1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,   // 45
			    1,      3,                      9,     11,                       // 61
			                    5,      7,                     13,     15,       // 65
	                    1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,   // 69
	                    1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16 }; // 85

	// Row starts.
	const int iA[NROWS + 1] = { 1, 2, 3, 7, 11, 12, 13, 17, 21, 25, 29, 45, 61, 65, 69, 85, 101 };

	// Values of the function and its derivatives.
	double x[NCOLS] = { f00, f10, f01, f11, Dx_f00, Dx_f10, Dx_f01, Dx_f11, Dy_f00, Dy_f10, Dy_f01, Dy_f11, Dxy_f00, Dxy_f10, Dxy_f01, Dxy_f11 };

	// Coefficients.
	double a[NROWS] = { 0.0 };

	// Perform CSR matrix multiplication.
	mkl_dcsrgemv(&T, &NROWS, A, iA, jA, x, a);

	// Now do dot product using step sizes obtained as follows:
	double d00 = 1.0;
	double d10 = dx;
	double d01 = dy;
	double d11 = d10 * d01;
	double d20 = dx * d10;
	double d21 = dx * d11;
	double d02 = dy * d01;
	double d12 = dy * d11;
	double d22 = d20 * d02;
	double d30 = dx * d20;
	double d31 = dx * d21;
	double d32 = dx * d22;
	double d03 = dy * d02;
	double d13 = dy * d12;
	double d23 = dy * d22;
	double d33 = d30 * d03;

	// Step sizes vector.
	double d[NROWS] = { 1.0, d10, d20, d30, d01, d11, d21, d31, d02, d12, d22, d32, d03, d13, d23, d33 };

	// Result.
	double p = 0.0;

	// Do dot product.
	p = cblas_ddot(NROWS, a, 1, d, 1);

	// Return interpolated value.
	return p;
#endif
}

double bicubic(const int i1, const int j1, const double di, const double dj, const double *u, const double *Dr_u, const double *Dz_u, const double *Drz_u)
{
	return bicubic_csr_interpolator(u[IDX(i1, j1)], u[IDX(i1, j1 + 1)], u[IDX(i1 + 1, j1)], u[IDX(i1 + 1, j1 + 1)],
	dr * Dr_u[IDX(i1, j1)], dr * Dr_u[IDX(i1, j1 + 1)], dr * Dr_u[IDX(i1 + 1, j1)], dr * Dr_u[IDX(i1 + 1, j1 + 1)],
	dz * Dz_u[IDX(i1, j1)], dz * Dz_u[IDX(i1, j1 + 1)], dz * Dz_u[IDX(i1 + 1, j1)], dz * Dz_u[IDX(i1 + 1, j1 + 1)],
	dr * dz * Drz_u[IDX(i1, j1)], dr * dz * Drz_u[IDX(i1, j1 + 1)], dr * dz * Drz_u[IDX(i1 + 1, j1)], dr * dz * Drz_u[IDX(i1 + 1, j1 + 1)],
	di, dj);
}

void ah_bicubic_interpolator(int *flag, const double rr, const double th, 
		double *i_a, double *i_b, double *i_h, 
		double *i_c, double *i_phi,
		double *i_Aa, double *i_Ab, 
		double *i_Ah, double *i_Ac, 
		double *i_K,
		double *i_Dr_a, double *i_Dr_b, 
		double *i_Dr_h, double *i_Dr_c, double *i_Dr_phi,
		double *i_Dz_a, double *i_Dz_b, 
		double *i_Dz_h, double *i_Dz_c, double *i_Dz_phi,
		const double *a, 
		const double *Dr_a,
		const double *Dz_a,
		const double *Drz_a,
		const double *b, 
		const double *Dr_b,
		const double *Dz_b,
		const double *Drz_b,
		const double *h,
		const double *Dr_h,
		const double *Dz_h,
		const double *Drz_h,
		const double *c, 
		const double *Dr_c,
		const double *Dz_c,
		const double *Drz_c,
		const double *phi,
		const double *Dr_phi,
		const double *Dz_phi,
		const double *Drz_phi,
		const double *Aa, 
		const double *Dr_Aa,
		const double *Dz_Aa,
		const double *Drz_Aa,
		const double *Ab, 
		const double *Dr_Ab,
		const double *Dz_Ab,
		const double *Drz_Ab,
		const double *Ah, 
		const double *Dr_Ah,
		const double *Dz_Ah,
		const double *Drz_Ah,
		const double *Ac, 
		const double *Dr_Ac,
		const double *Dz_Ac,
		const double *Drz_Ac,
		const double *K,
		const double *Dr_K,
		const double *Dz_K,
		const double *Drz_K,
		//const double *Dr_a,
		const double *Drr_a, 
		//const double *Drz_a,
		const double *Drrz_a,
		//const double *Dr_b,
		const double *Drr_b, 
		//const double *Drz_b,
		const double *Drrz_b,
		//const double *Dr_h,
		const double *Drr_h, 
		//const double *Drz_h,
		const double *Drrz_h,
		//const double *Dr_c, 
		const double *Drr_c,
		//const doubel *Drz_c,
		const double *Drrz_c,
		//const double *Dr_phi,
		const double *Drr_phi,
		//const double *Drz_phi,
		const double *Drrz_phi,
		//const double *Dz_a,
		//const double *Drz_a, 
		const double *Dzz_a,
		const double *Drzz_a,
		//const double *Dz_b,
		//const double *Drz_b, 
		const double *Dzz_b,
		const double *Drzz_b,
		//const double *Dz_h,
		//const double *Drz_h, 
		const double *Dzz_h,
		const double *Drzz_h,
		//const double *Dz_c, 
		//const double *Drz_c,
		const double *Dzz_c,
		const double *Drzz_c,
		//const double *Dz_phi,
		//const double *Drz_phi,
		const double *Dzz_phi,
		const double *Drzz_phi)
{
	// First calculate the r, z coordinates.
	double t_r  = rr * sin(th);
	double t_z  = rr * cos(th);

	// Now we must find where these coordinates are.
	// First find the floating point coordinate of each.
	double i = t_r / dr - 0.5 + ghost;
	double j = t_z / dz - 0.5 + ghost;

	// We need the floor and ceiling integers.
	int i1 = (int)floor(i);
	int j1 = (int)floor(j);
	int i2 = i1 + 1;
	int j2 = j1 + 1;

	// Check if interpolation is in-bounds.
	if ((i1 < 0) || (i1 >= NrTotal - 1))
	{
#ifdef VERBOSE
		fprintf(stderr, "INTERPOLATION OUT OF GRID RANGE! rr = %3.3E\n", rr);
#endif
		*flag = 1;
		return;
	}
	else if ((j1 < 0) || (j1 >= NzTotal - 1))
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

	// Get separation from floor integer.
	double di = i - (double)i1;
	double dj = j - (double)j1;

	// With these coordinates, we can fetch the values via bicubic interpolation.
	*i_a = bicubic(i1, j1, di, dj, a, Dr_a, Dz_a, Drz_a);
	*i_b = bicubic(i1, j1, di, dj, b, Dr_b, Dz_b, Drz_b);
	*i_h = bicubic(i1, j1, di, dj, h, Dr_h, Dz_h, Drz_h);
	*i_c = bicubic(i1, j1, di, dj, c, Dr_c, Dz_c, Drz_c);
	*i_phi = bicubic(i1, j1, di, dj, phi, Dr_phi, Dz_phi, Drz_phi);

	*i_Aa = bicubic(i1, j1, di, dj, Aa, Dr_Aa, Dz_Aa, Drz_Aa);
	*i_Ab = bicubic(i1, j1, di, dj, Ab, Dr_Ab, Dz_Ab, Drz_Ab);
	*i_Ah = bicubic(i1, j1, di, dj, Ah, Dr_Ah, Dz_Ah, Drz_Ah);
	*i_Ac = bicubic(i1, j1, di, dj, Ac, Dr_Ac, Dz_Ac, Drz_Ac);
	*i_K = bicubic(i1, j1, di, dj, K, Dr_K, Dz_K, Drz_K);

	*i_Dr_a = bicubic(i1, j1, di, dj, Dr_a, Drr_a, Drz_a, Drrz_a);
	*i_Dr_b = bicubic(i1, j1, di, dj, Dr_b, Drr_b, Drz_b, Drrz_b);
	*i_Dr_h = bicubic(i1, j1, di, dj, Dr_h, Drr_h, Drz_h, Drrz_h);
	*i_Dr_c = bicubic(i1, j1, di, dj, Dr_c, Drr_c, Drz_c, Drrz_c);
	*i_Dr_phi = bicubic(i1, j1, di, dj, Dr_phi, Drr_phi, Drz_phi, Drrz_phi);

	*i_Dz_a = bicubic(i1, j1, di, dj, Dz_a, Drz_a, Dzz_a, Drzz_a);
	*i_Dz_b = bicubic(i1, j1, di, dj, Dz_b, Drz_b, Dzz_b, Drzz_b);
	*i_Dz_h = bicubic(i1, j1, di, dj, Dz_h, Drz_h, Dzz_h, Drzz_h);
	*i_Dz_c = bicubic(i1, j1, di, dj, Dz_c, Drz_c, Dzz_c, Drzz_c);
	*i_Dz_phi = bicubic(i1, j1, di, dj, Dz_phi, Drz_phi, Dzz_phi, Drzz_phi);
}
