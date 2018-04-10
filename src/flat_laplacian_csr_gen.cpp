#include "tools.h"

// One-based indexing.
#define BASE 1

// Print CSR matrix for debug.
#undef DEBUG

// Nonzero calculator.
int nnz_flat_laplacian(const int NrInterior, const int NzInterior, const int order, const int robin)
{
	int nnz;

	int n_robin = robin + order;

	// Second order Laplacian.
	if (order == 2)
	{
		nnz = 5 * NrInterior * NzInterior
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	// Fourth-order Laplacian.
	else if (order == 4)
	{
		nnz = 9 * (NrInterior - 2) * (NzInterior - 2)
			+ (8 + 10) * (NrInterior + NzInterior - 4)
			+ (7 + 9 + 9 + 11)
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	return nnz;
}

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
	const int z_sym)
{
	// Constant numbers.
	const double third = 1.0 / 3.0;
	const double sixth = 1.0 / 6.0;
	const double twelfth = 1.0 / 12.0;

	// Number of nonzero elements.
	int nnz = nnz_flat_laplacian(NrInterior, NzInterior, order, robin);

	// Number of Robin elements.
	int n_robin = robin + order;

	// Number of elements we have filled in.
	int offset = 0;

	// Auxiliary variables.
	double r, z, rr2, ir;
	double robin1, robin2, robin3, robin4, robin5, robin6, robin7;
	int i, j, t_offset;

	// SECOND-ORDER FLAT LAPLACIAN.
	if (order == 2)
	{
		// Lower-left corner: diagonal symmetry.
		A.ia[IDX(0, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)(r_sym * z_sym);
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		f[IDX(0, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
#pragma omp parallel shared(A, f) private(offset)
		{
#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each j iteration fills 2 elements.
				offset = t_offset + 2 * (j - 1);
				A.ia[IDX(0, j)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)r_sym;
				A.ja[offset] = BASE + IDX(0, j);
				A.ja[offset + 1] = BASE + IDX(1, j);
				f[IDX(0, j)] = 0.0;
			}
		}

		// We have now filled:
		offset = 2 + 2 * NzInterior;

		// Upper-left corner: also axis symmetry.
		A.ia[IDX(0, NzInterior + 1)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)r_sym;
		A.ja[offset] = BASE + IDX(0, NzInterior + 1);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior + 1);
		f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Now come the interior points plus the top and bottom boundaries with
		// Robin and equatorial symmetry respectively.
#pragma omp parallel shared(A, f) private(offset, j, r, z, ir, rr2, robin1, robin2, robin3)
		{
#pragma omp for schedule(guided)
			for (i = 1; i < NrInterior + 1; i++)
			{
				// Each iteration of i loop will fill 5 * NzInterior + (2 + n_robin) values.
				offset = t_offset + (i - 1) * (5 * NzInterior + 2 + n_robin);


				// R coordinate.
				r = (double)i - 0.5;
				// Inverse.
				ir = 1.0 / r;

				// Do bottom boundary first with equatorial symmetry.
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Now loop over interior points.
				for (j = 1; j < NzInterior + 1; j++)
				{
					// Row begins at offset.
					A.ia[IDX(i, j)] = BASE + offset;
					// Values.
					A.a[offset] = 1.0 - 0.5 * ir;
					A.a[offset + 1] = 1.0;
					A.a[offset + 2] = dr * dz * s[IDX(i, j)] - 4.0;
					A.a[offset + 3] = 1.0;
					A.a[offset + 4] = 1.0 + 0.5 * ir;
					// Columns.
					A.ja[offset] = BASE + IDX(i - 1, j);
					A.ja[offset + 1] = BASE + IDX(i, j - 1);
					A.ja[offset + 2] = BASE + IDX(i, j);
					A.ja[offset + 3] = BASE + IDX(i, j + 1);
					A.ja[offset + 4] = BASE + IDX(i + 1, j);
					// Multiply RHS by dr * dz.
					f[IDX(i, j)] *= dr * dz;
					// Increase offset by five.
					offset += 5;
				}

				// Now fill the top boundary with Robin.
				j = NzInterior + 1;
				// Z coordinate.
				z = (double)j - 0.5;
				// Radial coordinate.
				rr2 = r * r + z * z;
				A.ia[IDX(i, NzInterior + 1)] = BASE + offset;
				switch (robin)
				{
				case 1:
					A.a[offset] = 0.5 * rr2 / z;
					A.a[offset + 1] = -2.0 * rr2 / z;
					A.a[offset + 2] = 1.0 + 1.5 * rr2 / z;
					A.ja[offset] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior + 1);
					break;
				case 2:
					robin1 = 0.5 * rr2 * rr2 / (z * z);
					robin2 = 0.25 * (4.0 * rr2 / z - r * r / sqrt(rr2));
					A.a[offset] = -robin1;
					A.a[offset + 1] = 4.0 * robin1 + robin2;
					A.a[offset + 2] = -5.0 * robin1 - 4.0 * robin2;
					A.a[offset + 3] = 2.0 * robin1 + 3.0 * robin2 + 1.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior + 1);
					break;
				case 3:
					robin3 = pow(rr2 / z, 3);
					robin2 = 3.0 * pow(rr2 / z, 2) * (3.0 - pow(r / z, 2));
					robin1 = 3.0 * (rr2 / z) * (6.0 + (-3.0 + rr2 / pow(z, 2)) * pow(r / z, 2));
					A.a[offset] = robin3 / 4.0;
					A.a[offset + 1] = -7.0 * robin3 / 6.0 - robin2 / 6.0;
					A.a[offset + 2] = 2.0 * robin3 + 2.0 * robin2 / 3.0 + robin1 / 12.0;
					A.a[offset + 3] = -3.0 * robin3 / 2.0 - 5.0 * robin2 / 6.0 - robin1 / 3.0;
					A.a[offset + 4] = 1.0 + 5.0 * robin3 / 12.0 + robin2 / 3.0 + robin1 / 4.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
					break;
				}
				// Also fill RHS term.
				f[IDX(i, NzInterior + 1)] = uInf;
			}
		}

		// At this point we have now filled:
		offset = 4 + 2 * NzInterior + 5 * NrInterior * NzInterior + (2 + n_robin) * NrInterior;

		// Lower-right corner: equatorial symmetry.
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Robin boundary.
		r = (double)NrInterior + 0.5;
#pragma omp parallel shared(A, f) private(offset, z, rr2, robin1, robin2, robin3)
		{
#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{

				z = (double)j - 0.5;
				rr2 = r * r + z * z;
				// Each iteration of the loop fills n_robin elements.
				offset = t_offset + n_robin * (j - 1);

				A.ia[IDX(NrInterior + 1, j)] = BASE + offset;
				switch (robin)
				{
				case 1:
					A.a[offset] = 0.5 * rr2 / r;
					A.a[offset + 1] = -2.0 * rr2 / r;
					A.a[offset + 2] = 1.0 + 1.5 * rr2 / r;
					A.ja[offset] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior + 1, j);
					break;
				case 2:
					robin1 = 0.5 * rr2 * rr2 / (r * r);
					robin2 = 0.25 * (4.0 * rr2 / r - z * z / sqrt(rr2));
					A.a[offset] = -robin1;
					A.a[offset + 1] = 4.0 * robin1 + robin2;
					A.a[offset + 2] = -5.0 * robin1 - 4.0 * robin2;
					A.a[offset + 3] = 2.0 * robin1 + 3.0 * robin2 + 1.0;
					A.ja[offset] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior + 1, j);
					break;
				case 3:
					robin3 = pow(rr2 / r, 3);
					robin2 = 3.0 * pow(rr2 / r, 2) * (3.0 - pow(z / r, 2));
					robin1 = 3.0 * (rr2 / r) * (6.0 + (-3.0 + rr2 / pow(r, 2)) * pow(z / r, 2));
					A.a[offset] = robin3 / 4.0;
					A.a[offset + 1] = -7.0 * robin3 / 6.0 - robin2 / 6.0;
					A.a[offset + 2] = 2.0 * robin3 + 2.0 * robin2 / 3.0 + robin1 / 12.0;
					A.a[offset + 3] = -3.0 * robin3 / 2.0 - 5.0 * robin2 / 6.0 - robin1 / 3.0;
					A.a[offset + 4] = 1.0 + 5.0 * robin3 / 12.0 + robin2 / 3.0 + robin1 / 4.0;
					A.ja[offset] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior + 1, j);
					break;
				}
				// Also fill RHS term.
				f[IDX(NrInterior + 1, j)] = uInf;
			}
		}

		// At this point, we have now filled:
		offset = 6 + (2 + n_robin) * (NzInterior + NrInterior) + 5 * NrInterior * NzInterior;

		// Upper-right corner: fill with Robin.
		r = (double)NrInterior + 0.5;
		z = (double)NzInterior + 0.5;
		rr2 = r * r + z * z;
		A.ia[IDX(NrInterior + 1, NzInterior + 1)] = BASE + offset;
		switch (robin)
		{
		case 1:
			A.a[offset] = 0.5 * sqrt(rr2 / 2.0);
			A.a[offset + 1] = -2.0 * sqrt(rr2 / 2.0);
			A.a[offset + 2] = 1.0 + 1.5 * sqrt(rr2 / 2.0);
			A.ja[offset] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 1] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 2] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
		case 2:
			robin1 = 0.25 * rr2;
			robin2 = 2.0 * sqrt(rr2 / 2.0);
			A.a[offset] = -robin1;
			A.a[offset + 1] = 4.0 * robin1 + robin2;
			A.a[offset + 2] = -5.0 * robin1 - 4.0 * robin2;
			A.a[offset + 3] = 2.0 * robin1 + 3.0 * robin2 + 1.0;
			A.ja[offset] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 2] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 3] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
		case 3:
			A.a[offset] = rr2 * sqrt(rr2 / 2.0) / 8.0;
			A.a[offset + 1] = -rr2 * (18.0 + 7.0 * sqrt(2.0 * rr2)) / 24.0;
			A.a[offset + 2] = (3.0 * sqrt(2.0 * rr2) + 12.0 * rr2 + 2.0 * rr2 * sqrt(2.0 * rr2)) / 4.0;
			A.a[offset + 3] = -3.0 * (8.0 * sqrt(2.0 * rr2) + 10.0 * rr2 + rr2 * sqrt(2.0 * rr2)) / 8.0;
			A.a[offset + 4] = 1.0 + (108.0 * sqrt(2.0 * rr2) + 72.0 * rr2 + 5.0 * rr2 * sqrt(2.0 * rr2)) / 48.0;
			A.ja[offset] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 3] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 4] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
		}
		offset += n_robin;
		// Also fill RHS term.
		f[IDX(NrInterior + 1, NzInterior + 1)] = uInf;

		// Assert fill-in.
		assert(nnz == offset);

		// Fill last element of row offsets.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}
	// FOURTH-ORDER FLAT LAPLACIAN.
	else if (order == 4)
	{
		// Set offset to zero.
		offset = 0;

		// Lower-left corner: diagonal symmetry.
		A.ia[IDX(0, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)(r_sym * z_sym);
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		f[IDX(0, 0)] = 0.0;
		offset += 2;

		// offset = 2.

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
#pragma omp parallel shared(A, f) private(offset)
		{
#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each j iteration fills 2 elements.
				offset = t_offset + 2 * (j - 1);
				A.ia[IDX(0, j)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)r_sym;
				A.ja[offset] = BASE + IDX(0, j);
				A.ja[offset + 1] = BASE + IDX(1, j);
				f[IDX(0, j)] = 0.0;
			}
		}

		// We have now filled:
		offset = 2 + 2 * NzInterior;

		// offset = 2 + NzInterior * 2.

		// Upper-left corner: also axis symmetry.
		A.ia[IDX(0, NzInterior + 1)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)r_sym;
		A.ja[offset] = BASE + IDX(0, NzInterior + 1);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior + 1);
		f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior.

		// Left-interior points with semi-one sided finite difference.
		i = 1;

		// Bottom boundary with equatorial symmetry.
		j = 0;
		A.ia[IDX(1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(1, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		f[IDX(1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2.

		// Semi-one sided bottom with extra ghost zone.
		//
		//    o
		//    |
		//    o
		//    |
		// o--x--o--o
		//    |
		//    o
		//
		// Thus:
		// 
		// u(i-1, j  ) : 16/12 - 8h/12r.
		// u(i  , j-1) : 16/12
		// u(i  , j  ) : s(i, j) * h**2 - 30/12 - 30/12
		// u(i  , j+1) : 16/12 + (z_sym)*(-1/12)
		// u(i  , j+2) : -1/12
		// u(i+1, j  ) : 16/12 + 8h/12r + (r_sym)*(-1/12 + 1h/12r)
		// u(i+2, j  ) : -1/12 - 1h/12r
		//
		j = 1;
		r = 0.5;
		ir = 2.0;

		A.ia[IDX(1, 1)] = BASE + offset;
		A.a[offset] = 4.0 * third - 2.0 * third * ir;
		A.a[offset + 1] = 4.0 * third;
		A.a[offset + 2] = dr * dz * s[IDX(1, 1)] - 5.0;
		A.a[offset + 3] = 4.0 * third + (double)z_sym * (-twelfth);
		A.a[offset + 4] = -twelfth;
		A.a[offset + 5] = 4.0 * third + 2.0 * third * ir + (double)r_sym * (-twelfth + twelfth * ir);
		A.a[offset + 6] = -twelfth - twelfth * ir;
		A.ja[offset] = BASE + IDX(0, 1);
		A.ja[offset + 1] = BASE + IDX(1, 0);
		A.ja[offset + 2] = BASE + IDX(1, 1);
		A.ja[offset + 3] = BASE + IDX(1, 2);
		A.ja[offset + 4] = BASE + IDX(1, 3);
		A.ja[offset + 5] = BASE + IDX(2, 1);
		A.ja[offset + 6] = BASE + IDX(3, 1);
		f[IDX(1, 1)] *= dr * dz;
		offset += 7;

		// offset = (2 + 2) +  2 * NzInterior 
		// + (2 + 7).

		// Set temporary offset.
		t_offset = offset;

		// Combination of semi-onesided and centered.
		//
		//    o
		//    |
		//    o
		//    |
		// o--x--o--o
		//    |
		//    o
		//    |
		//    o
		//
		// Thus:
		//
		// u(i-1, j  ) : 16/12 - 8h/12r
		// u(i  , j-2) : -1/12
		// u(i  , j-1) : 16/12
		// u(i  , j  ) : s(i, j) * h**2 - 30/12 - 30/12
		// u(i  , j+1) : 16/12
		// u(i  , j+2) : -1/12
		// u(i+1, j  ) : 16/12 + 8h/12r + (r_sym)*(-1/12 + 1h/12r)
		// u(i+2, j  ) : -1/12 - 1h/12r
		//
#pragma omp parallel shared(A, f) private(offset)
		{
#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Each iteration fills 8 elements.
				offset = t_offset + 8 * (j - 2);
				A.ia[IDX(1, j)] = BASE + offset;
				A.a[offset] = 4.0 * third - 2.0 * third * ir;
				A.a[offset + 1] = -twelfth;
				A.a[offset + 2] = 4.0 * third;
				A.a[offset + 3] = dr * dz * s[IDX(1, j)] - 5.0;
				A.a[offset + 4] = 4.0 * third;
				A.a[offset + 5] = -twelfth;
				A.a[offset + 6] = 4.0 * third + 2.0 * third * ir + (double)r_sym * (-twelfth + twelfth * ir);
				A.a[offset + 7] = -twelfth - twelfth * ir;;
				A.ja[offset] = BASE + IDX(0, j);
				A.ja[offset + 1] = BASE + IDX(1, j - 2);
				A.ja[offset + 2] = BASE + IDX(1, j - 1);
				A.ja[offset + 3] = BASE + IDX(1, j);
				A.ja[offset + 4] = BASE + IDX(1, j + 1);
				A.ja[offset + 5] = BASE + IDX(1, j + 2);
				A.ja[offset + 6] = BASE + IDX(2, j);
				A.ja[offset + 7] = BASE + IDX(3, j);
				f[IDX(1, j)] *= dr * dz;
			}
		}

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2).


		// Semi-one sided top.
		//
		//    o
		//    |
		// o--x--o--o
		//    |
		//    o
		//    |
		//    o
		//    |
		//    o
		//    |
		//    o
		//
		// Recall:
		// - Semi-one sided on z.
		// u''(x) = (10u(i + 1) - 15u(i) - 4u(i - 1) + 14u(i - 2) - 6u(i - 3) + u(i - 4)/12h**2
		// 
		// Thus:
		//
		// u(i-1, j  ) : 16/12 - 8h/12r
		// u(i  , j-4) : 1/12
		// u(i  , j-3) : -6/12
		// u(i  , j-2) : 14/12
		// u(i  , j-1) : -4/12
		// u(i  , j  ) : s(i, j) * h **2  - 30/12 - 15/12
		// u(i  , j+1) : 10/12
		// u(i+1, j  ) : 16/12 + 8h/12r + (r_sym)*(-1/12 + 1h/12r)
		// u(i+2, j  ) : -1/12 - 1h/12r
		//
		j = NzInterior;
		A.ia[IDX(1, NzInterior)] = BASE + offset;
		A.a[offset] = 4.0 * third - 2.0 * third * ir;
		A.a[offset + 1] = twelfth;
		A.a[offset + 2] = -0.5;
		A.a[offset + 3] = 7.0 * sixth;
		A.a[offset + 4] = -third;
		A.a[offset + 5] = dr * dz * s[IDX(i, j)] - 3.75;
		A.a[offset + 6] = 5.0 * sixth;
		A.a[offset + 7] = 4.0 * third + 2.0 * third * ir + (double)r_sym * (-twelfth + twelfth * ir);
		A.a[offset + 8] = -twelfth - twelfth * ir;
		A.ja[offset] = BASE + IDX(0, NzInterior);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior - 4);
		A.ja[offset + 2] = BASE + IDX(1, NzInterior - 3);
		A.ja[offset + 3] = BASE + IDX(1, NzInterior - 2);
		A.ja[offset + 4] = BASE + IDX(1, NzInterior - 1);
		A.ja[offset + 5] = BASE + IDX(1, NzInterior);
		A.ja[offset + 6] = BASE + IDX(1, NzInterior + 1);
		A.ja[offset + 7] = BASE + IDX(2, NzInterior);
		A.ja[offset + 8] = BASE + IDX(3, NzInterior);
		f[IDX(1, NzInterior)] *= dr * dz;
		offset += 9;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9.

		j = NzInterior + 1;
		z = (double)j - 0.5;
		rr2 = r * r + z * z;
		A.ia[IDX(1, NzInterior + 1)] = BASE + offset;
		switch (robin)
		{
		case 1:
			robin1 = rr2 / z;
			A.a[offset] = robin1 / 4.0;
			A.a[offset + 1] = -4.0 * robin1 / 3.0;
			A.a[offset + 2] = 3.0 * robin1;
			A.a[offset + 3] = -4.0 * robin1;
			A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
			break;
		case 2:
			robin2 = pow(rr2 / z, 2);
			robin1 = (rr2 / z) * (4.0 - pow(r / z, 2));
			A.a[offset] = -5.0 * robin2 / 12.0;
			A.a[offset + 1] = 61.0 * robin2 / 24.0 + robin1 / 8.0;
			A.a[offset + 2] = -13.0 * robin2 / 2.0 - 2.0 * robin1 / 3.0;
			A.a[offset + 3] = 107.0 * robin2 / 12.0 + 3.0 * robin1 / 2.0;
			A.a[offset + 4] = -77.0 * robin2 / 12.0 - 2.0 * robin1;
			A.a[offset + 5] = 1.0 + 15.0 * robin2 / 8.0 + 25.0 * robin1 / 24.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior + 1);
			break;
		case 3:
			robin3 = pow(rr2 / z, 3);
			robin2 = 3.0 * pow(rr2 / z, 2) * (3.0 - pow(r / z, 2));
			robin1 = 3.0 * (rr2 / z) * (6.0 + (-3.0 + rr2 / pow(z, 2)) * pow(r / z, 2));
			A.a[offset] = 5.0 * robin3 / 16.0;
			A.a[offset + 1] = -13.0 * robin3 / 6.0 - 5.0 * robin2 / 36.0;
			A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
			A.a[offset + 3] = -31.0 * robin3 / 3.0 - 13.0 * robin2 / 6.0 - 2.0 * robin1 / 9.0;
			A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + robin1 / 2.0;
			A.a[offset + 5] = -29.0 * robin3 / 6 - 77.0 * robin2 / 36.0 - 2.0 * robin1 / 3.0;
			A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 5.0 * robin2 / 8.0 + 25.0 * robin1 / 72.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 5);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior);
			A.ja[offset + 6] = BASE + IDX(i, NzInterior + 1);
			break;
#ifdef R4
		case 4:
			robin4 = pow(rr2 / z, 4);
			robin3 = 2.0 * pow(rr2 / z, 3) * (8.0 - 3.0 * pow(r / z, 2));
			robin2 = 3.0 * pow(rr2 / z, 2) * (24.0 - 12.0 * pow(r / z, 2) + 5.0 * pow(r / z, 4));
			robin1 = 3.0 * (rr2 / z) * (32.0 + 12.0 * (-2.0 + rr2 / pow(z, 2)) * pow(r / z, 2) - 5.0 * (rr2 / pow(z, 2)) * pow(r / z, 4));
			A.a[offset] = -7.0 * robin4 / 48.0;
			A.a[offset + 1] = 41.0 * robin4 / 36.0 + 5.0 * robin3 / 64.0;
			A.a[offset + 2] = -185.0 * robin4 / 48.0 - 13.0 * robin3 / 24.0 - 5.0 * robin2 / 144.0;
			A.a[offset + 3] = 22.0 * robin4 / 3.0 + 307.0 * robin3 / 192.0 + 61.0 * robin2 / 288.0 + robin1 / 96.0;
			A.a[offset + 4] = -1219.0 * robin4 / 144.0 - 31.0 * robin3 / 12.0 - 13.0 * robin2 / 24.0 - robin1 / 18.0;
			A.a[offset + 5] = 71.0 * robin4 / 12.0 + 461.0 * robin3 / 192.0 + 107.0 * robin2 / 144.0 + robin1 / 8.0;
			A.a[offset + 6] = -37.0 * robin4 / 16.0 - 29.0 * robin3 / 24.0 - 77.0 * robin2 / 144.0 - robin1 / 6.0;
			A.a[offset + 7] = 1.0 + 7.0 * robin4 / 18.0 + 49.0 * robin3 / 192.0 + 5.0 * robin2 / 32.0 + 25.0 * robin1 / 288.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 6);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 5);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 6] = BASE + IDX(i, NzInterior);
			A.ja[offset + 7] = BASE + IDX(i, NzInterior + 1);
			break;
#endif
#ifdef R7
		case 7:
			robin7 = pow(rr2 / z, 7);
			robin6 = 7.0 * pow(rr2 / z, 6) * (7.0 - 3.0 * pow(r / z, 2));
			robin5 = 42.0 * pow(rr2 / z, 5) * (21.0 - 15.0 * pow(r / z, 2) + 5.0 * pow(r / z, 4));
			robin4 = 210.0 * pow(rr2 / z, 4) * (35.0 - 30.0 * pow(r / z, 2) + 17.0 * pow(r / z, 4) - 6.0 * pow(r / z, 6));
			robin3 = 105.0 * pow(rr2 / z, 3) * (280.0 - 252.0 * pow(r / z, 2) + 45.0 * (rr2 / pow(z, 2)) * pow(r / z, 6) + 3.0 * (4.0 * (rr2 / pow(z, 2)) * pow(r / z, 2) + 35.0 * pow(r / z, 4)) + (51.0 * (rr2 / pow(z, 2)) * pow(r / z, 4) - 196.0 * pow(r / z, 6)));
			robin2 = 315.0 * rr2 * (168.0 - 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(z, 4) + pow(rr2, 3) / pow(z, 6)) * pow(r / z, 2) + 2.0 * (rr2 / pow(z, 2)) * (175.0 - 49.0 * rr2 / pow(z, 2) - 18.0 * pow(rr2, 2) / pow(z, 4)) * pow(r / z, 4) + 3.0 * (pow(rr2, 2) / pow(z, 4)) * (49.0 - 11.0 * rr2 / pow(z, 2)) * pow(r / z, 6));
			robin1 = 315.0 * (rr2 / z) * (112.0 + 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(z, 4) + pow(rr2, 3) / pow(z, 6)) * pow(r / z, 2) + 2.0 * (rr2 / pow(z, 2)) * (-175.0 + 49.0 * rr2 / pow(z, 2) + 18.0 * pow(rr2, 2) / pow(z, 4)) * pow(r / z, 4) + 3.0 * (pow(rr2, 2) / pow(z, 4)) * (-49.0 + 11.0 * rr2 / pow(z, 2)) * pow(r / z, 6));
			A.a[offset] = robin7 / 384.0;
			A.a[offset + 1] = -119.0 * robin7 / 4320.0 - robin6 / 560.0;
			A.a[offset + 2] = 757.0 * robin7 / 5760.0 + 347.0 * robin6 / 20160.0 + robin5 / 864.0;
			A.a[offset + 3] = -1877.0 * robin7 / 5040.0 - 373.0 * robin6 / 5040.0 - 61.0 * robin5 / 6048.0 - robin4 / 1440.0;
			A.a[offset + 4] = 1999.0 * robin7 / 2880.0 + 313.0 * robin6 / 1680.0 + 13.0 * robin5 / 336.0 + 41.0 * robin4 / 7560.0 + robin3 / 2688.0;
			A.a[offset + 5] = -8.0 * robin7 / 9.0 - 305.0 * robin6 / 1008.0 - 2581.0 * robin5 / 30240.0 - 37.0 * robin4 / 2016.0 - 13.0 * robin3 / 5040.0 - robin2 / 6048.0;
			A.a[offset + 6] = 2281.0 * robin7 / 2880.0 + 3313.0 * robin6 / 10080.0 + 179.0 * robin5 / 1512.0 + 11.0 * robin4 / 315.0 + 307.0 * robin3 / 40320.0 + 61.0 * robin2 / 60480.0 + robin1 / 20160.0;
			A.a[offset + 7] = -349.0 * robin7 / 720.0 - 401.0 * robin6 / 1680.0 - 71.0 * robin5 / 672.0 - 1219.0 * robin4 / 30240.0 - 31.0 * robin3 / 2520.0 - 13.0 * robin2 / 5040.0 - robin1 / 3780.0;
			A.a[offset + 8] = 1123.0 * robin7 / 5760.0 + 563.0 * robin6 / 5040.0 + 179.0 * robin5 / 3024.0 + 71.0 * robin4 / 2520.0 + 461.0 * robin3 / 40320.0 + 107.0 * robin2 / 30240.0 + robin1 / 1680.0;
			A.a[offset + 9] = -67.0 * robin7 / 1440.0 - 11.0 * robin6 / 360.0 - 115.0 * robin5 / 6048.0 - 37.0 * robin4 / 3360.0 - 29.0 * robin3 / 5040.0 - 11.0 * robin2 / 4320.0 - robin1 / 1260.0;
			A.a[offset + 10] = 1.0 + 121.0 * robin7 / 24192.0 + 5.0 * robin6 / 1344.0 + 3.0 * robin5 / 1120.0 + robin4 / 540.0 + 7.0 * robin3 / 5760.0 + robin2 / 1344.0 + 5.0 * robin1 / 12096.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 9);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 8);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 7);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 6);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior - 5);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 6] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 7] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 8] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 9] = BASE + IDX(i, NzInterior);
			A.ja[offset + 10] = BASE + IDX(i, NzInterior + 1);
			break;
#endif
		}
		f[IDX(1, NzInterior + 1)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin.

		// Set temporary offset.
		t_offset = offset;

#pragma omp parallel shared(A, f) private(j, offset, r, z, ir, rr2,\
		robin1, robin2, robin3, robin4, robin5, robin6, robin7)
		{
#pragma omp for schedule(guided)
			for (i = 2; i < NrInterior; i++)
			{
				// Each iteration fills 2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin points.
				offset = t_offset + (i - 2) * (9 * NzInterior + 2 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;

				// Inverse.
				ir = 1.0 / r;

				// Equatorial symmetry.
				j = 0;
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Semi-onesided and centered difference.
				//
				//       o
				//       |
				//       o
				//       |
				// o--o--x--o--o
				//       |
				//       o
				//
				// Thus:
				// 
				// u(i-2, j  ) : -1/12 + 1h/12r
				// u(i-1, j  ) : 16/12 - 8h/12r
				// u(i  , j-1) : 16/12
				// u(i  , j  ) : s(i,j) * h**2 -30/12 - 30/12
				// u(i  , j+1) : 16/12 + (z_sym)*(-1/12)
				// u(i  , j+2) : -1/12
				// u(i+1, j  ) : 16/12 + 8h/12r
				// u(i+2, j  ) : -1/12 - 1h/12r
				//
				j = 1;
				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = -twelfth + twelfth * ir;
				A.a[offset + 1] = 4.0 * third - 2.0 * third * ir;
				A.a[offset + 2] = 4.0 * third;
				A.a[offset + 3] = dr * dz * s[IDX(i, j)] - 5.0;
				A.a[offset + 4] = 4.0 * third + (double)z_sym * (-twelfth);
				A.a[offset + 5] = -twelfth;
				A.a[offset + 6] = 4.0 * third + 2.0 * third * ir;
				A.a[offset + 7] = -twelfth - twelfth * ir;
				A.ja[offset] = BASE + IDX(i - 2, 1);
				A.ja[offset + 1] = BASE + IDX(i - 1, 1);
				A.ja[offset + 2] = BASE + IDX(i, 0);
				A.ja[offset + 3] = BASE + IDX(i, 1);
				A.ja[offset + 4] = BASE + IDX(i, 2);
				A.ja[offset + 5] = BASE + IDX(i, 3);
				A.ja[offset + 6] = BASE + IDX(i + 1, 1);
				A.ja[offset + 7] = BASE + IDX(i + 2, 1);
				f[IDX(i, j)] *= dr * dz;
				offset += 8;

				// Centered interior points.
				//
				//       o
				//       |
				//       o
				//       |
				// o--o--x--o--o
				//       |
				//       o
				//       |
				//       o
				//
				// Recall:
				// - Centered on r:
				// u' (x) = (-u(i+2)+8u(i+1)-8u(i-1)+u(i-2))/12h
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				// - Centered on z:
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				//
				// Thus:
				//
				// u(i-2, j  ) : -1/12 + 1h/12r
				// u(i-1, j  ) : 16/12 - 8h/12r
				// u(i  , j-2) : -1/12
				// u(i  , j-1) : 16/12
				// u(i  , j  ) : s(i,j) * h**2 -30/12 - 30/12
				// u(i  , j+1) : 16/12
				// u(i  , j+2) : -1/12
				// u(i+1, j  ) : 16/12 + 8h/12r
				// u(i+2, j  ) : -1/12 - 1h/12r
				//
				for (j = 2; j < NzInterior; j++)
				{
					A.ia[IDX(i, j)] = BASE + offset;
					A.a[offset] = -twelfth + twelfth * ir;
					A.a[offset + 1] = 4.0 * third - 2.0 * third * ir;
					A.a[offset + 2] = -twelfth;
					A.a[offset + 3] = 4.0 * third;
					A.a[offset + 4] = dr * dz * s[IDX(i, j)] - 5.0;
					A.a[offset + 5] = 4.0 * third;
					A.a[offset + 6] = -twelfth;
					A.a[offset + 7] = 4.0 * third + 2.0 * third * ir;
					A.a[offset + 8] = -twelfth - twelfth * ir;
					A.ja[offset] = BASE + IDX(i - 2, j);
					A.ja[offset + 1] = BASE + IDX(i - 1, j);
					A.ja[offset + 2] = BASE + IDX(i, j - 2);
					A.ja[offset + 3] = BASE + IDX(i, j - 1);
					A.ja[offset + 4] = BASE + IDX(i, j);
					A.ja[offset + 5] = BASE + IDX(i, j + 1);
					A.ja[offset + 6] = BASE + IDX(i, j + 2);
					A.ja[offset + 7] = BASE + IDX(i + 1, j);
					A.ja[offset + 8] = BASE + IDX(i + 2, j);
					f[IDX(i, j)] *= dr * dz;
					offset += 9;
				}

				// Semi-onesided and centered difference.
				// 
				//       o
				//       |
				// o--o--x--o--o
				//       |
				//       o
				//       |
				//       o
				//       |
				//       o
				//       |
				//       o
				//
				// - Centered on r:
				// u' (x) = (-u(i+2)+8u(i+1)-8u(i-1)+u(i-2))/12h
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				// - Semi-onesided on z:
				// u''(x) = (10u(i + 1) - 15u(i) - 4u(i - 1) + 14u(i - 2) - 6u(i - 3) + u(i - 4)/12h**2
				//
				// Thus:
				//
				// u(i-2, j  ) : -1/12 + 1h/12r
				// u(i-1, j  ) : 16/12 - 8h/12r
				// u(i  , j-4) : 1/12
				// u(i  , j-3) : -6/12
				// u(i  , j-2) : 14/12
				// u(i  , j-1) : -4/12
				// u(i  , j  ) : s(i,j) * h**2 -30/12 -15/12
				// u(i  , j+1) : 10/12
				// u(i+1, j  ) : 16/12 + 8h/12r
				// u(i+2, j  ) : -1/12 - 1h/12r
				//
				j = NzInterior;
				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = -twelfth + twelfth * ir;
				A.a[offset + 1] = 4.0 * third - 2.0 * third * ir;
				A.a[offset + 2] = twelfth;
				A.a[offset + 3] = -0.5;
				A.a[offset + 4] = 7.0 * sixth;
				A.a[offset + 5] = -third;
				A.a[offset + 6] = dr * dz * s[IDX(i, j)] - 3.75;
				A.a[offset + 7] = 5.0 * sixth;
				A.a[offset + 8] = 4.0 * third + 2.0 * third * ir;
				A.a[offset + 9] = -twelfth - twelfth * ir;
				A.ja[offset] = BASE + IDX(i - 2, j);
				A.ja[offset + 1] = BASE + IDX(i - 1, j);
				A.ja[offset + 2] = BASE + IDX(i, j - 4);
				A.ja[offset + 3] = BASE + IDX(i, j - 3);
				A.ja[offset + 4] = BASE + IDX(i, j - 2);
				A.ja[offset + 5] = BASE + IDX(i, j - 1);
				A.ja[offset + 6] = BASE + IDX(i, j);
				A.ja[offset + 7] = BASE + IDX(i, j + 1);
				A.ja[offset + 8] = BASE + IDX(i + 1, j);
				A.ja[offset + 9] = BASE + IDX(i + 2, j);
				f[IDX(i, j)] *= dr * dz;
				offset += 10;

				j = NzInterior + 1;
				z = (double)j - 0.5;
				rr2 = r * r + z * z;
				A.ia[IDX(i, j)] = BASE + offset;
				switch (robin)
				{
				case 1:
					robin1 = rr2 / z;
					A.a[offset] = robin1 / 4.0;
					A.a[offset + 1] = -4.0 * robin1 / 3.0;
					A.a[offset + 2] = 3.0 * robin1;
					A.a[offset + 3] = -4.0 * robin1;
					A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
					break;
				case 2:
					robin2 = pow(rr2 / z, 2);
					robin1 = (rr2 / z) * (4.0 - pow(r / z, 2));
					A.a[offset] = -5.0 * robin2 / 12.0;
					A.a[offset + 1] = 61.0 * robin2 / 24.0 + robin1 / 8.0;
					A.a[offset + 2] = -13.0 * robin2 / 2.0 - 2.0 * robin1 / 3.0;
					A.a[offset + 3] = 107.0 * robin2 / 12.0 + 3.0 * robin1 / 2.0;
					A.a[offset + 4] = -77.0 * robin2 / 12.0 - 2.0 * robin1;
					A.a[offset + 5] = 1.0 + 15.0 * robin2 / 8.0 + 25.0 * robin1 / 24.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 4);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior);
					A.ja[offset + 5] = BASE + IDX(i, NzInterior + 1);
					break;
				case 3:
					robin3 = pow(rr2 / z, 3);
					robin2 = 3.0 * pow(rr2 / z, 2) * (3.0 - pow(r / z, 2));
					robin1 = 3.0 * (rr2 / z) * (6.0 + (-3.0 + rr2 / pow(z, 2)) * pow(r / z, 2));
					A.a[offset] = 5.0 * robin3 / 16.0;
					A.a[offset + 1] = -13.0 * robin3 / 6.0 - 5.0 * robin2 / 36.0;
					A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
					A.a[offset + 3] = -31.0 * robin3 / 3.0 - 13.0 * robin2 / 6.0 - 2.0 * robin1 / 9.0;
					A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + robin1 / 2.0;
					A.a[offset + 5] = -29.0 * robin3 / 6 - 77.0 * robin2 / 36.0 - 2.0 * robin1 / 3.0;
					A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 5.0 * robin2 / 8.0 + 25.0 * robin1 / 72.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 5);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 4);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 5] = BASE + IDX(i, NzInterior);
					A.ja[offset + 6] = BASE + IDX(i, NzInterior + 1);
					break;
#ifdef R4
				case 4:
					robin4 = pow(rr2 / z, 4);
					robin3 = 2.0 * pow(rr2 / z, 3) * (8.0 - 3.0 * pow(r / z, 2));
					robin2 = 3.0 * pow(rr2 / z, 2) * (24.0 - 12.0 * pow(r / z, 2) + 5.0 * pow(r / z, 4));
					robin1 = 3.0 * (rr2 / z) * (32.0 + 12.0 * (-2.0 + rr2 / pow(z, 2)) * pow(r / z, 2) - 5.0 * (rr2 / pow(z, 2)) * pow(r / z, 4));
					A.a[offset] = -7.0 * robin4 / 48.0;
					A.a[offset + 1] = 41.0 * robin4 / 36.0 + 5.0 * robin3 / 64.0;
					A.a[offset + 2] = -185.0 * robin4 / 48.0 - 13.0 * robin3 / 24.0 - 5.0 * robin2 / 144.0;
					A.a[offset + 3] = 22.0 * robin4 / 3.0 + 307.0 * robin3 / 192.0 + 61.0 * robin2 / 288.0 + robin1 / 96.0;
					A.a[offset + 4] = -1219.0 * robin4 / 144.0 - 31.0 * robin3 / 12.0 - 13.0 * robin2 / 24.0 - robin1 / 18.0;
					A.a[offset + 5] = 71.0 * robin4 / 12.0 + 461.0 * robin3 / 192.0 + 107.0 * robin2 / 144.0 + robin1 / 8.0;
					A.a[offset + 6] = -37.0 * robin4 / 16.0 - 29.0 * robin3 / 24.0 - 77.0 * robin2 / 144.0 - robin1 / 6.0;
					A.a[offset + 7] = 1.0 + 7.0 * robin4 / 18.0 + 49.0 * robin3 / 192.0 + 5.0 * robin2 / 32.0 + 25.0 * robin1 / 288.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 6);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 5);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 4);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 5] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 6] = BASE + IDX(i, NzInterior);
					A.ja[offset + 7] = BASE + IDX(i, NzInterior + 1);
					break;
#endif
#ifdef R7
				case 7:
					robin7 = pow(rr2 / z, 7);
					robin6 = 7.0 * pow(rr2 / z, 6) * (7.0 - 3.0 * pow(r / z, 2));
					robin5 = 42.0 * pow(rr2 / z, 5) * (21.0 - 15.0 * pow(r / z, 2) + 5.0 * pow(r / z, 4));
					robin4 = 210.0 * pow(rr2 / z, 4) * (35.0 - 30.0 * pow(r / z, 2) + 17.0 * pow(r / z, 4) - 6.0 * pow(r / z, 6));
					robin3 = 105.0 * pow(rr2 / z, 3) * (280.0 - 252.0 * pow(r / z, 2) + 45.0 * (rr2 / pow(z, 2)) * pow(r / z, 6) + 3.0 * (4.0 * (rr2 / pow(z, 2)) * pow(r / z, 2) + 35.0 * pow(r / z, 4)) + (51.0 * (rr2 / pow(z, 2)) * pow(r / z, 4) - 196.0 * pow(r / z, 6)));
					robin2 = 315.0 * rr2 * (168.0 - 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(z, 4) + pow(rr2, 3) / pow(z, 6)) * pow(r / z, 2) + 2.0 * (rr2 / pow(z, 2)) * (175.0 - 49.0 * rr2 / pow(z, 2) - 18.0 * pow(rr2, 2) / pow(z, 4)) * pow(r / z, 4) + 3.0 * (pow(rr2, 2) / pow(z, 4)) * (49.0 - 11.0 * rr2 / pow(z, 2)) * pow(r / z, 6));
					robin1 = 315.0 * (rr2 / z) * (112.0 + 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(z, 4) + pow(rr2, 3) / pow(z, 6)) * pow(r / z, 2) + 2.0 * (rr2 / pow(z, 2)) * (-175.0 + 49.0 * rr2 / pow(z, 2) + 18.0 * pow(rr2, 2) / pow(z, 4)) * pow(r / z, 4) + 3.0 * (pow(rr2, 2) / pow(z, 4)) * (-49.0 + 11.0 * rr2 / pow(z, 2)) * pow(r / z, 6));
					A.a[offset] = robin7 / 384.0;
					A.a[offset + 1] = -119.0 * robin7 / 4320.0 - robin6 / 560.0;
					A.a[offset + 2] = 757.0 * robin7 / 5760.0 + 347.0 * robin6 / 20160.0 + robin5 / 864.0;
					A.a[offset + 3] = -1877.0 * robin7 / 5040.0 - 373.0 * robin6 / 5040.0 - 61.0 * robin5 / 6048.0 - robin4 / 1440.0;
					A.a[offset + 4] = 1999.0 * robin7 / 2880.0 + 313.0 * robin6 / 1680.0 + 13.0 * robin5 / 336.0 + 41.0 * robin4 / 7560.0 + robin3 / 2688.0;
					A.a[offset + 5] = -8.0 * robin7 / 9.0 - 305.0 * robin6 / 1008.0 - 2581.0 * robin5 / 30240.0 - 37.0 * robin4 / 2016.0 - 13.0 * robin3 / 5040.0 - robin2 / 6048.0;
					A.a[offset + 6] = 2281.0 * robin7 / 2880.0 + 3313.0 * robin6 / 10080.0 + 179.0 * robin5 / 1512.0 + 11.0 * robin4 / 315.0 + 307.0 * robin3 / 40320.0 + 61.0 * robin2 / 60480.0 + robin1 / 20160.0;
					A.a[offset + 7] = -349.0 * robin7 / 720.0 - 401.0 * robin6 / 1680.0 - 71.0 * robin5 / 672.0 - 1219.0 * robin4 / 30240.0 - 31.0 * robin3 / 2520.0 - 13.0 * robin2 / 5040.0 - robin1 / 3780.0;
					A.a[offset + 8] = 1123.0 * robin7 / 5760.0 + 563.0 * robin6 / 5040.0 + 179.0 * robin5 / 3024.0 + 71.0 * robin4 / 2520.0 + 461.0 * robin3 / 40320.0 + 107.0 * robin2 / 30240.0 + robin1 / 1680.0;
					A.a[offset + 9] = -67.0 * robin7 / 1440.0 - 11.0 * robin6 / 360.0 - 115.0 * robin5 / 6048.0 - 37.0 * robin4 / 3360.0 - 29.0 * robin3 / 5040.0 - 11.0 * robin2 / 4320.0 - robin1 / 1260.0;
					A.a[offset + 10] = 1.0 + 121.0 * robin7 / 24192.0 + 5.0 * robin6 / 1344.0 + 3.0 * robin5 / 1120.0 + robin4 / 540.0 + 7.0 * robin3 / 5760.0 + robin2 / 1344.0 + 5.0 * robin1 / 12096.0;
					A.ja[offset] = BASE + IDX(i, NzInterior - 9);
					A.ja[offset + 1] = BASE + IDX(i, NzInterior - 8);
					A.ja[offset + 2] = BASE + IDX(i, NzInterior - 7);
					A.ja[offset + 3] = BASE + IDX(i, NzInterior - 6);
					A.ja[offset + 4] = BASE + IDX(i, NzInterior - 5);
					A.ja[offset + 5] = BASE + IDX(i, NzInterior - 4);
					A.ja[offset + 6] = BASE + IDX(i, NzInterior - 3);
					A.ja[offset + 7] = BASE + IDX(i, NzInterior - 2);
					A.ja[offset + 8] = BASE + IDX(i, NzInterior - 1);
					A.ja[offset + 9] = BASE + IDX(i, NzInterior);
					A.ja[offset + 10] = BASE + IDX(i, NzInterior + 1);
					break;
#endif
				}
				f[IDX(i, j)] = uInf;
			}
		}

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
			+ (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).

		// Interior bottom-right corner with equatorial symmetry.
		i = NrInterior;
		r = (double)i - 0.5;
		ir = 1.0 / r;

		j = 0;
		A.ia[IDX(NrInterior, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior, 1);
		f[IDX(NrInterior, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin)
		// + 2.

		// Semi-onesided and centered difference.
		//
		//             o
		//             |
		//             o
		//             |
		// o--o--o--o--x--o
		//             |
		//             o
		//
		// Recall:
		// - Semi-onesided on r:
		// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
		// u' (x) = (3u(i+1)+10u(i)-18u(i-1)+6u(i-2)-u(i-3))/12h
		// - Semi-onesided on z:
		// u''(x) = (-u(i+3)+4u(i+2)+6u(i+1)-20u(i)-11u(i-1))/12h**2.
		// 
		// Thus:
		// 
		// u(i-4, j  ) : 1/12
		// u(i-3, j  ) : -6/12 - 1h/12r
		// u(i-2, j  ) : 14/12 + 6h/12r
		// u(i-1, j  ) : -4/12 - 18h/12r
		// u(i  , j-1) : 16/12
		// u(i  , j  ) : s(i,j) * h**2 -15/12 - 30/12 +10h/12r
		// u(i  , j+1) : 16/12 + (z_sym)*(-1/12)
		// u(i  , j+2) : -1/12
		// u(i+1, j  ) : 10/12 + 3h/12r
		//
		j = 1;
		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = twelfth;
		A.a[offset + 1] = -0.5 - twelfth * ir;
		A.a[offset + 2] = 7.0 * sixth + 0.5 * ir;
		A.a[offset + 3] = -third - 1.5 * ir;
		A.a[offset + 4] = 4.0 * third;
		A.a[offset + 5] = dr * dz * s[IDX(i, j)] - 3.75 + 5.0 * sixth * ir;
		A.a[offset + 6] = 4.0 * third + (double)z_sym * (-twelfth);
		A.a[offset + 7] = -twelfth;
		A.a[offset + 8] = 5.0 * sixth + 0.25 * ir;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j);
		A.ja[offset + 2] = BASE + IDX(i - 2, j);
		A.ja[offset + 3] = BASE + IDX(i - 1, j);
		A.ja[offset + 4] = BASE + IDX(i, j - 1);
		A.ja[offset + 5] = BASE + IDX(i, j);
		A.ja[offset + 6] = BASE + IDX(i, j + 1);
		A.ja[offset + 7] = BASE + IDX(i, j + 2);
		A.ja[offset + 8] = BASE + IDX(i + 1, j);
		f[IDX(i, j)] *= dr * dz;
		offset += 9;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9.

		// Set temporary offset.
		t_offset = offset;

#pragma omp parallel shared(A, f) private(offset)
		{
#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Eac iteration fills 10 elements.
				offset = t_offset + (j - 2) * 10;
				// Semi-onesided and centered difference.
				//
				//             o
				//             |
				//             o
				//             |
				// o--o--o--o--x--o
				//             |
				//             o
				//             |
				//             o
				//
				// Recall:
				// - Semi-onesided on r:
				// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
				// u' (x) = (3u(i+1)+10u(i)-18u(i-1)+6u(i-2)-u(i-3))/12h
				// - Centered on z:
				// u''(x) = (-(u(i+2)+u(i-2))+16(u(i+1)+u(i-1))-30u(i))/12h**2
				// 
				// Thus:
				// 
				// u(i-4, j  ) : 1/12
				// u(i-3, j  ) : -6/12 - 1h/12r
				// u(i-2, j  ) : 14/12 + 6h/12r
				// u(i-1, j  ) : -4/12 - 18h/12r
				// u(i  , j-2) : -1/12
				// u(i  , j-1) : 16/12
				// u(i  , j  ) : s(i,j) * h**2 -15/12 + 10h/12r - 30/12
				// u(i  , j+1) : 16/12
				// u(i  , j+2) : -1/12
				// u(i+1, j  ) : 10/12 + 3h/12r
				//
				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = twelfth;
				A.a[offset + 1] = -0.5 - twelfth * ir;
				A.a[offset + 2] = 7.0 * sixth + 0.5 * ir;
				A.a[offset + 3] = -third - 1.5 * ir;
				A.a[offset + 4] = -twelfth;
				A.a[offset + 5] = 4.0 * third;
				A.a[offset + 6] = dr * dz * s[IDX(i, j)] - 3.75 + 5.0 * sixth * ir;
				A.a[offset + 7] = 4.0 * third;
				A.a[offset + 8] = -twelfth;
				A.a[offset + 9] = 5.0 * sixth + 0.25 * ir;
				A.ja[offset] = BASE + IDX(i - 4, j);
				A.ja[offset + 1] = BASE + IDX(i - 3, j);
				A.ja[offset + 2] = BASE + IDX(i - 2, j);
				A.ja[offset + 3] = BASE + IDX(i - 1, j);
				A.ja[offset + 4] = BASE + IDX(i, j - 2);
				A.ja[offset + 5] = BASE + IDX(i, j - 1);
				A.ja[offset + 6] = BASE + IDX(i, j);
				A.ja[offset + 7] = BASE + IDX(i, j + 1);
				A.ja[offset + 8] = BASE + IDX(i, j + 2);
				A.ja[offset + 9] = BASE + IDX(i + 1, j);
				f[IDX(i, j)] *= dr * dz;
			}
		}

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
			+ (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin)
			+ 2 + 9 + 10 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2)

		// Fill interior top-right corner.
		//
		//             o
		//             |
		// o--o--o--o--x--o
		//             |
		//             o
		//             |
		//             o
		//             |
		//             o
		//             |
		//             o
		//
		// Recall:
		// - Semi-onesided on r:
		// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
		// u' (x) = (3u(i+1)+10u(i)-18u(i-1)+6u(i-2)-u(i-3))/12h
		// - Semi-onesided on z:
		// u''(x) = (10u(i+1)-15u(i)-4u(i-1)+14u(i-2)-6u(i-3)+u(i-4)/12h**2
		//
		// Thus:
		//
		// u(i-4, j  ) : 1/12
		// u(i-3, j  ) : -6/12 - 1h/12r
		// u(i-2, j  ) : 14/12 + 6h/12r
		// u(i-1, j  ) : -4/12 - 18h/12r
		// u(i  , j-4) : 1/12
		// u(i  , j-3) : -6/12
		// u(i  , j-2) : 14/12
		// u(i  , j-1) : -4/12
		// u(i  , j  ) : s(i,j) * h**2 - 30/12 + 10h/12r
		// u(i  , j+1) : 10/12
		// u(i+1, j  ) : 10/12 + 3h/12r
		//
		j = NzInterior;
		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = twelfth;
		A.a[offset + 1] = -0.5 - twelfth * ir;
		A.a[offset + 2] = 7.0 * sixth + 0.5 * ir;
		A.a[offset + 3] = -third - 1.5 * ir;
		A.a[offset + 4] = twelfth;
		A.a[offset + 5] = -0.5;
		A.a[offset + 6] = 7.0 * sixth;
		A.a[offset + 7] = -third;
		A.a[offset + 8] = dr * dz * s[IDX(i, j)] - 2.5 + 5.0 * sixth * ir;
		A.a[offset + 9] = 5.0 * sixth;
		A.a[offset + 10] = 5.0 * sixth + 0.25 * ir;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j);
		A.ja[offset + 2] = BASE + IDX(i - 2, j);
		A.ja[offset + 3] = BASE + IDX(i - 1, j);
		A.ja[offset + 4] = BASE + IDX(i, j - 4);
		A.ja[offset + 5] = BASE + IDX(i, j - 3);
		A.ja[offset + 6] = BASE + IDX(i, j - 2);
		A.ja[offset + 7] = BASE + IDX(i, j - 1);
		A.ja[offset + 8] = BASE + IDX(i, j);
		A.ja[offset + 9] = BASE + IDX(i, j + 1);
		A.ja[offset + 10] = BASE + IDX(i + 1, j);
		f[IDX(i, j)] *= dr * dz;
		offset += 11;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11

		j = NzInterior + 1;
		z = (double)j - 0.5;
		rr2 = r * r + z * z;
		A.ia[IDX(i, j)] = BASE + offset;
		switch (robin)
		{
		case 1:
			robin1 = rr2 / z;
			A.a[offset] = robin1 / 4.0;
			A.a[offset + 1] = -4.0 * robin1 / 3.0;
			A.a[offset + 2] = 3.0 * robin1;
			A.a[offset + 3] = -4.0 * robin1;
			A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior + 1);
			break;
		case 2:
			robin2 = pow(rr2 / z, 2);
			robin1 = (rr2 / z) * (4.0 - pow(r / z, 2));
			A.a[offset] = -5.0 * robin2 / 12.0;
			A.a[offset + 1] = 61.0 * robin2 / 24.0 + robin1 / 8.0;
			A.a[offset + 2] = -13.0 * robin2 / 2.0 - 2.0 * robin1 / 3.0;
			A.a[offset + 3] = 107.0 * robin2 / 12.0 + 3.0 * robin1 / 2.0;
			A.a[offset + 4] = -77.0 * robin2 / 12.0 - 2.0 * robin1;
			A.a[offset + 5] = 1.0 + 15.0 * robin2 / 8.0 + 25.0 * robin1 / 24.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior + 1);
			break;
		case 3:
			robin3 = pow(rr2 / z, 3);
			robin2 = 3.0 * pow(rr2 / z, 2) * (3.0 - pow(r / z, 2));
			robin1 = 3.0 * (rr2 / z) * (6.0 + (-3.0 + rr2 / pow(z, 2)) * pow(r / z, 2));
			A.a[offset] = 5.0 * robin3 / 16.0;
			A.a[offset + 1] = -13.0 * robin3 / 6.0 - 5.0 * robin2 / 36.0;
			A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
			A.a[offset + 3] = -31.0 * robin3 / 3.0 - 13.0 * robin2 / 6.0 - 2.0 * robin1 / 9.0;
			A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + robin1 / 2.0;
			A.a[offset + 5] = -29.0 * robin3 / 6 - 77.0 * robin2 / 36.0 - 2.0 * robin1 / 3.0;
			A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 5.0 * robin2 / 8.0 + 25.0 * robin1 / 72.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 5);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior);
			A.ja[offset + 6] = BASE + IDX(i, NzInterior + 1);
			break;
#ifdef R4
		case 4:
			robin4 = pow(rr2 / z, 4);
			robin3 = 2.0 * pow(rr2 / z, 3) * (8.0 - 3.0 * pow(r / z, 2));
			robin2 = 3.0 * pow(rr2 / z, 2) * (24.0 - 12.0 * pow(r / z, 2) + 5.0 * pow(r / z, 4));
			robin1 = 3.0 * (rr2 / z) * (32.0 + 12.0 * (-2.0 + rr2 / pow(z, 2)) * pow(r / z, 2) - 5.0 * (rr2 / pow(z, 2)) * pow(r / z, 4));
			A.a[offset] = -7.0 * robin4 / 48.0;
			A.a[offset + 1] = 41.0 * robin4 / 36.0 + 5.0 * robin3 / 64.0;
			A.a[offset + 2] = -185.0 * robin4 / 48.0 - 13.0 * robin3 / 24.0 - 5.0 * robin2 / 144.0;
			A.a[offset + 3] = 22.0 * robin4 / 3.0 + 307.0 * robin3 / 192.0 + 61.0 * robin2 / 288.0 + robin1 / 96.0;
			A.a[offset + 4] = -1219.0 * robin4 / 144.0 - 31.0 * robin3 / 12.0 - 13.0 * robin2 / 24.0 - robin1 / 18.0;
			A.a[offset + 5] = 71.0 * robin4 / 12.0 + 461.0 * robin3 / 192.0 + 107.0 * robin2 / 144.0 + robin1 / 8.0;
			A.a[offset + 6] = -37.0 * robin4 / 16.0 - 29.0 * robin3 / 24.0 - 77.0 * robin2 / 144.0 - robin1 / 6.0;
			A.a[offset + 7] = 1.0 + 7.0 * robin4 / 18.0 + 49.0 * robin3 / 192.0 + 5.0 * robin2 / 32.0 + 25.0 * robin1 / 288.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 6);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 5);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 6] = BASE + IDX(i, NzInterior);
			A.ja[offset + 7] = BASE + IDX(i, NzInterior + 1);
			break;
#endif
#ifdef R7
		case 7:
			robin7 = pow(rr2 / z, 7);
			robin6 = 7.0 * pow(rr2 / z, 6) * (7.0 - 3.0 * pow(r / z, 2));
			robin5 = 42.0 * pow(rr2 / z, 5) * (21.0 - 15.0 * pow(r / z, 2) + 5.0 * pow(r / z, 4));
			robin4 = 210.0 * pow(rr2 / z, 4) * (35.0 - 30.0 * pow(r / z, 2) + 17.0 * pow(r / z, 4) - 6.0 * pow(r / z, 6));
			robin3 = 105.0 * pow(rr2 / z, 3) * (280.0 - 252.0 * pow(r / z, 2) + 45.0 * (rr2 / pow(z, 2)) * pow(r / z, 6) + 3.0 * (4.0 * (rr2 / pow(z, 2)) * pow(r / z, 2) + 35.0 * pow(r / z, 4)) + (51.0 * (rr2 / pow(z, 2)) * pow(r / z, 4) - 196.0 * pow(r / z, 6)));
			robin2 = 315.0 * rr2 * (168.0 - 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(z, 4) + pow(rr2, 3) / pow(z, 6)) * pow(r / z, 2) + 2.0 * (rr2 / pow(z, 2)) * (175.0 - 49.0 * rr2 / pow(z, 2) - 18.0 * pow(rr2, 2) / pow(z, 4)) * pow(r / z, 4) + 3.0 * (pow(rr2, 2) / pow(z, 4)) * (49.0 - 11.0 * rr2 / pow(z, 2)) * pow(r / z, 6));
			robin1 = 315.0 * (rr2 / z) * (112.0 + 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(z, 4) + pow(rr2, 3) / pow(z, 6)) * pow(r / z, 2) + 2.0 * (rr2 / pow(z, 2)) * (-175.0 + 49.0 * rr2 / pow(z, 2) + 18.0 * pow(rr2, 2) / pow(z, 4)) * pow(r / z, 4) + 3.0 * (pow(rr2, 2) / pow(z, 4)) * (-49.0 + 11.0 * rr2 / pow(z, 2)) * pow(r / z, 6));
			A.a[offset] = robin7 / 384.0;
			A.a[offset + 1] = -119.0 * robin7 / 4320.0 - robin6 / 560.0;
			A.a[offset + 2] = 757.0 * robin7 / 5760.0 + 347.0 * robin6 / 20160.0 + robin5 / 864.0;
			A.a[offset + 3] = -1877.0 * robin7 / 5040.0 - 373.0 * robin6 / 5040.0 - 61.0 * robin5 / 6048.0 - robin4 / 1440.0;
			A.a[offset + 4] = 1999.0 * robin7 / 2880.0 + 313.0 * robin6 / 1680.0 + 13.0 * robin5 / 336.0 + 41.0 * robin4 / 7560.0 + robin3 / 2688.0;
			A.a[offset + 5] = -8.0 * robin7 / 9.0 - 305.0 * robin6 / 1008.0 - 2581.0 * robin5 / 30240.0 - 37.0 * robin4 / 2016.0 - 13.0 * robin3 / 5040.0 - robin2 / 6048.0;
			A.a[offset + 6] = 2281.0 * robin7 / 2880.0 + 3313.0 * robin6 / 10080.0 + 179.0 * robin5 / 1512.0 + 11.0 * robin4 / 315.0 + 307.0 * robin3 / 40320.0 + 61.0 * robin2 / 60480.0 + robin1 / 20160.0;
			A.a[offset + 7] = -349.0 * robin7 / 720.0 - 401.0 * robin6 / 1680.0 - 71.0 * robin5 / 672.0 - 1219.0 * robin4 / 30240.0 - 31.0 * robin3 / 2520.0 - 13.0 * robin2 / 5040.0 - robin1 / 3780.0;
			A.a[offset + 8] = 1123.0 * robin7 / 5760.0 + 563.0 * robin6 / 5040.0 + 179.0 * robin5 / 3024.0 + 71.0 * robin4 / 2520.0 + 461.0 * robin3 / 40320.0 + 107.0 * robin2 / 30240.0 + robin1 / 1680.0;
			A.a[offset + 9] = -67.0 * robin7 / 1440.0 - 11.0 * robin6 / 360.0 - 115.0 * robin5 / 6048.0 - 37.0 * robin4 / 3360.0 - 29.0 * robin3 / 5040.0 - 11.0 * robin2 / 4320.0 - robin1 / 1260.0;
			A.a[offset + 10] = 1.0 + 121.0 * robin7 / 24192.0 + 5.0 * robin6 / 1344.0 + 3.0 * robin5 / 1120.0 + robin4 / 540.0 + 7.0 * robin3 / 5760.0 + robin2 / 1344.0 + 5.0 * robin1 / 12096.0;
			A.ja[offset] = BASE + IDX(i, NzInterior - 9);
			A.ja[offset + 1] = BASE + IDX(i, NzInterior - 8);
			A.ja[offset + 2] = BASE + IDX(i, NzInterior - 7);
			A.ja[offset + 3] = BASE + IDX(i, NzInterior - 6);
			A.ja[offset + 4] = BASE + IDX(i, NzInterior - 5);
			A.ja[offset + 5] = BASE + IDX(i, NzInterior - 4);
			A.ja[offset + 6] = BASE + IDX(i, NzInterior - 3);
			A.ja[offset + 7] = BASE + IDX(i, NzInterior - 2);
			A.ja[offset + 8] = BASE + IDX(i, NzInterior - 1);
			A.ja[offset + 9] = BASE + IDX(i, NzInterior);
			A.ja[offset + 10] = BASE + IDX(i, NzInterior + 1);
			break;
#endif
		}
		f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin

		// Do bottom-right corner.
		i = NrInterior + 1;
		r = (double)i - 0.5;

		j = 0;
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin

		// Set temporary offset.
		t_offset = offset;

#pragma omp parallel shared(A, f) private(offset, z, rr2,\
		robin1, robin2, robin3, robin4, robin5, robin6, robin7)
		{
#pragma omp for schedule(guided)
			for (j = 1; j < NzInterior + 1; j++)
			{
				// Each iteration fills n_robin elements.
				offset = t_offset + n_robin * (j - 1);

				// Z coordinate.
				z = (double)j - 0.5;
				rr2 = r * r + z * z;
				A.ia[IDX(i, j)] = BASE + offset;
				switch (robin)
				{
				case 1:
					robin1 = rr2 / r;
					A.a[offset] = robin1 / 4.0;
					A.a[offset + 1] = -4.0 * robin1 / 3.0;
					A.a[offset + 2] = 3.0 * robin1;
					A.a[offset + 3] = -4.0 * robin1;
					A.a[offset + 4] = 1.0 + 25.0 * robin1 / 12.0;
					A.ja[offset] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior + 1, j);
					break;
				case 2:
					robin2 = pow(rr2 / r, 2);
					robin1 = (rr2 / r) * (4.0 - pow(z / r, 2));
					A.a[offset] = -5.0 * robin2 / 12.0;
					A.a[offset + 1] = 61.0 * robin2 / 24.0 + robin1 / 8.0;
					A.a[offset + 2] = -13.0 * robin2 / 2.0 - 2.0 * robin1 / 3.0;
					A.a[offset + 3] = 107.0 * robin2 / 12.0 + 3.0 * robin1 / 2.0;
					A.a[offset + 4] = -77.0 * robin2 / 12.0 - 2.0 * robin1;
					A.a[offset + 5] = 1.0 + 15.0 * robin2 / 8.0 + 25.0 * robin1 / 24.0;
					A.ja[offset] = BASE + IDX(NrInterior - 4, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior, j);
					A.ja[offset + 5] = BASE + IDX(NrInterior + 1, j);
					break;
				case 3:
					robin3 = pow(rr2 / r, 3);
					robin2 = 3.0 * pow(rr2 / r, 2) * (3.0 - pow(z / r, 2));
					robin1 = 3.0 * (rr2 / r) * (6.0 + (-3.0 + rr2 / pow(r, 2)) * pow(z / r, 2));
					A.a[offset] = 5.0 * robin3 / 16.0;
					A.a[offset + 1] = -13.0 * robin3 / 6.0 - 5.0 * robin2 / 36.0;
					A.a[offset + 2] = 307.0 * robin3 / 48.0 + 61.0 * robin2 / 72.0 + robin1 / 24.0;
					A.a[offset + 3] = -31.0 * robin3 / 3.0 - 13.0 * robin2 / 6.0 - 2.0 * robin1 / 9.0;
					A.a[offset + 4] = 461.0 * robin3 / 48.0 + 107.0 * robin2 / 36.0 + robin1 / 2.0;
					A.a[offset + 5] = -29.0 * robin3 / 6 - 77.0 * robin2 / 36.0 - 2.0 * robin1 / 3.0;
					A.a[offset + 6] = 1.0 + 49.0 * robin3 / 48.0 + 5.0 * robin2 / 8.0 + 25.0 * robin1 / 72.0;
					A.ja[offset] = BASE + IDX(NrInterior - 5, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 4, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 5] = BASE + IDX(NrInterior, j);
					A.ja[offset + 6] = BASE + IDX(NrInterior + 1, j);
					break;
#ifdef R4
				case 4:
					robin4 = pow(rr2 / r, 4);
					robin3 = 2.0 * pow(rr2 / r, 3) * (8.0 - 3.0 * pow(z / r, 2));
					robin2 = 3.0 * pow(rr2 / r, 2) * (24.0 - 12.0 * pow(z / r, 2) + 5.0 * pow(z / r, 4));
					robin1 = 3.0 * (rr2 / r) * (32.0 + 12.0 * (-2.0 + rr2 / pow(r, 2)) * pow(z / r, 2) - 5.0 * (rr2 / pow(r, 2)) * pow(z / r, 4));
					A.a[offset] = -7.0 * robin4 / 48.0;
					A.a[offset + 1] = 41.0 * robin4 / 36.0 + 5.0 * robin3 / 64.0;
					A.a[offset + 2] = -185.0 * robin4 / 48.0 - 13.0 * robin3 / 24.0 - 5.0 * robin2 / 144.0;
					A.a[offset + 3] = 22.0 * robin4 / 3.0 + 307.0 * robin3 / 192.0 + 61.0 * robin2 / 288.0 + robin1 / 96.0;
					A.a[offset + 4] = -1219.0 * robin4 / 144.0 - 31.0 * robin3 / 12.0 - 13.0 * robin2 / 24.0 - robin1 / 18.0;
					A.a[offset + 5] = 71.0 * robin4 / 12.0 + 461.0 * robin3 / 192.0 + 107.0 * robin2 / 144.0 + robin1 / 8.0;
					A.a[offset + 6] = -37.0 * robin4 / 16.0 - 29.0 * robin3 / 24.0 - 77.0 * robin2 / 144.0 - robin1 / 6.0;
					A.a[offset + 7] = 1.0 + 7.0 * robin4 / 18.0 + 49.0 * robin3 / 192.0 + 5.0 * robin2 / 32.0 + 25.0 * robin1 / 288.0;
					A.ja[offset] = BASE + IDX(NrInterior - 6, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 5, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 4, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 5] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 6] = BASE + IDX(NrInterior, j);
					A.ja[offset + 7] = BASE + IDX(NrInterior + 1, j);
					break;
#endif
#ifdef R7
				case 7:
					robin7 = pow(rr2 / r, 7);
					robin6 = 7.0 * pow(rr2 / r, 6) * (7.0 - 3.0 * pow(z / r, 2));
					robin5 = 42.0 * pow(rr2 / r, 5) * (21.0 - 15.0 * pow(z / r, 2) + 5.0 * pow(z / r, 4));
					robin4 = 210.0 * pow(rr2 / r, 4) * (35.0 - 30.0 * pow(z / r, 2) + 17.0 * pow(z / r, 4) - 6.0 * pow(z / r, 6));
					robin3 = 105.0 * pow(rr2 / r, 3) * (280.0 - 252.0 * pow(z / r, 2) + 45.0 * (rr2 / pow(r, 2)) * pow(z / r, 6) + 3.0 * (4.0 * (rr2 / pow(r, 2)) * pow(z / r, 2) + 35.0 * pow(z / r, 4)) + (51.0 * (rr2 / pow(r, 2)) * pow(z / r, 4) - 196.0 * pow(z / r, 6)));
					robin2 = 315.0 * rr2 * (168.0 - 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(r, 4) + pow(rr2, 3) / pow(r, 6)) * pow(z / r, 2) + 2.0 * (rr2 / pow(r, 2)) * (175.0 - 49.0 * rr2 / pow(r, 2) - 18.0 * pow(rr2, 2) / pow(r, 4)) * pow(z / r, 4) + 3.0 * (pow(rr2, 2) / pow(r, 4)) * (49.0 - 11.0 * rr2 / pow(r, 2)) * pow(z / r, 6));
					robin1 = 315.0 * (rr2 / r) * (112.0 + 8.0 * (-21.0 + 14.0 * pow(rr2, 2) / pow(r, 4) + pow(rr2, 3) / pow(r, 6)) * pow(z / r, 2) + 2.0 * (rr2 / pow(r, 2)) * (-175.0 + 49.0 * rr2 / pow(r, 2) + 18.0 * pow(rr2, 2) / pow(r, 4)) * pow(z / r, 4) + 3.0 * (pow(rr2, 2) / pow(r, 4)) * (-49.0 + 11.0 * rr2 / pow(r, 2)) * pow(z / r, 6));
					A.a[offset] = robin7 / 384.0;
					A.a[offset + 1] = -119.0 * robin7 / 4320.0 - robin6 / 560.0;
					A.a[offset + 2] = 757.0 * robin7 / 5760.0 + 347.0 * robin6 / 20160.0 + robin5 / 864.0;
					A.a[offset + 3] = -1877.0 * robin7 / 5040.0 - 373.0 * robin6 / 5040.0 - 61.0 * robin5 / 6048.0 - robin4 / 1440.0;
					A.a[offset + 4] = 1999.0 * robin7 / 2880.0 + 313.0 * robin6 / 1680.0 + 13.0 * robin5 / 336.0 + 41.0 * robin4 / 7560.0 + robin3 / 2688.0;
					A.a[offset + 5] = -8.0 * robin7 / 9.0 - 305.0 * robin6 / 1008.0 - 2581.0 * robin5 / 30240.0 - 37.0 * robin4 / 2016.0 - 13.0 * robin3 / 5040.0 - robin2 / 6048.0;
					A.a[offset + 6] = 2281.0 * robin7 / 2880.0 + 3313.0 * robin6 / 10080.0 + 179.0 * robin5 / 1512.0 + 11.0 * robin4 / 315.0 + 307.0 * robin3 / 40320.0 + 61.0 * robin2 / 60480.0 + robin1 / 20160.0;
					A.a[offset + 7] = -349.0 * robin7 / 720.0 - 401.0 * robin6 / 1680.0 - 71.0 * robin5 / 672.0 - 1219.0 * robin4 / 30240.0 - 31.0 * robin3 / 2520.0 - 13.0 * robin2 / 5040.0 - robin1 / 3780.0;
					A.a[offset + 8] = 1123.0 * robin7 / 5760.0 + 563.0 * robin6 / 5040.0 + 179.0 * robin5 / 3024.0 + 71.0 * robin4 / 2520.0 + 461.0 * robin3 / 40320.0 + 107.0 * robin2 / 30240.0 + robin1 / 1680.0;
					A.a[offset + 9] = -67.0 * robin7 / 1440.0 - 11.0 * robin6 / 360.0 - 115.0 * robin5 / 6048.0 - 37.0 * robin4 / 3360.0 - 29.0 * robin3 / 5040.0 - 11.0 * robin2 / 4320.0 - robin1 / 1260.0;
					A.a[offset + 10] = 1.0 + 121.0 * robin7 / 24192.0 + 5.0 * robin6 / 1344.0 + 3.0 * robin5 / 1120.0 + robin4 / 540.0 + 7.0 * robin3 / 5760.0 + robin2 / 1344.0 + 5.0 * robin1 / 12096.0;
					A.ja[offset] = BASE + IDX(NrInterior - 9, j);
					A.ja[offset + 1] = BASE + IDX(NrInterior - 8, j);
					A.ja[offset + 2] = BASE + IDX(NrInterior - 7, j);
					A.ja[offset + 3] = BASE + IDX(NrInterior - 6, j);
					A.ja[offset + 4] = BASE + IDX(NrInterior - 5, j);
					A.ja[offset + 5] = BASE + IDX(NrInterior - 4, j);
					A.ja[offset + 6] = BASE + IDX(NrInterior - 3, j);
					A.ja[offset + 7] = BASE + IDX(NrInterior - 2, j);
					A.ja[offset + 8] = BASE + IDX(NrInterior - 1, j);
					A.ja[offset + 9] = BASE + IDX(NrInterior, j);
					A.ja[offset + 10] = BASE + IDX(NrInterior + 1, j);
					break;
#endif
				}
				f[IDX(i, j)] = uInf;
			}
		}

		offset = (2 + 2 + 2) + 2 * NzInterior
			+ 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
			+ (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin)
			+ 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin
			+ n_robin * NzInterior;

		// offset = (2 + 2 + 2) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin
		// + n_robin * NzInterior

		// Finally, fill top corner with Robin
		j = NzInterior + 1;
		z = (double)j - 0.5;
		rr2 = r * r + z * z;
		A.ia[IDX(i, j)] = BASE + offset;
		switch (robin)
		{
		case 1:
			A.a[offset] = sqrt(rr2 / 2.0) / 4.0;
			A.a[offset + 1] = -2.0 * sqrt(2.0 * rr2) / 3.0;
			A.a[offset + 2] = 3.0 * sqrt(rr2 / 2.0);
			A.a[offset + 3] = -2.0 * sqrt(2.0 * rr2);
			A.a[offset + 4] = 1.0 + 25.0 * sqrt(rr2 / 2.0) / 12.0;
			A.ja[offset] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 3] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 4] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
		case 2:
			A.a[offset] = -5.0 * rr2 / 24.0;
			A.a[offset + 1] = (12.0 * sqrt(2.0 * rr2) + 61.0 * rr2) / 48.0;
			A.a[offset + 2] = -4.0 * sqrt(2.0 * rr2) / 3.0 - 13.0 * rr2 / 4.0;
			A.a[offset + 3] = 3.0 * sqrt(2.0 * rr2) + 107.0 * rr2 / 24.0;
			A.a[offset + 4] = -4.0 * sqrt(2.0 * rr2) - 77.0 * rr2 / 24.0;
			A.a[offset + 5] = 1.0 + 5.0 * (20.0 * sqrt(2.0 * rr2) + 9.0 * rr2) / 48.0;
			A.ja[offset] = BASE + IDX(NrInterior - 4, NzInterior - 4);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 3] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 4] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 5] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
		case 3:
			A.a[offset] = 5.0 * rr2 * sqrt(rr2 / 2.0) / 32.0;
			A.a[offset + 1] = -rr2 * (15.0 + 13.0 * sqrt(2.0 * rr2)) / 24.0;
			A.a[offset + 2] = (72.0 * sqrt(2.0 * rr2) + 732.0 * rr2 + 307.0 * rr2 * sqrt(2.0 * rr2)) / 192.0;
			A.a[offset + 3] = -(24.0 * sqrt(2.0 * rr2) + 117.0 * rr2 + 31.0 * rr2 * sqrt(2.0 * rr2)) / 12.0;
			A.a[offset + 4] = 9.0 * sqrt(rr2 / 2.0) + 107.0 * rr2 / 8.0 + 461.0 * rr2 * sqrt(rr2 / 2.0) / 96.0;
			A.a[offset + 5] = -(144.0 * sqrt(2.0 * rr2) + 231.0 * rr2 + 29.0 * rr2 * sqrt(2.0 * rr2)) / 24.0;
			A.a[offset + 6] = 1.0 + (600.0 * sqrt(2.0 * rr2) + 540.0 * rr2 + 49.0 * rr2 * sqrt(2.0 * rr2)) / 192.0;
			A.ja[offset] = BASE + IDX(NrInterior - 5, NzInterior - 5);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 4, NzInterior - 4);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 3] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 4] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 5] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 6] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
#ifdef R4
		case 4:
			A.a[offset] = -7.0 * pow(rr2, 2) / 192.0;
			A.a[offset + 1] = rr2 * (45.0 * sqrt(2.0 * rr2) + 41.0 * rr2) / 144.0;
			A.a[offset + 2] = -rr2 * (240.0 + 416.0 * sqrt(2.0 * rr2) + 185.0 * rr2) / 192.0;
			A.a[offset + 3] = (24.0 * sqrt(2.0 * rr2) + 366.0 * rr2 + 307.0 * rr2 * sqrt(2.0 * rr2) + 88.0 * pow(rr2, 2)) / 48.0;
			A.a[offset + 4] = -8.0 * sqrt(2.0 * rr2) / 3.0 - 39.0 * rr2 / 2.0 - 31.0 * rr2 * sqrt(2.0 * rr2) / 3.0 - 1219.0 * pow(rr2, 2) / 576.0;
			A.a[offset + 5] = (288.0 * sqrt(2.0 * rr2) + 1284.0 * rr2 + 461.0 * rr2 * sqrt(2.0 * rr2) + 71.0 * pow(rr2, 2)) / 48.0;
			A.a[offset + 6] = -8.0 * sqrt(2.0 * rr2) - 77.0 * rr2 / 4.0 - 29.0 * rr2 * sqrt(rr2 / 2.0) / 3.0 - 37.0 * pow(rr2, 2) / 64.0;
			A.a[offset + 7] = 1.0 + (600.0 * sqrt(2.0 * rr2) + 810.0 * rr2 + 147.0 * rr2 * sqrt(2.0 * rr2) + 14.0 * pow(rr2, 2)) / 144.0;
			A.ja[offset] = BASE + IDX(NrInterior - 6, NzInterior - 6);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 5, NzInterior - 5);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 4, NzInterior - 4);
			A.ja[offset + 3] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 4] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 5] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 6] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 7] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
#endif
#ifdef R7
		case 7:
			A.a[offset] = pow(rr2, 3) * sqrt(rr2 / 2.0) / 3072.0;
			A.a[offset + 1] = -7.0 * pow(rr2, 3) * (108.0 + 17.0 * sqrt(2.0 * rr2)) / 69120.0;
			A.a[offset + 2] = pow(rr2, 2) * (11760.0 * sqrt(2.0 * rr2) + 9716.0 * rr2 + 757.0 * rr2 * sqrt(2.0 * rr2)) / 92160.0;
			A.a[offset + 3] = -pow(rr2, 2) * (102900.0 + 89670.0 * sqrt(2.0 * rr2) + 36554.0 * rr2 + 1877.0 * rr2 * sqrt(2.0 * rr2)) / 80640.0;
			A.a[offset + 4] = rr2 * (126000.0 * sqrt(2.0 * rr2) + 459200.0 * rr2 + 196560.0 * rr2 * sqrt(2.0 * rr2) + 52584.0 * pow(rr2, 2) + 1999.0 * pow(rr2, 2) * sqrt(2.0 * rr2)) / 46080.0;
			A.a[offset + 5] = -rr2 * (25200.0 + 109200.0 * sqrt(2.0 * rr2) + 194250.0 * rr2 + 54201.0 * rr2 * sqrt(2.0 * rr2) + 10675.0 * pow(rr2, 2) + 320.0 * pow(rr2, 2) * sqrt(2.0 * rr2)) / 5760.0;
			A.a[offset + 6] = (40320.0 * sqrt(2.0 * rr2) + 1229760.0 * rr2 + 2578800.0 * rr2 * sqrt(2.0 * rr2) + 2956800.0 * pow(rr2, 2) + 601440.0 * pow(rr2, 2) * sqrt(2.0 * rr2) + 92764.0 * pow(rr2, 3) + 2281.0 * pow(rr2, 3) * sqrt(2.0 * rr2)) / 46080.0;
			A.a[offset + 7] = -(53760.0 * sqrt(2.0 * rr2) + 786240.0 * rr2 + 1041600.0 * rr2 * sqrt(2.0 * rr2) + 853300.0 * pow(rr2, 2) + 134190.0 * pow(rr2, 2) * sqrt(2.0 * rr2) + 16842.0 * pow(rr2, 3) + 349.0 * pow(rr2, 3) * sqrt(2.0 * rr2)) / 11520.0;
			A.a[offset + 8] = (967680 * sqrt(2.0 * rr2) + 8628480.0 * rr2 + 7744800.0 * rr2 * sqrt(2.0 * rr2) + 4771200.0 * pow(rr2, 2) + 601440.0 * pow(rr2, 2) * sqrt(2.0 * rr2) + 63056.0 * pow(rr2, 3) + 1123.0 * pow(rr2, 3) * sqrt(2.0 * rr2)) / 92160.0;
			A.a[offset + 9] = -(322560.0 * sqrt(2.0 * rr2) + 1552320.0 * rr2 + 974400.0 * rr2 * sqrt(2.0 * rr2) + 466200.0 * pow(rr2, 2) + 48300.0 * pow(rr2, 2) * sqrt(2.0 * rr2) + 4312.0 * pow(rr2, 3) + 67.0 * pow(rr2, 3) * sqrt(2.0 * rr2)) / 23040.0;
			A.a[offset + 10] = 1.0 + 175.0 * sqrt(rr2 / 2.0) / 12.0 + 315.0 * rr2 / 16.0 + 1715.0 * rr2 * sqrt(rr2 / 2.0) / 96.0 + 245.0 * pow(rr2, 2) / 72.0 + 189.0 * pow(rr2, 2) * sqrt(rr2 / 2.0) / 320.0 + 35.0 * pow(rr2, 3) / 1536.0 + 121.0 * pow(rr2, 3) * sqrt(rr2 / 2.0) / 193536.0;
			A.ja[offset] = BASE + IDX(NrInterior - 9, NzInterior - 9);
			A.ja[offset + 1] = BASE + IDX(NrInterior - 8, NzInterior - 8);
			A.ja[offset + 2] = BASE + IDX(NrInterior - 7, NzInterior - 7);
			A.ja[offset + 3] = BASE + IDX(NrInterior - 6, NzInterior - 6);
			A.ja[offset + 4] = BASE + IDX(NrInterior - 5, NzInterior - 5);
			A.ja[offset + 5] = BASE + IDX(NrInterior - 4, NzInterior - 4);
			A.ja[offset + 6] = BASE + IDX(NrInterior - 3, NzInterior - 3);
			A.ja[offset + 7] = BASE + IDX(NrInterior - 2, NzInterior - 2);
			A.ja[offset + 8] = BASE + IDX(NrInterior - 1, NzInterior - 1);
			A.ja[offset + 9] = BASE + IDX(NrInterior, NzInterior);
			A.ja[offset + 10] = BASE + IDX(NrInterior + 1, NzInterior + 1);
			break;
#endif
		}
		f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2 + 2 + n_robin) + 2 * NzInterior 
		// + 2 + 7 + 8 * (NzInterior - 2) + 9 + n_robin
		// + (NrInterior - 2) * (2 + 8 + 9 * (NzInterior - 2) + 10 + n_robin).
		// + 2 + 9 + 10 * (NzInterior - 2) + 11 + n_robin
		// + n_robin * NzInterior

		// Fill last element.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}

#ifdef DEBUG
	csr_print(A, "A.asc", "iA.asc", "jA.asc");
#endif

	// All done.
	return;
}
