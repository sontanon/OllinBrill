#include "tools.h"
#include "param.h"

// One-based indexing.
#define BASE 1

// Print CSR matrix for debug.
#undef DEBUG

// Nonzero calculator.
int nnz_general_elliptic(const int NrInterior, const int NzInterior, const int order, const int robin)
{
	int nnz;

	int n_robin = robin + order;

	// Second order Laplacian.
	if (order == 2)
	{
		nnz = 9 * NrInterior * NzInterior
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	// Fourth-order Laplacian.
	else if (order == 4)
	{
		nnz = 17 * (NrInterior - 2) * (NzInterior - 2)
			+ (16 + 26) * (NrInterior + NzInterior - 4)
			+ (14 + 21 + 21 + 27)
			+ (2 + n_robin) * (NrInterior + NzInterior)
			+ (2 + 2 + 2 + n_robin);
	}
	return nnz;
}

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
	const int z_sym)
{
	// Number of Robin elements.
	int n_robin = robin + order;

	// Constant numbers.
	const double third = 1.0 / 3.0;
	const double sixth = 1.0 / 6.0;
	const double twelfth = 1.0 / 12.0;

	// Number of nonzero elements.
	int nnz = nnz_general_elliptic(NrInterior, NzInterior, order, robin);

	// Number of differing elements.

	// Number of elements we have filled in.
	int offset = 0;

	// Auxiliary variables.
	double r, z, rr2;
	int i, j, t_offset;
	double aux_a, aux_b, aux_c, aux_d, aux_e, aux_s;
	double robin1, robin2, robin3, robin4, robin5, robin6, robin7;

	// SECOND-ORDER GENERAL ELLIPTIC EQUATION.
	if (order == 2)
	{
		offset = 0;

		// Lower-left corner: diagonal symmetry.
		A.ia[IDX(0, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)(r_sym * z_sym);
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		ell_f[IDX(0, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
#pragma omp parallel shared(A, ell_f) private(offset)
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
				ell_f[IDX(0, j)] = 0.0;
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
		ell_f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Now come the interior points plus the top and bottom boundaries with
		// Robin and equatorial symmetry respectively.
#pragma omp parallel shared(A, ell_f) private(offset, j, r, z, rr2,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s,\
		robin1, robin2, robin3)
		{
#pragma omp for schedule(guided)
			for (i = 1; i < NrInterior + 1; i++)
			{
				// Each iteration of i loop will fill 9 * NzInterior + (2 + n_robin) values.
				offset = t_offset + (i - 1) * (9 * NzInterior + 2 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;

				// Do bottom boundary first with equatorial symmetry.
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				ell_f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Now loop over interior points.
				for (j = 1; j < NzInterior + 1; j++)
				{
					// Fetch values.
					aux_a = ell_a[IDX(i, j)];
					aux_b = 0.25 * ell_b[IDX(i, j)];
					aux_c = ell_c[IDX(i, j)];
					aux_d = 0.5 * dz * ell_d[IDX(i, j)];
					aux_e = 0.5 * dr * ell_e[IDX(i, j)];
					aux_s = dr * dz * ell_s[IDX(i, j)];
					// Row begins at offset.
					A.ia[IDX(i, j)] = BASE + offset;
					// Values.
					A.a[offset] = aux_b;
					A.a[offset + 1] = aux_a - aux_d;
					A.a[offset + 2] = -aux_b;
					A.a[offset + 3] = aux_c - aux_e;
					A.a[offset + 4] = aux_s - 2.0 * (aux_a + aux_c);
					A.a[offset + 5] = aux_c + aux_e;
					A.a[offset + 6] = -aux_b;
					A.a[offset + 7] = aux_a + aux_d;
					A.a[offset + 8] = aux_b;
					// Columns.
					A.ja[offset] = BASE + IDX(i - 1, j - 1);
					A.ja[offset + 1] = BASE + IDX(i - 1, j);
					A.ja[offset + 2] = BASE + IDX(i - 1, j + 1);
					A.ja[offset + 3] = BASE + IDX(i, j - 1);
					A.ja[offset + 4] = BASE + IDX(i, j);
					A.ja[offset + 5] = BASE + IDX(i, j + 1);
					A.ja[offset + 6] = BASE + IDX(i + 1, j - 1);
					A.ja[offset + 7] = BASE + IDX(i + 1, j);
					A.ja[offset + 8] = BASE + IDX(i + 1, j + 1);
					// Multiply RHS by dr * dz.
					ell_f[IDX(i, j)] *= dr * dz;
					// Increase offset.
					offset += 9;
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
				ell_f[IDX(i, NzInterior + 1)] = uInf;
			}
		}

		// At this point we have now filled:
		offset = 4 + 2 * NzInterior + 9 * NrInterior * NzInterior + (2 + n_robin) * NrInterior;

		// Lower-right corner: equatorial symmetry.
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		ell_f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// Set temporary offset.
		t_offset = offset;

		// Robin boundary.
		r = (double)NrInterior + 0.5;
#pragma omp parallel shared(A, ell_f) private(offset, z, rr2,\
		robin1, robin2, robin3)
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
				ell_f[IDX(NrInterior + 1, j)] = uInf;
			}
		}

		// At this point, we have now filled:
		offset = 6 + (2 + n_robin) * (NzInterior + NrInterior) + 9 * NrInterior * NzInterior;

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
		ell_f[IDX(NrInterior + 1, NzInterior + 1)] = uInf;

		// Assert fill-in.
		assert(nnz == offset);

		// Fill last element of row offsets.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}
	// FOURTH-ORDER GENERAL ELLIPTIC EQUATION
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
		ell_f[IDX(0, 0)] = 0.0;
		offset += 2;

		// offset = 2.

		// Set temporary offset.
		t_offset = offset;

		// Fill left-boundary using axis symmetry.
#pragma omp parallel shared(A, ell_f) private(offset)
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
				ell_f[IDX(0, j)] = 0.0;
			}
		}

#ifdef DEBUG
		printf("Stage 1.\n");
#endif

		// We have now filled:
		offset = 2 + 2 * NzInterior;

		// offset = 2 + NzInterior * 2.

		// Upper-left corner: also axis symmetry.
		A.ia[IDX(0, NzInterior + 1)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)r_sym;
		A.ja[offset] = BASE + IDX(0, NzInterior + 1);
		A.ja[offset + 1] = BASE + IDX(1, NzInterior + 1);
		ell_f[IDX(0, NzInterior + 1)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior.

		// Left-interior points with semi-one sided finite difference.
		i = 1;
		r = 0.5;

		// Bottom boundary with equatorial symmetry.
		j = 0;
		A.ia[IDX(1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(1, 0);
		A.ja[offset + 1] = BASE + IDX(1, 1);
		ell_f[IDX(1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior + 2.

		// Stencil for lower-left corner.
		//
		//     o   o   o
		//     |   | / 
		// o   o   o---o
		//   \ | / 
		// o---x---o---o
		//   / | \
		// o   o   o
		//
		// Thus:
		// 
		// u(i-1, j-1) : b/3
		// u(i-1, j  ) : 2(2a - d)/3
		// u(i-1, j+1) : -b/3
		// u(i  , j-1) : 2(2c - e)/3
		// u(i  , j  ) : -5(a + c)/2 + s
		// u(i  , j+1) : 2(2c + e)/3 + z_sym*((-c + e)/12)
		// u(i  , j+2) : -(c + e)/12
		// u(i+1, j-1) : -b/3
		// u(i+1, j  ) : 2(2a + d)/3 + r_sym*((-a + d)/12)
		// u(i+1, j+1) : b/3 + r_sym*z_sym*(-b/48)
		// u(i+1, j+2) : r_sym*(b/48)
		// u(i+2, j  ) : -(a + d)/12
		// u(i+2, j+1) : z_sym*(b/48)
		// u(i+2, j+2) : -b/48
		//
		j = 1;
		aux_a = ell_a[IDX(i, j)];
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)];
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(1, 1)] = BASE + offset;
		A.a[offset] = aux_b * third;
		A.a[offset + 1] = 2.0 * (2.0 * aux_a - aux_d) * third;
		A.a[offset + 2] = -aux_b * third;
		A.a[offset + 3] = 2.0 * (2.0 * aux_c - aux_e) * third;
		A.a[offset + 4] = -2.5 * (aux_a + aux_c) + aux_s;
		A.a[offset + 5] = 2.0 * (2.0 * aux_c + aux_e) * third + (double)z_sym * ((-aux_c + aux_e) * twelfth);
		A.a[offset + 6] = -(aux_c + aux_e) * twelfth;
		A.a[offset + 7] = -aux_b * third;
		A.a[offset + 8] = 2.0 * (2.0 * aux_a + aux_d) * third + (double)r_sym * ((-aux_a + aux_d) * twelfth);
		A.a[offset + 9] = aux_b * third + (double)(r_sym * z_sym) * (-aux_b / 48.0);
		A.a[offset + 10] = (double)r_sym * (aux_b / 48.0);
		A.a[offset + 11] = -(aux_a + aux_d) * twelfth;
		A.a[offset + 12] = (double)z_sym * (aux_b / 48.0);
		A.a[offset + 13] = -aux_b / 48.0;
		A.ja[offset] = BASE + IDX(0, 0);
		A.ja[offset + 1] = BASE + IDX(0, 1);
		A.ja[offset + 2] = BASE + IDX(0, 2);
		A.ja[offset + 3] = BASE + IDX(1, 0);
		A.ja[offset + 4] = BASE + IDX(1, 1);
		A.ja[offset + 5] = BASE + IDX(1, 2);
		A.ja[offset + 6] = BASE + IDX(1, 3);
		A.ja[offset + 7] = BASE + IDX(2, 0);
		A.ja[offset + 8] = BASE + IDX(2, 1);
		A.ja[offset + 9] = BASE + IDX(2, 2);
		A.ja[offset + 10] = BASE + IDX(2, 3);
		A.ja[offset + 11] = BASE + IDX(3, 1);
		A.ja[offset + 12] = BASE + IDX(3, 2);
		A.ja[offset + 13] = BASE + IDX(3, 3);
		ell_f[IDX(1, 1)] *= dr * dz;
		offset += 14;

		// offset = (2 + 2) +  2 * NzInterior + (2 + 14).

		// Set temporary offset.
		t_offset = offset;

		// Stencil for left strip.
		//
		//     o   o   o
		//     | /   / 
		// o   o   o    
		//   \ | / 
		// o---x---o---o
		//   / | \
		// o   o   o
		//     | \   \
		//     o   o   o
		//
		// Thus:
		// 
		// u(i-1, j-1) : b/3
		// u(i-1, j  ) : 2(2a - d)/3
		// u(i-1, j+1) : -b/3
		// u(i  , j-2) : (-c + e)/12
		// u(i  , j-1) : 2(2c - e)/3
		// u(i  , j  ) : -5(a + c)/2 + s
		// u(i  , j+1) : 2(2c + e)/3
		// u(i  , j+2) : -(c + e)/12
		// u(i+1, j-2) : r_sym*(-b/48)
		// u(i+1, j-1) : -b/3
		// u(i+1, j  ) : 2(2a + d)/3 + r_sym*((-a + d)/12)
		// u(i+1, j+1) : b/3
		// u(i+1, j+2) : r_sym*(b/48)
		// u(i+2, j-2) : b/48
		// u(i+2, j  ) : -(a + d)/12
		// u(i+2, j+2) : -b/48
		// 
#pragma omp parallel shared(A, ell_f) private(offset,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s)
		{
#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Each iterations fills 16 elements.
				offset = t_offset + 16 * (j - 2);
				aux_a = ell_a[IDX(i, j)];
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)];
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(1, j)] = BASE + offset;
				A.a[offset] = aux_b * third;
				A.a[offset + 1] = 2.0 * (2.0 * aux_a - aux_d) * third;
				A.a[offset + 2] = -aux_b * third;
				A.a[offset + 3] = (-aux_c + aux_e) * twelfth;
				A.a[offset + 4] = 2.0 * (2.0 * aux_c - aux_e) * third;
				A.a[offset + 5] = -2.5 * (aux_a + aux_c) + aux_s;
				A.a[offset + 6] = 2.0 * (2.0 * aux_c + aux_e) * third;
				A.a[offset + 7] = -(aux_c + aux_e) * twelfth;
				A.a[offset + 8] = (double)r_sym * (-aux_b / 48.0);
				A.a[offset + 9] = -aux_b * third;
				A.a[offset + 10] = 2.0 * (2.0 * aux_a + aux_d) * third + (double)r_sym * ((-aux_a + aux_d) * twelfth);
				A.a[offset + 11] = aux_b * third;
				A.a[offset + 12] = (double)r_sym * (aux_b / 48.0);
				A.a[offset + 13] = aux_b / 48.0;
				A.a[offset + 14] = -(aux_a + aux_d) * twelfth;
				A.a[offset + 15] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 1] = BASE + IDX(i - 1, j);
				A.ja[offset + 2] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 3] = BASE + IDX(i, j - 2);
				A.ja[offset + 4] = BASE + IDX(i, j - 1);
				A.ja[offset + 5] = BASE + IDX(i, j);
				A.ja[offset + 6] = BASE + IDX(i, j + 1);
				A.ja[offset + 7] = BASE + IDX(i, j + 2);
				A.ja[offset + 8] = BASE + IDX(i + 1, j - 2);
				A.ja[offset + 9] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 10] = BASE + IDX(i + 1, j);
				A.ja[offset + 11] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 12] = BASE + IDX(i + 1, j + 2);
				A.ja[offset + 13] = BASE + IDX(i + 2, j - 2);
				A.ja[offset + 14] = BASE + IDX(i + 2, j);
				A.ja[offset + 15] = BASE + IDX(i + 2, j + 2);
				ell_f[IDX(1, j)] *= dr * dz;
			}
		}

#ifdef DEBUG
		printf("Stage 2.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2).


		// Stencil for top left corner.
		//
		// o   o   o   o
		//   \ | /   /
		// o---x---o---o
		//   / | \   \
		// o   o   o   o
		//   / | \   \
		// o   o   o   o
		//   / | \   \
		// o   o   o   o
		//     |      
		//     o        
		//
		// Thus:
		// 
		// u(i-1, j-3) : b/18
		// u(i-1, j-2) : -b/3
		// u(i-1, j-1) : b
		// u(i-1, j  ) : (12a - 5b - 6d)/9
		// u(i-1, j+1) : -b/6
		// u(i  , j-4) : c/12
		// u(i  , j-3) : -(6c + e)/12
		// u(i  , j-2) : (7c + 3e)/6
		// u(i  , j-1) : -(2c + 9e)/6
		// u(i  , j  ) : -5a/2 - 5c/4 + 5e/6 + s
		// u(i  , j+1) : 5c/6 + e/4
		// u(i+1, j-3) : -b/18 + r_sym*(-b/144)
		// u(i+1, j-2) : b/3 + r_sym*(b/24)
		// u(i+1, j-1) : -b + r_sym*(-b/8)
		// u(i+1, j  ) : (12a + 5b + 6d)/9 + r_sym*(-6a + 5b + 6d)/72
		// u(i+1, j+1) : b/6 + r_sym*(b/48)
		// u(i+2, j-3) : b/144
		// u(i+2, j-2) : -b/24
		// u(i+2, j-1) : b/8
		// u(i+2, j  ) : -(6a + 5b + 6d)/72
		// u(i+2, j+1) : -b/48
		//
		j = NzInterior;
		aux_a = ell_a[IDX(i, j)];
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)];
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(1, j)] = BASE + offset;
		A.a[offset] = aux_b / 18.0;
		A.a[offset + 1] = -aux_b * third;
		A.a[offset + 2] = aux_b;
		A.a[offset + 3] = (12.0 * aux_a - 5.0 * aux_b - 6.0 * aux_d) / 9.0;
		A.a[offset + 4] = -aux_b * sixth;
		A.a[offset + 5] = aux_c * twelfth;
		A.a[offset + 6] = -(6.0 * aux_c + aux_e) * twelfth;
		A.a[offset + 7] = (7.0 * aux_c + 3.0 * aux_e) * sixth;
		A.a[offset + 8] = -(2.0 * aux_c + 9.0 * aux_e) * sixth;
		A.a[offset + 9] = -2.5 * aux_a - 1.25 * aux_c + 5.0 * sixth * aux_e + aux_s;
		A.a[offset + 10] = 5.0 * aux_c * sixth + aux_e * 0.25;
		A.a[offset + 11] = -aux_b / 18.0 + (double)r_sym * (-aux_b / 144.0);
		A.a[offset + 12] = aux_b * third + (double)r_sym * (aux_b / 24.0);
		A.a[offset + 13] = -aux_b + (double)r_sym * (-aux_b * 0.125);
		A.a[offset + 14] = (12.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 9.0 + (double)r_sym * ((-6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0);
		A.a[offset + 15] = aux_b * sixth + (double)r_sym * (aux_b / 48.0);
		A.a[offset + 16] = aux_b / 144.0;
		A.a[offset + 17] = -aux_b / 24.0;
		A.a[offset + 18] = aux_b * 0.125;
		A.a[offset + 19] = -(6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
		A.a[offset + 20] = -aux_b / 48.0;
		A.ja[offset] = BASE + IDX(i - 1, j - 3);
		A.ja[offset + 1] = BASE + IDX(i - 1, j - 2);
		A.ja[offset + 2] = BASE + IDX(i - 1, j - 1);
		A.ja[offset + 3] = BASE + IDX(i - 1, j);
		A.ja[offset + 4] = BASE + IDX(i - 1, j + 1);
		A.ja[offset + 5] = BASE + IDX(i, j - 4);
		A.ja[offset + 6] = BASE + IDX(i, j - 3);
		A.ja[offset + 7] = BASE + IDX(i, j - 2);
		A.ja[offset + 8] = BASE + IDX(i, j - 1);
		A.ja[offset + 9] = BASE + IDX(i, j);
		A.ja[offset + 10] = BASE + IDX(i, j + 1);
		A.ja[offset + 11] = BASE + IDX(i + 1, j - 3);
		A.ja[offset + 12] = BASE + IDX(i + 1, j - 2);
		A.ja[offset + 13] = BASE + IDX(i + 1, j - 1);
		A.ja[offset + 14] = BASE + IDX(i + 1, j);
		A.ja[offset + 15] = BASE + IDX(i + 1, j + 1);
		A.ja[offset + 16] = BASE + IDX(i + 2, j - 3);
		A.ja[offset + 17] = BASE + IDX(i + 2, j - 2);
		A.ja[offset + 18] = BASE + IDX(i + 2, j - 1);
		A.ja[offset + 19] = BASE + IDX(i + 2, j);
		A.ja[offset + 20] = BASE + IDX(i + 2, j + 1);
		ell_f[IDX(1, j)] *= dr * dz;
		offset += 21;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21.
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
		ell_f[IDX(1, NzInterior + 1)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin.

		// Set temporary offset.
		t_offset = offset;

#pragma omp parallel shared(A, ell_f) private(j, offset, r, z, rr2,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s,\
		robin1, robin2, robin3, robin4, robin5, robin6, robin7)
		{
#pragma omp for schedule(guided)
			for (i = 2; i < NrInterior; i++)
			{
				// Each iteration fills 2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin points.
				offset = t_offset + (i - 2) * (17 * NzInterior + 10 + n_robin);

				// R coordinate.
				r = (double)i - 0.5;

				// Equatorial symmetry.
				j = 0;
				A.ia[IDX(i, 0)] = BASE + offset;
				A.a[offset] = 1.0;
				A.a[offset + 1] = -(double)z_sym;
				A.ja[offset] = BASE + IDX(i, 0);
				A.ja[offset + 1] = BASE + IDX(i, 1);
				ell_f[IDX(i, 0)] = 0.0;
				offset += 2;

				// Stencil for lower strip.
				//
				// o       o       o
				//   \     |     /
				// o   o   o   o   o
				//   \   \ | /   /
				// o---o---x---o---o
				//       / | \
				//     o   o   o
				//
				// Thus:
				// 
				// u(i-2, j  ) : (-a + d)/12
				// u(i-2, j+1) : z_sym*(-b/48)
				// u(i-2, j+2) : b/48
				// u(i-1, j-1) : b/3
				// u(i-1, j  ) : 2(2a - d)/3
				// u(i-1, j+1) : -b/3
				// u(i  , j-1) : 2(2c - e)/3
				// u(i  , j  ) : -5(a + c)/2 + s
				// u(i  , j+1) : 2(2c + e)/3 + z_sym*((-c + e)/12)
				// u(i  , j+2) : -(c + e)/12
				// u(i+1, j-1) : -b/3
				// u(i+1, j  ) : 2(2a + d)/3
				// u(i+1, j+1) : b/3
				// u(i+2, j  ) : -(a + d)/12
				// u(i+2, j+1) : z_sym*(b/48)
				// u(i+2, j+2) : -b/48
				//
				j = 1;
				aux_a = ell_a[IDX(i, j)];
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)];
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = (-aux_a + aux_d) * twelfth;
				A.a[offset + 1] = (double)z_sym * (-aux_b / 48.0);
				A.a[offset + 2] = aux_b / 48.0;
				A.a[offset + 3] = aux_b * third;
				A.a[offset + 4] = 2.0 * (2.0 * aux_a - aux_d) * third;
				A.a[offset + 5] = -aux_b * third;
				A.a[offset + 6] = 2.0 * (2.0 * aux_c - aux_e) * third;
				A.a[offset + 7] = -2.5 * (aux_a + aux_c) + aux_s;
				A.a[offset + 8] = 2.0 * (2.0 * aux_c + aux_e) * third + (double)z_sym * ((-aux_c + aux_e) * twelfth);
				A.a[offset + 9] = -(aux_c + aux_e) * twelfth;
				A.a[offset + 10] = -aux_b * third;
				A.a[offset + 11] = 2.0 * (2.0 * aux_a + aux_d) * third;
				A.a[offset + 12] = aux_b * third;
				A.a[offset + 13] = -(aux_a + aux_d) * twelfth;
				A.a[offset + 14] = (double)z_sym * (aux_b / 48.0);
				A.a[offset + 15] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 2, j);
				A.ja[offset + 1] = BASE + IDX(i - 2, j + 1);
				A.ja[offset + 2] = BASE + IDX(i - 2, j + 2);
				A.ja[offset + 3] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 4] = BASE + IDX(i - 1, j);
				A.ja[offset + 5] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 6] = BASE + IDX(i, j - 1);
				A.ja[offset + 7] = BASE + IDX(i, j);
				A.ja[offset + 8] = BASE + IDX(i, j + 1);
				A.ja[offset + 9] = BASE + IDX(i, j + 2);
				A.ja[offset + 10] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 11] = BASE + IDX(i + 1, j);
				A.ja[offset + 12] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 13] = BASE + IDX(i + 2, j);
				A.ja[offset + 14] = BASE + IDX(i + 2, j + 1);
				A.ja[offset + 15] = BASE + IDX(i + 2, j + 2);
				ell_f[IDX(i, j)] *= dr * dz;
				offset += 16;

				// Centered interior points.
				//
				// o       o       o
				//   \     |     /
				//     o   o   o
				//       \ | /
				// o---o---x---o---o
				//       / | \
				//     o   o   o
				//   /     |     \
				// o       o       o
				//
				// Thus:
				// 
				// u(i-2, j-2) : -b/48
				// u(i-2, j  ) : (-a + d)/12
				// u(i-2, j+2) : b/48
				// u(i-1, j-1) : b/3
				// u(i-1, j  ) : 2(2a - d)/3
				// u(i-1, j+1) : -b/3
				// u(i  , j-2) : (-c + e)/12
				// u(i  , j-1) : 2(2c - e)/3
				// u(i  , j  ) : -5(a + c)/2 + s
				// u(i  , j+1) : 2(2c + e)/3
				// u(i  , j+2) : -(c + e)/12
				// u(i+1, j-1) : -b/3
				// u(i+1, j  ) : 2(2a + d)/3
				// u(i+1, j+1) : b/3
				// u(i+2, j-2) : b/48
				// u(i+2, j  ) : -(a + d)/12
				// u(i+2, j+2) : -b/48
				//
				for (j = 2; j < NzInterior; j++)
				{
					aux_a = ell_a[IDX(i, j)];
					aux_b = ell_b[IDX(i, j)];
					aux_c = ell_c[IDX(i, j)];
					aux_d = dz * ell_d[IDX(i, j)];
					aux_e = dr * ell_e[IDX(i, j)];
					aux_s = dr * dz * ell_s[IDX(i, j)];

					A.ia[IDX(i, j)] = BASE + offset;
					A.a[offset] = -aux_b / 48.0;
					A.a[offset + 1] = (-aux_a + aux_d) * twelfth;
					A.a[offset + 2] = aux_b / 48.0;
					A.a[offset + 3] = aux_b * third;
					A.a[offset + 4] = 2.0 * (2.0 * aux_a - aux_d) * third;
					A.a[offset + 5] = -aux_b * third;
					A.a[offset + 6] = (-aux_c + aux_e) * twelfth;
					A.a[offset + 7] = 2.0 * (2.0 * aux_c - aux_e) * third;
					A.a[offset + 8] = -2.5 * (aux_a + aux_c) + aux_s;
					A.a[offset + 9] = 2.0 * (2.0 * aux_c + aux_e) * third;
					A.a[offset + 10] = -(aux_c + aux_e) * twelfth;
					A.a[offset + 11] = -aux_b * third;
					A.a[offset + 12] = 2.0 * (2.0 * aux_a + aux_d) * third;
					A.a[offset + 13] = aux_b * third;
					A.a[offset + 14] = aux_b / 48.0;
					A.a[offset + 15] = -(aux_a + aux_d) * twelfth;
					A.a[offset + 16] = -aux_b / 48.0;
					A.ja[offset] = BASE + IDX(i - 2, j - 2);
					A.ja[offset + 1] = BASE + IDX(i - 2, j);
					A.ja[offset + 2] = BASE + IDX(i - 2, j + 2);
					A.ja[offset + 3] = BASE + IDX(i - 1, j - 1);
					A.ja[offset + 4] = BASE + IDX(i - 1, j);
					A.ja[offset + 5] = BASE + IDX(i - 1, j + 1);
					A.ja[offset + 6] = BASE + IDX(i, j - 2);
					A.ja[offset + 7] = BASE + IDX(i, j - 1);
					A.ja[offset + 8] = BASE + IDX(i, j);
					A.ja[offset + 9] = BASE + IDX(i, j + 1);
					A.ja[offset + 10] = BASE + IDX(i, j + 2);
					A.ja[offset + 11] = BASE + IDX(i + 1, j - 1);
					A.ja[offset + 12] = BASE + IDX(i + 1, j);
					A.ja[offset + 13] = BASE + IDX(i + 1, j + 1);
					A.ja[offset + 14] = BASE + IDX(i + 2, j - 2);
					A.ja[offset + 15] = BASE + IDX(i + 2, j);
					A.ja[offset + 16] = BASE + IDX(i + 2, j + 2);
					ell_f[IDX(i, j)] *= dr * dz;
					offset += 17;
				}

				// Stencil for left strip.
				//
				// o   o   o   o   o
				//   \   \ | /   /
				// o---o---x---o---o
				//   /   / | \   \
				// o   o   o   o   o
				//   /   / | \   \
				// o   o   o   o   o
				//   /   / | \   \
				// o   o   o   o   o
				//         |      
				//         o        
				//
				// Thus:
				// 
				// u(i-2, j-3) : -b/144
				// u(i-2, j-2) : b/24
				// u(i-2, j-1) : -b/8
				// u(i-2, j  ) : (-6a + 5b + 6d)/72
				// u(i-2, j+1) : b/48
				// u(i-1, j-3) : b/18
				// u(i-1, j-2) : -b/3
				// u(i-1, j-1) : b
				// u(i-1, j  ) : (12a - 5b - 6d)/9
				// u(i-1, j+1) : -b/6
				// u(i  , j-4) : c/12
				// u(i  , j-3) : -(6c + e)/12
				// u(i  , j-2) : (7c + 3e)/6
				// u(i  , j-1) : -(2c + 9e)/6
				// u(i  , j  ) : -5a/2 -5c/4 + 5e/6 + s
				// u(i  , j+1) : 5c/6 + e/4
				// u(i+1, j-3) : -b/18
				// u(i+1, j-2) : b/3
				// u(i+1, j-1) : -b
				// u(i+1, j  ) : (12a + 5b + 6d)/9
				// u(i+1, j+1) : b/6
				// u(i+2, j-3) : b/144
				// u(i+2, j-2) : -b/24
				// u(i+2, j-1) : b/8
				// u(i+2, j  ) : -(6a + 5b + 6d)/72
				// u(i+2, j+1) : -b/48
				//
				j = NzInterior;
				aux_a = ell_a[IDX(i, j)];
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)];
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = -aux_b / 144.0;
				A.a[offset + 1] = aux_b / 24.0;
				A.a[offset + 2] = -aux_b * 0.125;
				A.a[offset + 3] = (-6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
				A.a[offset + 4] = aux_b / 48.0;
				A.a[offset + 5] = aux_b / 18.0;
				A.a[offset + 6] = -aux_b * third;
				A.a[offset + 7] = aux_b;
				A.a[offset + 8] = (12.0 * aux_a - 5.0 * aux_b - 6.0 * aux_d) / 9.0;
				A.a[offset + 9] = -aux_b * sixth;
				A.a[offset + 10] = aux_c * twelfth;
				A.a[offset + 11] = -(6.0 * aux_c + aux_e) * twelfth;
				A.a[offset + 12] = (7.0 * aux_c + 3.0 * aux_e) * sixth;
				A.a[offset + 13] = -(2.0 * aux_c + 9.0 * aux_e) * sixth;
				A.a[offset + 14] = -2.5 * aux_a - 1.25 * aux_c + 5.0 * aux_e * sixth + aux_s;
				A.a[offset + 15] = 5.0 * aux_c * sixth + aux_e * 0.25;
				A.a[offset + 16] = -aux_b / 18.0;
				A.a[offset + 17] = aux_b * third;
				A.a[offset + 18] = -aux_b;
				A.a[offset + 19] = (12.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 9.0;
				A.a[offset + 20] = aux_b * sixth;
				A.a[offset + 21] = aux_b / 144.0;
				A.a[offset + 22] = -aux_b / 24.0;
				A.a[offset + 23] = aux_b * 0.125;
				A.a[offset + 24] = -(6.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
				A.a[offset + 25] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 2, j - 3);
				A.ja[offset + 1] = BASE + IDX(i - 2, j - 2);
				A.ja[offset + 2] = BASE + IDX(i - 2, j - 1);
				A.ja[offset + 3] = BASE + IDX(i - 2, j);
				A.ja[offset + 4] = BASE + IDX(i - 2, j + 1);
				A.ja[offset + 5] = BASE + IDX(i - 1, j - 3);
				A.ja[offset + 6] = BASE + IDX(i - 1, j - 2);
				A.ja[offset + 7] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 8] = BASE + IDX(i - 1, j);
				A.ja[offset + 9] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 10] = BASE + IDX(i, j - 4);
				A.ja[offset + 11] = BASE + IDX(i, j - 3);
				A.ja[offset + 12] = BASE + IDX(i, j - 2);
				A.ja[offset + 13] = BASE + IDX(i, j - 1);
				A.ja[offset + 14] = BASE + IDX(i, j);
				A.ja[offset + 15] = BASE + IDX(i, j + 1);
				A.ja[offset + 16] = BASE + IDX(i + 1, j - 3);
				A.ja[offset + 17] = BASE + IDX(i + 1, j - 2);
				A.ja[offset + 18] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 19] = BASE + IDX(i + 1, j);
				A.ja[offset + 20] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 21] = BASE + IDX(i + 2, j - 3);
				A.ja[offset + 22] = BASE + IDX(i + 2, j - 2);
				A.ja[offset + 23] = BASE + IDX(i + 2, j - 1);
				A.ja[offset + 24] = BASE + IDX(i + 2, j);
				A.ja[offset + 25] = BASE + IDX(i + 2, j + 1);
				ell_f[IDX(i, j)] *= dr * dz;
				offset += 26;

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
				ell_f[IDX(i, j)] = uInf;
			}
		}

#ifdef DEBUG
		printf("Stage 3.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior
			+ 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
			+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin);

		// offset = (2 + 2) + 2 * NzInterior 
		// + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin.
		// + (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin).

		// Interior bottom-right corner with equatorial symmetry.
		i = NrInterior;
		r = (double)i - 0.5;

		j = 0;
		A.ia[IDX(NrInterior, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior, 1);
		ell_f[IDX(NrInterior, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2.


		// Right strip lower corner stencil.
		//
		//     o   o   o   o   o
		//       \   \   \ | /
		//     o   o   o   o   o
		//       \   \   \ | /
		// o---o---o---o---x---o
		//       /   /   / | \
		//     o   o   o   o   o
		// 
		// Thus:
		// 
		// u(i-4, j  ) : a/12
		// u(i-3, j-1) : b/18
		// u(i-3, j  ) : -(6a + d)/12
		// u(i-3, j+1) : -b/18 + z_sym*(-b/144)
		// u(i-3, j+2) : b/144
		// u(i-2, j-1) : -b/3
		// u(i-2, j  ) : (7a + 3d)/6
		// u(i-2, j+1) : b/3 + z_sym*(b/24)
		// u(i-2, j+2) : -b/24
		// u(i-1, j-1) : b
		// u(i-1, j  ) : -(2a + 9d)/6
		// u(i-1, j+1) : -b + z_sym*(-b/8)
		// u(i-1, j+2) : b/8
		// u(i  , j-1) : (-5b + 12c - 6e)/9
		// u(i  , j  ) : -5a/4 - 5c/2 + 5d/6 + s
		// u(i  , j+1) : (5b + 12c + 6e)/9 + z_sym*((5b - 6c + 6e)/72)
		// u(i  , j+2) : -(5b + 6c + 6e)/72
		// u(i+1, j-1) : -b/6
		// u(i+1, j  ) : 5a/6 + d/4
		// u(i+1, j+1) : b/6 + z_sym*(b/48)
		// u(i+1, j+2) : -b/48
		//
		j = 1;
		aux_a = ell_a[IDX(i, j)];
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)];
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = aux_a * twelfth;
		A.a[offset + 1] = aux_b / 18.0;
		A.a[offset + 2] = -(6.0 * aux_a + aux_d) * twelfth;
		A.a[offset + 3] = -aux_b / 18.0 + (double)z_sym * (-aux_b / 144.0);
		A.a[offset + 4] = aux_b / 144.0;
		A.a[offset + 5] = -aux_b * third;
		A.a[offset + 6] = (7.0 * aux_a + 3.0 * aux_d) * sixth;
		A.a[offset + 7] = aux_b * third + (double)z_sym * (aux_b / 24.0);
		A.a[offset + 8] = -aux_b / 24.0;
		A.a[offset + 9] = aux_b;
		A.a[offset + 10] = -(2.0 * aux_a + 9.0 * aux_d) * sixth;
		A.a[offset + 11] = -aux_b + (double)z_sym * (-aux_b * 0.125);
		A.a[offset + 12] = aux_b * 0.125;
		A.a[offset + 13] = (-5.0 * aux_b + 12.0 * aux_c - 6.0 * aux_e) / 9.0;
		A.a[offset + 14] = -1.25 * aux_a - 2.5 * aux_c + 5.0 * aux_d * sixth + aux_s;
		A.a[offset + 15] = (5.0 * aux_b + 12.0 * aux_c + 6.0 * aux_e) / 9.0 + (double)z_sym * ((5.0 * aux_b - 6.0 * aux_c + 6.0 * aux_e) / 72.0);
		A.a[offset + 16] = -(5.0 * aux_b + 6.0 * aux_c + 6.0 * aux_e) / 72.0;
		A.a[offset + 17] = -aux_b * sixth;
		A.a[offset + 18] = 5.0 * aux_a * sixth + aux_d * 0.25;
		A.a[offset + 19] = aux_b * sixth + (double)z_sym * (aux_b / 48.0);
		A.a[offset + 20] = -aux_b / 48.0;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j - 1);
		A.ja[offset + 2] = BASE + IDX(i - 3, j);
		A.ja[offset + 3] = BASE + IDX(i - 3, j + 1);
		A.ja[offset + 4] = BASE + IDX(i - 3, j + 2);
		A.ja[offset + 5] = BASE + IDX(i - 2, j - 1);
		A.ja[offset + 6] = BASE + IDX(i - 2, j);
		A.ja[offset + 7] = BASE + IDX(i - 2, j + 1);
		A.ja[offset + 8] = BASE + IDX(i - 2, j + 2);
		A.ja[offset + 9] = BASE + IDX(i - 1, j - 1);
		A.ja[offset + 10] = BASE + IDX(i - 1, j);
		A.ja[offset + 11] = BASE + IDX(i - 1, j + 1);
		A.ja[offset + 12] = BASE + IDX(i - 1, j + 2);
		A.ja[offset + 13] = BASE + IDX(i, j - 1);
		A.ja[offset + 14] = BASE + IDX(i, j);
		A.ja[offset + 15] = BASE + IDX(i, j + 1);
		A.ja[offset + 16] = BASE + IDX(i, j + 2);
		A.ja[offset + 17] = BASE + IDX(i + 1, j - 1);
		A.ja[offset + 18] = BASE + IDX(i + 1, j);
		A.ja[offset + 19] = BASE + IDX(i + 1, j + 1);
		A.ja[offset + 20] = BASE + IDX(i + 1, j + 2);
		ell_f[IDX(i, j)] *= dr * dz;
		offset += 21;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21.

		// Set temporary offset.
		t_offset = offset;

#pragma omp parallel shared(A, ell_f) private(offset,\
		aux_a, aux_b, aux_c, aux_d, aux_e, aux_s)
		{
#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Eac iteration fills 26 elements.
				offset = t_offset + (j - 2) * 26;

				// Right strip stencil.
				//
				//     o   o   o   o   o
				//       \   \   \ | /
				//     o   o   o   o   o
				//       \   \   \ | /
				// o---o---o---o---x---o
				//       /   /   / | \
				//     o   o   o   o   o
				//       /   /   /   \
				//     o   o   o   o   o
				// 
				// Thus:
				// 
				// u(i-4, j  ) : a/12
				// u(i-3, j-2) : -b/144
				// u(i-3, j-1) : b/18
				// u(i-3, j  ) : -(6a + d)/12
				// u(i-3, j+1) : -b/18
				// u(i-3, j+2) : b/144
				// u(i-2, j-2) : b/24
				// u(i-2, j-1) : -b/3
				// u(i-2, j  ) : (7a + 3d)/6
				// u(i-2, j+1) : b/3
				// u(i-2, j+2) : -b/24
				// u(i-1, j-2) : -b/8
				// u(i-1, j-1) : b
				// u(i-1, j  ) : -(2a + 9d)/6
				// u(i-1, j+1) : -b
				// u(i-1, j+2) : b/8
				// u(i  , j-2) : (5b - 6c + 6e)/72
				// u(i  , j-1) : (-5b + 12c - 6e)/9
				// u(i  , j  ) : -5a/4 - 5c/2 + 5d/6 + s
				// u(i  , j+1) : (5b + 12c + 6e)/9
				// u(i  , j+2) : -(5b + 6c + 6e)/72
				// u(i+1, j-2) : b/48
				// u(i+1, j-1) : -b/6
				// u(i+1, j  ) : 5a/6 + d/4
				// u(i+1, j+1) : b/6
				// u(i+1, j+2) : -b/48
				//
				aux_a = ell_a[IDX(i, j)];
				aux_b = ell_b[IDX(i, j)];
				aux_c = ell_c[IDX(i, j)];
				aux_d = dz * ell_d[IDX(i, j)];
				aux_e = dr * ell_e[IDX(i, j)];
				aux_s = dr * dz * ell_s[IDX(i, j)];

				A.ia[IDX(i, j)] = BASE + offset;
				A.a[offset] = aux_a * twelfth;
				A.a[offset + 1] = -aux_b / 144.0;
				A.a[offset + 2] = aux_b / 18.0;
				A.a[offset + 3] = -(6.0 * aux_a + aux_d) * twelfth;
				A.a[offset + 4] = -aux_b / 18.0;
				A.a[offset + 5] = aux_b / 144.0;
				A.a[offset + 6] = aux_b / 24.0;
				A.a[offset + 7] = -aux_b * third;
				A.a[offset + 8] = (7.0 * aux_a + 3.0 * aux_d) * sixth;
				A.a[offset + 9] = aux_b * third;
				A.a[offset + 10] = -aux_b / 24.0;
				A.a[offset + 11] = -aux_b * 0.125;
				A.a[offset + 12] = aux_b;
				A.a[offset + 13] = -(2.0 * aux_a + 9.0 * aux_d) * sixth;
				A.a[offset + 14] = -aux_b;
				A.a[offset + 15] = aux_b * 0.125;
				A.a[offset + 16] = (5.0 * aux_b - 6.0 * aux_c + 6.0 * aux_e) / 72.0;
				A.a[offset + 17] = (-5.0 * aux_b + 12.0 * aux_c - 6.0 * aux_e) / 9.0;
				A.a[offset + 18] = -1.25 * aux_a - 2.5 * aux_c + 5.0 * aux_d * sixth + aux_s;
				A.a[offset + 19] = (5.0 * aux_b + 12.0 * aux_c + 6.0 * aux_e) / 9.0;
				A.a[offset + 20] = -(5.0 * aux_b + 6.0 * aux_c + 6.0 * aux_e) / 72.0;
				A.a[offset + 21] = aux_b / 48.0;
				A.a[offset + 22] = -aux_b * sixth;
				A.a[offset + 23] = 5.0 * aux_a * sixth + aux_d * 0.25;
				A.a[offset + 24] = aux_b * sixth;
				A.a[offset + 25] = -aux_b / 48.0;
				A.ja[offset] = BASE + IDX(i - 4, j);
				A.ja[offset + 1] = BASE + IDX(i - 3, j - 2);
				A.ja[offset + 2] = BASE + IDX(i - 3, j - 1);
				A.ja[offset + 3] = BASE + IDX(i - 3, j);
				A.ja[offset + 4] = BASE + IDX(i - 3, j + 1);
				A.ja[offset + 5] = BASE + IDX(i - 3, j + 2);
				A.ja[offset + 6] = BASE + IDX(i - 2, j - 2);
				A.ja[offset + 7] = BASE + IDX(i - 2, j - 1);
				A.ja[offset + 8] = BASE + IDX(i - 2, j);
				A.ja[offset + 9] = BASE + IDX(i - 2, j + 1);
				A.ja[offset + 10] = BASE + IDX(i - 2, j + 2);
				A.ja[offset + 11] = BASE + IDX(i - 1, j - 2);
				A.ja[offset + 12] = BASE + IDX(i - 1, j - 1);
				A.ja[offset + 13] = BASE + IDX(i - 1, j);
				A.ja[offset + 14] = BASE + IDX(i - 1, j + 1);
				A.ja[offset + 15] = BASE + IDX(i - 1, j + 2);
				A.ja[offset + 16] = BASE + IDX(i, j - 2);
				A.ja[offset + 17] = BASE + IDX(i, j - 1);
				A.ja[offset + 18] = BASE + IDX(i, j);
				A.ja[offset + 19] = BASE + IDX(i, j + 1);
				A.ja[offset + 20] = BASE + IDX(i, j + 2);
				A.ja[offset + 21] = BASE + IDX(i + 1, j - 2);
				A.ja[offset + 22] = BASE + IDX(i + 1, j - 1);
				A.ja[offset + 23] = BASE + IDX(i + 1, j);
				A.ja[offset + 24] = BASE + IDX(i + 1, j + 1);
				A.ja[offset + 25] = BASE + IDX(i + 1, j + 2);
				ell_f[IDX(i, j)] *= dr * dz;
			}
		}


#ifdef DEBUG
		printf("Stage 4.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
			+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin)
			+ 2 + 21 + 26 * (NzInterior - 2);

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2).

		// Fill interior top-right corner.
		//
		//     o   o   o   o   o
		//       \   \   \ | /
		// o---o---o---o---x---o
		//       \  \    / | \
		//     o   o   o   o   o
		//       \   /   \ | \
		//     o   o   o   o   o
		//       /   \   \ | \
		//     o   o   o   o   o
		//                 |
		//                 o   
		//
		// Thus:
		//
		// u(i-4, j  ) : a/12
		// u(i-3, j-3) : b/144
		// u(i-3, j-2) : -b/24
		// u(i-3, j-1) : b/8
		// u(i-3, j  ) : -(36a + 5b + 6d)/72
		// u(i-3, j+1) : -b/48
		// u(i-2, j-3) : -b/24
		// u(i-2, j-2) : b/4
		// u(i-2, j-1) : -3b/4
		// u(i-2, j  ) : (14a + 5b + 6d)/12
		// u(i-2, j+1) : b/8
		// u(i-1, j-3) : b/8
		// u(i-1, j-2) : -3b/4
		// u(i-1, j-1) : 9b/4
		// u(i-1, j  ) : -a/3 - 5b/4 - 3d/2
		// u(i-1, j+1) : -3b/8
		// u(i  , j-4) : c/12
		// u(i  , j-3) : -(36c + 5b + 6e)/72
		// u(i  , j-2) : (14c + 5b + 6e)/12
		// u(i  , j-1) : -c/3 - 5b/4 - 3e/2
		// u(i  , j  ) : -5(9a - 5b + 9c - 6d - 6e)/36 + s
		// u(i  , j+1) : (5b + 20c + 6e)/24
		// u(i+1, j-3) : -b/48
		// u(i+1, j-2) : b/8
		// u(i+1, j-1) : -3b/8
		// u(i+1, j  ) : (20a + 5b + 6d)/24
		// u(i+1, j+1) : b/16
		//
		j = NzInterior;
		aux_a = ell_a[IDX(i, j)];
		aux_b = ell_b[IDX(i, j)];
		aux_c = ell_c[IDX(i, j)];
		aux_d = dz * ell_d[IDX(i, j)];
		aux_e = dr * ell_e[IDX(i, j)];
		aux_s = dr * dz * ell_s[IDX(i, j)];

		A.ia[IDX(i, j)] = BASE + offset;
		A.a[offset] = aux_a * twelfth;
		A.a[offset + 1] = aux_b / 144.0;
		A.a[offset + 2] = -aux_b / 24.0;
		A.a[offset + 3] = aux_b * 0.125;
		A.a[offset + 4] = -(36.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) / 72.0;
		A.a[offset + 5] = -aux_b / 48.0;
		A.a[offset + 6] = -aux_b / 24.0;
		A.a[offset + 7] = aux_b * 0.25;
		A.a[offset + 8] = -aux_b * 0.75;
		A.a[offset + 9] = (14.0 * aux_a + 5.0 * aux_b + 6.0 * aux_d) * twelfth;
		A.a[offset + 10] = aux_b * 0.125;
		A.a[offset + 11] = aux_b * 0.125;
		A.a[offset + 12] = -aux_b * 0.75;
		A.a[offset + 13] = 9.0 * aux_b * 0.25;
		A.a[offset + 14] = -aux_a * third - 5.0 * aux_b * 0.25 - 1.5 * aux_d;
		A.a[offset + 15] = -3.0 * aux_b * 0.125;
		A.a[offset + 16] = aux_c * twelfth;
		A.a[offset + 17] = -(36.0 * aux_c + 5.0 * aux_b + 6.0 * aux_e) / 72.0;
		A.a[offset + 18] = (14.0 * aux_c + 5.0 * aux_b + 6.0 * aux_e) * twelfth;
		A.a[offset + 19] = -aux_c * third - 5.0 * aux_b * 0.25 - 1.5 * aux_e;
		A.a[offset + 20] = -5.0 * (9.0 * aux_a - 5.0 * aux_b + 9.0 * aux_c - 6.0 * aux_d - 6.0 * aux_e) / 36.0 + aux_s;
		A.a[offset + 21] = (5.0 * aux_b + 20.0 * aux_c + 6.0 * aux_e) / 24.0;
		A.a[offset + 22] = -aux_b / 48.0;
		A.a[offset + 23] = aux_b * 0.125;
		A.a[offset + 24] = -3.0 * aux_b * 0.125;
		A.a[offset + 25] = (5.0 * aux_b + 20.0 * aux_a + 6.0 * aux_d) / 24.0;
		A.a[offset + 26] = aux_b * 0.0625;
		A.ja[offset] = BASE + IDX(i - 4, j);
		A.ja[offset + 1] = BASE + IDX(i - 3, j - 3);
		A.ja[offset + 2] = BASE + IDX(i - 3, j - 2);
		A.ja[offset + 3] = BASE + IDX(i - 3, j - 1);
		A.ja[offset + 4] = BASE + IDX(i - 3, j);
		A.ja[offset + 5] = BASE + IDX(i - 3, j + 1);
		A.ja[offset + 6] = BASE + IDX(i - 2, j - 3);
		A.ja[offset + 7] = BASE + IDX(i - 2, j - 2);
		A.ja[offset + 8] = BASE + IDX(i - 2, j - 1);
		A.ja[offset + 9] = BASE + IDX(i - 2, j);
		A.ja[offset + 10] = BASE + IDX(i - 2, j + 1);
		A.ja[offset + 11] = BASE + IDX(i - 1, j - 3);
		A.ja[offset + 12] = BASE + IDX(i - 1, j - 2);
		A.ja[offset + 13] = BASE + IDX(i - 1, j - 1);
		A.ja[offset + 14] = BASE + IDX(i - 1, j);
		A.ja[offset + 15] = BASE + IDX(i - 1, j + 1);
		A.ja[offset + 16] = BASE + IDX(i, j - 4);
		A.ja[offset + 17] = BASE + IDX(i, j - 3);
		A.ja[offset + 18] = BASE + IDX(i, j - 2);
		A.ja[offset + 19] = BASE + IDX(i, j - 1);
		A.ja[offset + 20] = BASE + IDX(i, j);
		A.ja[offset + 21] = BASE + IDX(i, j + 1);
		A.ja[offset + 22] = BASE + IDX(i + 1, j - 3);
		A.ja[offset + 23] = BASE + IDX(i + 1, j - 2);
		A.ja[offset + 24] = BASE + IDX(i + 1, j - 1);
		A.ja[offset + 25] = BASE + IDX(i + 1, j);
		A.ja[offset + 26] = BASE + IDX(i + 1, j + 1);
		ell_f[IDX(i, j)] *= dr * dz;
		offset += 27;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27

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
		ell_f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.

		// Do bottom-right corner.
		i = NrInterior + 1;
		r = (double)i - 0.5;

		j = 0;
		A.ia[IDX(NrInterior + 1, 0)] = BASE + offset;
		A.a[offset] = 1.0;
		A.a[offset + 1] = -(double)z_sym;
		A.ja[offset] = BASE + IDX(NrInterior + 1, 0);
		A.ja[offset + 1] = BASE + IDX(NrInterior + 1, 1);
		ell_f[IDX(NrInterior + 1, 0)] = 0.0;
		offset += 2;

		// offset = (2 + 2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.

		// Set temporary offset.
		t_offset = offset;

#pragma omp parallel shared(A, ell_f) private(offset, z, rr2,\
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
				ell_f[IDX(i, j)] = uInf;
			}
		}

#ifdef DEBUG
		printf("Stage 5.\n");
#endif

		// At this point we have filled:
		offset = (2 + 2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
			+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin)
			+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin
			+ n_robin * NzInterior;

		// offset = (2 + 2 + 2) + 2 * NzInterior + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// 	+ (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 +  n_robin)
		// 	+ 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.
		// 	+ n_robin * NzInterior.

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
		ell_f[IDX(i, j)] = uInf;
		offset += n_robin;

		// offset = (2 + 2 + 2 + n_robin) + 2 * NzInterior 
		// + 2 + 14 + 16 * (NzInterior - 2) + 21 + n_robin
		// + (NrInterior - 2) * (2 + 16 + 17 * (NzInterior - 2) + 26 + n_robin)
		// + 2 + 21 + 26 * (NzInterior - 2) + 27 + n_robin.
		// + n_robin * NzInterior.

		// Assert fill-in.
		assert(offset == nnz);

		// Fill last element.
		A.ia[IDX(NrInterior + 1, NzInterior + 1) + 1] = BASE + nnz;
	}

#ifdef DEBUG
	csr_print(A, "ge_A.asc", "ge_iA.asc", "ge_jA.asc");
	writeSingleFile(ell_f, "ge_f.asc");
#endif
	// All done.
	return;
}
