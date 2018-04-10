#include "tools.h"

const double half = 0.5;
const double twelfth = 1.0 / 12.0;

// R 2ND ORDER DERIVATIVES.
double r_2nd_interior_plus(const double *u, const int i, const int j, const double idr)
{
	return  half * idr * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i - 1, j)] + u[IDX(i - 2, j)]);
}
double r_2nd_interior_minus(const double *u, const int i, const int j, const double idr)
{
	return -half * idr * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i + 1, j)] + u[IDX(i + 2, j)]);
}
double r_2nd_1edge_plus(const double *u, const int i, const int j, const double idr)
{
	return  half * idr * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i - 1, j)] + u[IDX(i - 2, j)]);
}
double r_2nd_1edge_minus(const double *u, const int i, const int j, const double idr)
{
	return half * idr * (u[IDX(i + 1, j)] - u[IDX(i - 1, j)]);
}
double r_2nd_edge(const double *u, const int i, const int j, const double idr)
{
	return  half * idr * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i - 1, j)] + u[IDX(i - 2, j)]);
}

// R 4TH ORDER DERIVATIVES.
double r_4th_interior_plus(const double *u, const int i, const int j, const double idr)
{
	return twelfth * idr * (3.0 * u[IDX(i + 1, j)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i - 1, j)] + 6.0 * u[IDX(i - 2, j)] - u[IDX(i - 3, j)]);
}
double r_4th_interior_minus(const double *u, const int i, const int j, const double idr)
{
	return -twelfth * idr * (3.0 * u[IDX(i - 1, j)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i + 1, j)] + 6.0 * u[IDX(i + 2, j)] - u[IDX(i + 3, j)]);
}
double r_4th_1edge_plus(const double *u, const int i, const int j, const double idr)
{
	return twelfth * idr * (3.0 * u[IDX(i + 1, j)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i - 1, j)] + 6.0 * u[IDX(i - 2, j)] - u[IDX(i - 3, j)]);
}
double r_4th_1edge_minus(const double *u, const int i, const int j, const double idr)
{
	return twelfth * idr * (8.0 * (u[IDX(i + 1, j)] - u[IDX(i - 1, j)]) - (u[IDX(i + 2, j)] - u[IDX(i - 2, j)]));
}
double r_4th_2edge(const double *u, const int i, const int j, const double idr)
{
	return twelfth * idr * (3.0 * u[IDX(i + 1, j)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i - 1, j)] + 6.0 * u[IDX(i - 2, j)] - u[IDX(i - 3, j)]);
}
double r_4th_edge(const double *u, const int i, const int j, const double idr)
{
	return twelfth * idr * (25.0 * u[IDX(i, j)] - 48.0 * u[IDX(i - 1, j)] + 36.0 * u[IDX(i - 2, j)] - 16.0 * u[IDX(i - 3, j)] + 3.0 * u[IDX(i - 4, j)]);
}


// Z 2ND ORDER DERIVATIVES.
double z_2nd_interior_plus(const double *u, const int i, const int j, const double idz)
{
	return  half * idz * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i, j - 1)] + u[IDX(i, j - 2)]);
}
double z_2nd_interior_minus(const double *u, const int i, const int j, const double idz)
{
	return -half * idz * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i, j + 1)] + u[IDX(i, j + 2)]);
}
double z_2nd_1edge_plus(const double *u, const int i, const int j, const double idz)
{
	return  half * idz * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i, j - 1)] + u[IDX(i, j - 2)]);
}
double z_2nd_1edge_minus(const double *u, const int i, const int j, const double idz)
{
	return half * idz * (u[IDX(i, j + 1)] - u[IDX(i, j - 1)]);
}
double z_2nd_edge(const double *u, const int i, const int j, const double idz)
{
	return  half * idz * (3.0 * u[IDX(i, j)] - 4.0 * u[IDX(i, j - 1)] + u[IDX(i, j - 2)]);
}

// Z 4TH ORDER DERIVATIVES.
double z_4th_interior_plus(const double *u, const int i, const int j, const double idz)
{
	return twelfth * idz * (3.0 * u[IDX(i, j + 1)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i, j - 1)] + 6.0 * u[IDX(i, j - 2)] - u[IDX(i, j - 3)]);
}
double z_4th_interior_minus(const double *u, const int i, const int j, const double idz)
{
	return -twelfth * idz * (3.0 * u[IDX(i, j - 1)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i, j + 1)] + 6.0 * u[IDX(i, j + 2)] - u[IDX(i, j + 3)]);
}
double z_4th_1edge_plus(const double *u, const int i, const int j, const double idz)
{
	return twelfth * idz * (3.0 * u[IDX(i, j + 1)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i, j - 1)] + 6.0 * u[IDX(i, j - 2)] - u[IDX(i, j - 3)]);
}
double z_4th_1edge_minus(const double *u, const int i, const int j, const double idz)
{
	return twelfth * idz * (8.0 * (u[IDX(i, j + 1)] - u[IDX(i, j - 1)]) - (u[IDX(i, j + 2)] - u[IDX(i, j - 2)]));
}
double z_4th_2edge(const double *u, const int i, const int j, const double idz)
{
	return twelfth * idz * (3.0 * u[IDX(i, j + 1)] + 10.0 * u[IDX(i, j)] - 18.0 * u[IDX(i, j - 1)] + 6.0 * u[IDX(i, j - 2)] - u[IDX(i, j - 3)]);
}
double z_4th_edge(const double *u, const int i, const int j, const double idz)
{
	return twelfth * idz * (25.0 * u[IDX(i, j)] - 48.0 * u[IDX(i, j - 1)] + 36.0 * u[IDX(i, j - 2)] - 16.0 * u[IDX(i, j - 3)] + 3.0 * u[IDX(i, j - 4)]);
}

#include "param.h"
#include "arrays.h"

void shift_diffadvr(void)
{
	int i, j, k;

	double idr = 1.0 / dr;

	// Second order advective derivatives.
	if (strcmp(order, "two") == 0)
	{
		// Loop over Z.
		#pragma omp parallel shared(DAr_alpha, DAr_phi, DAr_a, DAr_b, DAr_h, DAr_c, DAr_lambda,\
		DAr_K, DAr_A_a, DAr_A_b, DAr_A_h, DAr_A_c, DAr_A_lambda, DAr_beta_r, DAr_beta_z,\
		DAr_dtbeta_r, DAr_dtbeta_z, DAr_Deltar, DAr_Deltaz) private(i)
		{
			#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				for (i = ghost; i < NrTotal - 2; i++)
				{
					if (-beta_r[IDX(i, j)] >= 0.0)
					{
						DAr_alpha[IDX(i, j)] = r_2nd_interior_plus(alpha, i, j, idr);
						DAr_phi[IDX(i, j)] = r_2nd_interior_plus(phi, i, j, idr);
						DAr_a[IDX(i, j)] = r_2nd_interior_plus(a, i, j, idr);
						DAr_b[IDX(i, j)] = r_2nd_interior_plus(b, i, j, idr);
						DAr_h[IDX(i, j)] = r_2nd_interior_plus(h, i, j, idr);
						DAr_c[IDX(i, j)] = r_2nd_interior_plus(c, i, j, idr);
						DAr_lambda[IDX(i, j)] = r_2nd_interior_plus(lambda, i, j, idr);
						DAr_K[IDX(i, j)] = r_2nd_interior_plus(K, i, j, idr);
						DAr_A_a[IDX(i, j)] = r_2nd_interior_plus(A_a, i, j, idr);
						DAr_A_b[IDX(i, j)] = r_2nd_interior_plus(A_b, i, j, idr);
						DAr_A_h[IDX(i, j)] = r_2nd_interior_plus(A_h, i, j, idr);
						DAr_A_c[IDX(i, j)] = r_2nd_interior_plus(A_c, i, j, idr);
						DAr_A_lambda[IDX(i, j)] = r_2nd_interior_plus(A_lambda, i, j, idr);
						DAr_beta_r[IDX(i, j)] = r_2nd_interior_plus(beta_r, i, j, idr);
						DAr_beta_z[IDX(i, j)] = r_2nd_interior_plus(beta_z, i, j, idr);
						DAr_dtbeta_r[IDX(i, j)] = r_2nd_interior_plus(dtbeta_r, i, j, idr);
						DAr_dtbeta_z[IDX(i, j)] = r_2nd_interior_plus(dtbeta_z, i, j, idr);
						DAr_Deltar[IDX(i, j)] = r_2nd_interior_plus(Deltar, i, j, idr);
						DAr_Deltaz[IDX(i, j)] = r_2nd_interior_plus(Deltaz, i, j, idr);
					}
					else
					{
						DAr_alpha[IDX(i, j)] = r_2nd_interior_minus(alpha, i, j, idr);
						DAr_phi[IDX(i, j)] = r_2nd_interior_minus(phi, i, j, idr);
						DAr_a[IDX(i, j)] = r_2nd_interior_minus(a, i, j, idr);
						DAr_b[IDX(i, j)] = r_2nd_interior_minus(b, i, j, idr);
						DAr_h[IDX(i, j)] = r_2nd_interior_minus(h, i, j, idr);
						DAr_c[IDX(i, j)] = r_2nd_interior_minus(c, i, j, idr);
						DAr_lambda[IDX(i, j)] = r_2nd_interior_minus(lambda, i, j, idr);
						DAr_K[IDX(i, j)] = r_2nd_interior_minus(K, i, j, idr);
						DAr_A_a[IDX(i, j)] = r_2nd_interior_minus(A_a, i, j, idr);
						DAr_A_b[IDX(i, j)] = r_2nd_interior_minus(A_b, i, j, idr);
						DAr_A_h[IDX(i, j)] = r_2nd_interior_minus(A_h, i, j, idr);
						DAr_A_c[IDX(i, j)] = r_2nd_interior_minus(A_c, i, j, idr);
						DAr_A_lambda[IDX(i, j)] = r_2nd_interior_minus(A_lambda, i, j, idr);
						DAr_beta_r[IDX(i, j)] = r_2nd_interior_minus(beta_r, i, j, idr);
						DAr_beta_z[IDX(i, j)] = r_2nd_interior_minus(beta_z, i, j, idr);
						DAr_dtbeta_r[IDX(i, j)] = r_2nd_interior_minus(dtbeta_r, i, j, idr);
						DAr_dtbeta_z[IDX(i, j)] = r_2nd_interior_minus(dtbeta_z, i, j, idr);
						DAr_Deltar[IDX(i, j)] = r_2nd_interior_minus(Deltar, i, j, idr);
						DAr_Deltaz[IDX(i, j)] = r_2nd_interior_minus(Deltaz, i, j, idr);
					}
				}

				// Second-to-last point.
				i = NrTotal - 2;
				if (-beta_r[IDX(i, j)] >= 0.0)
				{
					DAr_alpha[IDX(i, j)] = r_2nd_1edge_plus(alpha, i, j, idr);
					DAr_phi[IDX(i, j)] = r_2nd_1edge_plus(phi, i, j, idr);
					DAr_a[IDX(i, j)] = r_2nd_1edge_plus(a, i, j, idr);
					DAr_b[IDX(i, j)] = r_2nd_1edge_plus(b, i, j, idr);
					DAr_h[IDX(i, j)] = r_2nd_1edge_plus(h, i, j, idr);
					DAr_c[IDX(i, j)] = r_2nd_1edge_plus(c, i, j, idr);
					DAr_lambda[IDX(i, j)] = r_2nd_1edge_plus(lambda, i, j, idr);
					DAr_K[IDX(i, j)] = r_2nd_1edge_plus(K, i, j, idr);
					DAr_A_a[IDX(i, j)] = r_2nd_1edge_plus(A_a, i, j, idr);
					DAr_A_b[IDX(i, j)] = r_2nd_1edge_plus(A_b, i, j, idr);
					DAr_A_h[IDX(i, j)] = r_2nd_1edge_plus(A_h, i, j, idr);
					DAr_A_c[IDX(i, j)] = r_2nd_1edge_plus(A_c, i, j, idr);
					DAr_A_lambda[IDX(i, j)] = r_2nd_1edge_plus(A_lambda, i, j, idr);
					DAr_beta_r[IDX(i, j)] = r_2nd_1edge_plus(beta_r, i, j, idr);
					DAr_beta_z[IDX(i, j)] = r_2nd_1edge_plus(beta_z, i, j, idr);
					DAr_dtbeta_r[IDX(i, j)] = r_2nd_1edge_plus(dtbeta_r, i, j, idr);
					DAr_dtbeta_z[IDX(i, j)] = r_2nd_1edge_plus(dtbeta_z, i, j, idr);
					DAr_Deltar[IDX(i, j)] = r_2nd_1edge_plus(Deltar, i, j, idr);
					DAr_Deltaz[IDX(i, j)] = r_2nd_1edge_plus(Deltaz, i, j, idr);
				}
				else
				{
					DAr_alpha[IDX(i, j)] = r_2nd_1edge_minus(alpha, i, j, idr);
					DAr_phi[IDX(i, j)] = r_2nd_1edge_minus(phi, i, j, idr);
					DAr_a[IDX(i, j)] = r_2nd_1edge_minus(a, i, j, idr);
					DAr_b[IDX(i, j)] = r_2nd_1edge_minus(b, i, j, idr);
					DAr_h[IDX(i, j)] = r_2nd_1edge_minus(h, i, j, idr);
					DAr_c[IDX(i, j)] = r_2nd_1edge_minus(c, i, j, idr);
					DAr_lambda[IDX(i, j)] = r_2nd_1edge_minus(lambda, i, j, idr);
					DAr_K[IDX(i, j)] = r_2nd_1edge_minus(K, i, j, idr);
					DAr_A_a[IDX(i, j)] = r_2nd_1edge_minus(A_a, i, j, idr);
					DAr_A_b[IDX(i, j)] = r_2nd_1edge_minus(A_b, i, j, idr);
					DAr_A_h[IDX(i, j)] = r_2nd_1edge_minus(A_h, i, j, idr);
					DAr_A_c[IDX(i, j)] = r_2nd_1edge_minus(A_c, i, j, idr);
					DAr_A_lambda[IDX(i, j)] = r_2nd_1edge_minus(A_lambda, i, j, idr);
					DAr_beta_r[IDX(i, j)] = r_2nd_1edge_minus(beta_r, i, j, idr);
					DAr_beta_z[IDX(i, j)] = r_2nd_1edge_minus(beta_z, i, j, idr);
					DAr_dtbeta_r[IDX(i, j)] = r_2nd_1edge_minus(dtbeta_r, i, j, idr);
					DAr_dtbeta_z[IDX(i, j)] = r_2nd_1edge_minus(dtbeta_z, i, j, idr);
					DAr_Deltar[IDX(i, j)] = r_2nd_1edge_minus(Deltar, i, j, idr);
					DAr_Deltaz[IDX(i, j)] = r_2nd_1edge_minus(Deltaz, i, j, idr);
				}

				// Last point.
				i = NrTotal - 1;
				DAr_alpha[IDX(i, j)] = r_2nd_edge(alpha, i, j, idr);
				DAr_phi[IDX(i, j)] = r_2nd_edge(phi, i, j, idr);
				DAr_a[IDX(i, j)] = r_2nd_edge(a, i, j, idr);
				DAr_b[IDX(i, j)] = r_2nd_edge(b, i, j, idr);
				DAr_h[IDX(i, j)] = r_2nd_edge(h, i, j, idr);
				DAr_c[IDX(i, j)] = r_2nd_edge(c, i, j, idr);
				DAr_lambda[IDX(i, j)] = r_2nd_edge(lambda, i, j, idr);
				DAr_K[IDX(i, j)] = r_2nd_edge(K, i, j, idr);
				DAr_A_a[IDX(i, j)] = r_2nd_edge(A_a, i, j, idr);
				DAr_A_b[IDX(i, j)] = r_2nd_edge(A_b, i, j, idr);
				DAr_A_h[IDX(i, j)] = r_2nd_edge(A_h, i, j, idr);
				DAr_A_c[IDX(i, j)] = r_2nd_edge(A_c, i, j, idr);
				DAr_A_lambda[IDX(i, j)] = r_2nd_edge(A_lambda, i, j, idr);
				DAr_beta_r[IDX(i, j)] = r_2nd_edge(beta_r, i, j, idr);
				DAr_beta_z[IDX(i, j)] = r_2nd_edge(beta_z, i, j, idr);
				DAr_dtbeta_r[IDX(i, j)] = r_2nd_edge(dtbeta_r, i, j, idr);
				DAr_dtbeta_z[IDX(i, j)] = r_2nd_edge(dtbeta_z, i, j, idr);
				DAr_Deltar[IDX(i, j)] = r_2nd_edge(Deltar, i, j, idr);
				DAr_Deltaz[IDX(i, j)] = r_2nd_edge(Deltaz, i, j, idr);
			}
		}
	}
	else if (strcmp(order, "four") == 0)
	{
		// Loop over Z.
		#pragma omp parallel shared(DAr_alpha, DAr_phi, DAr_a, DAr_b, DAr_h, DAr_c, DAr_lambda,\
		DAr_K, DAr_A_a, DAr_A_b, DAr_A_h, DAr_A_c, DAr_A_lambda, DAr_beta_r, DAr_beta_z,\
		DAr_dtbeta_r, DAr_dtbeta_z, DAr_Deltar, DAr_Deltaz) private(i)
		{
			#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				for (i = ghost; i < NrTotal - 3; i++)
				{
					if (-beta_r[IDX(i, j)] >= 0.0)
					{
						DAr_alpha[IDX(i, j)] = r_4th_interior_plus(alpha, i, j, idr);
						DAr_phi[IDX(i, j)] = r_4th_interior_plus(phi, i, j, idr);
						DAr_a[IDX(i, j)] = r_4th_interior_plus(a, i, j, idr);
						DAr_b[IDX(i, j)] = r_4th_interior_plus(b, i, j, idr);
						DAr_h[IDX(i, j)] = r_4th_interior_plus(h, i, j, idr);
						DAr_c[IDX(i, j)] = r_4th_interior_plus(c, i, j, idr);
						DAr_lambda[IDX(i, j)] = r_4th_interior_plus(lambda, i, j, idr);
						DAr_K[IDX(i, j)] = r_4th_interior_plus(K, i, j, idr);
						DAr_A_a[IDX(i, j)] = r_4th_interior_plus(A_a, i, j, idr);
						DAr_A_b[IDX(i, j)] = r_4th_interior_plus(A_b, i, j, idr);
						DAr_A_h[IDX(i, j)] = r_4th_interior_plus(A_h, i, j, idr);
						DAr_A_c[IDX(i, j)] = r_4th_interior_plus(A_c, i, j, idr);
						DAr_A_lambda[IDX(i, j)] = r_4th_interior_plus(A_lambda, i, j, idr);
						DAr_beta_r[IDX(i, j)] = r_4th_interior_plus(beta_r, i, j, idr);
						DAr_beta_z[IDX(i, j)] = r_4th_interior_plus(beta_z, i, j, idr);
						DAr_dtbeta_r[IDX(i, j)] = r_4th_interior_plus(dtbeta_r, i, j, idr);
						DAr_dtbeta_z[IDX(i, j)] = r_4th_interior_plus(dtbeta_z, i, j, idr);
						DAr_Deltar[IDX(i, j)] = r_4th_interior_plus(Deltar, i, j, idr);
						DAr_Deltaz[IDX(i, j)] = r_4th_interior_plus(Deltaz, i, j, idr);
					}
					else
					{
						DAr_alpha[IDX(i, j)] = r_4th_interior_minus(alpha, i, j, idr);
						DAr_phi[IDX(i, j)] = r_4th_interior_minus(phi, i, j, idr);
						DAr_a[IDX(i, j)] = r_4th_interior_minus(a, i, j, idr);
						DAr_b[IDX(i, j)] = r_4th_interior_minus(b, i, j, idr);
						DAr_h[IDX(i, j)] = r_4th_interior_minus(h, i, j, idr);
						DAr_c[IDX(i, j)] = r_4th_interior_minus(c, i, j, idr);
						DAr_lambda[IDX(i, j)] = r_4th_interior_minus(lambda, i, j, idr);
						DAr_K[IDX(i, j)] = r_4th_interior_minus(K, i, j, idr);
						DAr_A_a[IDX(i, j)] = r_4th_interior_minus(A_a, i, j, idr);
						DAr_A_b[IDX(i, j)] = r_4th_interior_minus(A_b, i, j, idr);
						DAr_A_h[IDX(i, j)] = r_4th_interior_minus(A_h, i, j, idr);
						DAr_A_c[IDX(i, j)] = r_4th_interior_minus(A_c, i, j, idr);
						DAr_A_lambda[IDX(i, j)] = r_4th_interior_minus(A_lambda, i, j, idr);
						DAr_beta_r[IDX(i, j)] = r_4th_interior_minus(beta_r, i, j, idr);
						DAr_beta_z[IDX(i, j)] = r_4th_interior_minus(beta_z, i, j, idr);
						DAr_dtbeta_r[IDX(i, j)] = r_4th_interior_minus(dtbeta_r, i, j, idr);
						DAr_dtbeta_z[IDX(i, j)] = r_4th_interior_minus(dtbeta_z, i, j, idr);
						DAr_Deltar[IDX(i, j)] = r_4th_interior_minus(Deltar, i, j, idr);
						DAr_Deltaz[IDX(i, j)] = r_4th_interior_minus(Deltaz, i, j, idr);
					}
				}

				// Third-to-last point.
				i = NrTotal - 3;
				if (-beta_r[IDX(i, j)] >= 0.0)
				{
					DAr_alpha[IDX(i, j)] = r_4th_1edge_plus(alpha, i, j, idr);
					DAr_phi[IDX(i, j)] = r_4th_1edge_plus(phi, i, j, idr);
					DAr_a[IDX(i, j)] = r_4th_1edge_plus(a, i, j, idr);
					DAr_b[IDX(i, j)] = r_4th_1edge_plus(b, i, j, idr);
					DAr_h[IDX(i, j)] = r_4th_1edge_plus(h, i, j, idr);
					DAr_c[IDX(i, j)] = r_4th_1edge_plus(c, i, j, idr);
					DAr_lambda[IDX(i, j)] = r_4th_1edge_plus(lambda, i, j, idr);
					DAr_K[IDX(i, j)] = r_4th_1edge_plus(K, i, j, idr);
					DAr_A_a[IDX(i, j)] = r_4th_1edge_plus(A_a, i, j, idr);
					DAr_A_b[IDX(i, j)] = r_4th_1edge_plus(A_b, i, j, idr);
					DAr_A_h[IDX(i, j)] = r_4th_1edge_plus(A_h, i, j, idr);
					DAr_A_c[IDX(i, j)] = r_4th_1edge_plus(A_c, i, j, idr);
					DAr_A_lambda[IDX(i, j)] = r_4th_1edge_plus(A_lambda, i, j, idr);
					DAr_beta_r[IDX(i, j)] = r_4th_1edge_plus(beta_r, i, j, idr);
					DAr_beta_z[IDX(i, j)] = r_4th_1edge_plus(beta_z, i, j, idr);
					DAr_dtbeta_r[IDX(i, j)] = r_4th_1edge_plus(dtbeta_r, i, j, idr);
					DAr_dtbeta_z[IDX(i, j)] = r_4th_1edge_plus(dtbeta_z, i, j, idr);
					DAr_Deltar[IDX(i, j)] = r_4th_1edge_plus(Deltar, i, j, idr);
					DAr_Deltaz[IDX(i, j)] = r_4th_1edge_plus(Deltaz, i, j, idr);
				}
				else
				{
					DAr_alpha[IDX(i, j)] = r_4th_1edge_minus(alpha, i, j, idr);
					DAr_phi[IDX(i, j)] = r_4th_1edge_minus(phi, i, j, idr);
					DAr_a[IDX(i, j)] = r_4th_1edge_minus(a, i, j, idr);
					DAr_b[IDX(i, j)] = r_4th_1edge_minus(b, i, j, idr);
					DAr_h[IDX(i, j)] = r_4th_1edge_minus(h, i, j, idr);
					DAr_c[IDX(i, j)] = r_4th_1edge_minus(c, i, j, idr);
					DAr_lambda[IDX(i, j)] = r_4th_1edge_minus(lambda, i, j, idr);
					DAr_K[IDX(i, j)] = r_4th_1edge_minus(K, i, j, idr);
					DAr_A_a[IDX(i, j)] = r_4th_1edge_minus(A_a, i, j, idr);
					DAr_A_b[IDX(i, j)] = r_4th_1edge_minus(A_b, i, j, idr);
					DAr_A_h[IDX(i, j)] = r_4th_1edge_minus(A_h, i, j, idr);
					DAr_A_c[IDX(i, j)] = r_4th_1edge_minus(A_c, i, j, idr);
					DAr_A_lambda[IDX(i, j)] = r_4th_1edge_minus(A_lambda, i, j, idr);
					DAr_beta_r[IDX(i, j)] = r_4th_1edge_minus(beta_r, i, j, idr);
					DAr_beta_z[IDX(i, j)] = r_4th_1edge_minus(beta_z, i, j, idr);
					DAr_dtbeta_r[IDX(i, j)] = r_4th_1edge_minus(dtbeta_r, i, j, idr);
					DAr_dtbeta_z[IDX(i, j)] = r_4th_1edge_minus(dtbeta_z, i, j, idr);
					DAr_Deltar[IDX(i, j)] = r_4th_1edge_minus(Deltar, i, j, idr);
					DAr_Deltaz[IDX(i, j)] = r_4th_1edge_minus(Deltaz, i, j, idr);
				}

				// Second to last point.
				i = NrTotal - 2;
				DAr_alpha[IDX(i, j)] = r_4th_2edge(alpha, i, j, idr);
				DAr_phi[IDX(i, j)] = r_4th_2edge(phi, i, j, idr);
				DAr_a[IDX(i, j)] = r_4th_2edge(a, i, j, idr);
				DAr_b[IDX(i, j)] = r_4th_2edge(b, i, j, idr);
				DAr_h[IDX(i, j)] = r_4th_2edge(h, i, j, idr);
				DAr_c[IDX(i, j)] = r_4th_2edge(c, i, j, idr);
				DAr_lambda[IDX(i, j)] = r_4th_2edge(lambda, i, j, idr);
				DAr_K[IDX(i, j)] = r_4th_2edge(K, i, j, idr);
				DAr_A_a[IDX(i, j)] = r_4th_2edge(A_a, i, j, idr);
				DAr_A_b[IDX(i, j)] = r_4th_2edge(A_b, i, j, idr);
				DAr_A_h[IDX(i, j)] = r_4th_2edge(A_h, i, j, idr);
				DAr_A_c[IDX(i, j)] = r_4th_2edge(A_c, i, j, idr);
				DAr_A_lambda[IDX(i, j)] = r_4th_2edge(A_lambda, i, j, idr);
				DAr_beta_r[IDX(i, j)] = r_4th_2edge(beta_r, i, j, idr);
				DAr_beta_z[IDX(i, j)] = r_4th_2edge(beta_z, i, j, idr);
				DAr_dtbeta_r[IDX(i, j)] = r_4th_2edge(dtbeta_r, i, j, idr);
				DAr_dtbeta_z[IDX(i, j)] = r_4th_2edge(dtbeta_z, i, j, idr);
				DAr_Deltar[IDX(i, j)] = r_4th_2edge(Deltar, i, j, idr);
				DAr_Deltaz[IDX(i, j)] = r_4th_2edge(Deltaz, i, j, idr);

				// Last point.
				i = NrTotal - 1;
				DAr_alpha[IDX(i, j)] = r_4th_edge(alpha, i, j, idr);
				DAr_phi[IDX(i, j)] = r_4th_edge(phi, i, j, idr);
				DAr_a[IDX(i, j)] = r_4th_edge(a, i, j, idr);
				DAr_b[IDX(i, j)] = r_4th_edge(b, i, j, idr);
				DAr_h[IDX(i, j)] = r_4th_edge(h, i, j, idr);
				DAr_c[IDX(i, j)] = r_4th_edge(c, i, j, idr);
				DAr_lambda[IDX(i, j)] = r_4th_edge(lambda, i, j, idr);
				DAr_K[IDX(i, j)] = r_4th_edge(K, i, j, idr);
				DAr_A_a[IDX(i, j)] = r_4th_edge(A_a, i, j, idr);
				DAr_A_b[IDX(i, j)] = r_4th_edge(A_b, i, j, idr);
				DAr_A_h[IDX(i, j)] = r_4th_edge(A_h, i, j, idr);
				DAr_A_c[IDX(i, j)] = r_4th_edge(A_c, i, j, idr);
				DAr_A_lambda[IDX(i, j)] = r_4th_edge(A_lambda, i, j, idr);
				DAr_beta_r[IDX(i, j)] = r_4th_edge(beta_r, i, j, idr);
				DAr_beta_z[IDX(i, j)] = r_4th_edge(beta_z, i, j, idr);
				DAr_dtbeta_r[IDX(i, j)] = r_4th_edge(dtbeta_r, i, j, idr);
				DAr_dtbeta_z[IDX(i, j)] = r_4th_edge(dtbeta_z, i, j, idr);
				DAr_Deltar[IDX(i, j)] = r_4th_edge(Deltar, i, j, idr);
				DAr_Deltaz[IDX(i, j)] = r_4th_edge(Deltaz, i, j, idr);
			}
		}
	}

	// Axis symmetries.
	for (k = 0; k < ghost; k++)
	{
		#pragma omp parallel shared(DAr_alpha, DAr_phi, DAr_a, DAr_b, DAr_h, DAr_c, DAr_lambda,\
		DAr_K, DAr_A_a, DAr_A_b, DAr_A_h, DAr_A_c, DAr_A_lambda, DAr_beta_r, DAr_beta_z,\
		DAr_dtbeta_r, DAr_dtbeta_z, DAr_Deltar, DAr_Deltaz)
		{
			#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				DAr_alpha[IDX(ghost - 1 - k, j)] = -DAr_alpha[IDX(ghost + k, j)];
				DAr_phi[IDX(ghost - 1 - k, j)] = -DAr_phi[IDX(ghost + k, j)];
				DAr_a[IDX(ghost - 1 - k, j)] = -DAr_a[IDX(ghost + k, j)];
				DAr_b[IDX(ghost - 1 - k, j)] = -DAr_b[IDX(ghost + k, j)];
				DAr_h[IDX(ghost - 1 - k, j)] = -DAr_h[IDX(ghost + k, j)];
				DAr_c[IDX(ghost - 1 - k, j)] = -DAr_c[IDX(ghost + k, j)];
				DAr_lambda[IDX(ghost - 1 - k, j)] = -DAr_lambda[IDX(ghost + k, j)];
				DAr_K[IDX(ghost - 1 - k, j)] = -DAr_K[IDX(ghost + k, j)];
				DAr_A_a[IDX(ghost - 1 - k, j)] = -DAr_A_a[IDX(ghost + k, j)];
				DAr_A_b[IDX(ghost - 1 - k, j)] = -DAr_A_b[IDX(ghost + k, j)];
				DAr_A_h[IDX(ghost - 1 - k, j)] = -DAr_A_h[IDX(ghost + k, j)];
				DAr_A_c[IDX(ghost - 1 - k, j)] = -DAr_A_c[IDX(ghost + k, j)];
				DAr_A_lambda[IDX(ghost - 1 - k, j)] = -DAr_A_lambda[IDX(ghost + k, j)];
				DAr_beta_r[IDX(ghost - 1 - k, j)] = DAr_beta_r[IDX(ghost + k, j)];
				DAr_beta_z[IDX(ghost - 1 - k, j)] = -DAr_beta_z[IDX(ghost + k, j)];
				DAr_dtbeta_r[IDX(ghost - 1 - k, j)] = DAr_dtbeta_r[IDX(ghost + k, j)];
				DAr_dtbeta_z[IDX(ghost - 1 - k, j)] = -DAr_dtbeta_z[IDX(ghost + k, j)];
				DAr_Deltar[IDX(ghost - 1 - k, j)] = DAr_Deltar[IDX(ghost + k, j)];
				DAr_Deltaz[IDX(ghost - 1 - k, j)] = -DAr_Deltaz[IDX(ghost + k, j)];
			}
		}
	}

	return;
}

void shift_diffadvz(void)
{
	int i, j, k;

	double idz = 1.0 / dz;

	// Second order advective derivatives.
	if (strcmp(order, "two") == 0)
	{
		// Loop over Z.
		#pragma omp parallel shared(DAz_alpha, DAz_phi, DAz_a, DAz_b, DAz_h, DAz_c, DAz_lambda,\
		DAz_K, DAz_A_a, DAz_A_b, DAz_A_h, DAz_A_c, DAz_A_lambda, DAz_beta_r, DAz_beta_z,\
		DAz_dtbeta_r, DAz_dtbeta_z, DAz_Deltar, DAz_Deltaz) private(j)
		{
			#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				for (j = ghost; j < NzTotal - 2; j++)
				{
					if (-beta_r[IDX(i, j)] >= 0.0)
					{
						DAz_alpha[IDX(i, j)] = z_2nd_interior_plus(alpha, i, j, idz);
						DAz_phi[IDX(i, j)] = z_2nd_interior_plus(phi, i, j, idz);
						DAz_a[IDX(i, j)] = z_2nd_interior_plus(a, i, j, idz);
						DAz_b[IDX(i, j)] = z_2nd_interior_plus(b, i, j, idz);
						DAz_h[IDX(i, j)] = z_2nd_interior_plus(h, i, j, idz);
						DAz_c[IDX(i, j)] = z_2nd_interior_plus(c, i, j, idz);
						DAz_lambda[IDX(i, j)] = z_2nd_interior_plus(lambda, i, j, idz);
						DAz_K[IDX(i, j)] = z_2nd_interior_plus(K, i, j, idz);
						DAz_A_a[IDX(i, j)] = z_2nd_interior_plus(A_a, i, j, idz);
						DAz_A_b[IDX(i, j)] = z_2nd_interior_plus(A_b, i, j, idz);
						DAz_A_h[IDX(i, j)] = z_2nd_interior_plus(A_h, i, j, idz);
						DAz_A_c[IDX(i, j)] = z_2nd_interior_plus(A_c, i, j, idz);
						DAz_A_lambda[IDX(i, j)] = z_2nd_interior_plus(A_lambda, i, j, idz);
						DAz_beta_r[IDX(i, j)] = z_2nd_interior_plus(beta_r, i, j, idz);
						DAz_beta_z[IDX(i, j)] = z_2nd_interior_plus(beta_z, i, j, idz);
						DAz_dtbeta_r[IDX(i, j)] = z_2nd_interior_plus(dtbeta_r, i, j, idz);
						DAz_dtbeta_z[IDX(i, j)] = z_2nd_interior_plus(dtbeta_z, i, j, idz);
						DAz_Deltar[IDX(i, j)] = z_2nd_interior_plus(Deltar, i, j, idz);
						DAz_Deltaz[IDX(i, j)] = z_2nd_interior_plus(Deltaz, i, j, idz);
					}
					else
					{
						DAz_alpha[IDX(i, j)] = z_2nd_interior_minus(alpha, i, j, idz);
						DAz_phi[IDX(i, j)] = z_2nd_interior_minus(phi, i, j, idz);
						DAz_a[IDX(i, j)] = z_2nd_interior_minus(a, i, j, idz);
						DAz_b[IDX(i, j)] = z_2nd_interior_minus(b, i, j, idz);
						DAz_h[IDX(i, j)] = z_2nd_interior_minus(h, i, j, idz);
						DAz_c[IDX(i, j)] = z_2nd_interior_minus(c, i, j, idz);
						DAz_lambda[IDX(i, j)] = z_2nd_interior_minus(lambda, i, j, idz);
						DAz_K[IDX(i, j)] = z_2nd_interior_minus(K, i, j, idz);
						DAz_A_a[IDX(i, j)] = z_2nd_interior_minus(A_a, i, j, idz);
						DAz_A_b[IDX(i, j)] = z_2nd_interior_minus(A_b, i, j, idz);
						DAz_A_h[IDX(i, j)] = z_2nd_interior_minus(A_h, i, j, idz);
						DAz_A_c[IDX(i, j)] = z_2nd_interior_minus(A_c, i, j, idz);
						DAz_A_lambda[IDX(i, j)] = z_2nd_interior_minus(A_lambda, i, j, idz);
						DAz_beta_r[IDX(i, j)] = z_2nd_interior_minus(beta_r, i, j, idz);
						DAz_beta_z[IDX(i, j)] = z_2nd_interior_minus(beta_z, i, j, idz);
						DAz_dtbeta_r[IDX(i, j)] = z_2nd_interior_minus(dtbeta_r, i, j, idz);
						DAz_dtbeta_z[IDX(i, j)] = z_2nd_interior_minus(dtbeta_z, i, j, idz);
						DAz_Deltar[IDX(i, j)] = z_2nd_interior_minus(Deltar, i, j, idz);
						DAz_Deltaz[IDX(i, j)] = z_2nd_interior_minus(Deltaz, i, j, idz);
					}
				}

				// Second-to-last point.
				j = NzTotal - 2;
				if (-beta_r[IDX(i, j)] >= 0.0)
				{
					DAz_alpha[IDX(i, j)] = z_2nd_1edge_plus(alpha, i, j, idz);
					DAz_phi[IDX(i, j)] = z_2nd_1edge_plus(phi, i, j, idz);
					DAz_a[IDX(i, j)] = z_2nd_1edge_plus(a, i, j, idz);
					DAz_b[IDX(i, j)] = z_2nd_1edge_plus(b, i, j, idz);
					DAz_h[IDX(i, j)] = z_2nd_1edge_plus(h, i, j, idz);
					DAz_c[IDX(i, j)] = z_2nd_1edge_plus(c, i, j, idz);
					DAz_lambda[IDX(i, j)] = z_2nd_1edge_plus(lambda, i, j, idz);
					DAz_K[IDX(i, j)] = z_2nd_1edge_plus(K, i, j, idz);
					DAz_A_a[IDX(i, j)] = z_2nd_1edge_plus(A_a, i, j, idz);
					DAz_A_b[IDX(i, j)] = z_2nd_1edge_plus(A_b, i, j, idz);
					DAz_A_h[IDX(i, j)] = z_2nd_1edge_plus(A_h, i, j, idz);
					DAz_A_c[IDX(i, j)] = z_2nd_1edge_plus(A_c, i, j, idz);
					DAz_A_lambda[IDX(i, j)] = z_2nd_1edge_plus(A_lambda, i, j, idz);
					DAz_beta_r[IDX(i, j)] = z_2nd_1edge_plus(beta_r, i, j, idz);
					DAz_beta_z[IDX(i, j)] = z_2nd_1edge_plus(beta_z, i, j, idz);
					DAz_dtbeta_r[IDX(i, j)] = z_2nd_1edge_plus(dtbeta_r, i, j, idz);
					DAz_dtbeta_z[IDX(i, j)] = z_2nd_1edge_plus(dtbeta_z, i, j, idz);
					DAz_Deltar[IDX(i, j)] = z_2nd_1edge_plus(Deltar, i, j, idz);
					DAz_Deltaz[IDX(i, j)] = z_2nd_1edge_plus(Deltaz, i, j, idz);
				}
				else
				{
					DAz_alpha[IDX(i, j)] = z_2nd_1edge_minus(alpha, i, j, idz);
					DAz_phi[IDX(i, j)] = z_2nd_1edge_minus(phi, i, j, idz);
					DAz_a[IDX(i, j)] = z_2nd_1edge_minus(a, i, j, idz);
					DAz_b[IDX(i, j)] = z_2nd_1edge_minus(b, i, j, idz);
					DAz_h[IDX(i, j)] = z_2nd_1edge_minus(h, i, j, idz);
					DAz_c[IDX(i, j)] = z_2nd_1edge_minus(c, i, j, idz);
					DAz_lambda[IDX(i, j)] = z_2nd_1edge_minus(lambda, i, j, idz);
					DAz_K[IDX(i, j)] = z_2nd_1edge_minus(K, i, j, idz);
					DAz_A_a[IDX(i, j)] = z_2nd_1edge_minus(A_a, i, j, idz);
					DAz_A_b[IDX(i, j)] = z_2nd_1edge_minus(A_b, i, j, idz);
					DAz_A_h[IDX(i, j)] = z_2nd_1edge_minus(A_h, i, j, idz);
					DAz_A_c[IDX(i, j)] = z_2nd_1edge_minus(A_c, i, j, idz);
					DAz_A_lambda[IDX(i, j)] = z_2nd_1edge_minus(A_lambda, i, j, idz);
					DAz_beta_r[IDX(i, j)] = z_2nd_1edge_minus(beta_r, i, j, idz);
					DAz_beta_z[IDX(i, j)] = z_2nd_1edge_minus(beta_z, i, j, idz);
					DAz_dtbeta_r[IDX(i, j)] = z_2nd_1edge_minus(dtbeta_r, i, j, idz);
					DAz_dtbeta_z[IDX(i, j)] = z_2nd_1edge_minus(dtbeta_z, i, j, idz);
					DAz_Deltar[IDX(i, j)] = z_2nd_1edge_minus(Deltar, i, j, idz);
					DAz_Deltaz[IDX(i, j)] = z_2nd_1edge_minus(Deltaz, i, j, idz);
				}

				// Last point.
				j = NrTotal - 1;
				DAz_alpha[IDX(i, j)] = z_2nd_edge(alpha, i, j, idz);
				DAz_phi[IDX(i, j)] = z_2nd_edge(phi, i, j, idz);
				DAz_a[IDX(i, j)] = z_2nd_edge(a, i, j, idz);
				DAz_b[IDX(i, j)] = z_2nd_edge(b, i, j, idz);
				DAz_h[IDX(i, j)] = z_2nd_edge(h, i, j, idz);
				DAz_c[IDX(i, j)] = z_2nd_edge(c, i, j, idz);
				DAz_lambda[IDX(i, j)] = z_2nd_edge(lambda, i, j, idz);
				DAz_K[IDX(i, j)] = z_2nd_edge(K, i, j, idz);
				DAz_A_a[IDX(i, j)] = z_2nd_edge(A_a, i, j, idz);
				DAz_A_b[IDX(i, j)] = z_2nd_edge(A_b, i, j, idz);
				DAz_A_h[IDX(i, j)] = z_2nd_edge(A_h, i, j, idz);
				DAz_A_c[IDX(i, j)] = z_2nd_edge(A_c, i, j, idz);
				DAz_A_lambda[IDX(i, j)] = z_2nd_edge(A_lambda, i, j, idz);
				DAz_beta_r[IDX(i, j)] = z_2nd_edge(beta_r, i, j, idz);
				DAz_beta_z[IDX(i, j)] = z_2nd_edge(beta_z, i, j, idz);
				DAz_dtbeta_r[IDX(i, j)] = z_2nd_edge(dtbeta_r, i, j, idz);
				DAz_dtbeta_z[IDX(i, j)] = z_2nd_edge(dtbeta_z, i, j, idz);
				DAz_Deltar[IDX(i, j)] = z_2nd_edge(Deltar, i, j, idz);
				DAz_Deltaz[IDX(i, j)] = z_2nd_edge(Deltaz, i, j, idz);
			}
		}
	}
	else if (strcmp(order, "four") == 0)
	{
		// Loop over Z.
		#pragma omp parallel shared(DAz_alpha, DAz_phi, DAz_a, DAz_b, DAz_h, DAz_c, DAz_lambda,\
		DAz_K, DAz_A_a, DAz_A_b, DAz_A_h, DAz_A_c, DAz_A_lambda, DAz_beta_r, DAz_beta_z,\
		DAz_dtbeta_r, DAz_dtbeta_z, DAz_Deltar, DAz_Deltaz) private(j)
		{
			#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				for (j = ghost; j < NzTotal - 3; j++)
				{
					if (-beta_r[IDX(i, j)] >= 0.0)
					{
						DAz_alpha[IDX(i, j)] = z_4th_interior_plus(alpha, i, j, idz);
						DAz_phi[IDX(i, j)] = z_4th_interior_plus(phi, i, j, idz);
						DAz_a[IDX(i, j)] = z_4th_interior_plus(a, i, j, idz);
						DAz_b[IDX(i, j)] = z_4th_interior_plus(b, i, j, idz);
						DAz_h[IDX(i, j)] = z_4th_interior_plus(h, i, j, idz);
						DAz_c[IDX(i, j)] = z_4th_interior_plus(c, i, j, idz);
						DAz_lambda[IDX(i, j)] = z_4th_interior_plus(lambda, i, j, idz);
						DAz_K[IDX(i, j)] = z_4th_interior_plus(K, i, j, idz);
						DAz_A_a[IDX(i, j)] = z_4th_interior_plus(A_a, i, j, idz);
						DAz_A_b[IDX(i, j)] = z_4th_interior_plus(A_b, i, j, idz);
						DAz_A_h[IDX(i, j)] = z_4th_interior_plus(A_h, i, j, idz);
						DAz_A_c[IDX(i, j)] = z_4th_interior_plus(A_c, i, j, idz);
						DAz_A_lambda[IDX(i, j)] = z_4th_interior_plus(A_lambda, i, j, idz);
						DAz_beta_r[IDX(i, j)] = z_4th_interior_plus(beta_r, i, j, idz);
						DAz_beta_z[IDX(i, j)] = z_4th_interior_plus(beta_z, i, j, idz);
						DAz_dtbeta_r[IDX(i, j)] = z_4th_interior_plus(dtbeta_r, i, j, idz);
						DAz_dtbeta_z[IDX(i, j)] = z_4th_interior_plus(dtbeta_z, i, j, idz);
						DAz_Deltar[IDX(i, j)] = z_4th_interior_plus(Deltar, i, j, idz);
						DAz_Deltaz[IDX(i, j)] = z_4th_interior_plus(Deltaz, i, j, idz);
					}
					else
					{
						DAz_alpha[IDX(i, j)] = z_4th_interior_minus(alpha, i, j, idz);
						DAz_phi[IDX(i, j)] = z_4th_interior_minus(phi, i, j, idz);
						DAz_a[IDX(i, j)] = z_4th_interior_minus(a, i, j, idz);
						DAz_b[IDX(i, j)] = z_4th_interior_minus(b, i, j, idz);
						DAz_h[IDX(i, j)] = z_4th_interior_minus(h, i, j, idz);
						DAz_c[IDX(i, j)] = z_4th_interior_minus(c, i, j, idz);
						DAz_lambda[IDX(i, j)] = z_4th_interior_minus(lambda, i, j, idz);
						DAz_K[IDX(i, j)] = z_4th_interior_minus(K, i, j, idz);
						DAz_A_a[IDX(i, j)] = z_4th_interior_minus(A_a, i, j, idz);
						DAz_A_b[IDX(i, j)] = z_4th_interior_minus(A_b, i, j, idz);
						DAz_A_h[IDX(i, j)] = z_4th_interior_minus(A_h, i, j, idz);
						DAz_A_c[IDX(i, j)] = z_4th_interior_minus(A_c, i, j, idz);
						DAz_A_lambda[IDX(i, j)] = z_4th_interior_minus(A_lambda, i, j, idz);
						DAz_beta_r[IDX(i, j)] = z_4th_interior_minus(beta_r, i, j, idz);
						DAz_beta_z[IDX(i, j)] = z_4th_interior_minus(beta_z, i, j, idz);
						DAz_dtbeta_r[IDX(i, j)] = z_4th_interior_minus(dtbeta_r, i, j, idz);
						DAz_dtbeta_z[IDX(i, j)] = z_4th_interior_minus(dtbeta_z, i, j, idz);
						DAz_Deltar[IDX(i, j)] = z_4th_interior_minus(Deltar, i, j, idz);
						DAz_Deltaz[IDX(i, j)] = z_4th_interior_minus(Deltaz, i, j, idz);
					}
				}

				// Third-to-last point.
				j = NrTotal - 3;
				if (-beta_r[IDX(i, j)] >= 0.0)
				{
					DAz_alpha[IDX(i, j)] = z_4th_1edge_plus(alpha, i, j, idz);
					DAz_phi[IDX(i, j)] = z_4th_1edge_plus(phi, i, j, idz);
					DAz_a[IDX(i, j)] = z_4th_1edge_plus(a, i, j, idz);
					DAz_b[IDX(i, j)] = z_4th_1edge_plus(b, i, j, idz);
					DAz_h[IDX(i, j)] = z_4th_1edge_plus(h, i, j, idz);
					DAz_c[IDX(i, j)] = z_4th_1edge_plus(c, i, j, idz);
					DAz_lambda[IDX(i, j)] = z_4th_1edge_plus(lambda, i, j, idz);
					DAz_K[IDX(i, j)] = z_4th_1edge_plus(K, i, j, idz);
					DAz_A_a[IDX(i, j)] = z_4th_1edge_plus(A_a, i, j, idz);
					DAz_A_b[IDX(i, j)] = z_4th_1edge_plus(A_b, i, j, idz);
					DAz_A_h[IDX(i, j)] = z_4th_1edge_plus(A_h, i, j, idz);
					DAz_A_c[IDX(i, j)] = z_4th_1edge_plus(A_c, i, j, idz);
					DAz_A_lambda[IDX(i, j)] = z_4th_1edge_plus(A_lambda, i, j, idz);
					DAz_beta_r[IDX(i, j)] = z_4th_1edge_plus(beta_r, i, j, idz);
					DAz_beta_z[IDX(i, j)] = z_4th_1edge_plus(beta_z, i, j, idz);
					DAz_dtbeta_r[IDX(i, j)] = z_4th_1edge_plus(dtbeta_r, i, j, idz);
					DAz_dtbeta_z[IDX(i, j)] = z_4th_1edge_plus(dtbeta_z, i, j, idz);
					DAz_Deltar[IDX(i, j)] = z_4th_1edge_plus(Deltar, i, j, idz);
					DAz_Deltaz[IDX(i, j)] = z_4th_1edge_plus(Deltaz, i, j, idz);
				}
				else
				{
					DAz_alpha[IDX(i, j)] = z_4th_1edge_minus(alpha, i, j, idz);
					DAz_phi[IDX(i, j)] = z_4th_1edge_minus(phi, i, j, idz);
					DAz_a[IDX(i, j)] = z_4th_1edge_minus(a, i, j, idz);
					DAz_b[IDX(i, j)] = z_4th_1edge_minus(b, i, j, idz);
					DAz_h[IDX(i, j)] = z_4th_1edge_minus(h, i, j, idz);
					DAz_c[IDX(i, j)] = z_4th_1edge_minus(c, i, j, idz);
					DAz_lambda[IDX(i, j)] = z_4th_1edge_minus(lambda, i, j, idz);
					DAz_K[IDX(i, j)] = z_4th_1edge_minus(K, i, j, idz);
					DAz_A_a[IDX(i, j)] = z_4th_1edge_minus(A_a, i, j, idz);
					DAz_A_b[IDX(i, j)] = z_4th_1edge_minus(A_b, i, j, idz);
					DAz_A_h[IDX(i, j)] = z_4th_1edge_minus(A_h, i, j, idz);
					DAz_A_c[IDX(i, j)] = z_4th_1edge_minus(A_c, i, j, idz);
					DAz_A_lambda[IDX(i, j)] = z_4th_1edge_minus(A_lambda, i, j, idz);
					DAz_beta_r[IDX(i, j)] = z_4th_1edge_minus(beta_r, i, j, idz);
					DAz_beta_z[IDX(i, j)] = z_4th_1edge_minus(beta_z, i, j, idz);
					DAz_dtbeta_r[IDX(i, j)] = z_4th_1edge_minus(dtbeta_r, i, j, idz);
					DAz_dtbeta_z[IDX(i, j)] = z_4th_1edge_minus(dtbeta_z, i, j, idz);
					DAz_Deltar[IDX(i, j)] = z_4th_1edge_minus(Deltar, i, j, idz);
					DAz_Deltaz[IDX(i, j)] = z_4th_1edge_minus(Deltaz, i, j, idz);
				}

				// Second to last point.
				j = NrTotal - 2;
				DAz_alpha[IDX(i, j)] = z_4th_2edge(alpha, i, j, idz);
				DAz_phi[IDX(i, j)] = z_4th_2edge(phi, i, j, idz);
				DAz_a[IDX(i, j)] = z_4th_2edge(a, i, j, idz);
				DAz_b[IDX(i, j)] = z_4th_2edge(b, i, j, idz);
				DAz_h[IDX(i, j)] = z_4th_2edge(h, i, j, idz);
				DAz_c[IDX(i, j)] = z_4th_2edge(c, i, j, idz);
				DAz_lambda[IDX(i, j)] = z_4th_2edge(lambda, i, j, idz);
				DAz_K[IDX(i, j)] = z_4th_2edge(K, i, j, idz);
				DAz_A_a[IDX(i, j)] = z_4th_2edge(A_a, i, j, idz);
				DAz_A_b[IDX(i, j)] = z_4th_2edge(A_b, i, j, idz);
				DAz_A_h[IDX(i, j)] = z_4th_2edge(A_h, i, j, idz);
				DAz_A_c[IDX(i, j)] = z_4th_2edge(A_c, i, j, idz);
				DAz_A_lambda[IDX(i, j)] = z_4th_2edge(A_lambda, i, j, idz);
				DAz_beta_r[IDX(i, j)] = z_4th_2edge(beta_r, i, j, idz);
				DAz_beta_z[IDX(i, j)] = z_4th_2edge(beta_z, i, j, idz);
				DAz_dtbeta_r[IDX(i, j)] = z_4th_2edge(dtbeta_r, i, j, idz);
				DAz_dtbeta_z[IDX(i, j)] = z_4th_2edge(dtbeta_z, i, j, idz);
				DAz_Deltar[IDX(i, j)] = z_4th_2edge(Deltar, i, j, idz);
				DAz_Deltaz[IDX(i, j)] = z_4th_2edge(Deltaz, i, j, idz);

				// Last point.
				j = NrTotal - 1;
				DAz_alpha[IDX(i, j)] = z_4th_edge(alpha, i, j, idz);
				DAz_phi[IDX(i, j)] = z_4th_edge(phi, i, j, idz);
				DAz_a[IDX(i, j)] = z_4th_edge(a, i, j, idz);
				DAz_b[IDX(i, j)] = z_4th_edge(b, i, j, idz);
				DAz_h[IDX(i, j)] = z_4th_edge(h, i, j, idz);
				DAz_c[IDX(i, j)] = z_4th_edge(c, i, j, idz);
				DAz_lambda[IDX(i, j)] = z_4th_edge(lambda, i, j, idz);
				DAz_K[IDX(i, j)] = z_4th_edge(K, i, j, idz);
				DAz_A_a[IDX(i, j)] = z_4th_edge(A_a, i, j, idz);
				DAz_A_b[IDX(i, j)] = z_4th_edge(A_b, i, j, idz);
				DAz_A_h[IDX(i, j)] = z_4th_edge(A_h, i, j, idz);
				DAz_A_c[IDX(i, j)] = z_4th_edge(A_c, i, j, idz);
				DAz_A_lambda[IDX(i, j)] = z_4th_edge(A_lambda, i, j, idz);
				DAz_beta_r[IDX(i, j)] = z_4th_edge(beta_r, i, j, idz);
				DAz_beta_z[IDX(i, j)] = z_4th_edge(beta_z, i, j, idz);
				DAz_dtbeta_r[IDX(i, j)] = z_4th_edge(dtbeta_r, i, j, idz);
				DAz_dtbeta_z[IDX(i, j)] = z_4th_edge(dtbeta_z, i, j, idz);
				DAz_Deltar[IDX(i, j)] = z_4th_edge(Deltar, i, j, idz);
				DAz_Deltaz[IDX(i, j)] = z_4th_edge(Deltaz, i, j, idz);
			}
		}
	}

	// Axis symmetries.
	for (k = 0; k < ghost; k++)
	{
		#pragma omp parallel shared(DAz_alpha, DAz_phi, DAz_a, DAz_b, DAz_h, DAz_c, DAz_lambda,\
		DAz_K, DAz_A_a, DAz_A_b, DAz_A_h, DAz_A_c, DAz_A_lambda, DAz_beta_r, DAz_beta_z,\
		DAz_dtbeta_r, DAz_dtbeta_z, DAz_Deltar, DAz_Deltaz)
		{
			#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				DAz_alpha[IDX(i, ghost - 1 - k)] = -DAz_alpha[IDX(i, ghost + k)];
				DAz_phi[IDX(i, ghost - 1 - k)] = -DAz_phi[IDX(i, ghost + k)];
				DAz_a[IDX(i, ghost - 1 - k)] = -DAz_a[IDX(i, ghost + k)];
				DAz_b[IDX(i, ghost - 1 - k)] = -DAz_b[IDX(i, ghost + k)];
				DAz_h[IDX(i, ghost - 1 - k)] = -DAz_h[IDX(i, ghost + k)];
				DAz_c[IDX(i, ghost - 1 - k)] = DAz_c[IDX(i, ghost + k)];
				DAz_lambda[IDX(i, ghost - 1 - k)] = -DAz_lambda[IDX(i, ghost + k)];
				DAz_K[IDX(i, ghost - 1 - k)] = -DAz_K[IDX(i, ghost + k)];
				DAz_A_a[IDX(i, ghost - 1 - k)] = -DAz_A_a[IDX(i, ghost + k)];
				DAz_A_b[IDX(i, ghost - 1 - k)] = -DAz_A_b[IDX(i, ghost + k)];
				DAz_A_h[IDX(i, ghost - 1 - k)] = -DAz_A_h[IDX(i, ghost + k)];
				DAz_A_c[IDX(i, ghost - 1 - k)] = DAz_A_c[IDX(i, ghost + k)];
				DAz_A_lambda[IDX(i, ghost - 1 - k)] = -DAz_A_lambda[IDX(i, ghost + k)];
				DAz_beta_r[IDX(i, ghost - 1 - k)] = -DAz_beta_r[IDX(i, ghost + k)];
				DAz_beta_z[IDX(i, ghost - 1 - k)] = DAz_beta_z[IDX(i, ghost + k)];
				DAz_dtbeta_r[IDX(i, ghost - 1 - k)] = -DAz_dtbeta_r[IDX(i, ghost + k)];
				DAz_dtbeta_z[IDX(i, ghost - 1 - k)] = DAz_dtbeta_z[IDX(i, ghost + k)];
				DAz_Deltar[IDX(i, ghost - 1 - k)] = -DAz_Deltar[IDX(i, ghost + k)];
				DAz_Deltaz[IDX(i, ghost - 1 - k)] = DAz_Deltaz[IDX(i, ghost + k)];
			}
		}
	}

	return;
}
