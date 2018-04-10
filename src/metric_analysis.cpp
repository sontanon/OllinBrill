#include "tools.h"

void spherical_metric(double *g_rr, double *g_thth, double *g_phph, double *g_rth,
	double *ig_rr, double *ig_thth, double *ig_phph, double *ig_rth,
	const double r, const double z, const double a, const double b, const double h, const double c)
{
	double rr2 = r * r + z * z;
	double rr = sqrt(rr2);
	double r2 = r * r;

	double cos_th = z / rr;
	double sin_th = r / rr;

	double cos_th2 = cos_th * cos_th;
	double sin_th2 = sin_th * sin_th;

	double sin_2th = 2.0 * cos_th * sin_th;
	double cos_2th = cos_th2 - sin_th2;

	*g_rr = b * cos_th2 + a * sin_th2 + c * r * sin_2th;
	*g_thth = (a * cos_th2 - 2.0 * c * r * cos_th * sin_th + b * sin_th2) * rr2;
	*g_phph = h * r2;
	*g_rth = (c * r * cos_2th + (a - b) * cos_th * sin_th) * rr;

	double det = -(*g_rth) * (*g_rth) + (*g_rr) * (*g_thth);

	*ig_rr = *g_thth / det;
	*ig_thth = *g_rr / det;
	*ig_rth = -(*g_rth) / det;
	*ig_phph = 1.0 / (*g_phph);
}

#include "param.h"
#include "arrays.h"

FILE *f_c_rr;
FILE *f_c_th;
FILE *f_c_ph;

void metric_analysis(void)
{
	int k = 0;

	const double courant = 1.0 / sqrt(3.0);

#pragma omp parallel shared(phys_a, phys_b, phys_h, phys_c, phys_lambda, \
	g_rr, g_thth, g_phph, g_rth, ig_rr, ig_thth, ig_phph, ig_rth, c_rr, c_th, c_ph,\
	courant_rr, courant_th, courant_ph)
	{
#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Calculate physical metric.
			phys_a[k] = psi4[k] * a[k];
			phys_b[k] = psi4[k] * b[k];
			phys_h[k] = psi4[k] * h[k];
			phys_c[k] = psi4[k] * c[k];
			phys_lambda[k] = psi4[k] * lambda[k];

			// Calculate spherical metric.
			spherical_metric(g_rr + k, g_thth + k, g_phph + k, g_rth + k,
				ig_rr + k, ig_thth + k, ig_phph + k, ig_rth + k,
				r[k], z[k], phys_a[k], phys_b[k], phys_h[k], phys_c[k]);

			// Calculate light speed.
			/*
			c_rr[k] = alpha[k] * sqrt(ig_rr[k]);
			c_th[k] = alpha[k] * sqrt(ig_thth[k]);
			c_ph[k] = alpha[k] * sqrt(ig_phph[k]);
			*/
			c_rr[k] = alpha[k] / sqrt(fabs(g_rr[k]));
			c_th[k] = alpha[k] / sqrt(fabs(g_thth[k]));
			c_ph[k] = alpha[k] / sqrt(fabs(g_phph[k]));

			// Check Courant condition.
			/*
			courant_rr[k] = courant / fabs(c_rr[k]);
			courant_th[k] = courant / fabs(c_th[k]);
			courant_ph[k] = courant / fabs(c_ph[k]);
			*/
			courant_rr[k] = courant * sqrt(fabs(g_rr[k])) / alpha[k];
			courant_th[k] = courant * sqrt(fabs(g_thth[k])) / alpha[k];
			courant_ph[k] = courant * sqrt(fabs(g_phph[k])) / alpha[k];
		}
	}

	// Check Courant condition.
	double dtfac_rr = fabs(courant_rr[cblas_idamin(DIM, courant_rr, 1)]);
	double dtfac_th = fabs(courant_th[cblas_idamin(DIM, courant_th, 1)]);
	double dtfac_ph = fabs(courant_ph[cblas_idamin(DIM, courant_ph, 1)]);

	// Open file at zero time.
	if (t == 0.0)
	{
		f_c_rr = fopen("c_rr.asc", "w");
		f_c_th = fopen("c_th.asc", "w");
		f_c_ph = fopen("c_ph.asc", "w");
	}

	// Select power of two Courant factor.
	if (1.0 <= dtfac_rr)
	{
		dtfac_new = 0.5;
		CFL_n_new = 1;
	}
	else if (0.5 <= dtfac_rr)
	{
		dtfac_new = 0.25;
		CFL_n_new = 2;
	}
	else if (0.25 <= dtfac_rr)
	{
		dtfac_new = 0.125;
		CFL_n_new = 3;
	}
	else if (0.125 <= dtfac_rr)
	{
		dtfac_new = 0.0625;
		CFL_n_new = 4;
	}
	else
	{
		printf("\nERROR: COURANT VIOLATION PROBABLY TAKING PLACE: dtfac = %3.3E.\n", dtfac_rr);
		dtfac_new = 0.0625;
		CFL_n_new = 4;
	}

	// Write to files.
	fprintf(f_c_rr, "%9.18E\t%9.18E\n", t, dtfac_rr);
	fprintf(f_c_th, "%9.18E\t%9.18E\n", t, dtfac_th);
	fprintf(f_c_ph, "%9.18E\t%9.18E\n", t, dtfac_ph);

	// Close files at final time.
	if (t >= FinalTime)
	{
		fclose(f_c_rr);
		fclose(f_c_th);
		fclose(f_c_ph);
	}
}
