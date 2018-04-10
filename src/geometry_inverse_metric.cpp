#include <math.h>

void calculate_inverse_metric(double *g_a, double *g_b, double *g_h, double *g_c, double *g_lambda, double *hdet,
	const double r, const double a, const double b, const double h, const double c, const double lambda)
{
	*g_a = b * h;
	*g_b = a * h;
	*g_h = 1. / h;
	*g_c = -c * h;

	*g_lambda = c * c - b * lambda;

	*hdet = 1.;
	
	return;
}
