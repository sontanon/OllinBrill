#include <math.h>

#define FT 0.5

void calculate_A2(double *A2_a, double *A2_b, double *A2_h, double *A2_c, double *A2_lambda, double *A2,
	const double r, const double A_a, const double A_b, const double A_h, const double A_c, const double A_lambda,
	const double g_a, const double g_b, const double g_h, const double g_c, const double g_lambda)
{
	*A2_a = A_a * A_a * g_a + r * r * A_c * (A_c * g_b + 2. * A_a * g_c);
	*A2_b = A_b * A_b * g_b + r * r * A_c * (A_c * g_a + 2. * A_b * g_c);
	*A2_c = A_a * A_c * g_a + A_b * A_c * g_b + A_a * A_b * g_c + r * r * A_c * A_c * g_c;
	*A2_h = A_h * A_h * g_h;
	*A2_lambda = A_c * (A_c * g_b + 2.0 * A_a * g_c)
		// Lambda regularization.
		+ FT * (A_a * A_a * g_lambda + A_lambda * (A_a + A_h) * g_h)
		+ (1.0 - FT) * (A_h * A_h * g_lambda + A_lambda * (A_a + A_h) * g_a);

	// Trace.
	*A2 = g_a * *A2_a + g_b * *A2_b + g_h * *A2_h + 2. * r * r * g_c * *A2_c;

	return;
}
