void calculate_upA(double *iA_a, double *iA_b, double *iA_h, double *iA_c,
	const double r, const double g_a, const double g_b, const double g_h, const double g_c,
	const double A_a, const double A_b, const double A_h, const double A_c)
{
	*iA_a = A_a * g_a * g_a + g_c * (2. * A_c * g_a + A_b * g_c) * r * r;
	*iA_b = A_b * g_b * g_b + g_c * (2. * A_c * g_b + A_a * g_c) * r * r;
	*iA_h = A_h * g_h * g_h;
	*iA_c = A_a * g_a * g_c + A_b * g_b * g_c + (g_a * g_b + r * r * g_c * g_c) * A_c;

	return;
}