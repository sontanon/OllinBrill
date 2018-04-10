#include "tools.h"

void calculate_maximal_slicing_elliptic_coefficients(double *ell_a, double *ell_b, double *ell_c, double *ell_d, double *ell_e, double *ell_s, double *ell_f,
	const double r, const double a, const double b, const double h, const double c,
	const double dRa, const double dRb, const double dRh, const double dRc,
	const double dZa, const double dZb, const double dZh, const double dZc,
	const double dRphi, const double dZphi, const double psi4, const double R)
{
	double aux = a * b - r * r * c * c;

	*ell_a = r * (b / aux) / psi4;
	*ell_b = r * (2.0 * r * (-c / aux)) / psi4;
	*ell_c = r * (a / aux) / psi4;
	*ell_d = (a*(h*h)*(2 * (b*b) + c*dZb*(r*r) - b*r*(dRb + 2 * dZc*r)) + r*(2 * dRb*h - b*b*dRa*(h*h) - c*(dZh + 4 * dZphi*h)*r + b*(dRh + h*(4 * dRphi + c*h*r*(dZa + 2 * dRc*r))))) /
		(2. * psi4);
	*ell_e = (r*(a*dZh + 2 * (dZa + 2 * a*dZphi)*h - a*(h*h)*(a*dZb + b*(dZa + 2 * dRc*r)) + c*(-(dRh*r) - 2 * h*(1 + 2 * dRphi*r) + h*h*(b*dRa*r + a*(-2 * b + r*(dRb + 2 * dZc*r)))))) /
		(2. * psi4);
	*ell_s = -r * R;
	*ell_f = 0.0;
	return;
}
