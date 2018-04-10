#include <math.h>

// This subroutine calculates the second covariant derivative of the lapse alpha for a given axisymmetric
// conformal metric.
//
// ___ ___          2           ___ k
// \ / \ / alpha = d   alpha - |       d alpha
//    i   j         ij              ij  k
//
// This subroutine also calculates the Laplacian by taking the trace of the previous tensor.
void calculate_D2alpha(double *D2alpha_a, double *D2alpha_b, double *D2alpha_h, double *D2alpha_c,
	double *D2alpha_lambda, double *LaplaAlpha,
	const double r, const double dRalpha, const double dZalpha,
	const double dRRalpha, const double dZZalpha, const double dRZalpha, const double ddRalpha,
	const double ga, const double gb, const double gh, const double gc, const double glambda,
	const double dRa, const double dRb, const double dRh, const double dRc,
	const double dZa, const double dZb, const double dZh, const double dZc,
	const double c, const double h, const double lambda, const double dZlambda)
{
	*D2alpha_a = dRRalpha - (dZalpha*(dRa*gc*r + gb*(2 * c - dZa + 2 * dRc*r))) / 2. - (dRalpha*(dRa*ga + gc*r*(2 * c - dZa + 2 * dRc*r))) / 2.;
	*D2alpha_b = dZZalpha + (dRalpha*(-(dZb*gc*r) + ga*(dRb - 2 * dZc*r))) / 2. + (dZalpha*(-(dZb*gb) + gc*r*(dRb - 2 * dZc*r))) / 2.;
	*D2alpha_c = -((dRalpha*dRb + dZa*dZalpha)*gc) / 2.
		// Regularization.
		- ga * dZa * (dRalpha / r) / 2. - gb * dZalpha * (dRb / r) / 2. + (dRZalpha / r);
	*D2alpha_h = ((dRalpha*dRh*ga + dZalpha*dZh*gb + gc*(dRalpha*dZh*r + dZalpha*(2 * h + dRh*r))) / 2.
		// Regularization.
		+ h * ga * (dRalpha / r));

	// Lambda term.
	*D2alpha_lambda = (-gc * dRalpha * dRc
		- gc * dZalpha * (dRa / r) / 2.
		- gb * dZalpha * (dRc / r)
		- gc * dZalpha * (dRh / r) / 2.
		+ (-c * gc + gc * dZa / 2. - gc * dZh / 2. - ga * (dRa / r) / 2. - ga * (dRh / r) / 2.) * (dRalpha / r)
		+ (ddRalpha / r) - glambda * h * (dRalpha / r)
		+ gb * dZalpha * dZlambda / 2.
		+ gc * lambda * dZalpha);

	// Laplacian.
	*LaplaAlpha = ga * *D2alpha_a + gb * *D2alpha_b + gh * *D2alpha_h + 2.0 * r * r * gc * *D2alpha_c;
	return;
}

// This subroutine takes the second conformal covariant derivative of the lapse and 
// completes the terms required for the physical connection in accordance to
//
// ___ ___         _^_ _^_                                  ^   ^kl
// \ / \ / alpha = \ / \ /  alpha - 4 (d phi) (d alpha) + 2 g   g   (d alpha) (d phi)
//    i   j           i   j             (i       j)          ij       l         k
//

void calculate_phys_D2alpha(double *D2alpha_a, double *D2alpha_b, double *D2alpha_h,
	double *D2alpha_c, double *D2alpha_lambda, double *LaplaAlpha,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double ga, const double gb, const double gh, const double gc,
	const double dRalpha, const double dZalpha, const double dRphi, const double dZphi,
	const double psi4)
{
	double t_a, t_b, t_c, t_lambda;
	double aux;

	t_a = dRalpha * dRphi;
	t_b = dZalpha * dZphi;
	t_c = 0.5 * (dZalpha * (dRphi / r) + dZphi * (dRalpha / r));
	// t_h = 0.0;
	t_lambda = (dRalpha / r) * (dRphi / r);

	aux = ga * t_a + gb * t_b + 2.0 * r * r * gc * t_c;

	*D2alpha_a += 2.0 * (-2.0 * t_a + a * aux);
	*D2alpha_b += 2.0 * (-2.0 * t_b + b * aux);
	*D2alpha_h += 2.0 * h * aux;
	*D2alpha_c += 2.0 * (-2.0 * t_c + c * aux);
	*D2alpha_lambda += 2.0 * (-2.0 * t_lambda + lambda * aux);

	// Physical trace.
	*LaplaAlpha = (ga * *D2alpha_a + gb * *D2alpha_b + gh * *D2alpha_h + 2.0 * r * r * gc * *D2alpha_c) / psi4;

	return;
}
