#include <math.h>

// Fine-tuning parameter for lambda regularization where there are two ways to make the regularization.
#define FT1 0.5
#define FT2 0.5
#define FT3 0.5

// This subroutine calculates the Ricci tensor and scalar for a given metric in axisymmetry and zero angular momentum.
//
// Output parameters are the Ricci tensor components given by Ra, Rb, Rh, Rc and Rlambda (for regularization).
// The Ricci scalar is R (capital R).
//
// Input parameters are the radial grid r, the metric a, b, h, c, lambda; its inverse ga, gb, ...; 
// first and second derivatives for the metric; the Delta vector and its first derivatives, including
// the regularization derivative.
void calculate_ricci(double *Ra, double *Rb, double *Rh, double *Rc, double *Rlambda, double *R,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double ga, const double gb, const double gh, const double gc, const double glambda,
	const double dRa, const double dRb, const double dRh, const double dRc, const double dRlambda,
	const double dZa, const double dZb, const double dZh, const double dZc, const double dZlambda,
	const double dRRa, const double dRRb, const double dRRh, const double dRRc, const double dRRlambda,
	const double dZZa, const double dZZb, const double dZZh, const double dZZc, const double dZZlambda,
	const double dRZa, const double dRZb, const double dRZh, const double dRZc, const double dRZlambda,
	const double deltaR, const double deltaZ,
	const double dRdeltaR, const double dRdeltaZ,
	const double dZdeltaR, const double dZdeltaZ,
	const double dDRdeltaR)
{
	// Tensor components.
	*Ra = (2 * deltaR*dRa + 4 * a*dRdeltaR + 2 * deltaZ*dZa - 2 * dRRa*ga + 3 * (dRa*dRa)*(ga*ga) - 2 * dZZa*gb + 4 * (c*c)*ga*gb +
		2 * (dZa*dZa)*ga*gb - dRb * dRb*(gb*gb) - dRh * dRh*(gh*gh) + 4 * gh*lambda +
		4 * (dRb*dZc*(gb*gb) - dRZa * gc + dRa * dZa*ga*gc + dRb * dZa*gb*gc + dRh * (gh*gh)*lambda)*r +
		(4 * (dRc*dRc)*ga*gb + 8 * dRa*dRc*ga*gc - 2 * (-(dRa*dRb) + dZa * dZa)*(gc*gc))*(r*r) +
		4 * gc*(2 * dRc*dZc*gb + (2 * dRc*dZa + dRa * dZc)*gc)*(r*r*r) +
		4 * c*r*(dRdeltaZ + 2 * dRc*ga*gb + 2 * dZa*(gc*gc)*r + 2 * gc*(dRa*ga + dZc * gb*r))) / 4.
		// Regularization.
		+ (-0.5 * gh * (dRa / r));
	*Rb = (2 * deltaR*dRb + 2 * deltaZ*dZb + 4 * b*dZdeltaZ - 2 * dZZb*gb + 3 * (dZb*dZb)*(gb*gb) - dZh * dZh*(gh*gh) +
		4 * (-dRZb + dRb * dZb*gb)*gc*r + 2 * gc*(4 * dZb*dZc*gb + (-(dRb*dRb) + dZa * dZb)*gc)*(r*r) +
		4 * (dRc*dZb + 2 * dRb*dZc)*(gc*gc)*(r*r*r) - dZa * (ga*ga)*(dZa - 4 * dRc*r) +
		4 * c*(dZa*(ga*ga) + dZh * (gh*gh) + dZdeltaR * r + gc * (2 * dZc*ga + dZb * gc)*(r*r)) +
		2 * ga*(-dRRb + gb * (dRb*dRb + 2 * (dZc*dZc)*(r*r)) + 2 * gc*(dRb*dZa*r + 2 * dRc*dZc*(r*r*r)))) / 4.
		// Regularization.
		+ (-0.5 * gh * (dRb / r));
	*Rc = (2 * deltaR*dRc + 2 * deltaZ*dZc + 2 * dRa*dRc*(ga*ga) - 2 * dZZc*gb + 2 * dZb*dZc*(gb*gb) - 4 * dZc*gc + 4 * (c*c)*ga*gc +
		2 * (dRb*dRb)*gb*gc + 2 * dZa*dZb*gb*gc + 2 * dZh*(gh*gh)*lambda +
		gc * (-4 * dRZc + 4 * dRc*dZb*gb + (dRb*dZa + 3 * dRa*dZb)*gc)*r +
		2 * gc*(2 * (dZc*dZc)*gb + (dRb*dRc + dZa * dZc)*gc)*(r*r) + 4 * dRc*dZc*(gc*gc)*(r*r*r) +
		2 * c*(dRdeltaR + dZdeltaZ + 2 * dZb*gb*gc + 4 * dRc*ga*gc*r + gc * gc*r*(dRb + 2 * dZc*r)) +
		2 * ga*(-dRRc + 2 * (dRb*dRc + dZa * dZc)*gb + gc * (dRa*dRb + dZa * dZa + 2 * dRa*dZc*r + 2 * (dRc*dRc)*(r*r)))) / 4.
		// Regularization.
		+ 0.5 * c * (deltaR / r) + 0.5 * a * (dZdeltaR / r)
		+ 0.5 * c * ga * ga * (dRa / r) + 0.25 * ga * ga * dZa * (dRa / r)
		+ c * ga * gb * (dRb / r) - 0.5 * ga * gb * dZa * (dRb / r)
		+ 0.25 * gb * gb * dZb * (dRb / r) - ga * (dRc / r) - 0.5 * gh * (dRc / r)
		+ 0.5 * b * (dRdeltaZ / r) + 0.5 * c * gh * gh * (dRh / r)
		- 0.25 * gh * gh * dZh * (dRh / r);
	*Rh = (deltaR*dRh + deltaZ * dZh - dRRh * ga - dZZh * gb + 2 * (c*c)*gb*gh - 2 * dRZh*gc*r + 4 * c*gc*gh*lambda*(r*r) +
		gh * (dRh*dRh*ga + dZh * dZh*gb - 2 * lambda + 2 * dRh*dZh*gc*r + 2 * ga*(lambda*lambda)*(r*r))) / 2.
		// Regularization.
		+ h * (deltaR / r) - 0.5 * gh *(dRh / r);
	*Rlambda = 2 * c*gc*(dZc*gb + dZa * gc - gh * lambda) +
		ga * (dRc*(dRc*gb + 2 * dRa*gc) - gh * (lambda*lambda)) +
		(gc*(4 * dRc*dZc*gb*r + gc * (-(dZa*dZa) + 4 * dRc*dZa*r + dRa * (dRb + 2 * dZc*r)))) / 2.
		// Regularization.
		- (-2 * (deltaR / r) *((dRa / r) - (dRh / r)) - 3 * ((dRa / r)*(dRa / r))*(ga*ga) + (dRb / r) * (dRb / r) *(gb*gb) +
		2 * ((dRh / r)*(dRh / r))*ga*gh + (dRh / r) * (dRh / r)*(gh*gh)) / 4. 
		+ ((dRb / r)*dZc*(gb*gb) - (dRZa / r) * gc + (dRZh / r) * gc + (dRa / r) * dZa*ga*gc + (dRb / r) * dZa*gb*gc +
		c * ((dRdeltaZ / r) + 2 * ga*((dRc / r)*gb + (dRa / r) * gc)) - (dRh / r) * dZh*gc*gh + (dRh / r) * (gh*gh)*lambda)
		// Lambda Regularization.
		+ c*c*gb*glambda + 0.5 * deltaZ * dZlambda - 0.5 * gb * dZZlambda
		+ FT1 * (lambda * (deltaR / r) + a * (dDRdeltaR / r))
		+ (1.0 - FT1) * (lambda * dRdeltaR + h * (dDRdeltaR / r)) 
		+ FT2 * (0.5*gb*glambda *dZa * dZa + 0.5 * gb * gh * dZlambda * (dZa + dZh))
		+ (1.0 - FT2) * (0.5*gb*glambda*dZh*dZh + 0.5 *ga*gb * dZlambda * (dZa + dZh))
		+ FT3 * (0.5 * glambda * (dRRh - dRRa) - 0.5 * gh *(dRRlambda + 5.0 * (dRlambda / r)))
		+ (1.0 - FT3) * (0.5 * glambda * (r * dRlambda - 2.0 * lambda) - 0.5 * ga * (dRRlambda + 5.0 * (dRlambda / r)));

	// Scalar is simply the trace of the tensor components.
	*R = *Ra * ga + *Rb * gb + *Rh * gh + 2. * r * r * gc * *Rc;

	return;
}

// This subroutine corrects a previously calculated Ricci tensor and scalar by adding terms asociated to the conformal factor phi.
//
// Output parameters are the Ricci tensor and scalar as in the above subroutine.
//
// Input parameters are the radial grid; the **conformal** metric and its first derivatives; the inverse **conformal** metric;
// the conformal factor's first and second derivatives and regularization derivative; finally, psi4 = exp(4*phi) as the conformal factor.
void calculate_phys_ricci(double *Ra, double *Rb, double *Rh, double *Rc, double *Rlambda, double *R,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double ga, const double gb, const double gh, const double gc, const double glambda,
	const double dRa, const double dRb, const double dRh, const double dRc, const double dRlambda,
	const double dZa, const double dZb, const double dZh, const double dZc, const double dZlambda,
	const double dRphi, const double dZphi, const double dRRphi, const double dZZphi, const double dRZphi,
	const double ddRphi, const double psi4)
{
	double aux;

	double t_R_a, t_R_b, t_R_h, t_R_c;
	// Calculate terms related to the second covariant derivative of the conformal factor.
	//
	//    _^_ _^_
	// -2 \ / \ /  phi
	//      i   j
	//
	t_R_a = -2.0 * (dRRphi - (dZphi*(dRa*gc*r + gb*(2 * c - dZa + 2 * dRc*r))) / 2. - (dRphi*(dRa*ga + gc*r*(2 * c - dZa + 2 * dRc*r))) / 2.);
	t_R_b = -2.0 * (dZZphi + (dRphi*(-(dZb*gc*r) + ga*(dRb - 2 * dZc*r))) / 2. + (dZphi*(-(dZb*gb) + gc*r*(dRb - 2 * dZc*r))) / 2.);
	t_R_c = -2.0 * (-((dRphi*dRb + dZa*dZphi)*gc) / 2.
		// Regularization.
		- ga * dZa * (dRphi / r) / 2. - gb * dZphi * (dRb / r) / 2. + (dRZphi / r));
	t_R_h = -2.0 * ((dRphi*dRh*ga + dZphi*dZh*gb + gc*(dRphi*dZh*r + dZphi*(2 * h + dRh*r))) / 2.
		// Regularization.
		+ h * ga * (dRphi / r));

	// Add these terms.
	*Ra += t_R_a;
	*Rb += t_R_b;
	*Rh += t_R_h;
	*Rc += t_R_c;

	// Lambda term.
	*Rlambda += -2.0 * (-gc * dRphi * dRc
		- gc * dZphi * (dRa / r) / 2.
		- gb * dZphi * (dRc / r)
		- gc * dZphi * (dRh / r) / 2.
		+ (-c * gc + gc * dZa / 2. - gc * dZh / 2. - ga * (dRa / r) / 2. - ga * (dRh / r) / 2.) * (dRphi / r)
		+ (ddRphi / r) - glambda * h * (dRphi / r)
		+ gb * dZphi * dZlambda / 2.
		+ gc * lambda * dZphi);

	// Correct terms by taking the trace.
	//
	//    ^   _^_k _^_
	// -2 g   \ /  \ /  phi
	//     ij        k
	//
	aux = ga * t_R_a + gb * t_R_b + gh * t_R_h + 2. * r * r * gc * t_R_c;
	*Ra += aux * a;
	*Rb += aux * b;
	*Rh += aux * h;
	*Rc += aux * c;
	*Rlambda += aux * lambda;

	// Add other terms.
	//   _^_      _^_          ^   _^_k    _^_
	// 4 \ /  phi \ /  phi - 4 g   \ / phi \ / phi
	//     i        j           ij           k
	//
	aux = ga * dRphi * dRphi + gb * dZphi * dZphi + 2. * r * gc * dZphi * dRphi;
	*Ra += 4. * (dRphi * dRphi - a * aux);
	*Rb += 4. * (dZphi * dZphi - b * aux);
	*Rh += -4. * h * aux;
	*Rc += 4. * (dZphi * (dRphi / r) - c * aux);
	*Rlambda += 4.0 * ((dRphi / r) * (dRphi / r) - lambda * aux);

	// Now we can trace the Ricci tensor to get the Ricci scalar.
	// Remember that this is a physical trace. Thus, we divide by psi4.
	*R = (*Ra * ga + *Rb * gb + *Rh * gh + 2. * r * r * gc * *Rc) / psi4;

	return;
}
