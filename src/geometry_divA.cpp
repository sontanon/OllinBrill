#include <math.h>

#define FT 0.5

void calculate_divA(double *divA_r, double *divA_z,
	const double r, const double a, const double b, const double h, const double c,
	const double Aa, const double Ab, const double Ah, const double Ac, const double Alambda,
	const double ga, const double gb, const double gh, const double gc, const double glambda,
	const double dRa, const double dRb, const double dRh, const double dRc,
	const double dZa, const double dZb, const double dZh, const double dZc,
	const double dRAa, const double dRAb, const double dRAh, const double dRAc,
	const double dZAa, const double dZAb, const double dZAh, const double dZAc)
{
	*divA_r = ((-4 * Ab*dRh - Ab * dRb*ga*(gb*gb) - Ah * dRh*ga*(gh*gh) - 2 * dRAb*h + 4 * a*Ab*b*dRh*h +
		2 * (b*b)*dRAa*(h*h) + 2 * a*b*dRAb*(h*h) +
		(-(Ah*dZh*gc*(gh*gh)) + Ab * gb*
		(ga*(2 * dZc*gb + 3 * dZa*gc) + gc * (2 * dZb*gb + dZh * gh)))*r -
		4 * a*Ab*c*dZh*h*r + (-2 * dZAc + 2 * Ab*(gc*gc)*gh)*h*r + 4 * Ab*(c*c)*(h*h)*r -
		2 * Ab*c*dZa*(h*h)*r - 2 * b*c*dZAa*(h*h)*r - 2 * a*c*dZAb*(h*h)*r +
		4 * a*b*dZAc*(h*h)*r - 2 * a*Ab*dZc*(h*h)*r +
		Ab * (gc*gc)*(2 * dRa*ga + 3 * dRb*gb + dRh * gh)*(r*r) - 4 * b*c*dRAc*(h*h)*(r*r) +
		4 * Ab*c*dRc*(h*h)*(r*r) + 4 * Ab*c*(gc*gc*gc)*(r*r*r) +
		Ab * (gc*gc)*(2 * dZc*gb - dZa * gc)*(r*r*r) + 4 * Ab*dRc*(gc*gc*gc)*(r*r*r*r) +
		Aa * (2 * dRa*(ga*ga*ga) + 4 * (b*b)*dRh*h + 4 * b*dRb*(h*h) - 2 * b*dZc*(h*h)*r +
			2 * c*(2 * (ga*ga)*gc - h * (2 * b*dZh + dZb * h))*r + dZb * (gc*gc*gc)*(r*r*r) +
			ga * ga*(dRb*gb + dRh * gh + 2 * gc*r*(dZa + 2 * dRc*r)) +
			ga * gc*r*(dZb*gb + dZh * gh + gc * r*(dRb + 4 * dZc*r))) +
		Ac * r*(-4 * dZh + 8 * a*b*dZh*h + 2 * b*dZa*(h*h) + 2 * a*dZb*(h*h) -
			4 * b*dRc*(h*h)*r + gc * gc*(3 * dZb*gb + dZh * gh)*(r*r) +
			2 * dRb*(gc*gc*gc)*(r*r*r) + 2 * dZc*(gc*gc*gc)*(r*r*r*r) +
			ga * ga*(3 * dZa*gb + 4 * dRa*gc*r) +
			4 * c*(2 * ga*(gc*gc)*(r*r) - h * (2 * b*(h + dRh * r) + h * r*(dRb - dZc * r))) +
			ga * (gb*(dZb*gb + dZh * gh) + gc * gc*(r*r)*(dZa + 8 * dRc*r) +
				2 * gc*(gh*(2 * h + dRh * r) + gb * r*(dRb + 3 * dZc*r))))) / 2.
		// Regularization.
		+ FT * ga * r * (glambda * Aa + gh * Alambda)
		+ (1.0 - FT) * ga * r * (ga * Alambda + glambda * Ah));

	*divA_z = ((-2 * h*((a*Ab + Aa * b)*c*h + Ac * (1 - 2 * a*b*h)) -
		4 * dRh*((a*Ab + Aa * b)*c*h + Ac * (1 - 2 * a*b*h))*r +
		4 * dZh*h*(a*a*Ab - 2 * a*Ac*c*(r*r) + Aa * (c*c)*(r*r)) +
		dZh * gh*(Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r)) +
		(dZa*ga + dRb * gc*r)*(Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r)) +
		(Aa*(ga*ga) + gc * (2 * Ac*ga + Ab * gc)*(r*r))*
		(dRa*gc*r + gb * (2 * c - dZa + 2 * dRc*r)) -
		Ah * (gh*gh)*(dZh*gb + gc * (2 * h + dRh * r)) +
		2 * (Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r))*
		(dZb*gb + gc * r*(-dRb + 2 * dZc*r)) +
		2 * (h*h)*(a*a*dZAb + c * (-2 * Ac*dZa + c * dZAa + 2 * Aa*dZc)*(r*r) +
			2 * a*(Ab*dZa - (c*dZAc + Ac * dZc)*(r*r))) +
		gh * (2 * h + dRh * r)*((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r))) +
		3 * r*(dRb*gb + dZa * gc*r)*((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r))) +
		r * (dRa*ga + gc * r*(2 * c - dZa + 2 * dRc*r))*
		((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r))) +
		2 * (h*h)*r*((Ac*b - Ab * c)*dRa + a * (-(c*dRAb) + b * dRAc + Ac * dRb - Ab * dRc) +
			dRc * (-(Aa*b) + Ac * c*(r*r)) +
			c * (-(b*dRAa) - Aa * dRb + c * dRAc*(r*r) + Ac * r*(2 * c + dRc * r)))) / 2.);

	return;
}
