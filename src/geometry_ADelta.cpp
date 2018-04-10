void calculate_ADelta(double *ADelta_r, double *ADelta_z,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double ga, const double gb, const double gh, const double gc,
	const double Aa, const double Ab, const double Ah, const double Ac,
	const double dRa, const double dRb, const double dRh, const double dRc,
	const double dZa, const double dZb, const double dZh, const double dZc)
{
	*ADelta_r = ((Ah*(2 * c - dZh)*gc*(gh*gh)*r + Ah * ga*(gh*gh)*(-dRh + 2 * lambda*r) +
		dZb * gc*r*(Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r)) +
		ga * (-dRb + 2 * dZc*r)*(Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r)) +
		dRa * ga*(Aa*(ga*ga) + gc * (2 * Ac*ga + Ab * gc)*(r*r)) +
		gc * r*(2 * c - dZa + 2 * dRc*r)*(Aa*(ga*ga) + gc * (2 * Ac*ga + Ab * gc)*(r*r)) +
		2 * dZa*ga*r*((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r))) +
		2 * dRb*gc*(r*r)*((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r)))) / 2.);

	*ADelta_z = ((Ah*(2 * c - dZh)*gb*(gh*gh) + Ah * gc*(gh*gh)*r*(-dRh + 2 * lambda*r) +
		dZb * gb*(Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r)) +
		gc * r*(-dRb + 2 * dZc*r)*(Ab*(gb*gb) + gc * (2 * Ac*gb + Aa * gc)*(r*r)) +
		dRa * gc*r*(Aa*(ga*ga) + gc * (2 * Ac*ga + Ab * gc)*(r*r)) +
		gb * (2 * c - dZa + 2 * dRc*r)*(Aa*(ga*ga) + gc * (2 * Ac*ga + Ab * gc)*(r*r)) +
		2 * dRb*gb*r*((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r))) +
		2 * dZa*gc*(r*r)*((Aa*ga + Ab * gb)*gc + Ac * (ga*gb + gc * gc*(r*r)))) / 2.);

	return;
}
