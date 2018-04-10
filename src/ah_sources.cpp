#include <math.h>

void ah_sources(double *sH, double *sdH,
               const double th, const double rr, const double dH,
               const double a, const double b, const double h, const double c, const double phi,
               const double Aa, const double Ab, const double Ah, const double Ac, const double K,
               const double Dr_a, const double Dr_b, const double Dr_h, const double Dr_c, const double Dr_phi,
               const double Dz_a, const double Dz_b, const double Dz_h, const double Dz_c, const double Dz_phi)
{
    // Common constants.
    double sin_th = sin(th);
    double cos_th = cos(th);
    double rr2 = pow(rr,2);
    double psi4 = exp(4*phi);

    // Calculate physical extrinsic curvature from
    // its traceless conformal counterpart.
    const double third = 1./3.;
    double Ka = psi4*(Aa + third*a*K);
    double Kb = psi4*(Ab + third*b*K);
    double Kh = psi4*(Ah + third*h*K);
    double Kc = psi4*(Ac + third*c*K);

    // Auxiliary inverse determinant.
    double ihdet = 1./((a * b - pow(rr * sin_th * c, 2)) * h);

    // Source for H is its derivative.
    *sH = dH;

    // Calculate u2 and u.
    double u2 = (h*ihdet*((b*pow(dH,2) + a*rr2)*pow(cos_th,2) + sin_th*(2*c*dH*rr2*cos(2*th) + (a*pow(dH,2) + b*rr2)*sin_th + 2*rr*cos_th*((a - b)*dH + c*(pow(dH,2) - rr2)*sin_th))))/(psi4*rr2);

    double u = sqrt(u2);

    // Term proportional to dH/sin_th.
    double aux = -psi4*u2*cos_th*b;

    // Source for H's derivative.
    // Add normalizable terms first.
    *sdH = -(psi4*u2*(h*ihdet*cos_th*(rr*(b*dH*(Dr_h + 4*Dr_phi*h) + 2*c*h*rr - a*(Dz_h + 4*Dz_phi*h)*rr) + c*rr2*(-(dH*(Dz_h + 4*Dz_phi*h)) + (Dr_h + 4*Dr_phi*h)*rr)*sin_th) +
        rr*(-2*b*pow(h,2)*ihdet + 2*Kh*rr*u - h*ihdet*(-2*c*dH*h + a*dH*(Dz_h + 4*Dz_phi*h) + b*(Dr_h + 4*Dr_phi*h)*rr)*sin_th + c*h*ihdet*rr*(dH*(Dr_h + 4*Dr_phi*h) + (Dz_h + 4*Dz_phi*h)*rr)*pow(sin_th,2))))/
   (2.*pow(h,2)*ihdet);

    // Add remaining terms.
    *sdH += -(4*(-2*psi4*rr2*u2*((a - b)*cos_th + c*rr*cos(2*th))*sin_th + h*ihdet*(b*dH*pow(cos_th,2) + sin_th*(c*rr2*cos(2*th) + a*dH*sin_th) + rr*cos_th*sin_th*(a - b + 2*c*dH*sin_th))*
          (2*a*rr*pow(cos_th,2) + 2*sin_th*((a - b)*dH*cos_th + c*dH*rr*cos(2*th) + rr*(b - 2*c*rr*cos_th)*sin_th)))*
       (rr*u*((Ka - Kb)*cos_th + Kc*rr*cos(2*th))*sin_th + (h*ihdet*rr*(-(b*(Dz_a + 4*a*Dz_phi)*pow(sin_th,3)) + c*(Dr_b + 4*b*Dr_phi)*rr*pow(sin_th,4) +
              pow(cos_th,3)*(a*(Dr_b + 4*b*Dr_phi) - c*(Dz_a + 4*a*Dz_phi)*rr*sin_th) +
              cos_th*pow(sin_th,2)*(b*Dr_a - a*Dr_b + b*Dr_b + 4*pow(b,2)*Dr_phi - c*Dr_c*rr2 - 4*pow(c,2)*Dr_phi*rr2 + c*(Dr_c + 4*c*Dr_phi)*rr2*cos(2*th) -
                 (2*pow(c,2) + 2*b*Dz_c - c*(2*Dz_a + Dz_b + 8*a*Dz_phi - 4*b*Dz_phi))*rr*sin_th) -
              pow(cos_th,2)*sin_th*(-(b*Dz_a) + a*(-2*c + Dz_a + Dz_b) + 4*pow(a,2)*Dz_phi + (c*(Dr_a + 2*Dr_b) - 2*a*Dr_c)*rr*sin_th - 2*c*(Dz_c + 4*c*Dz_phi)*rr2*pow(sin_th,2)) +
              (a - 2*b)*c*Dr_phi*rr*pow(sin(2*th),2)))/2. - (dH*(2*pow(sin_th,2) + a*(Dr_b + 4*b*Dr_phi)*h*ihdet*rr*pow(sin_th,3) - c*(Dz_a + 4*a*Dz_phi)*h*ihdet*rr2*pow(sin_th,4) -
              h*ihdet*rr*pow(cos_th,3)*(-(b*(Dz_a + 4*a*Dz_phi)) + c*(Dr_b + 4*b*Dr_phi)*rr*sin_th) +
              h*ihdet*rr*cos_th*pow(sin_th,2)*(-2*a*c + a*Dz_a - b*Dz_a + a*Dz_b + 4*pow(a,2)*Dz_phi - c*Dz_c*rr2 - 4*pow(c,2)*Dz_phi*rr2 + c*(Dz_c + 4*c*Dz_phi)*rr2*cos(2*th) +
                 (-2*a*Dr_c + c*(Dr_a + 2*Dr_b - 4*a*Dr_phi + 8*b*Dr_phi))*rr*sin_th) + pow(cos_th,2)*
               (2 + (-(a*Dr_b) + b*(Dr_a + Dr_b) + 4*pow(b,2)*Dr_phi)*h*ihdet*rr*sin_th - (2*pow(c,2) - c*(2*Dz_a + Dz_b) + 2*b*Dz_c)*h*ihdet*rr2*pow(sin_th,2) -
                 2*c*(Dr_c + 4*c*Dr_phi)*h*ihdet*pow(rr,3)*pow(sin_th,3)) + (2*a - b)*c*Dz_phi*h*ihdet*rr2*pow(sin(2*th),2)))/(2.*rr)) +
      2*(psi4*rr2*u2*(b*pow(cos_th,2) + (a + 2*c*rr*cos_th)*pow(sin_th,2)) -
         h*ihdet*pow(b*dH*pow(cos_th,2) + sin_th*(c*rr2*cos(2*th) + a*dH*sin_th) + rr*cos_th*sin_th*(a - b + 2*c*dH*sin_th),2))*
       (-2*pow(sin_th,2) - b*(Dr_b + 4*b*Dr_phi)*h*ihdet*rr*pow(sin_th,3) + (-(c*Dz_b) + 2*b*Dz_c + 4*b*c*Dz_phi)*h*ihdet*rr2*pow(sin_th,4) +
         h*ihdet*rr*cos_th*pow(sin_th,2)*(-2*b*Dz_a + a*Dz_b - 4*a*b*Dz_phi - c*Dz_c*rr2 - 4*pow(c,2)*Dz_phi*rr2 + c*(Dz_c + 4*c*Dz_phi)*rr2*cos(2*th) + 3*c*(Dr_b + 4*b*Dr_phi)*rr*sin_th) +
         h*ihdet*rr*pow(cos_th,3)*(-(a*(-2*c + Dz_a + 4*a*Dz_phi)) + (-(c*Dr_a) + 2*a*Dr_c + 4*a*c*Dr_phi)*rr*sin_th) + 2*rr*u*(Ka*pow(cos_th,2) + Kb*pow(sin_th,2) - 2*Kc*rr*cos_th*pow(sin_th,2)) -
         pow(cos_th,2)*(2 + (-(b*Dr_a) + 2*a*Dr_b + 4*a*b*Dr_phi)*h*ihdet*rr*sin_th + c*(2*c - 3*Dz_a)*h*ihdet*rr2*pow(sin_th,2) + 2*c*(Dr_c + 4*c*Dr_phi)*h*ihdet*pow(rr,3)*pow(sin_th,3)) +
         3*a*c*Dz_phi*h*ihdet*rr2*pow(sin(2*th),2) + dH*h*ihdet*(a*(Dz_b + 4*b*Dz_phi)*pow(sin_th,3) + c*(Dr_b + 4*b*Dr_phi)*rr*pow(sin_th,4) - 2*c*(Dz_c + 4*c*Dz_phi)*rr2*pow(sin_th,5) +
            pow(cos_th,2)*sin_th*(2*b*Dz_a - 4*pow(a,2)*Dz_phi + a*(2*c - Dz_a + 8*b*Dz_phi) - (c*(Dr_a + 2*Dr_b) - 2*a*Dr_c)*rr*sin_th) +
            cos_th*pow(sin_th,2)*(-((2*a - b)*(Dr_b + 4*b*Dr_phi)) + (-2*b*Dz_c + c*(2*Dz_a + Dz_b + 8*a*Dz_phi - 4*b*Dz_phi))*rr*sin_th) -
            pow(cos_th,3)*(b*(Dr_a + 4*a*Dr_phi) - c*(2*c - Dz_a - 4*a*Dz_phi)*rr*sin_th - 2*c*Dr_c*rr2*pow(sin_th,2)) + c*Dr_phi*rr*(a - 2*b + 2*c*rr*cos_th)*pow(sin(2*th),2))) +
      pow(rr,3)*(4*psi4*u2*(a*pow(cos_th,2) + (b - 2*c*rr*cos_th)*pow(sin_th,2)) -
         (h*ihdet*pow(2*a*rr*pow(cos_th,2) + 2*sin_th*((a - b)*dH*cos_th + c*dH*rr*cos(2*th) + rr*(b - 2*c*rr*cos_th)*sin_th),2))/rr2)*
       (u*(Kb*pow(cos_th,2) + (Ka + 2*Kc*rr*cos_th)*pow(sin_th,2)) + (dH*h*ihdet*
            (-(a*(-2*c + Dz_a + 4*a*Dz_phi)*pow(sin_th,3)) + (-(c*Dr_a) + 2*a*Dr_c + 4*a*c*Dr_phi)*rr*pow(sin_th,4) + pow(cos_th,3)*(b*(Dr_b + 4*b*Dr_phi) + (-2*b*Dz_c + c*(Dz_b - 4*b*Dz_phi))*rr*sin_th) +
              cos_th*pow(sin_th,2)*(-(b*(Dr_a - 4*a*Dr_phi)) + c*(2*c - 3*(Dz_a + 4*a*Dz_phi))*rr*sin_th + 2*c*(Dr_c + 4*c*Dr_phi)*rr2*pow(sin_th,2)) +
              pow(cos_th,2)*sin_th*(a*Dz_b - 2*b*(Dz_a + 2*a*Dz_phi) + 3*c*Dr_b*rr*sin_th - 2*c*(Dz_c + 4*c*Dz_phi)*rr2*pow(sin_th,2)) + a*Dr_b*sin_th*sin(2*th) + 3*b*c*Dr_phi*rr*pow(sin(2*th),2)))/(2.*rr)
           - (h*ihdet*(-(b*(Dr_a + 4*a*Dr_phi)*pow(sin_th,3)) + c*(2*c - Dz_a - 4*a*Dz_phi)*rr*pow(sin_th,4) + 2*c*(Dr_c + 4*c*Dr_phi)*rr2*pow(sin_th,5) +
              cos_th*pow(sin_th,2)*(-2*a*c + a*Dz_a - 2*b*Dz_a + 4*pow(a,2)*Dz_phi - 8*a*b*Dz_phi + 2*c*Dz_a*rr*cos_th + (-2*a*Dr_c + c*(Dr_a + 2*Dr_b - 4*a*Dr_phi + 8*b*Dr_phi))*rr*sin_th) -
              pow(cos_th,2)*sin_th*((2*a - b)*(Dr_b + 4*b*Dr_phi) - (c*Dz_b - 2*b*Dz_c)*rr*sin_th) - pow(cos_th,3)*(a*(Dz_b + 4*b*Dz_phi) + c*(Dr_b + 4*b*Dr_phi)*rr*sin_th - 2*c*Dz_c*rr2*pow(sin_th,2)) +
              c*Dz_phi*rr*(2*a - b + 2*c*rr*cos_th)*pow(sin(2*th),2)))/2.))/(4.*rr);

    // Regularization.
    // If th = 0, then dH/sin_th = sdH as th -> 0.
    if (th == 0.0)
        *sdH *= 1./(1. - aux);
    // Else we can treat it normally.
    else
        *sdH += aux * (dH/sin_th);

    // All done.
    return;
}
