
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

#define SMALL (1.e-9)


void UAv_Kerr_test(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
  const CCTK_REAL dxsq = CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0);
  const CCTK_REAL dysq = CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzsq = CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2);
  */
 
  /* let's create arrays for the F1, F2, F0, phi0, W functions */
  const CCTK_INT N_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *F1, *F2, *F0, *phi0, *W;
  CCTK_REAL *dW_dr, *dW_dth, *d2W_dth2, *d2W_drth;

  F1   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  F2   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  F0   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  phi0 = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  W    = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));

  dW_dr    = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  dW_dth   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  d2W_dth2 = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  d2W_drth = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));


  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        // const CCTK_REAL rho  = sqrt(rho2);

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL RR  = sqrt(RR2);

        const CCTK_REAL r  = RR * (1 + 0.25 * rH/RR) * (1 + 0.25 * rH/RR);

        const CCTK_REAL costh  = z1/RR;
        const CCTK_REAL costh2 = costh*costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth  = sqrt(sinth2);

        /* note that there are divisions by RR in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        F1[ind] = (1.0/2.0)*log(pow(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)), 2) + 256*ct*(ct - rH)*pow(z1, 2)/(pow(RR, 4)*pow(4 + rH/RR, 4)));

        F2[ind] = -F1[ind] + (1.0/2.0)*log(pow(pow(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)), 2) + 256*ct*(ct - rH)/(pow(RR, 2)*pow(4 + rH/RR, 4)), 2) - 256*ct*pow(4*RR - rH, 2)*(ct - rH)*rho2/(pow(RR, 4)*pow(4 + rH/RR, 4)*pow(4*RR + rH, 2)));

        F0[ind] = -F2[ind];

        W[ind] = 4096*sqrt(ct*(ct - rH))*(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)))*(2*ct - rH)*exp(-2*F1[ind] - 2*F2[ind])/(pow(RR, 3)*pow(4 + rH/RR, 6));


        dW_dr[ind] = sqrt(ct*(ct - rH))*(2*ct - rH)*(4*pow(ct, 4) + 2*pow(ct, 3)*r*pow(sinth, 2) - 16*pow(ct, 3)*r - pow(ct, 3)*rH*pow(sinth, 2) - pow(ct, 2)*pow(r, 2)*pow(sinth, 2) + 20*pow(ct, 2)*pow(r, 2) - 2*pow(ct, 2)*r*rH*pow(sinth, 2) + 4*pow(ct, 2)*r*rH + pow(ct, 2)*pow(rH, 2)*pow(sinth, 2) - pow(ct, 2)*pow(rH, 2) - 12*ct*pow(r, 3) + ct*pow(r, 2)*rH*pow(sinth, 2) - 2*ct*pow(r, 2)*rH + 3*pow(r, 4))/pow(-4*pow(ct, 4) + 8*pow(ct, 3)*r + 4*pow(ct, 3)*rH + pow(ct, 2)*pow(r, 2)*pow(sinth, 2) - 8*pow(ct, 2)*pow(r, 2) - pow(ct, 2)*r*rH*pow(sinth, 2) - 4*pow(ct, 2)*r*rH - pow(ct, 2)*pow(rH, 2) + 4*ct*pow(r, 3) - ct*pow(r, 2)*rH*pow(sinth, 2) + 2*ct*pow(r, 2)*rH + ct*r*pow(rH, 2)*pow(sinth, 2) - pow(r, 4), 2);

        dW_dth[ind] = 2*ct*r*sqrt(ct*(ct - rH))*(ct - r)*(ct - rH)*(2*ct - rH)*(r - rH)*sinth*costh/pow(-4*pow(ct, 4) + 8*pow(ct, 3)*r + 4*pow(ct, 3)*rH + pow(ct, 2)*pow(r, 2)*pow(sinth, 2) - 8*pow(ct, 2)*pow(r, 2) - pow(ct, 2)*r*rH*pow(sinth, 2) - 4*pow(ct, 2)*r*rH - pow(ct, 2)*pow(rH, 2) + 4*ct*pow(r, 3) - ct*pow(r, 2)*rH*pow(sinth, 2) + 2*ct*pow(r, 2)*rH + ct*r*pow(rH, 2)*pow(sinth, 2) - pow(r, 4), 2);

        d2W_dth2[ind] = -2*ct*r*sqrt(ct*(ct - rH))*(ct - r)*(ct - rH)*(2*ct - rH)*(r - rH)*(-4*pow(ct, 4)*pow(sinth, 2) + 4*pow(ct, 4)*pow(costh, 2) + 8*pow(ct, 3)*r*pow(sinth, 2) - 8*pow(ct, 3)*r*pow(costh, 2) + 4*pow(ct, 3)*rH*pow(sinth, 2) - 4*pow(ct, 3)*rH*pow(costh, 2) + pow(ct, 2)*pow(r, 2)*pow(sinth, 4) + 3*pow(ct, 2)*pow(r, 2)*pow(sinth, 2)*pow(costh, 2) - 8*pow(ct, 2)*pow(r, 2)*pow(sinth, 2) + 8*pow(ct, 2)*pow(r, 2)*pow(costh, 2) - pow(ct, 2)*r*rH*pow(sinth, 4) - 3*pow(ct, 2)*r*rH*pow(sinth, 2)*pow(costh, 2) - 4*pow(ct, 2)*r*rH*pow(sinth, 2) + 4*pow(ct, 2)*r*rH*pow(costh, 2) - pow(ct, 2)*pow(rH, 2)*pow(sinth, 2) + pow(ct, 2)*pow(rH, 2)*pow(costh, 2) + 4*ct*pow(r, 3)*pow(sinth, 2) - 4*ct*pow(r, 3)*pow(costh, 2) - ct*pow(r, 2)*rH*pow(sinth, 4) - 3*ct*pow(r, 2)*rH*pow(sinth, 2)*pow(costh, 2) + 2*ct*pow(r, 2)*rH*pow(sinth, 2) - 2*ct*pow(r, 2)*rH*pow(costh, 2) + ct*r*pow(rH, 2)*pow(sinth, 4) + 3*ct*r*pow(rH, 2)*pow(sinth, 2)*pow(costh, 2) - pow(r, 4)*pow(sinth, 2) + pow(r, 4)*pow(costh, 2))/pow(-4*pow(ct, 4) + 8*pow(ct, 3)*r + 4*pow(ct, 3)*rH + pow(ct, 2)*pow(r, 2)*pow(sinth, 2) - 8*pow(ct, 2)*pow(r, 2) - pow(ct, 2)*r*rH*pow(sinth, 2) - 4*pow(ct, 2)*r*rH - pow(ct, 2)*pow(rH, 2) + 4*ct*pow(r, 3) - ct*pow(r, 2)*rH*pow(sinth, 2) + 2*ct*pow(r, 2)*rH + ct*r*pow(rH, 2)*pow(sinth, 2) - pow(r, 4), 3);

        d2W_drth[ind] = -2*ct*sqrt(ct*(ct - rH))*(ct - rH)*(2*ct - rH)*(8*pow(ct, 5)*r - 4*pow(ct, 5)*rH - 12*pow(ct, 4)*pow(r, 2) - 8*pow(ct, 4)*r*rH + 4*pow(ct, 4)*pow(rH, 2) + 2*pow(ct, 3)*pow(r, 3)*pow(sinth, 2) - 8*pow(ct, 3)*pow(r, 3) - 3*pow(ct, 3)*pow(r, 2)*rH*pow(sinth, 2) + 36*pow(ct, 3)*pow(r, 2)*rH + pow(ct, 3)*r*pow(rH, 2)*pow(sinth, 2) - 2*pow(ct, 3)*r*pow(rH, 2) - pow(ct, 3)*pow(rH, 3) - pow(ct, 2)*pow(r, 4)*pow(sinth, 2) + 24*pow(ct, 2)*pow(r, 4) - pow(ct, 2)*pow(r, 3)*rH*pow(sinth, 2) - 36*pow(ct, 2)*pow(r, 3)*rH + 3*pow(ct, 2)*pow(r, 2)*pow(rH, 2)*pow(sinth, 2) - 9*pow(ct, 2)*pow(r, 2)*pow(rH, 2) - pow(ct, 2)*r*pow(rH, 3)*pow(sinth, 2) + 2*pow(ct, 2)*r*pow(rH, 3) - 18*ct*pow(r, 5) + ct*pow(r, 4)*rH*pow(sinth, 2) + 21*ct*pow(r, 4)*rH - ct*pow(r, 3)*pow(rH, 2)*pow(sinth, 2) + 4*ct*pow(r, 3)*pow(rH, 2) + 5*pow(r, 6) - 6*pow(r, 5)*rH)*sinth*costh/pow(-4*pow(ct, 4) + 8*pow(ct, 3)*r + 4*pow(ct, 3)*rH + pow(ct, 2)*pow(r, 2)*pow(sinth, 2) - 8*pow(ct, 2)*pow(r, 2) - pow(ct, 2)*r*rH*pow(sinth, 2) - 4*pow(ct, 2)*r*rH - pow(ct, 2)*pow(rH, 2) + 4*ct*pow(r, 3) - ct*pow(r, 2)*rH*pow(sinth, 2) + 2*ct*pow(r, 2)*rH + ct*r*pow(rH, 2)*pow(sinth, 2) - pow(r, 4), 3);


        phi0[ind] = 0.;

      }
    }
  }

  /* printf("F1 = %g\n", F1[0]); */
  /* printf("F2 = %g\n", F2[0]); */
  /* printf("F0 = %g\n", F0[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W = %g\n", W[0]); */

  const CCTK_REAL tt = cctk_time;

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL RR  = sqrt(RR2);

        /*
        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        const CCTK_REAL rho  = sqrt(rho2);
        */

        const CCTK_REAL costh  = z1/RR;
        const CCTK_REAL costh2 = costh*costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth  = sqrt(sinth2);

        /*
        const CCTK_REAL R_x = x1/RR;   // dR/dx
        const CCTK_REAL R_y = y1/RR;   // dR/dy
        const CCTK_REAL R_z = z1/RR;   // dR/dz
        */

        const CCTK_REAL sinth2ph_x = -y1/RR2; // sin(th)^2 dphi/dx
        const CCTK_REAL sinth2ph_y =  x1/RR2; // sin(th)^2 dphi/dy

        const CCTK_REAL R2sinth2ph_x = -y1;  // R^2 sin(th)^2 dphi/dx
        const CCTK_REAL R2sinth2ph_y =  x1;  // R^2 sin(th)^2 dphi/dy

        const CCTK_REAL Rsinthth_x  = z1*x1/RR2; // R sin(th) dth/dx
        const CCTK_REAL Rsinthth_y  = z1*y1/RR2; // R sin(th) dth/dy
        const CCTK_REAL Rsinthth_z  = -sinth2;   // R sin(th) dth/dz

        const CCTK_REAL ph = atan2(y1, x1);

        const CCTK_REAL cosph  = cos(ph);
        const CCTK_REAL sinph  = sin(ph);

        const CCTK_REAL cosmph = cos(mm*ph);
        const CCTK_REAL sinmph = sin(mm*ph);

        /* note the division by RR in the following. divisions by zero should be
           avoided by choosing a non-zero value for z0 (for instance) */
        const CCTK_REAL aux  = 1. + 0.25 * rH/RR;
        const CCTK_REAL aux2 = aux  * aux;
        const CCTK_REAL aux4 = aux2 * aux2;
        const CCTK_REAL aux5 = aux4 * aux;
        const CCTK_REAL aux6 = aux4 * aux2;
        const CCTK_REAL psi4 = exp(2. * F1[ind]) * aux4;
        const CCTK_REAL psi2 = sqrt(psi4);
        const CCTK_REAL psi1 = sqrt(psi2);


        const CCTK_REAL h_rho2 = exp(2. * (F2[ind] - F1[ind])) - 1.;

        // 3-metric
        gxx[ind] = psi4 * (1. + h_rho2 * sinph * sinph);
        gxy[ind] = -psi4 * h_rho2 * sinph * cosph;
        gxz[ind] = 0;
        gyy[ind] = psi4 * (1. + h_rho2 * cosph * cosph);
        gyz[ind] = 0;
        gzz[ind] = psi4;

        /*
          KRph/(R sin(th)^2)  = - 1/2 exp(2F2-F0) (1 + rH/(4R))^6 R dW/dr
          Kthph/(R sin(th))^3 = - 1/2 exp(2F2-F0) (1 + rH/(4R))^5 dW/dth / sin(th) 1/(R - rH/4)
        */

        // KRph/(R sin(th)^2)
        const CCTK_REAL KRph_o_Rsinth2 = -0.5 * exp(2. * F2[ind] - F0[ind]) * aux6 * RR * dW_dr[ind];

        const CCTK_REAL den = RR - 0.25 * rH;
        const CCTK_REAL eps = den/RR; // epsilon = 1 - rh/(4R) : small parameter close to the horizon
        const CCTK_REAL eps_o_1meps = eps/(1 - eps);


        // dW/dth / sin(th) 1/(R - rH/4)
        CCTK_REAL dWdth_o_sinth_den;

        // Axis and horizon regularizations
        //    if at the rho = 0 axis, no need to regularize the division by sin(th), Kij=0 
        //    (also valid in the vicinity of the horizon straight on the axis).
        //    We just need to avoid dividing by sinth=0 at machine precision, otherwise for now it seems the normal operation is enough.
        //    The threshold is chosen such that rho/z1 (~ rho/RR = sinth) < 1e-8, because it enters as a square in the computation of RR.
        if (fabs(sinth) < 1e-8)
          dWdth_o_sinth_den = 0.;
        //    if at R ~ rH/4 we need to regularize the division by R - rH/4. With f(R) = dWdth
        //    f(R) ~ f(rH/4)  +  df_dR(R=rH/4) * (R-rH/4)  +  1/2 * d2f_dR2(R=rH/4) * (R-rH/4)^2  +  1/6 * d3f_dR3(R=rH/4) * (R-rH/4)^3
        //         ~     0    +       0                    +  ... !=0
        //    f(R) / (R-rH/4) ~ [4R/rH * (1-4R/rH) - (4R/rH)^2 * (1-4R/rH)^2] * df_dr(R=rH/4)
        //                    ~ [...] * df_dr(R)      FURTHER ASSUMING df_dr - df_dr(rH/4) ~ (r(R)-rH)*d2f_dr2(rH/4) should be small 
        //                                            (at worst like the order 3 above?)
        //    [...] = eps/(1-eps) - (eps/(1-eps))^2
        else if (fabs(eps) < 1e-4)
          dWdth_o_sinth_den = (eps_o_1meps - eps_o_1meps*eps_o_1meps) * d2W_drth[ind] / sinth;
        else
          dWdth_o_sinth_den = dW_dth[ind] / (den * sinth);

        // Kthph/(R sin(th))^3
        const CCTK_REAL Kthph_o_R3sinth3 = -0.5 * exp(2. * F2[ind] - F0[ind]) * aux5 * dWdth_o_sinth_den;


        // extrinsic curvature
        kxx[ind] = 2.*KRph_o_Rsinth2 *  x1 * sinth2ph_x                     +  2.*Kthph_o_R3sinth3 *  Rsinthth_x * R2sinth2ph_x;
        kxy[ind] =    KRph_o_Rsinth2 * (x1 * sinth2ph_y + y1 * sinth2ph_x)  +     Kthph_o_R3sinth3 * (Rsinthth_x * R2sinth2ph_y + Rsinthth_y * R2sinth2ph_x);
        kxz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_x   +     Kthph_o_R3sinth3 *                              Rsinthth_z * R2sinth2ph_x;
        kyy[ind] = 2.*KRph_o_Rsinth2 *  y1 * sinth2ph_y                     +  2.*Kthph_o_R3sinth3 *  Rsinthth_y * R2sinth2ph_y;
        kyz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_y   +     Kthph_o_R3sinth3 *                              Rsinthth_z * R2sinth2ph_y;
        kzz[ind] = 0.;


        CCTK_REAL alph = exp(F0[ind]) * (RR - 0.25*rH) / (RR + 0.25*rH);
        if (alph < SMALL)
          alph = SMALL;

        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "Kerr_test"))
          alp[ind] = alph;

        // shift
        if (CCTK_EQUALS(initial_shift, "Kerr_test")) {
          betax[ind] =  W[ind] * y1;
          betay[ind] = -W[ind] * x1;
          betaz[ind] =  0.;
        }

        const CCTK_REAL phi0_l = phi0[ind];
        const CCTK_REAL omega = mm * OmegaH;

        // scalar fields
        // Not taking care about the regularization here because there's actually no field
        // and alpha can't be 0 because of test above
        phi1[ind]  = phi0_l * (cos(omega * tt) * cosmph + sin(omega * tt) * sinmph);
        phi2[ind]  = phi0_l * (cos(omega * tt) * sinmph - sin(omega * tt) * cosmph);
        Kphi1[ind] = 0.5 * mm * (W[ind] - OmegaH) / alph * phi2[ind];
        Kphi2[ind] = 0.5 * mm * (OmegaH - W[ind]) / alph * phi1[ind];

      } /* for i */
    }   /* for j */
  }     /* for k */

  free(F1); free(F2); free(F0); free(phi0); free(W);
  free(dW_dr); free(dW_dth); free(d2W_dth2); free(d2W_drth);

  return;
}
