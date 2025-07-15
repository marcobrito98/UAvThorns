
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

void UAv_IDBHProcaHair_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL []);


void UAv_IDProcaBS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  // WARNING: rotating stars not being proof-tested yet (output of warning in ParamCheck)
  // TODO: Remove when proof-tested.
  // What needs to be checked if mostly the behavior/regularity of W and K_ij
  
  // TODO: Is there a way to write the Proca fields in a regular way 
  //       (without 1/rr in particular, but also dealing with the axis...)
  //       like for the scalar case?


  /*
  const CCTK_REAL dxsq = CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0);
  const CCTK_REAL dysq = CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzsq = CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2);
  */

  CCTK_INT NF;      // NF will be the actual size of the arrays
  CCTK_INT NX;      // NX will be the number of X points
  CCTK_INT Ntheta;  // Ntheta will be the number of theta points

  CCTK_REAL *Xtmp, *thtmp, *F1_in, *F2_in, *F0_in, *Wbar_in, *H1_in, *H2_in, *H3_in, *V_in;
  Xtmp     = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  thtmp    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F1_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F2_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F0_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  Wbar_in  = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  H1_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  H2_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  H3_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  V_in     = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));

  // we get the data from the input file
  UAv_IDBHProcaHair_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, Wbar_in, H1_in, H2_in, H3_in, V_in);

  Ntheta = NF/NX;

  CCTK_VInfo(CCTK_THORNSTRING, "NX     = %d", NX);
  CCTK_VInfo(CCTK_THORNSTRING, "Ntheta = %d", Ntheta);
  CCTK_VInfo(CCTK_THORNSTRING, "NF     = %d", NF);

  // Minimum number of theta points required by the FD stencils, and the j index assignment below
  const CCTK_INT min_theta_pts = 5;
  if (Ntheta<min_theta_pts){
        CCTK_VERROR ("The initial data file doesn't have enough points in theta to be consistent with the implementation.\n "
        "Ntheta = %d. min_theta_points = %d.",
        Ntheta, min_theta_pts);
  }

  // now we create arrays with the X and theta coordinates
  CCTK_REAL X[NX], theta[Ntheta];
  for (int i = 0; i < NX; i++) {
    X[i]     = Xtmp[i];
    /* printf("X[%3d] = %lf\n", i, X[i]); */
  }
  for (int i = 0; i < Ntheta; i++) {
    theta[i] = thtmp[i*NX];
    /* printf("theta[%3d] = %lf\n", i, theta[i]); */
  }

  // the spacing in each coordinate is
  const CCTK_REAL dX     = (X[NX-1] - X[0])/(NX-1);
  const CCTK_REAL dtheta = (theta[Ntheta-1] - theta[0])/(Ntheta-1);

  /* printf("dX     = %e\n", dX); */
  /* printf("dtheta = %e\n", dtheta); */

  // make sure spacing is uniform in the provided grid
  for (int i = 1; i < NX; i++) {
    if (fabs(X[i] - X[i - 1] - dX) > SMALL) {
      printf("i = %d\n", i);
      CCTK_WARN(0, "X grid is not uniformly spaced. Aborting.");
    }
  }
  for (int j = 1; j < Ntheta; j++) {
    if (fabs(theta[j] - theta[j - 1] - dtheta) > SMALL)
      CCTK_WARN(0, "theta grid is not uniformly spaced. Aborting.");
  }


  // To take care properly of z=0 symmetry (i.e. theta <-> pi-theta) in interpolation
  // we need to extend the arrays to z<0 values.
  // For convenience, we keep Ntheta as the number of points in the input half-space.

  NF = NX * (2*Ntheta - 1);

  CCTK_REAL *F1_extd, *F2_extd, *F0_extd, *H2_extd, *H3_extd, *V_extd;
  F1_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  F2_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  F0_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  H2_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  H3_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  V_extd     = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // We'll use A_r (H1r) rather than the input H1_in
  // Notation here: A ~ H1r dr + ... = H1_in/r dr + ...
  CCTK_REAL *H1r_extd;
  H1r_extd          = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // now we need to take the derivatives of the Wbar function
  // Then we convert to W and store the values

  CCTK_REAL *W_extd, *dW_dr_extd, *dW_dth_extd;
  W_extd        = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // Same for H3 and V
  
  CCTK_REAL *dH3_dr_extd, *dH3_dth_extd;
  dH3_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH3_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  
  CCTK_REAL *dV_dr_extd, *dV_dth_extd;
  dV_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dV_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));


  // // Some auxi file for debug
  // FILE* debugfile = fopen ("testdebug.txt", "w");
  // if (debugfile == NULL) {
  //   CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
  //   "Unable to open file %s\n", "testdebug.txt");
  // } else {
  //   CCTK_VInfo(CCTK_THORNSTRING, "Write test file %s", "testdebug.txt");
  // }

  const CCTK_REAL oodX       = 1. / dX;
  const CCTK_REAL oodX12     = 1. / (12. * dX);
  const CCTK_REAL oodXsq12   = oodX * oodX12;
  const CCTK_REAL oodth12    = 1. / (12. * dtheta);

  // First loop on z>=0 half-space (i.e. input values of 0 <= theta <= pi/2)

  for (int jj = 0; jj < Ntheta; jj++) {
    for (int i = 0; i < NX; i++) {

      /* Theta symmetries
        WARNING: Be careful with symmetries, they may not be valid for every quantity!
        
        Instead of spelling out the cases, use symmetry directly
        theta = 0    (relevant for jm.): theta <-> -theta    i.e we need   max (theta, -theta)     = |theta|
        theta = pi/2 (relevant for jp.): theta <-> pi-theta  i.e we need   min (theta, pi-theta)   = pi/2 - |pi/2 - theta|

        /!\ Make sure that the rationale here is consistent with min_theta_points defined above.
      */
      const CCTK_INT j   = jj;
      const CCTK_INT jm1 = abs(jj-1);
      const CCTK_INT jm2 = abs(jj-2);
      const CCTK_INT jm3 = abs(jj-3);
      const CCTK_INT jm4 = abs(jj-4);
      const CCTK_INT jp1 = Ntheta - 1 - abs(Ntheta - 1 - (jj+1));
      const CCTK_INT jp2 = Ntheta - 1 - abs(Ntheta - 1 - (jj+2));
      const CCTK_INT jp3 = Ntheta - 1 - abs(Ntheta - 1 - (jj+3));
      const CCTK_INT jp4 = Ntheta - 1 - abs(Ntheta - 1 - (jj+4));



      const CCTK_INT ind    = i + j*NX;

      const CCTK_INT indim1 = i-1 + j*NX;
      const CCTK_INT indip1 = i+1 + j*NX;
      const CCTK_INT indim2 = i-2 + j*NX;
      const CCTK_INT indip2 = i+2 + j*NX;
      const CCTK_INT indip3 = i+3 + j*NX;
      const CCTK_INT indip4 = i+4 + j*NX;
      const CCTK_INT indip5 = i+5 + j*NX;

      const CCTK_INT indjm1 = i + jm1*NX;
      const CCTK_INT indjm2 = i + jm2*NX;
      const CCTK_INT indjm3 = i + jm3*NX;
      const CCTK_INT indjm4 = i + jm4*NX;
      const CCTK_INT indjp1 = i + jp1*NX;
      const CCTK_INT indjp2 = i + jp2*NX;
      const CCTK_INT indjp3 = i + jp3*NX;
      const CCTK_INT indjp4 = i + jp4*NX;


      // Just copy input values of ansatz functions
      F1_extd[ind]   = F1_in[ind];
      F2_extd[ind]   = F2_in[ind];
      F0_extd[ind]   = F0_in[ind];
      H2_extd[ind]   = H2_in[ind];
      H3_extd[ind]   = H3_in[ind];
      V_extd[ind]    = V_in[ind];


      const CCTK_REAL lX = X[i];
      /* const CCTK_REAL lth = theta[j]; */
      /* printf("X[%3d] = %lf\n", i, lX); */


      // Theta derivatives
      // /!\ Make sure that the stencils here are consistent with min_theta_points defined above.


      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL Wbar_th = (-Wbar_in[indjp2] + 8 * Wbar_in[indjp1] - 8 * Wbar_in[indjm1] + Wbar_in[indjm2]) *
      oodth12;
      
      // WARNING/TODO (rotating stars): Do we need to be careful with theta derivatives, like for V? Depending on m?
      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL H3_th = (-H3_in[indjp2] + 8 * H3_in[indjp1] - 8 * H3_in[indjm1] + H3_in[indjm2]) *
        oodth12;

      // Symmetries of V on the axis and/or the equator can vary (theta = 0, pi/2 resp.).
      // In particular, it can occur that dV/dth != 0, which can't be captured by centered finite differences and theta symmetry.
      // Since different systems have different symmetries, we resort to non-symmetric stencils in any case
      CCTK_REAL V_th;
      if (jj==0) {
        // 1st derivative with 4th order accuracy (forward stencils)
        V_th = (- 25 * V_in[ind]    + 48 * V_in[indjp1] - 36 * V_in[indjp2] + 16 * V_in[indjp3] - 3 * V_in[indjp4]) *
          oodth12;
      } else if (jj==1) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (-  3 * V_in[indjm1] - 10 * V_in[ind]    + 18 * V_in[indjp1] -  6 * V_in[indjp2] +     V_in[indjp3]) * 
          oodth12;
      } else if (jj==Ntheta-2) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (   3 * V_in[indjp1] + 10 * V_in[ind]    - 18 * V_in[indjm1] +  6 * V_in[indjm2] -     V_in[indjm3]) * 
          oodth12;
      } else if (jj==Ntheta-1) {
        // 1st derivative with 4th order accuracy (backward stencils)
        V_th = (  25 * V_in[ind]    - 48 * V_in[indjm1] + 36 * V_in[indjm2] - 16 * V_in[indjm3] + 3 * V_in[indjm4]) *
          oodth12;
      } else {
        // 1st derivative with 4th order accuracy (centered stencils)
        V_th = (-V_in[indjp2] + 8 * V_in[indjp1] - 8 * V_in[indjm1] + V_in[indjm2]) *
          oodth12;
      }

      CCTK_REAL Wbar_X, H3_X, V_X;
      CCTK_REAL Wbar_XX = 0.; // Used for r=0 (i==0), if Wbar_r_power == 2.
      CCTK_REAL H1_X = 0.;    // Used for r=0 (i==0), due to H1_in/r.


      /*
      Regarding finite differencing orders: for Scalar BS, plotting W and dW_dr, there were small discontinuities near r=0
      during tests with the previous 2nd order accuracy for i==0 and i==1.
      Those vanish when moving to 4th order accuracy.

      For i==NX-1 and i==NX-2, we keep 2nd order for now. The issue is not appearing as clearly,
      and they represent points which are physically far, so maybe better to keep the computation more local.
      */

      if (i == 0) {
        /* For the Boson Star, there's no issue, dWbar/dX != 0 at X==0, and x and r coordinates coincide. */

        // 1st derivative with 4th order accuracy (forward stencils)
        Wbar_X =(- 25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] + 16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) * oodX12;

        // 1st derivative with 4th order accuracy (forward stencils)
        H3_X =(- 25 * H3_in[ind] + 48 * H3_in[indip1] - 36 * H3_in[indip2] + 16 * H3_in[indip3] - 3 * H3_in[indip4]) * oodX12;
        
        // 1st derivative with 4th order accuracy (forward stencils)
        V_X =(- 25 * V_in[ind] + 48 * V_in[indip1] - 36 * V_in[indip2] + 16 * V_in[indip3] - 3 * V_in[indip4]) * oodX12;


        // Special care at r=0

        // H1_X required for H1r computed from H1_in/r
        // 1st derivative with 4th order accuracy (forward stencils)
        H1_X =(- 25 * H1_in[ind] + 48 * H1_in[indip1] - 36 * H1_in[indip2] + 16 * H1_in[indip3] - 3 * H1_in[indip4]) * oodX12;


        if (Wbar_r_power == 2) {
          // If Wbar = r^2 * W, to compute W(r=0), we need to compute Wbar_XX.
          // 2nd derivative with 4th order accuracy (forward stencils)
          Wbar_XX = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] + 214 * Wbar_in[indip2] 
                    - 156 * Wbar_in[indip3] + 61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) * oodXsq12;
        }

      } else if (i == 1 ) {
        // 1st derivative, 4th order accuracy
        Wbar_X = (- 3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] - 6 * Wbar_in[indip2] + Wbar_in[indip3]) * oodX12;
        
        // 1st derivative, 4th order accuracy
        H3_X = (- 3 * H3_in[indim1] - 10 * H3_in[ind] + 18 * H3_in[indip1] - 6 * H3_in[indip2] + H3_in[indip3]) * oodX12;
        
        // 1st derivative, 4th order accuracy
        V_X = (- 3 * V_in[indim1] - 10 * V_in[ind] + 18 * V_in[indip1] - 6 * V_in[indip2] + V_in[indip3]) * oodX12;

      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4*Wbar_in[indim1] + 3*Wbar_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        H3_X = (H3_in[indim2] - 4*H3_in[indim1] + 3*H3_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        V_X = (V_in[indim2] - 4*V_in[indim1] + 3*V_in[ind]) * 0.5 * oodX;

      } else if (i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;
        
        // 1st derivative with 2nd order accuracy (central stencils)
        H3_X = (-H3_in[indim1] + H3_in[indip1]) * 0.5 * oodX;
        
        // 1st derivative with 2nd order accuracy (central stencils)
        V_X = (-V_in[indim1] + V_in[indip1]) * 0.5 * oodX;

      } else {
        // 4th order accurate stencils
        Wbar_X    = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
        
        // 4th order accurate stencils
        H3_X    = (-H3_in[indip2] + 8 * H3_in[indip1] - 8 * H3_in[indim1] + H3_in[indim2]) * oodX12;
        
        // 4th order accurate stencils
        V_X    = (-V_in[indip2] + 8 * V_in[indip1] - 8 * V_in[indim1] + V_in[indim2]) * oodX12;
      
      }

      // From the X coordinate used in the input files to the r coordinate (coincides with x for the Boson Star, rH=0).
      // We also do the conversion from Wbar to W here, and H1_in to H1r, to tackle r = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  r == 0
      if (i == 0) {

        // At r=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_extd[ind]    = 0.; 
        dW_dth_extd[ind]   = 0.;
        
        // At X==0 (rr==0), dXdr = 1/C0
        dH3_dr_extd[ind]       = H3_X / C0;
        dV_dr_extd[ind]        =  V_X / C0;

        // For W we need more care depending on the power
        switch (Wbar_r_power)
        {
        case 0:   // Wbar = W
          W_extd[ind]        = Wbar_in[ind];
          break;
        
        case 1:   // Wbar = r * W
          /*
          dWbar/dr = W + r * dW/dr
                   = W + 0 * 0     at r=0
          
          dWbar/dr = dWbar/dX * dX/dr
          dX/dr = C/(C+r)^2 = 1/C  at r=0
          */
          W_extd[ind]        = Wbar_X / C0;
          break;
        
        case 2:   // Wbar = r^2 * W
          /*
          dWbar/dr   = 2r * W + r^2 * dW/dr
          d2Wbar/dr2 = 2  * W + 4r  * dW/dr + r^2 * d2W/dr2
                     = 2  * W + 0 + 0        at r=0

          d2Wbar/dr2 = d2Wbar/dX2 * (dX/dr)^2 + dWbar/dX * d2X/dr2
          dX/dr   =   C/(C+r)^2 =  1/C      at r=0
          d2X/dr2 = -2C/(C+r)^3 = -2/C^2    at r=0

          W (r=0) = 1/C^2 * [1/2 * d2Wbar/dX^2 (X=0)  -  dWbar/dX (X=0)]
          */
          W_extd[ind]        = (0.5 * Wbar_XX - Wbar_X)/(C0*C0);
          break;
        
        default:  // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }

        // From H1_in/r to H1r (idem Wbar above, case 1)
        H1r_extd[ind] = H1_X / C0;
      }

      // We need to be careful at X == 1 (r == infty) for radial derivatives (coordinate change is singular)
      else if (i == NX - 1) {
        // W -> 0 for r -> infty
        W_extd[ind]        = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction paper) also gives:
        dW_dr_extd[ind]    = 0.; 
        dW_dth_extd[ind]   = 0.; 


        dH3_dr_extd[ind]    = 0.;
        dV_dr_extd[ind]    = 0.;

        H1r_extd[ind] = 0.; // A_r = 0 at infinity

      } else {

        const CCTK_REAL rr = C0*lX/(1. - lX);

        // corresponding derivatives
        // const CCTK_REAL dXdr = 1./(C0 + rr) - rr/((C0 + rr)*(C0 + rr));
        const CCTK_REAL dXdr = C0/((C0 + rr)*(C0 + rr));

        const CCTK_REAL Wbar_r = dXdr * Wbar_X;

        dH3_dr_extd[ind]       = dXdr * H3_X;
        dV_dr_extd[ind]        = dXdr * V_X;
        
        // Now translate from Wbar to W
        switch (Wbar_r_power) // We could put a generic power for the computation here I guess...
        {
        case 0:   // Wbar = W
          W_extd[ind]        = Wbar_in[ind];
          dW_dr_extd[ind]    = Wbar_r;
          dW_dth_extd[ind]   = Wbar_th;
          break;
        
        case 1:   // Wbar = r * W
          W_extd[ind]        = Wbar_in[ind] / rr;
          dW_dr_extd[ind]    = (Wbar_r - W_extd[ind]) / rr; // dW/dr  =  1/r * dWbar/dr - Wbar / r^2  =  (dWbar/dr - W) / r
          dW_dth_extd[ind]   = Wbar_th / rr;
          break;
        
        case 2: ; // Wbar = r^2 * W
          // empty statement after case to prevent compilation error on some gcc versions...
          const CCTK_REAL rr2 = rr*rr;
          W_extd[ind]        = Wbar_in[ind] / rr2;
          dW_dr_extd[ind]    = Wbar_r / rr2 - 2 * W_extd[ind] / rr; // dW/dr  =  1/r^2 * dWbar/dr - 2 * Wbar / r^3  =  1/r^2 * dWbar/dr - 2 * W / r
          dW_dth_extd[ind]   = Wbar_th / rr2;
          break;
        
        default:  // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }

        // From H1_in/r to H1r
        H1r_extd[ind] = H1_in[ind] / rr;
        
      } // if/else i==...
      
      dH3_dth_extd[ind]     = H3_th;
      dV_dth_extd[ind]      = V_th;
    
      // fprintf (debugfile, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f ", 
      //             W_extd[ind], dW_dr_extd[ind], dW_dth_extd[ind],
      //             H3_extd[ind], dH3_dr_extd[ind], dH3_dth_extd[ind],
      //             V_extd[ind], dV_dr_extd[ind], dV_dth_extd[ind],
      //             H1r_extd[ind]);
    } // for i
    // fprintf (debugfile, "\n");
  } // for jj

  // fclose(debugfile);


  // Second loop on z<0 half-space (completion by symmetry)

  // Even parity: F1, F2, F0, W; r derivatives of even functions; theta derivatives of odd functions
  // Odd parity:  r derivatives of odd functions; theta derivatives of even functions
  // A even (default): H1, H3 and V are even, and H2 is odd.
  // A odd           : H1, H3 and V are  odd, and H2 is even.

  const CCTK_INT H1_z_sign = Amu_z_sym_is_odd ? -1 : +1;
  const CCTK_INT H2_z_sign = Amu_z_sym_is_odd ? +1 : -1;
  const CCTK_INT H3_z_sign = Amu_z_sym_is_odd ? -1 : +1;
  const CCTK_INT  V_z_sign = Amu_z_sym_is_odd ? -1 : +1;

  for (int jj = 1; jj < Ntheta; jj++) { // don't repeat theta == pi/2
    for (int i = 0; i < NX; i++) {

      // j or jsym == Ntheta - 1  is theta == pi/2
      const CCTK_INT j    = Ntheta - 1 + jj;
      const CCTK_INT jsym = Ntheta - 1 - jj; 

      const CCTK_INT ind    = i + j   *NX;
      const CCTK_INT indsym = i + jsym*NX;

      // Even
      F1_extd[ind]       = F1_extd[indsym];
      F2_extd[ind]       = F2_extd[indsym];
      F0_extd[ind]       = F0_extd[indsym];
      
      W_extd[ind]        = W_extd[indsym];
      dW_dr_extd[ind]    = dW_dr_extd[indsym];
      
      // Odd
      dW_dth_extd[ind]      = - dW_dth_extd[indsym];
      
      // Vector potential
      
      H1r_extd[ind]      = H1_z_sign * H1r_extd[indsym];
      H2_extd[ind]       = H2_z_sign *  H2_extd[indsym];
      H3_extd[ind]       = H3_z_sign *  H3_extd[indsym];
      V_extd[ind]        =  V_z_sign *   V_extd[indsym];
      
      dH3_dr_extd[ind]   = H3_z_sign * dH3_dr_extd[indsym];
      dV_dr_extd[ind]    =  V_z_sign *  dV_dr_extd[indsym];

      dH3_dth_extd[ind]  = - H3_z_sign * dH3_dth_extd[indsym];
      dV_dth_extd[ind]   = -  V_z_sign *  dV_dth_extd[indsym];
      

      } // for i
  } // for jj


  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *X_g, *theta_g;
  X_g     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  theta_g = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL rr2 = x1*x1 + y1*y1 + z1*z1;

        CCTK_REAL rr  = sqrt(rr2);
        /* For the Boson Star, x, r and R coordinates coincide (rH=0). */
	/* note that there are divisions by rr in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */
        
        // From r to the X radial coordinate (used in input files)
        const CCTK_REAL lX = rr / (C0 + rr);

        const CCTK_REAL ltheta = rr < 1e-16 ? 0 : acos( z1/rr );    // There should be at most one point in the grid with rr~0. Not sure about the threshold.

        X_g[ind]     = lX;
        theta_g[ind] = ltheta;
      }
    }
  }

  /* now for the interpolation */

  const CCTK_INT N_dims  = 2;   // 2-D interpolation

  const CCTK_INT N_input_arrays  = 14;
  const CCTK_INT N_output_arrays = 14;

  /* origin and stride of the input coordinates. with this Cactus reconstructs
     the whole X and theta array. */
  CCTK_REAL origin[N_dims];
  CCTK_REAL delta [N_dims];
  origin[0] = X[0];  origin[1] = theta[0];
  delta[0]  = dX;    delta[1]  = dtheta;

  /* points onto which we want to interpolate, ie, the grid points themselves in
     (X, theta) coordinates (computed above) */
  const void *interp_coords[N_dims];
  interp_coords[0] = (const void *) X_g;
  interp_coords[1] = (const void *) theta_g;


  /* input arrays */
  const void *input_arrays[N_input_arrays];
  CCTK_INT input_array_type_codes[N_input_arrays];
  CCTK_INT input_array_dims[N_dims];
  input_array_dims[0] = NX;
  input_array_dims[1] = 2*Ntheta-1;

  input_array_type_codes[0] = CCTK_VARIABLE_REAL;
  input_array_type_codes[1] = CCTK_VARIABLE_REAL;
  input_array_type_codes[2] = CCTK_VARIABLE_REAL;
  input_array_type_codes[3] = CCTK_VARIABLE_REAL;
  input_array_type_codes[4] = CCTK_VARIABLE_REAL;
  input_array_type_codes[5] = CCTK_VARIABLE_REAL;
  input_array_type_codes[6] = CCTK_VARIABLE_REAL;
  input_array_type_codes[7] = CCTK_VARIABLE_REAL;
  input_array_type_codes[8] = CCTK_VARIABLE_REAL;
  input_array_type_codes[9] = CCTK_VARIABLE_REAL;
  input_array_type_codes[10]= CCTK_VARIABLE_REAL;
  input_array_type_codes[11]= CCTK_VARIABLE_REAL;
  input_array_type_codes[12]= CCTK_VARIABLE_REAL;
  input_array_type_codes[13]= CCTK_VARIABLE_REAL;

  /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
  input_arrays[0] = (const void *) F1_extd;
  input_arrays[1] = (const void *) F2_extd;
  input_arrays[2] = (const void *) F0_extd;
  input_arrays[3] = (const void *) W_extd;
  input_arrays[4] = (const void *) dW_dr_extd;
  input_arrays[5] = (const void *) dW_dth_extd;
  input_arrays[6] = (const void *) H1r_extd;
  input_arrays[7] = (const void *) H2_extd;
  input_arrays[8] = (const void *) H3_extd;
  input_arrays[9] = (const void *) V_extd;
  input_arrays[10]= (const void *) dH3_dr_extd;
  input_arrays[11]= (const void *) dH3_dth_extd;
  input_arrays[12]= (const void *) dV_dr_extd;
  input_arrays[13]= (const void *) dV_dth_extd;

  /* output arrays */
  void *output_arrays[N_output_arrays];
  CCTK_INT output_array_type_codes[N_output_arrays];
  CCTK_REAL *F1, *F2, *F0, *W, *H1r, *H2, *H3, *V;
  CCTK_REAL *dW_dr, *dW_dth;
  CCTK_REAL *dH3_dr, *dH3_dth;
  CCTK_REAL *dV_dr, *dV_dth;

  F1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F0          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  W           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dr       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dth      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H1r         = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H3          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  V           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dr      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dth     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dr       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dth      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  output_array_type_codes[0] = CCTK_VARIABLE_REAL;
  output_array_type_codes[1] = CCTK_VARIABLE_REAL;
  output_array_type_codes[2] = CCTK_VARIABLE_REAL;
  output_array_type_codes[3] = CCTK_VARIABLE_REAL;
  output_array_type_codes[4] = CCTK_VARIABLE_REAL;
  output_array_type_codes[5] = CCTK_VARIABLE_REAL;
  output_array_type_codes[6] = CCTK_VARIABLE_REAL;
  output_array_type_codes[7] = CCTK_VARIABLE_REAL;
  output_array_type_codes[8] = CCTK_VARIABLE_REAL;
  output_array_type_codes[9] = CCTK_VARIABLE_REAL;
  output_array_type_codes[10]= CCTK_VARIABLE_REAL;
  output_array_type_codes[11]= CCTK_VARIABLE_REAL;
  output_array_type_codes[12]= CCTK_VARIABLE_REAL;
  output_array_type_codes[13]= CCTK_VARIABLE_REAL;

  output_arrays[0] = (void *) F1;
  output_arrays[1] = (void *) F2;
  output_arrays[2] = (void *) F0;
  output_arrays[3] = (void *) W;
  output_arrays[4] = (void *) dW_dr;
  output_arrays[5] = (void *) dW_dth;
  output_arrays[6] = (void *) H1r;
  output_arrays[7] = (void *) H2;
  output_arrays[8] = (void *) H3;
  output_arrays[9] = (void *) V;
  output_arrays[10]= (void *) dH3_dr;
  output_arrays[11]= (void *) dH3_dth;
  output_arrays[12]= (void *) dV_dr;
  output_arrays[13]= (void *) dV_dth;


  /* handle and settings for the interpolation routine */
  int operator_handle, param_table_handle;
  operator_handle    = CCTK_InterpHandle("Lagrange polynomial interpolation");
  param_table_handle = Util_TableCreateFromString("order=4 boundary_extrapolation_tolerance={0.1 1.0 0.05 0.05}");

  CCTK_INFO("Interpolating result...");

  /* do the actual interpolation, and check for error returns */
  int status = CCTK_InterpLocalUniform(N_dims, operator_handle,
                                       param_table_handle,
                                       origin, delta,
                                       N_interp_points,
                                       CCTK_VARIABLE_REAL,
                                       interp_coords,
                                       N_input_arrays, input_array_dims,
                                       input_array_type_codes,
                                       input_arrays,
                                       N_output_arrays, output_array_type_codes,
                                       output_arrays);
  if (status < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
    "interpolation screwed up!");
  }

  free(X_g); free(theta_g);
  free(Xtmp); free(thtmp);
  free(F1_in); free(F2_in); free(F0_in); free(Wbar_in);
  free(H1_in); free(H2_in); free(H3_in); free(V_in);
  
  free(F1_extd); free(F2_extd); free(F0_extd); free(W_extd);
  free(H1r_extd); free(H2_extd); free(H3_extd); free(V_extd);
  free(dW_dr_extd); free(dW_dth_extd);
  free(dH3_dr_extd); free(dH3_dth_extd);
  free(dV_dr_extd); free(dV_dth_extd);


  /* printf("F1 = %g\n", F1[0]); */
  /* printf("F2 = %g\n", F2[0]); */
  /* printf("F0 = %g\n", F0[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W = %g\n", W[0]); */


  /* now we finally write the metric and all 3+1 quantities. first we write the
     3-metric and extrinsic curvature, then Proca fields, then lapse and shift */
  /* For the Boson Star, in order to avoid unneeded regularizations and divisions,
      we express K_ij in terms of dW/drho and dW/dz.
      For points close to the axis and the origin, K_ij = 0.
  */

  const CCTK_REAL tt = cctk_time;

  const CCTK_REAL coswt = cos(omega_BS * tt);
  const CCTK_REAL sinwt = sin(omega_BS * tt);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        // For the Boson Star, r = R, no coordinate change needed.
        const CCTK_REAL rr2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL rr  = sqrt(rr2);
	/* note that there are divisions by rr in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        const CCTK_REAL rho  = sqrt(rho2);

        const CCTK_REAL costh  = z1/rr;
        const CCTK_REAL costh2 = costh*costh;
        /*
          For some grid points actually on the axis, it occurred that costh = 1-1e-16, resulting in sinth ~ 1.5e-8 instead of 0.
          Thus we force it in that case. 
          Even if there is a legit grid point such that theta ~ a few 1e-8, it should mean RR >> rho and the axis treatment should be fine.
        */
        CCTK_REAL sinth, sinth2;
        if (1-costh2 < 1e-15) {
          sinth2 = 0.;
          sinth  = 0.;
        } else {
          sinth2 = 1. - costh2;
          sinth  = sqrt(sinth2);
        }
        

        const CCTK_REAL ph = atan2(y1, x1);
        // If x1=y1=0, should return 0? The other metric functions should vanish anyway to make sure that this doesn't matter,
        // but can this lead to nan depending on the C implementation?

        const CCTK_REAL cosph  = cos(ph);
        const CCTK_REAL sinph  = sin(ph);

        const CCTK_REAL cosmph = cos(mm*ph);
        const CCTK_REAL sinmph = sin(mm*ph);

        const CCTK_REAL psi4 = exp(2. * F1[ind]);
        const CCTK_REAL psi2 = sqrt(psi4);
        const CCTK_REAL psi1 = sqrt(psi2);

        const CCTK_REAL h_rho2 = exp(2. * (F2[ind] - F1[ind])) - 1.;

        // add non-axisymmetric perturbation on conformal factor
        // NOTE: the perturbation is only taken into account for the 3-metric grid functions (not extrinsic curvature, lapse, ...)
        const CCTK_REAL argpert_cf = (rr - R0pert_conf_fac)/Sigmapert_conf_fac;
        const CCTK_REAL pert_cf = 1. + Apert_conf_fac * (x1*x1 - y1*y1)*mu*mu * exp( -0.5*argpert_cf*argpert_cf );

        const CCTK_REAL conf_fac = psi4 * pert_cf;

        // 3-metric
        gxx[ind] = conf_fac * (1. + h_rho2 * sinph * sinph);
        gxy[ind] = -conf_fac * h_rho2 * sinph * cosph;
        gxz[ind] = 0;
        gyy[ind] = conf_fac * (1. + h_rho2 * cosph * cosph);
        gyz[ind] = 0;
        gzz[ind] = conf_fac;

        /*
          d/drho = rho/r * d/dr  +    z/r^2 * d/dth
          d/dz   =   z/r * d/dr  -  rho/r^2 * d/dth

          Kxx = 0.5 * 2xy/rho        * exp(2F2-F0) * dW/drho   = 0.5 * rho * sin(2phi) * exp(2F2-F0) * dW/drho
          Kyy = - Kxx
          Kzz = 0
          Kxy =-0.5 * (x^2-y^2)/rho  * exp(2F2-F0) * dW/drho   = 0.5 * rho * cos(2phi) * exp(2F2-F0) * dW/drho
          Kxz = 0.5 * y * exp(2F2-F0) * dW/dz
          Kyz =-0.5 * x * exp(2F2-F0) * dW/dz
        */

        /*
          Close to the axis and the origin, Kij = 0.
          The "coordinate" part of the expressions above behave like rho (or r).
          Let's first consider a threshold of rho < 1e-8. The sphere r < 1e-8 is included in this cylinder.
          In this case, we just set d/drho and d/dz = 0 as proxies.
        */

        CCTK_REAL dW_drho, dW_dz;
        const CCTK_REAL exp_auxi = exp(2. * F2[ind] - F0[ind]);

        if (rho < 1e-8) {
          dW_drho = 0.;
          dW_dz   = 0.;
        }
        else {
          dW_drho = rho/rr * dW_dr[ind]  +   z1/rr2 * dW_dth[ind];
          dW_dz   =  z1/rr * dW_dr[ind]  -  rho/rr2 * dW_dth[ind];
        }

        // extrinsic curvature
        kxx[ind] =  0.5 * rho * sin(2*ph) * exp_auxi * dW_drho;
        kxy[ind] = -0.5 * rho * cos(2*ph) * exp_auxi * dW_drho;
        kxz[ind] =  0.5 *  y1 * exp_auxi * dW_dz;
        kyy[ind] = -kxx[ind];
        kyz[ind] = -0.5 *  x1 * exp_auxi * dW_dz;
        kzz[ind] =  0.;

          

        // lapse value (field initialization below)
        // No lapse regularization needed for the BS, the lapse is non-zero
        const CCTK_REAL alph = exp(F0[ind]);


        // let's add a perturbation to the Proca field as well
        // NOTE: the perturbation is added directed to every instance of e^{i m \varphi}, hence its derivatives are not taken into account
        // TODO (?): Design perturbation more generically as ~ cos((m+1)\varphi)
        const CCTK_REAL argpert_Proca = (rr - R0pert_Proca)/Sigmapert_Proca;
        const CCTK_REAL pert_Proca = 1. + Apert_Proca * (x1*x1 - y1*y1)*mu*mu * exp( -0.5*argpert_Proca*argpert_Proca );


        // ----- Proca fields -----

        // TODO: check what happens with divisions by rr and sinth, can we work around them?

        // Real and imaginay part of the harmonic dependence: exp[i(m\varphi - \omega t)]
        const CCTK_REAL harm_re = (coswt * cosmph + sinwt * sinmph) * pert_Proca;
        const CCTK_REAL harm_im = (coswt * sinmph - sinwt * cosmph) * pert_Proca;

        // No need to change the radial component, R and r coincide
        // A_x
        A1x[ind] = x1/rr * H1r[ind] * harm_re + costh*cosph/rr * H2[ind] * harm_re + sinph/rr * H3[ind] * harm_im;
        A2x[ind] = x1/rr * H1r[ind] * harm_im + costh*cosph/rr * H2[ind] * harm_im - sinph/rr * H3[ind] * harm_re;
        
        // A_y
        A1y[ind] = y1/rr * H1r[ind] * harm_re + costh*sinph/rr * H2[ind] * harm_re - cosph/rr * H3[ind] * harm_im;
        A2y[ind] = y1/rr * H1r[ind] * harm_im + costh*sinph/rr * H2[ind] * harm_im + cosph/rr * H3[ind] * harm_re;
        
        // A_z
        A1z[ind] = (z1/rr * H1r[ind] - sinth/rr * H2[ind]) * harm_re;
        A2z[ind] = (z1/rr * H1r[ind] - sinth/rr * H2[ind]) * harm_im;

        // A_\phi
        /*
          A_\phi = -n^\mu A_\mu = - (A_t + W*A_ph)/alpha
                = -i * e^{i (m ph - w t)} * (V + W H3 sinth) / alpha
        */
        Aphi1[ind] = (V[ind] + W[ind] * sinth * H3[ind]) / alph * harm_im;
        Aphi2[ind] =-(V[ind] + W[ind] * sinth * H3[ind]) / alph * harm_re;

        // ----- Electric fields -----
        
        // First, E_i in (R, th, ph) coordinates
        CCTK_REAL E1d_r, E2d_r, E1d_th, E2d_th, E1d_ph_o_sinth, E2d_ph_o_sinth;

        // E_r
        /*
          E_r = i * e^{i(m phi - w t)} / alpha * [- (m*W - w) H1r + dV/dr + W sinth dH3/dr]
        */

        E1d_r = -(- (mm * W[ind] - omega_BS) * H1r[ind] + dV_dr[ind] + W[ind] * sinth * dH3_dr[ind]) / alph * harm_im;
        E2d_r =  (- (mm * W[ind] - omega_BS) * H1r[ind] + dV_dr[ind] + W[ind] * sinth * dH3_dr[ind]) / alph * harm_re;

        // E_th
        /*
          E_th = i * e^{i(m phi - w t)} / alpha * [- (m*W - w) * H2 + dV/dth + W * d(H3 * sinth)/dth]
        */
        E1d_th = -(- (mm * W[ind] - omega_BS) * H2[ind] + dV_dth[ind] + W[ind] * (sinth * dH3_dth[ind] + costh * H3[ind])) / alph * harm_im;
        E2d_th =  (- (mm * W[ind] - omega_BS) * H2[ind] + dV_dth[ind] + W[ind] * (sinth * dH3_dth[ind] + costh * H3[ind])) / alph * harm_re;

        // E_ph / sinth
        /*
          E_ph = - (m * V + w * H3 * sinth) / alpha * e^{i(m phi - w t)}

          We include here the division by sinth which arises when computing E^ph.
          
          on the axis sinth=0, we need to regularize the division by sinth
          We have V(theta=0,pi) = 0, so with l'Hôpital's rule
          V / sinth ~ \pm dV/dth    for theta = 0, pi resp.
        */
        if (fabs(sinth) < 1e-8) {
          // TODO: see how to deal with this if now we allow rr==0
          const CCTK_INT zsign = (costh>=0) ? 1 : -1; // costh==0 shouldn't happen on the axis for a grid point, this would mean rr==0 too...

          E1d_ph_o_sinth = - (mm * zsign * dV_dth[ind] + omega_BS * H3[ind]) / alph * harm_re;
          E2d_ph_o_sinth = - (mm * zsign * dV_dth[ind] + omega_BS * H3[ind]) / alph * harm_im;
        
        } else {
          E1d_ph_o_sinth = - (mm * V[ind] / sinth + omega_BS * H3[ind]) / alph * harm_re;
          E2d_ph_o_sinth = - (mm * V[ind] / sinth + omega_BS * H3[ind]) / alph * harm_im;
        }


        // E^i components
        // Spherical auxiliaries
        
        // E^r/r = e^{-2*F1} * E_r / r
        const CCTK_REAL E1u_r_o_r = exp(-2*F1[ind]) * E1d_r / rr;
        const CCTK_REAL E2u_r_o_r = exp(-2*F1[ind]) * E2d_r / rr;

        // E^th = e^{-2*F1} / r^2 * E_th
        const CCTK_REAL E1u_th = exp(-2*F1[ind]) / rr2 * E1d_th;
        const CCTK_REAL E2u_th = exp(-2*F1[ind]) / rr2 * E2d_th;

        // E^ph = e^{-2*F2} / r^2 / sinth^2 * E_ph
        // We compute r * sinth * E^ph. The other division by sinth is managed with E_ph above 
        const CCTK_REAL rsinthE1u_ph = exp(-2*F2[ind]) / rr * E1d_ph_o_sinth;
        const CCTK_REAL rsinthE2u_ph = exp(-2*F2[ind]) / rr * E2d_ph_o_sinth;


        // Finally Cartesian components
        // E^x
        E1x[ind] = x1 * E1u_r_o_r + z1 * cosph * E1u_th - sinph * rsinthE1u_ph;
        E2x[ind] = x1 * E2u_r_o_r + z1 * cosph * E2u_th - sinph * rsinthE2u_ph;

        // E^y
        E1y[ind] = y1 * E1u_r_o_r + z1 * sinph * E1u_th + cosph * rsinthE1u_ph;
        E2y[ind] = y1 * E2u_r_o_r + z1 * sinph * E2u_th + cosph * rsinthE2u_ph;

        // E^z
        E1z[ind] = z1 * E1u_r_o_r - rho * E1u_th;
        E2z[ind] = z1 * E2u_r_o_r - rho * E2u_th;




        // zero-initialize constraint damping variable Zeta
        Zeta1[ind] = 0;
        Zeta2[ind] = 0;


        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "ProcaBS")) {
          alp[ind] = alph;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // shift
        if (CCTK_EQUALS(initial_shift, "ProcaBS")) {
          betax[ind] =  W[ind] * y1;
          betay[ind] = -W[ind] * x1;
          betaz[ind] =  0.;
        }

      } /* for i */
    }   /* for j */
  }     /* for k */

  free(F1); free(F2); free(F0); free(W);
  free(H1r); free(H2); free(H3); free(V);
  free(dW_dr); free(dW_dth);
  free(dH3_dr); free(dH3_dth);
  free(dV_dr); free(dV_dth);

  return;
}
