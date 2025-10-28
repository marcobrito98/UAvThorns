
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

void check_nan_or_inf(const char* var_name, double value);


void UAv_IDProcaBSBH(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  // WARNING: rotating stars not being proof-tested yet (output of warning in ParamCheck)
  // TODO: Remove when proof-tested.
  // What needs to be checked if mostly the behavior/regularity of W_1 and K_ij
  
  // TODO: Is there a way to write the Proca fields in a regular way 
  //       (without 1/rr_1 in particular, but also dealing with the axis...)
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


  // To take care properly of z1_1=0 symmetry (i.e. theta <-> pi-theta) in interpolation
  // we need to extend the arrays to z1_1<0 values.
  // For convenience, we keep Ntheta as the number of points in the input half-space.

  NF = NX * (2*Ntheta - 1);

  CCTK_REAL *F1_extd, *F2_extd, *F0_extd, *H2_extd, *H3_extd, *V_extd;
  F1_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  F2_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  F0_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  H2_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  H3_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  V_extd     = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // We'll use A_r (H1r_1) rather than the input H1_in
  // Notation here: A ~ H1r_1 dr + ... = H1_in/r dr + ...
  CCTK_REAL *H1r_extd;
  H1r_extd          = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // now we need to take the derivatives of the Wbar function
  // Then we convert to W_1 and store the values

  CCTK_REAL *W_extd, *dW_dr_extd, *dW_dth_extd;
  W_extd        = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // Same for H3_1 and V_1
  
  CCTK_REAL *dH3_dr_extd, *dH3_dth_extd;
  dH3_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH3_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  
  CCTK_REAL *dV_dr_extd, *dV_dth_extd, *dH1_dr_extd, *dH1_dth_extd, *dH2_dr_extd, *dH2_dth_extd;
  dV_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dV_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH1_dr_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH1_dth_extd  = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH2_dr_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH2_dth_extd  = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // New: same for F0_1, F1_1, F2_1
  CCTK_REAL *dF0_dr_extd, *dF1_dr_extd, *dF2_dr_extd, *dF0_dth_extd, *dF1_dth_extd, *dF2_dth_extd;
  dF0_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dF1_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dF2_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dF0_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dF1_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dF2_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));


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

  // First loop on z1_1>=0 half-space (i.e. input values of 0 <= theta <= pi/2)

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
      
      // WARNING/TODO (rotating stars): Do we need to be careful with theta derivatives, like for V_1? Depending on m?
      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL H3_th = (-H3_in[indjp2] + 8 * H3_in[indjp1] - 8 * H3_in[indjm1] + H3_in[indjm2]) *
        oodth12;

      const CCTK_REAL F0_th = (-F0_in[indjp2] + 8 * F0_in[indjp1] - 8 * F0_in[indjm1] + F0_in[indjm2]) *
        oodth12;
      const CCTK_REAL F1_th = (-F1_in[indjp2] + 8 * F1_in[indjp1] - 8 * F1_in[indjm1] + F1_in[indjm2]) *
        oodth12;
      const CCTK_REAL F2_th = (-F2_in[indjp2] + 8 * F2_in[indjp1] - 8 * F2_in[indjm1] + F2_in[indjm2]) *
        oodth12;

      // Symmetries of V_1 on the axis and/or the equator can vary (theta = 0, pi/2 resp.).
      // In particular, it can occur that dV/dth != 0, which can't be captured by centered finite differences and theta symmetry.
      // Since different systems have different symmetries, we resort to non-symmetric stencils in any case
      CCTK_REAL V_th, H1_th, H2_th;
      if (jj==0) {
        // 1st derivative with 4th order accuracy (forward stencils)
        V_th = (- 25 * V_in[ind]    + 48 * V_in[indjp1] - 36 * V_in[indjp2] + 16 * V_in[indjp3] - 3 * V_in[indjp4]) *
          oodth12;
        H1_th = (- 25 * H1_in[ind]    + 48 * H1_in[indjp1] - 36 * H1_in[indjp2] + 16 * H1_in[indjp3] - 3 * H1_in[indjp4]) *
          oodth12;
        H2_th = (- 25 * H2_in[ind]    + 48 * H2_in[indjp1] - 36 * H2_in[indjp2] + 16 * H2_in[indjp3] - 3 * H2_in[indjp4]) *
          oodth12;
      } else if (jj==1) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (-  3 * V_in[indjm1] - 10 * V_in[ind]    + 18 * V_in[indjp1] -  6 * V_in[indjp2] +     V_in[indjp3]) * 
          oodth12;
        H1_th = (-  3 * H1_in[indjm1] - 10 * H1_in[ind]    + 18 * H1_in[indjp1] -  6 * H1_in[indjp2] +     H1_in[indjp3]) * 
          oodth12;
        H2_th = (-  3 * H2_in[indjm1] - 10 * H2_in[ind]    + 18 * H2_in[indjp1] -  6 * H2_in[indjp2] +     H2_in[indjp3]) * 
          oodth12;
      } else if (jj==Ntheta-2) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (   3 * V_in[indjp1] + 10 * V_in[ind]    - 18 * V_in[indjm1] +  6 * V_in[indjm2] -     V_in[indjm3]) * 
          oodth12;
        H1_th = (   3 * H1_in[indjp1] + 10 * H1_in[ind]    - 18 * H1_in[indjm1] +  6 * H1_in[indjm2] -     H1_in[indjm3]) * 
          oodth12;
        H2_th = (   3 * H2_in[indjp1] + 10 * H2_in[ind]    - 18 * H2_in[indjm1] +  6 * H2_in[indjm2] -     H2_in[indjm3]) * 
          oodth12;
      } else if (jj==Ntheta-1) {
        // 1st derivative with 4th order accuracy (backward stencils)
        V_th = (  25 * V_in[ind]    - 48 * V_in[indjm1] + 36 * V_in[indjm2] - 16 * V_in[indjm3] + 3 * V_in[indjm4]) *
          oodth12;
        H1_th = (  25 * H1_in[ind]    - 48 * H1_in[indjm1] + 36 * H1_in[indjm2] - 16 * H1_in[indjm3] + 3 * H1_in[indjm4]) *
          oodth12;
        H2_th = (  25 * H2_in[ind]    - 48 * H2_in[indjm1] + 36 * H2_in[indjm2] - 16 * H2_in[indjm3] + 3 * H2_in[indjm4]) *
          oodth12;
      } else {
        // 1st derivative with 4th order accuracy (centered stencils)
        V_th = (-V_in[indjp2] + 8 * V_in[indjp1] - 8 * V_in[indjm1] + V_in[indjm2]) *
          oodth12;
        H1_th = (-H1_in[indjp2] + 8 * H1_in[indjp1] - 8 * H1_in[indjm1] + H1_in[indjm2]) *
          oodth12;
        H2_th = (-H2_in[indjp2] + 8 * H2_in[indjp1] - 8 * H2_in[indjm1] + H2_in[indjm2]) *
          oodth12;
      }

      CCTK_REAL Wbar_X, H3_X, V_X,F0_X, F1_X, F2_X,H2_X; // radial derivatives
      CCTK_REAL Wbar_XX = 0.; // Used for r=0 (i==0), if Wbar_r_power == 2.
      CCTK_REAL H1_X = 0.;    // Used for r=0 (i==0), due to H1_in/r.

      /*
      Regarding finite differencing orders: for Scalar BS, plotting W_1 and dW_dr_1, there were small discontinuities near r=0
      during tests with the previous 2nd order accuracy for i==0 and i==1.
      Those vanish when moving to 4th order accuracy.

      For i==NX-1 and i==NX-2, we keep 2nd order for now. The issue is not appearing as clearly,
      and they represent points which are physically far, so maybe better to keep the computation more local.
      */

      if (i == 0) {
        /* For the Boson Star, there's no issue, dWbar/dX != 0 at X==0, and x*gamma and r coordinates coincide. */

        // 1st derivative with 4th order accuracy (forward stencils)
        Wbar_X =(- 25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] + 16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) * oodX12;

        // 1st derivative with 4th order accuracy (forward stencils)
        H3_X =(- 25 * H3_in[ind] + 48 * H3_in[indip1] - 36 * H3_in[indip2] + 16 * H3_in[indip3] - 3 * H3_in[indip4]) * oodX12;
        
        // 1st derivative with 4th order accuracy (forward stencils)
        V_X =(- 25 * V_in[ind] + 48 * V_in[indip1] - 36 * V_in[indip2] + 16 * V_in[indip3] - 3 * V_in[indip4]) * oodX12;

        H2_X =(- 25 * H2_in[ind] + 48 * H2_in[indip1] - 36 * H2_in[indip2] + 16 * H2_in[indip3] - 3 * H2_in[indip4]) * oodX12;

        // 1st derivative with 4th order accuracy (forward stencils)
        F0_X =(- 25 * F0_in[ind] + 48 * F0_in[indip1] - 36 * F0_in[indip2] + 16 * F0_in[indip3] - 3 * F0_in[indip4]) * oodX12;
        // 1st derivative with 4th order accuracy (forward stencils)
        F1_X =(- 25 * F1_in[ind] + 48 * F1_in[indip1] - 36 * F1_in[indip2] + 16 * F1_in[indip3] - 3 * F1_in[indip4]) * oodX12;
        // 1st derivative with 4th order accuracy (forward stencils)
        F2_X =(- 25 * F2_in[ind] + 48 * F2_in[indip1] - 36 * F2_in[indip2] + 16 * F2_in[indip3] - 3 * F2_in[indip4]) * oodX12;



        // Special care at r=0

        // H1_X required for H1r_1 computed from H1_in/r
        // 1st derivative with 4th order accuracy (forward stencils)
        H1_X =(- 25 * H1_in[ind] + 48 * H1_in[indip1] - 36 * H1_in[indip2] + 16 * H1_in[indip3] - 3 * H1_in[indip4]) * oodX12;


        if (Wbar_r_power == 2) {
          // If Wbar = r^2 * W_1, to compute W_1(r=0), we need to compute Wbar_XX.
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

        H2_X = (- 3 * H2_in[indim1] - 10 * H2_in[ind] + 18 * H2_in[indip1] - 6 * H2_in[indip2] + H2_in[indip3]) * oodX12;

        H1_X = (- 3 * H1_in[indim1] - 10 * H1_in[ind] + 18 * H1_in[indip1] - 6 * H1_in[indip2] + H1_in[indip3]) * oodX12;

        // 1st derivative, 4th order accuracy
        F0_X = (- 3 * F0_in[indim1] - 10 * F0_in[ind] + 18 * F0_in[indip1] - 6 * F0_in[indip2] + F0_in[indip3]) * oodX12;
        // 1st derivative, 4th order accuracy
        F1_X = (- 3 * F1_in[indim1] - 10 * F1_in[ind] + 18 * F1_in[indip1] - 6 * F1_in[indip2] + F1_in[indip3]) * oodX12;
        // 1st derivative, 4th order accuracy
        F2_X = (- 3 * F2_in[indim1] - 10 * F2_in[ind] + 18 * F2_in[indip1] - 6 * F2_in[indip2] + F2_in[indip3]) * oodX12;

      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4*Wbar_in[indim1] + 3*Wbar_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        H3_X = (H3_in[indim2] - 4*H3_in[indim1] + 3*H3_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        V_X = (V_in[indim2] - 4*V_in[indim1] + 3*V_in[ind]) * 0.5 * oodX;

        H2_X = (H2_in[indim2] - 4*H2_in[indim1] + 3*H2_in[ind]) * 0.5 * oodX;

        H1_X = (H1_in[indim2] - 4*H1_in[indim1] + 3*H1_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        F0_X = (F0_in[indim2] - 4*F0_in[indim1] + 3*F0_in[ind]) * 0.5 * oodX;
        // 1st derivative with 2nd order accuracy (backward stencils)
        F1_X = (F1_in[indim2] - 4*F1_in[indim1] + 3*F1_in[ind]) * 0.5 * oodX;
        // 1st derivative with 2nd order accuracy (backward stencils)
        F2_X = (F2_in[indim2] - 4*F2_in[indim1] + 3*F2_in[ind]) * 0.5 * oodX;

      } else if (i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;
        
        // 1st derivative with 2nd order accuracy (central stencils)
        H3_X = (-H3_in[indim1] + H3_in[indip1]) * 0.5 * oodX;
        
        // 1st derivative with 2nd order accuracy (central stencils)
        V_X = (-V_in[indim1] + V_in[indip1]) * 0.5 * oodX;

        H2_X = (-H2_in[indim1] + H2_in[indip1]) * 0.5 * oodX;

        H1_X = (-H1_in[indim1] + H1_in[indip1]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (central stencils)
        F0_X = (-F0_in[indim1] + F0_in[indip1]) * 0.5 * oodX;
        // 1st derivative with 2nd order accuracy (central stencils)
        F1_X = (-F1_in[indim1] + F1_in[indip1]) * 0.5 * oodX;
        // 1st derivative with 2nd order accuracy (central stencils)
        F2_X = (-F2_in[indim1] + F2_in[indip1]) * 0.5 * oodX;

      } else {
        // 4th order accurate stencils
        Wbar_X    = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
        
        // 4th order accurate stencils
        H3_X    = (-H3_in[indip2] + 8 * H3_in[indip1] - 8 * H3_in[indim1] + H3_in[indim2]) * oodX12;
        
        // 4th order accurate stencils
        V_X    = (-V_in[indip2] + 8 * V_in[indip1] - 8 * V_in[indim1] + V_in[indim2]) * oodX12;

        H2_X    = (-H2_in[indip2] + 8 * H2_in[indip1] - 8 * H2_in[indim1] + H2_in[indim2]) * oodX12;

        H1_X    = (-H1_in[indip2] + 8 * H1_in[indip1] - 8 * H1_in[indim1] + H1_in[indim2]) * oodX12;

        // 4th order accurate stencils
        F0_X    = (-F0_in[indip2] + 8 * F0_in[indip1] - 8 * F0_in[indim1] + F0_in[indim2]) * oodX12;
        // 4th order accurate stencils
        F1_X    = (-F1_in[indip2] + 8 * F1_in[indip1] - 8 * F1_in[indim1] + F1_in[indim2]) * oodX12;
        // 4th order accurate stencils
        F2_X    = (-F2_in[indip2] + 8 * F2_in[indip1] - 8 * F2_in[indim1] + F2_in[indim2]) * oodX12;

      
      }

      // From the X coordinate used in the input files to the r coordinate (coincides with x1 for the Boson Star, rH=0).
      // We also do the conversion from Wbar to W_1 here, and H1_in to H1r_1, to tackle r = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  r == 0
      if (i == 0) {

        // At r=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_extd[ind]    = 0.; 
        dW_dth_extd[ind]   = 0.;
        
        // At X==0 (rr_1==0), dXdr = 1/C0
        dH3_dr_extd[ind]       = H3_X / C0;
        dV_dr_extd[ind]        =  V_X / C0;
        dH1_dr_extd[ind]       = H1_X / C0;
        dH2_dr_extd[ind]       = H2_X / C0;
        dF0_dr_extd[ind]       = F0_X / C0;
        dF1_dr_extd[ind]       = F1_X / C0;
        dF2_dr_extd[ind]       = F2_X / C0;

        // For W_1 we need more care depending on the power
        switch (Wbar_r_power)
        {
        case 0:   // Wbar = W_1
          W_extd[ind]        = Wbar_in[ind];
          break;
        
        case 1:   // Wbar = r * W_1
          /*
          dWbar/dr = W_1 + r * dW/dr
                   = W_1 + 0 * 0     at r=0
          
          dWbar/dr = dWbar/dX * dX/dr
          dX/dr = C/(C+r)^2 = 1/C  at r=0
          */
          W_extd[ind]        = Wbar_X / C0;
          break;
        
        case 2:   // Wbar = r^2 * W_1
          /*
          dWbar/dr   = 2r * W_1 + r^2 * dW/dr
          d2Wbar/dr2 = 2  * W_1 + 4r  * dW/dr + r^2 * d2W/dr2
                     = 2  * W_1 + 0 + 0        at r=0

          d2Wbar/dr2 = d2Wbar/dX2 * (dX/dr)^2 + dWbar/dX * d2X/dr2
          dX/dr   =   C/(C+r)^2 =  1/C      at r=0
          d2X/dr2 = -2C/(C+r)^3 = -2/C^2    at r=0

          W_1 (r=0) = 1/C^2 * [1/2 * d2Wbar/dX^2 (X=0)  -  dWbar/dX (X=0)]
          */
          W_extd[ind]        = (0.5 * Wbar_XX - Wbar_X)/(C0*C0);
          break;
        
        default:  // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }

        // From H1_in/r to H1r_1 (idem Wbar above, case 1)
        H1r_extd[ind] = H1_X / C0;
      }

      // We need to be careful at X == 1 (r == infty) for radial derivatives (coordinate change is singular)
      else if (i == NX - 1) {
        // W_1 -> 0 for r -> infty
        W_extd[ind]        = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction paper) also gives:
        dW_dr_extd[ind]    = 0.; 
        dW_dth_extd[ind]   = 0.; 


        dH3_dr_extd[ind]    = 0.;
        dV_dr_extd[ind]    = 0.;
        dH1_dr_extd[ind]   = 0.;
        dH2_dr_extd[ind]   = 0.;
        dF0_dr_extd[ind]   = 0.;
        dF1_dr_extd[ind]   = 0.;
        dF2_dr_extd[ind]   = 0.;

        H1r_extd[ind] = 0.; // A_r = 0 at infinity

      } else {

        const CCTK_REAL rr_1 = C0*lX/(1. - lX);

        // corresponding derivatives
        // const CCTK_REAL dXdr = 1./(C0 + rr_1) - rr_1/((C0 + rr_1)*(C0 + rr_1));
        const CCTK_REAL dXdr = C0/((C0 + rr_1)*(C0 + rr_1));

        const CCTK_REAL Wbar_r = dXdr * Wbar_X;

        dH3_dr_extd[ind]       = dXdr * H3_X;
        dV_dr_extd[ind]        = dXdr * V_X;
        dH1_dr_extd[ind]       = dXdr * H1_X;
        dH2_dr_extd[ind]       = dXdr * H2_X;
        dF0_dr_extd[ind]       = dXdr * F0_X;
        dF1_dr_extd[ind]       = dXdr * F1_X;
        dF2_dr_extd[ind]       = dXdr * F2_X;
        
        // Now translate from Wbar to W_1
        switch (Wbar_r_power) // We could put a generic power for the computation here I guess...
        {
        case 0:   // Wbar = W_1
          W_extd[ind]        = Wbar_in[ind];
          dW_dr_extd[ind]    = Wbar_r;
          dW_dth_extd[ind]   = Wbar_th;
          break;
        
        case 1:   // Wbar = r * W_1
          W_extd[ind]        = Wbar_in[ind] / rr_1;
          dW_dr_extd[ind]    = (Wbar_r - W_extd[ind]) / rr_1; // dW/dr  =  1/r * dWbar/dr - Wbar / r^2  =  (dWbar/dr - W_1) / r
          dW_dth_extd[ind]   = Wbar_th / rr_1;
          break;
        
        case 2: ; // Wbar = r^2 * W_1
          // empty statement after case to prevent compilation error on some gcc versions...
          const CCTK_REAL rr2_1 = rr_1*rr_1;
          W_extd[ind]        = Wbar_in[ind] / rr2_1;
          dW_dr_extd[ind]    = Wbar_r / rr2_1 - 2 * W_extd[ind] / rr_1; // dW/dr  =  1/r^2 * dWbar/dr - 2 * Wbar / r^3  =  1/r^2 * dWbar/dr - 2 * W_1 / r
          dW_dth_extd[ind]   = Wbar_th / rr2_1;
          break;
        
        default:  // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }

        // From H1_in/r to H1r
        H1r_extd[ind] = H1_in[ind] / rr_1;
        
      } // if/else i==...
      
      dH3_dth_extd[ind]     = H3_th;
      dV_dth_extd[ind]      = V_th;
      dH1_dth_extd[ind]     = H1_th;
      dH2_dth_extd[ind]     = H2_th;
      dF0_dth_extd[ind]     = F0_th;
      dF1_dth_extd[ind]     = F1_th;
      dF2_dth_extd[ind]     = F2_th;
    
      // fprintf (debugfile, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f ", 
      //             W_extd[ind], dW_dr_extd[ind], dW_dth_extd[ind],
      //             H3_extd[ind], dH3_dr_extd[ind], dH3_dth_extd[ind],
      //             V_extd[ind], dV_dr_extd[ind], dV_dth_extd[ind],
      //             H1r_extd[ind]);
    } // for i
    // fprintf (debugfile, "\n");
  } // for jj

  // fclose(debugfile);


  // Second loop on z1_1<0 half-space (completion by symmetry)

  // Even parity: F1_1, F2_1, F0_1, W_1; r derivatives of even functions; theta derivatives of odd functions
  // Odd parity:  r derivatives of odd functions; theta derivatives of even functions
  // A even (default): H1, H3_1 and V_1 are even, and H2_1 is odd.
  // A odd           : H1, H3_1 and V_1 are  odd, and H2_1 is even.

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
      dF0_dr_extd[ind]   = dF0_dr_extd[indsym];
      dF1_dr_extd[ind]   = dF1_dr_extd[indsym];
      dF2_dr_extd[ind]   = dF2_dr_extd[indsym];
      
      // Odd
      dW_dth_extd[ind]      = - dW_dth_extd[indsym];
      dF0_dth_extd[ind]     = - dF0_dth_extd[indsym];
      dF1_dth_extd[ind]     = - dF1_dth_extd[indsym];
      dF2_dth_extd[ind]     = - dF2_dth_extd[indsym];
      
      // Vector potential
      
      H1r_extd[ind]      = H1_z_sign * H1r_extd[indsym];
      H2_extd[ind]       = H2_z_sign *  H2_extd[indsym];
      H3_extd[ind]       = H3_z_sign *  H3_extd[indsym];
      V_extd[ind]        =  V_z_sign *   V_extd[indsym];
      
      dH3_dr_extd[ind]   = H3_z_sign * dH3_dr_extd[indsym];
      dV_dr_extd[ind]    =  V_z_sign *  dV_dr_extd[indsym];
      dH1_dr_extd[ind]   = H1_z_sign * dH1_dr_extd[indsym];
      dH2_dr_extd[ind]   = H2_z_sign * dH2_dr_extd[indsym];
      

      dH3_dth_extd[ind]  = - H3_z_sign * dH3_dth_extd[indsym];
      dV_dth_extd[ind]   = -  V_z_sign *  dV_dth_extd[indsym];
      dH1_dth_extd[ind]  = - H1_z_sign * dH1_dth_extd[indsym];
      dH2_dth_extd[ind]  = - H2_z_sign * dH2_dth_extd[indsym];
      

      } // for i
  } // for jj


  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *X_g, *theta_g;
  X_g     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  theta_g = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  const CCTK_REAL bs_v2 = bs_v * bs_v;
  const CCTK_REAL gamma2 = 1. / (1. - bs_v2);
  const CCTK_REAL gamma = sqrt(gamma2);

  for (int k = 0; k < cctk_lsh[2]; ++k) { //tenho de por aqui os gammas? O resultado parece melhor
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1_1  = x[ind] - x0;
        const CCTK_REAL y1_1  = y[ind] - y0;
        const CCTK_REAL z1_1  = z[ind] - z0;

        const CCTK_REAL rr2_1 = x1_1*x1_1*gamma2 + y1_1*y1_1 + z1_1*z1_1;

        CCTK_REAL rr_1  = sqrt(rr2_1);
        /* For the Boson Star, x, r and R coordinates coincide (rH=0). */
	/* note that there are divisions by rr_1 in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */
        
        // From r to the X radial coordinate (used in input files)
        const CCTK_REAL lX = rr_1 / (C0 + rr_1);

        const CCTK_REAL ltheta = rr_1 < 1e-16 ? 0 : acos( z1_1/rr_1 );    // There should be at most one point in the grid with rr_1~0. Not sure about the threshold.

        X_g[ind]     = lX;
        theta_g[ind] = ltheta;
      }
    }
  }

  /* now for the interpolation */

  const CCTK_INT N_dims  = 2;   // 2-D interpolation

  const CCTK_INT N_input_arrays  = 24;
  const CCTK_INT N_output_arrays = 24;

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
  input_array_type_codes[14]= CCTK_VARIABLE_REAL;
  input_array_type_codes[15]= CCTK_VARIABLE_REAL;
  input_array_type_codes[16]= CCTK_VARIABLE_REAL;
  input_array_type_codes[17]= CCTK_VARIABLE_REAL;
  input_array_type_codes[18]= CCTK_VARIABLE_REAL;
  input_array_type_codes[19]= CCTK_VARIABLE_REAL;
  input_array_type_codes[20]= CCTK_VARIABLE_REAL;
  input_array_type_codes[21]= CCTK_VARIABLE_REAL;
  input_array_type_codes[22]= CCTK_VARIABLE_REAL;
  input_array_type_codes[23]= CCTK_VARIABLE_REAL;

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
  input_arrays[14]= (const void *) dF0_dr_extd;
  input_arrays[15]= (const void *) dF1_dr_extd;
  input_arrays[16]= (const void *) dF2_dr_extd;
  input_arrays[17]= (const void *) dF0_dth_extd;
  input_arrays[18]= (const void *) dF1_dth_extd;
  input_arrays[19]= (const void *) dF2_dth_extd;
  input_arrays[20]= (const void *) dH1_dr_extd;
  input_arrays[21]= (const void *) dH1_dth_extd;
  input_arrays[22]= (const void *) dH2_dr_extd;
  input_arrays[23]= (const void *) dH2_dth_extd;


  /* output arrays */
  void *output_arrays[N_output_arrays];
  CCTK_INT output_array_type_codes[N_output_arrays];
  CCTK_REAL *F1_1, *F2_1, *F0_1, *W_1, *H1r_1, *H2_1, *H3_1, *V_1;
  CCTK_REAL *dW_dr_1, *dW_dth_1;
  CCTK_REAL *dH3_dr_1, *dH3_dth_1;
  CCTK_REAL *dV_dr_1, *dV_dth_1, *dH1_dr_1, *dH1_dth_1, *dH2_dr_1, *dH2_dth_1;
  CCTK_REAL *dF0_dr_1, *dF1_dr_1, *dF2_dr_1, *dF0_dth_1, *dF1_dth_1, *dF2_dth_1;

  F1_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F2_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F0_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  W_1           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dr_1       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dth_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H1r_1         = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H2_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H3_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  V_1           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dr_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dth_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dr_1       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dth_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH1_dr_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH1_dth_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dF0_dr_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dF1_dr_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dF2_dr_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dF0_dth_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dF1_dth_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dF2_dth_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH2_dr_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH2_dth_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  

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
  output_array_type_codes[14]= CCTK_VARIABLE_REAL;
  output_array_type_codes[15]= CCTK_VARIABLE_REAL;
  output_array_type_codes[16]= CCTK_VARIABLE_REAL;
  output_array_type_codes[17]= CCTK_VARIABLE_REAL;
  output_array_type_codes[18]= CCTK_VARIABLE_REAL;
  output_array_type_codes[19]= CCTK_VARIABLE_REAL;
  output_array_type_codes[20]= CCTK_VARIABLE_REAL;
  output_array_type_codes[21]= CCTK_VARIABLE_REAL;
  output_array_type_codes[22]= CCTK_VARIABLE_REAL;
  output_array_type_codes[23]= CCTK_VARIABLE_REAL;


  output_arrays[0] = (void *) F1_1;
  output_arrays[1] = (void *) F2_1;
  output_arrays[2] = (void *) F0_1;
  output_arrays[3] = (void *) W_1;
  output_arrays[4] = (void *) dW_dr_1;
  output_arrays[5] = (void *) dW_dth_1;
  output_arrays[6] = (void *) H1r_1;
  output_arrays[7] = (void *) H2_1;
  output_arrays[8] = (void *) H3_1;
  output_arrays[9] = (void *) V_1;
  output_arrays[10]= (void *) dH3_dr_1;
  output_arrays[11]= (void *) dH3_dth_1;
  output_arrays[12]= (void *) dV_dr_1;
  output_arrays[13]= (void *) dV_dth_1;
  output_arrays[14]= (void *) dF0_dr_1;
  output_arrays[15]= (void *) dF1_dr_1;
  output_arrays[16]= (void *) dF2_dr_1;
  output_arrays[17]= (void *) dF0_dth_1;
  output_arrays[18]= (void *) dF1_dth_1;
  output_arrays[19]= (void *) dF2_dth_1;
  output_arrays[20]= (void *) dH1_dr_1;
  output_arrays[21]= (void *) dH1_dth_1;
  output_arrays[22]= (void *) dH2_dr_1;
  output_arrays[23]= (void *) dH2_dth_1;


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
  free(dH1_dr_extd); free(dH1_dth_extd);
  free(dF0_dr_extd); free(dF1_dr_extd); free(dF2_dr_extd);
  free(dF0_dth_extd); free(dF1_dth_extd); free(dF2_dth_extd);


  /* printf("F1_1 = %g\n", F1_1[0]); */
  /* printf("F2_1 = %g\n", F2_1[0]); */
  /* printf("F0_1 = %g\n", F0_1[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W_1 = %g\n", W_1[0]); */


  /* now we finally write the metric and all 3+1 quantities. first we write the
     3-metric and extrinsic curvature, then Proca fields, then lapse and shift */
  /* For the Boson Star, in order to avoid unneeded regularizations and divisions,
      we express K_ij in terms of dW/drho and dW/dz.
      For points close to the axis and the origin, K_ij = 0.
  */

  const CCTK_REAL tt = cctk_time;

  const CCTK_REAL coswt = cos(omega_BS * tt*gamma);
  const CCTK_REAL sinwt = sin(omega_BS * tt*gamma);


  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        //Proca Star A

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1_1  = x[ind] - x0;
        const CCTK_REAL y1_1  = y[ind] - y0;
        const CCTK_REAL z1_1  = z[ind] - z0;

        // For the Boson Star, r = R, no coordinate change needed.
        CCTK_REAL rr2_1 = x1_1*x1_1*gamma2 + y1_1*y1_1 + z1_1*z1_1;
        // if( rr2_1 < pow( eps_r, 2 ) ) {
        // rr2_1 = pow( eps_r, 2 );
        // }
        const CCTK_REAL rr_1  = sqrt(rr2_1);
	      /* note that there are divisions by rr_1 in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        CCTK_REAL rho2_1 = x1_1*x1_1*gamma2 + y1_1*y1_1;
        // if( rho2_1 < pow( eps_r, 2 ) ){
        // rho2_1 = pow( eps_r, 2 );
        // }
        const CCTK_REAL rho_1  = sqrt(rho2_1);

        

        const CCTK_REAL costh_1  = z1_1/rr_1;
        const CCTK_REAL costh2_1 = costh_1*costh_1;
        /*
          For some grid points actually on the axis, it occurred that costh_1 = 1-1e-16, resulting in sinth_1 ~ 1.5e-8 instead of 0.
          Thus we force it in that case. 
          Even if there is a legit grid point such that theta ~ a few 1e-8, it should mean RR >> rho_1 and the axis treatment should be fine.
        */
        CCTK_REAL sinth_1, sinth2_1;
        if (1-costh2_1 < 1e-15) {
          sinth2_1 = 0.;
          sinth_1  = 0.;
        } else {
          sinth2_1 = 1. - costh2_1;
          sinth_1  = sqrt(sinth2_1);
        }

        const CCTK_REAL d_sinph_dx = -x1_1 * gamma * y1_1 / (rho2_1 * sqrt(rho2_1));
        const CCTK_REAL d_sinph_dy =  (x1_1*gamma * x1_1*gamma) / (rho2_1 * sqrt(rho2_1));
        const CCTK_REAL d_sinph_dz = 0;


        const CCTK_REAL d_cosph_dx = (y1_1 * y1_1) / (rho2_1 * sqrt(rho2_1));
        const CCTK_REAL d_cosph_dy = -(x1_1*gamma * y1_1) / (rho2_1 * sqrt(rho2_1));
        const CCTK_REAL d_cosph_dz = 0;
        

        const CCTK_REAL ph_1 = atan2(y1_1, x1_1*gamma);
        // If x1_1=y1_1=0, should return 0? The other metric functions should vanish anyway to make sure that this doesn't matter,
        // but can this lead to nan depending on the C implementation?

        const CCTK_REAL cosph  = cos(ph_1);
        const CCTK_REAL sinph  = sin(ph_1);

        const CCTK_REAL cosmph = cos(mm*ph_1);
        const CCTK_REAL sinmph = sin(mm*ph_1);

        // const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        // const CCTK_REAL psi2_1 = sqrt(psi4_1);
        // const CCTK_REAL psi1_1 = sqrt(psi2_1);

        const CCTK_REAL h_rho2_1 = exp(2. * (F2_1[ind] - F1_1[ind])) - 1.;

        const CCTK_REAL R_x_1    = x1_1*gamma/rr_1 ;
        const CCTK_REAL R_y_1    = y1_1/rr_1 ;
        const CCTK_REAL R_z_1    = z1_1/rr_1 ;

        const CCTK_REAL th_x_1   = costh_1*R_x_1/rho_1 ;
        const CCTK_REAL th_y_1   = costh_1*R_y_1/rho_1 ;
        const CCTK_REAL th_z_1   = -rho_1/rr2_1;

        const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        const CCTK_REAL psi2_1 = sqrt(psi4_1);
        const CCTK_REAL psi1_1 = sqrt(psi2_1);


        const CCTK_REAL dF1_1_dx = dF1_dr_1[ind]*R_x_1 + dF1_dth_1[ind]*th_x_1;
        const CCTK_REAL dF1_1_dy = dF1_dr_1[ind]*R_y_1 + dF1_dth_1[ind]*th_y_1;
        const CCTK_REAL dF1_1_dz = dF1_dr_1[ind]*R_z_1 + dF1_dth_1[ind]*th_z_1;
        const CCTK_REAL dF2_1_dx = dF2_dr_1[ind]*R_x_1 + dF2_dth_1[ind]*th_x_1;
        const CCTK_REAL dF2_1_dy = dF2_dr_1[ind]*R_y_1 + dF2_dth_1[ind]*th_y_1;
        const CCTK_REAL dF2_1_dz = dF2_dr_1[ind]*R_z_1 + dF2_dth_1[ind]*th_z_1;
        const CCTK_REAL dF0_1_dx = dF0_dr_1[ind]*R_x_1 + dF0_dth_1[ind]*th_x_1;
        const CCTK_REAL dF0_1_dy = dF0_dr_1[ind]*R_y_1 + dF0_dth_1[ind]*th_y_1;
        const CCTK_REAL dF0_1_dz = dF0_dr_1[ind]*R_z_1 + dF0_dth_1[ind]*th_z_1;


        CCTK_REAL G[4][4]; // temporary storage for the 4-metric
        CCTK_REAL G3_inv[4][4]; // temporary storage for the inverse of the 3-metric
        CCTK_REAL Gb[4][4]; // temporary storage for the boosted metric
        CCTK_REAL gammaA_inv[4][4]; // temporary storage for the inverse of the boosted 3-metric

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            G[a][b] = 0.0;
            gammaA_inv[a][b] = 0.0;
            G3_inv[a][b] = 0.0;
            Gb[a][b] = 0.0;
          }
        }

        G[0][0] = - exp(2. * F0_1[ind]);
        G[1][1] = psi4_1 * (1. + h_rho2_1* sinph * sinph);
        G[1][2] = -psi4_1 * h_rho2_1 * sinph * cosph;
        G[2][1] = G[1][2];
        G[2][2] = psi4_1 * (1. + h_rho2_1 * cosph * cosph);
        G[3][3] = psi4_1;



        G3_inv[1][1] =  (pow(x1_1*gamma,2)/exp(2. * F1_1[ind]) + pow(y1_1,2)/exp(2. * \
                        F2_1[ind]))/rho2_1;
        G3_inv[1][2] = ((exp(-2. * F1_1[ind]) - \
                       exp(-2 * F2_1[ind]))*x1_1*gamma*y1_1)/rho2_1;
        G3_inv[1][3] =  0.0;
        G3_inv[2][1] = G3_inv[1][2];
        G3_inv[2][2] =  (pow(x1_1*gamma,2)/exp(2. * F2_1[ind]) + pow(y1_1,2)/exp(2. * \
                        F1_1[ind]))/rho2_1;
        G3_inv[2][3] = 0;
        G3_inv[3][1] = G3_inv[1][3];
        G3_inv[3][2] = G3_inv[2][3];
        G3_inv[3][3] = exp(-2. * F1_1[ind]);


        // Derivatives of the metric functions
        CCTK_REAL dG[4][4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              dG[a][b][c] = 0.0;
            }
          }
        }
        // dG[a][b][c] = dG_ab/dx^c

        dG[1][1][1] = (2*exp(2. * F1_1[ind])*x1_1*gamma*(pow(y1_1,2) + x1_1*gamma*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF1_1_dx) + 2*exp(2. * \
                      F2_1[ind])*pow(y1_1,2)*(-x1_1*gamma + (pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dx))/pow(pow(x1_1*gamma,2) \
                      + pow(y1_1,2),2);

        dG[1][1][2] = (2*exp(2. * F1_1[ind])*pow(x1_1*gamma,2)*(-y1_1 + (pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dy) + 2*exp(2. * \
                      F2_1[ind])*y1_1*(pow(x1_1*gamma,2) + y1_1*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dy))/pow(pow(x1_1*gamma,2) \
                      + pow(y1_1,2),2);

        dG[1][1][3] = (2*(exp(2. * \
                      F1_1[ind])*pow(x1_1*gamma,2)*dF1_1_dz + exp(2. \
                      * F2_1[ind])*pow(y1_1,2)*dF2_1_dz))/(pow(\
                      x1_1*gamma,2) + pow(y1_1,2));

        dG[1][2][1] = (exp(2. * F1_1[ind])*y1_1*(-pow(x1_1*gamma,2) + pow(y1_1,2) + \
                      2*x1_1*gamma*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF1_1_dx) - exp(2. * \
                      F2_1[ind])*y1_1*(-pow(x1_1*gamma,2) + pow(y1_1,2) + 2*x1_1*gamma*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dx))/pow(pow(x1_1*gamma,2) \
                      + pow(y1_1,2),2);


        dG[1][2][2] = (exp(2. * F1_1[ind])*x1_1*gamma*(pow(x1_1*gamma,2) - pow(y1_1,2) + \
                      2*y1_1*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF1_1_dy) - exp(2. * \
                      F2_1[ind])*x1_1*gamma*(pow(x1_1*gamma,2) - pow(y1_1,2) + 2*y1_1*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dy))/pow(pow(x1_1*gamma,2) \
                      + pow(y1_1,2),2);

        dG[1][2][3] = (2*x1_1*gamma*y1_1*(exp(2. * F1_1[ind])*dF1_1_dz - exp(2. * F2_1[ind])*dF2_1_dz))/rho2_1;

        dG[2][2][1] = (2*exp(2. * F1_1[ind])*pow(y1_1,2)*(-x1_1*gamma + (pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF1_1_dx) + 2*exp(2. * \
                      F2_1[ind])*x1_1*gamma*(pow(y1_1,2) + x1_1*gamma*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dx))/pow(pow(x1_1*gamma,2) \
                      + pow(y1_1,2),2);

        dG[2][2][2] = (2*exp(2. * F1_1[ind])*y1_1*(pow(x1_1*gamma,2) + y1_1*(pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF1_1_dy) + 2*exp(2. * \
                      F2_1[ind])*pow(x1_1*gamma,2)*(-y1_1 + (pow(x1_1*gamma,2) + \
                      pow(y1_1,2))*dF2_1_dy))/pow(pow(x1_1*gamma,2) \
                      + pow(y1_1,2),2);
        dG[2][2][3] = (2*(exp(2. * \
                      F1_1[ind])*pow(y1_1,2)*dF1_1_dz + exp(2. \
                      * F2_1[ind])*pow(x1_1*gamma,2)*dF2_1_dz))/(pow(\
                      x1_1*gamma,2) + pow(y1_1,2)); //dG23_dx^i = 0

        dG[3][3][1] = 2*exp(2. * F1_1[ind])*dF1_1_dx;
        dG[3][3][2] = 2*exp(2. * F1_1[ind])*dF1_1_dy;
        dG[3][3][3] = 2*exp(2. * F1_1[ind])*dF1_1_dz;


        dG[0][0][1] = -2*exp(2. * F0_1[ind])*dF0_1_dx;
        dG[0][0][2] = -2*exp(2. * F0_1[ind])*dF0_1_dy;
        dG[0][0][3] = -2*exp(2. * F0_1[ind])*dF0_1_dz;


        //symmetries
        dG[2][1][1] = dG[1][2][1];
        dG[2][1][2] = dG[1][2][2];
        dG[2][1][3] = dG[1][2][3];
        dG[3][1][1] = dG[1][3][1];
        dG[3][1][2] = dG[1][3][2];
        dG[3][1][3] = dG[1][3][3];
        dG[3][2][1] = dG[2][3][1];
        dG[3][2][2] = dG[2][3][2];
        dG[3][2][3] = dG[2][3][3];

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              if (isnan(dG[a][b][c]) || isinf(dG[a][b][c])) {
                fprintf(stderr, "Error: dG[%d][%d][%d] is nan or inf at grid point (%lf,%lf,%lf)\n", a, b, c, x1_1, y1_1, z1_1);
              }
            }
          }
        }

        CCTK_REAL invLambda[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            invLambda[a][b] = 0.0;
          }
        }
        invLambda[0][0] = gamma;
        invLambda[0][1] = gamma*bs_v;
        invLambda[1][0] = gamma*bs_v;
        invLambda[1][1] = gamma;
        invLambda[2][2] = 1.;
        invLambda[3][3] = 1.;


        CCTK_REAL Lambda[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            Lambda[a][b] = 0.0;
          }
        }
        Lambda[0][0] = gamma;
        Lambda[0][1] = -gamma*bs_v;
        Lambda[1][0] = -gamma*bs_v;
        Lambda[1][1] = gamma;
        Lambda[2][2] = 1.;
        Lambda[3][3] = 1.;

        // Boosted metric

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
              CCTK_REAL sum = 0.0;
              for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                  sum += invLambda[mu][a]*invLambda[nu][b]*G[mu][nu];
              Gb[a][b] = sum;
          }
        }


        gammaA_inv[1][1] = (exp(2*(F0_1[ind]+F2_1[ind]))*pow(x1_1*gamma,2) + \
                           exp(2*(F0_1[ind]+F1_1[ind]))*pow(y1_1,2) - \
                           pow(bs_v,2)*exp(2*(F1_1[ind]+F2_1[ind]))*(rho2_1))/(pow(-1 + \
                           pow(bs_v,2),2)*exp(2*(F0_1[ind]+F1_1[ind]+F2_1[ind]))*(rho2_1)*pow(gamma,2));
        gammaA_inv[1][2] = ((exp(2. * F1_1[ind]) - exp(2. * F2_1[ind]))*x1_1*gamma*y1_1)/((-1 + \
                           pow(bs_v,2))*exp(2*(F1_1[ind]+F2_1[ind]))*(rho2_1)*gamma);
        gammaA_inv[1][3] = 0;
        gammaA_inv[2][1] = gammaA_inv[1][2];
        gammaA_inv[2][2] = (pow(x1_1*gamma,2)/exp(2. * F2_1[ind]) + pow(y1_1,2)/exp(2. * \
                           F1_1[ind]))/(rho2_1);
        gammaA_inv[2][3] = 0;
        gammaA_inv[3][1] = gammaA_inv[1][3];
        gammaA_inv[3][2] = gammaA_inv[2][3];
        gammaA_inv[3][3] = exp(-2. * F1_1[ind]);


        // Check for NaN or Inf in gammaA_inv
        for (int a = 1; a < 4; ++a) {
          for (int b = 1; b < 4; ++b) {
            if (isnan(gammaA_inv[a][b]) || isinf(gammaA_inv[a][b])) {
              fprintf(stderr, "Error: gammaA_inv[%d][%d] is nan or inf at grid point (%lf,%lf,%lf)\n", a, b, x1_1, y1_1, z1_1);
            }
          }
        }

        CCTK_REAL dGb[4][4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              dGb[a][b][c] = 0.0;
            }
          }
        }
        
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              CCTK_REAL sum = 0.0;
              for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                  for (int lam = 0; lam < 4; ++lam)
                    sum += invLambda[mu][a]*invLambda[nu][b]*invLambda[lam][c]*dG[mu][nu][lam];
              dGb[a][b][c] = sum;
            }
          }
        }
        // Check for NaN or Inf in dGb
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              if (isnan(dGb[a][b][c]) || isinf(dGb[a][b][c])) {
                fprintf(stderr, "Error: dGb[%d][%d][%d] is nan or inf at grid point (%lf,%lf,%lf)\n", a, b, c, x1_1, y1_1, z1_1);
              }
            }
          }
        }




        // Now we compute the 3+1 quantities
        
        // Shift
        CCTK_REAL beta1[4],betaup1[4];
        beta1[0] = 0;
        beta1[1] = Gb[0][1];
        beta1[2] = Gb[0][2];
        beta1[3] = Gb[0][3];
        betaup1[0] = 0;
        betaup1[1] = gammaA_inv[1][1]*beta1[1] + gammaA_inv[1][2]*beta1[2] + gammaA_inv[1][3]*beta1[3];
        betaup1[2] = gammaA_inv[2][1]*beta1[1] + gammaA_inv[2][2]*beta1[2] + gammaA_inv[2][3]*beta1[3];
        betaup1[3] = gammaA_inv[3][1]*beta1[1] + gammaA_inv[3][2]*beta1[2] + gammaA_inv[3][3]*beta1[3];

        // Lapse
        const CCTK_REAL alpha1 = sqrt(-Gb[0][0] + betaup1[1]*beta1[1] + betaup1[2]*beta1[2] + betaup1[3]*beta1[3]);

        // Check for NaN in beta1 and betaup1
        for (int idx = 0; idx < 4; ++idx) {
          if (isnan(beta1[idx]) || isinf(beta1[idx])) {
            fprintf(stderr, "Error: beta1[%d] is NaN at grid point (%lf,%lf,%lf)\n", idx, x1_1, y1_1, z1_1);
          }
          if (isnan(betaup1[idx]) || isinf(betaup1[idx])) {
            fprintf(stderr, "Error: betaup1[%d] is NaN at grid point (%lf,%lf,%lf)\n", idx, x1_1, y1_1, z1_1);
          }
        }

        
        check_nan_or_inf("1/alpha1", 1/alpha1);

        CCTK_REAL K_A[4][4]; // extrinsic curvature
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            K_A[a][b] = 0.0;
          }
        } //K_0\mu might not be zero but irrelevant for what i want to compute

        for (int a = 1; a < 4; ++a) { //since non rotating K=0 as a bypass
          for (int b = 1; b < 4; ++b) {
            CCTK_REAL sum1 = 0.0;
            CCTK_REAL sum2 = 0.0;
            CCTK_REAL sum3 = 0.0;
            for (int c = 1; c < 4; ++c) {
              sum1 += betaup1[c]*dGb[a][b][c];
              sum2 += betaup1[c]*dGb[b][c][a];
              sum3 += betaup1[c]*dGb[a][c][b];
            }
            K_A[a][b] = -1 / (2. * alpha1) * (dGb[a][b][0] - sum1 - (dGb[0][b][a] - sum2) - (dGb[0][a][b] - sum3));
          }
        }

        for (int a = 1; a < 4; ++a) {
          for (int b = 1; b < 4; ++b) {
            if (isnan(K_A[a][b]) || isinf(K_A[a][b])) {
              fprintf(stderr, "Error: K_{%d,%d} is nan at grid point (%lf,%lf,%lf)\n",a,b, x1_1, y1_1, z1_1);
            }
          }
        }

    //Black Hole B

      CCTK_REAL alpha0,psi1_2;

      // if (CCTK_EQUALS(bh_spin_direction, "z1_1")) { 

        CCTK_REAL x1_2  = x[ind] - x0_2;
        CCTK_REAL y1_2  = y[ind] - y0_2;
        CCTK_REAL z1_2  = z[ind] - z0_2;

        // const CCTK_REAL bh_v2 = bs_v * bs_v;
        const CCTK_REAL bh_spin2 = bh_spin*bh_spin;
        // const CCTK_REAL gamma2 = 1. / (1. - bh_v2);
        // const CCTK_REAL gamma = sqrt(gamma2);
        // const CCTK_REAL rr2_2 = gamma2*x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;
        // const CCTK_REAL rr_2  = sqrt(rr2_2);
        CCTK_REAL rr2_2 = x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;
        if( rr2_2 < pow( eps_r, 2 ) ) {
        rr2_2 = pow( eps_r, 2 );
        }
        const CCTK_REAL rr_2 = sqrt(rr2_2);
        // const CCTK_REAL rho2_2 = gamma2*x1_2*x1_2 + y1_2*y1_2;
        // const CCTK_REAL rho_2  = sqrt(rho2_2);
        CCTK_REAL rho2_2 = x1_2*x1_2 + y1_2*y1_2;
        if( rho2_2 < pow( eps_r, 2 ) ){
        rho2_2 = pow( eps_r, 2 );
        }
        const CCTK_REAL rho_2  = sqrt(rho2_2);
        

        const CCTK_REAL theta_2 = acos(z1_2/rr_2);

        const CCTK_REAL deltakerr2_2 = bh_mass*bh_mass - bh_spin2 ;
        const CCTK_REAL deltakerr  = sqrt(deltakerr2_2) ;

        const CCTK_REAL costh_2  = z1_2/rr_2 ;
        const CCTK_REAL costh2_2 = costh_2*costh_2 ;
        const CCTK_REAL sinth2_2 = 1. - costh2_2 ;
        const CCTK_REAL sinth_2  = sqrt(sinth2_2) ;

        // const CCTK_REAL R_x    = gamma*x1_2/rr_2 ;
        const CCTK_REAL R_x    = x1_2/rr_2 ;
        const CCTK_REAL R_y    = y1_2/rr_2 ;
        const CCTK_REAL R_z    = z1_2/rr_2 ;

        // const CCTK_REAL x_R    = gamma*x1_2/rr_2 ;
        const CCTK_REAL x_R    = x1_2/rr_2 ;
        const CCTK_REAL y_R    = y1_2/rr_2 ;
        const CCTK_REAL z_R    = z1_2/rr_2 ;

        const CCTK_REAL sinth2ph_x = -y1_2/rr2_2 ;
        // const CCTK_REAL sinth2ph_y =  gamma*x1_2/rr2_2 ;
        const CCTK_REAL sinth2ph_y =  x1_2/rr2_2 ;


        // const CCTK_REAL sinthth_x  = z1_2*gamma*x1_2/(rr_2*rr2_2) ;
        const CCTK_REAL sinthth_x  = z1_2*x1_2/(rr_2*rr2_2) ; 
        const CCTK_REAL sinthth_y  = z1_2*y1_2/(rr_2*rr2_2) ; 
        const CCTK_REAL sinthth_z  = -sinth2_2/rr_2 ; 

        // const CCTK_REAL sinthx_th  = gamma*x1_2 * costh_2 ;
        const CCTK_REAL sinthx_th  = x1_2 * costh_2 ;
        const CCTK_REAL sinthy_th  = y1_2 * costh_2 ;
        const CCTK_REAL sinthz_th  = -rr_2 * sinth2_2 ;


        const CCTK_REAL rBL    = rr_2 + bh_mass + 0.25*deltakerr2_2 / rr_2 ;   // Boyer-Lindquist coordinate r

        const CCTK_REAL RRrBL  = rr2_2 + rr_2*bh_mass + 0.25*deltakerr2_2 ;

        const CCTK_REAL rho2kerr   = rBL*rBL + bh_spin2 * costh2_2 ;
        const CCTK_REAL rhokerr    = sqrt(rho2kerr) ;

        const CCTK_REAL sigma  = (2.*bh_mass*rBL)/rho2kerr;
        const CCTK_REAL hh     = (1 + sigma) / (RRrBL*RRrBL + rr2_2*bh_spin*bh_spin * costh2_2) ;

        const CCTK_REAL psi4_2 = rho2kerr / rr2_2 ;
        const CCTK_REAL psi2_2 = sqrt(psi4_2) ;
        psi1_2 = sqrt(psi2_2) ;
        

        // non-axisymmetric perturbation.
        /* pert = 1. + AA * (x1_2*x1_2 - y1_2*y1_2)/(bh_mass*bh_mass) * exp( -2.*rr2_2/deltakerr2_2 ) ; */
        
        alpha0  = (rr_2 + 0.5*deltakerr)*(rr_2 - 0.5*deltakerr) / rr_2 * \
                 1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2_2)) ;
        const CCTK_REAL alpha02 = alpha0*alpha0 ;

       

        // add non-axisymmetric perturbation on conformal factor
        // NOTE: the perturbation is only taken into account for the 3-metric grid functions (not extrinsic curvature, lapse, ...)
        const CCTK_REAL argpert_cf = (rr_1 - R0pert_conf_fac)/Sigmapert_conf_fac;
        const CCTK_REAL pert_cf = 1. + Apert_conf_fac * (x1_1*x1_1*gamma2 - y1_1*y1_1)*mu*mu * exp( -0.5*argpert_cf*argpert_cf );

        const CCTK_REAL conf_fac = psi4_1 * pert_cf;

     

        CCTK_REAL gammaB[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            gammaB[a][b] = 0.0;
          }
        }
        gammaB[1][1] = psi4_2*(1+bh_spin2*hh*y1_2*y1_2);
        gammaB[1][2] = -psi4_2*bh_spin2*hh*y1_2*x1_2;
        gammaB[1][3] = 0;
        gammaB[2][1] = gammaB[1][2];
        gammaB[2][2] = psi4_2 * ( 1. + bh_spin2 * hh * x1_2*x1_2);
        gammaB[2][3] = 0;
        gammaB[3][1] = gammaB[1][3];
        gammaB[3][2] = gammaB[2][3];
        gammaB[3][3] = psi4_2; 

        CCTK_REAL gammaB_inv[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            gammaB_inv[a][b] = 0.0;
          }
        }
   
        gammaB_inv[1][1] = (1+bh_spin2*hh*x1_2*x1_2)/(psi4_2*(1+bh_spin2*hh*rho2_2));
        gammaB_inv[1][2] = bh_spin2*hh*x1_2*y1_2/(psi4_2*(1+bh_spin2*hh*rho2_2));
        gammaB_inv[1][3] = 0;
        gammaB_inv[2][1] = bh_spin2*hh*x1_2*y1_2/(psi4_2*(1+bh_spin2*hh*rho2_2));
        gammaB_inv[2][2] = (1+bh_spin2*hh*y1_2*y1_2)/(psi4_2*(1+bh_spin2*hh*rho2_2));
        gammaB_inv[2][3] = 0;
        gammaB_inv[3][1] = 0;
        gammaB_inv[3][2] = 0;
        gammaB_inv[3][3] = 1./psi4_2;

        

        // 3-metric
        gxx[ind] = psi4_2*(1+bh_spin2*hh*y1_2*y1_2) + Gb[1][1] - 1;
        gxy[ind] = -psi4_2*bh_spin2*hh*y1_2*x1_2 + Gb[1][2];
        gxz[ind] = 0 + Gb[1][3];
        gyy[ind] = psi4_2 * ( 1. + bh_spin2 * hh * x1_2*x1_2) + Gb[2][2] - 1;
        gyz[ind] = 0 + Gb[2][3];
        gzz[ind] = psi4_2 + Gb[3][3] - 1;


        check_nan_or_inf("gxx",gxx[ind]);
        check_nan_or_inf("gxy",gxy[ind]);
        check_nan_or_inf("gxz",gxz[ind]);
        check_nan_or_inf("gyy",gyy[ind]);
        check_nan_or_inf("gyz",gyz[ind]);
        check_nan_or_inf("gzz",gzz[ind]);


        const CCTK_REAL HF     = - bh_spin2*bh_spin * alpha0 * sigma/rhokerr * costh_2  ;  // we are dividing by sinth2_2
        const CCTK_REAL Athph  = HF / rr_2 ;                                        // we are dividing by sinth_2

        const CCTK_REAL aux    =  rho2kerr * (rBL*rBL - bh_spin2) + 2.*rBL*rBL * (rBL*rBL + bh_spin2);

        const CCTK_REAL HE     = bh_spin*bh_mass * aux / (rhokerr*rhokerr*rhokerr) * 
                 1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2_2)) ;

        const CCTK_REAL ARph   = HE / rr2_2 ;                                       // we are dividing by sinth2_2


        const CCTK_REAL Axx = 2.*ARph *  R_x * sinth2ph_x                     +  2.*Athph *  sinthth_x * sinth2ph_x ;
        const CCTK_REAL Axy =    ARph * (R_x * sinth2ph_y + R_y * sinth2ph_x) +     Athph * (sinthth_x * sinth2ph_y + sinthth_y * sinth2ph_x) ;
        const CCTK_REAL Axz =    ARph *                     R_z * sinth2ph_x  +     Athph *                           sinthth_z * sinth2ph_x  ; 
        const CCTK_REAL Ayy = 2.*ARph *  R_y * sinth2ph_y                     +  2.*Athph *  sinthth_y * sinth2ph_y ;
        const CCTK_REAL Ayz =    ARph *                     R_z * sinth2ph_y  +     Athph *                           sinthth_z * sinth2ph_y  ;

        CCTK_REAL dW_drho, dW_dz;
        const CCTK_REAL exp_auxi = exp(2. * F2_1[ind] - F0_1[ind]);

        if (rho_1 < 1e-8) {
          dW_drho = 0.;
          dW_dz   = 0.;
        }
        else {
          dW_drho = rho_1/rr_1 * dW_dr_1[ind]  +   z1_1/rr2_1 * dW_dth_1[ind];
          dW_dz   =  z1_1/rr_1 * dW_dr_1[ind]  -  rho_1/rr2_1 * dW_dth_1[ind];
        }

        CCTK_REAL gamma_final[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            gamma_final[a][b] = 0.0;
          }
        }
        gamma_final[1][1] = gxx[ind];
        gamma_final[1][2] = gxy[ind];
        gamma_final[1][3] = gxz[ind];
        gamma_final[2][1] = gxy[ind];
        gamma_final[2][2] = gyy[ind];
        gamma_final[2][3] = gyz[ind];
        gamma_final[3][1] = gxz[ind];
        gamma_final[3][2] = gyz[ind];
        gamma_final[3][3] = gzz[ind];


        CCTK_REAL K_B[4][4]; // extrinsic curvature
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            K_B[a][b] = 0.0;
          }
        } //K_0\mu might not be zero but irrelevant for what i want to compute

        K_B[1][1] =  Axx / psi2_2;
        K_B[1][2] =  Axy / psi2_2;
        K_B[1][3] =  Axz / psi2_2;
        K_B[2][1] =  Axy / psi2_2;
        K_B[2][2] =  Ayy / psi2_2;
        K_B[2][3] =  Ayz / psi2_2;
        K_B[3][1] =  Axz / psi2_2;
        K_B[3][2] =  Ayz / psi2_2;
        K_B[3][3] =  0;


        CCTK_REAL Kfinal[4][4]; // extrinsic curvature
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            Kfinal[a][b] = 0.0;
          }
        }

        for (int i = 1; i < 4; ++i) {
          for (int j = 1; j < 4; ++j) {
            CCTK_REAL sum1 = 0.0;
            CCTK_REAL sum2 = 0.0;
            for (int m = 1; m < 4; ++m) {
              for (int n = 1; n < 4; ++n) {
                sum1 += gamma_final[m][i] * (K_A[j][n] * gammaA_inv[n][m] + K_B[j][n] * gammaB_inv[n][m]);
                sum2 += gamma_final[m][j] * (K_A[i][n] * gammaA_inv[n][m] + K_B[i][n] * gammaB_inv[n][m]);
              }
            }
          Kfinal[i][j] = 0.5 * (sum1 + sum2);
          }
        }



        kxx[ind] = Kfinal[1][1];
        kxy[ind] = Kfinal[1][2];
        kxz[ind] = Kfinal[1][3];
        kyy[ind] = Kfinal[2][2];
        kyz[ind] = Kfinal[2][3];
        kzz[ind] = Kfinal[3][3];

        check_nan_or_inf("kxx",kxx[ind]);
        check_nan_or_inf("kxy",kxy[ind]);
        check_nan_or_inf("kxz",kxz[ind]);
        check_nan_or_inf("kyy",kyy[ind]);
        check_nan_or_inf("kyz",kyz[ind]);
        check_nan_or_inf("kzz",kzz[ind]);  

          
    //} // end if bh_spin_direction == "z1_1"


    // if (CCTK_EQUALS(bh_spin_direction, "y1_1")) { // rotation applied (x',y',z') = (x,z,-y)

    //     CCTK_REAL x1_2  = x[ind] - x0_2;
    //     CCTK_REAL y1_2  = y[ind] - y0_2;
    //     CCTK_REAL z1_2  = z[ind] - z0_2;

    //     // const CCTK_REAL bh_v2 = bs_v * bs_v;
    //     const CCTK_REAL bh_spin2 = bh_spin*bh_spin;
    //     // const CCTK_REAL gamma2 = 1. / (1. - bh_v2);
    //     // const CCTK_REAL gamma = sqrt(gamma2);
    //     CCTK_REAL rr2_2 = x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;
    //     if( rr2_2 < pow( eps_r, 2 ) ) {
    //     rr2_2 = pow( eps_r, 2 );
    //     }
    //     const CCTK_REAL rr_2  = sqrt(rr2_2);
    //     // const CCTK_REAL rho2_2 = gamma2*x1_2*x1_2 + y1_2*y1_2;
    //     // const CCTK_REAL rho_2  = sqrt(rho2_2);
    //     CCTK_REAL rho2_2 = x1_2*x1_2 + z1_2*z1_2;
    //     if( rho2_2 < pow( eps_r, 2 ) ){
    //     rho2_2 = pow( eps_r, 2 );
    //     }
    //     const CCTK_REAL rho_2  = sqrt(rho2_2);
        

    //     const CCTK_REAL theta_2 = acos(y1_2/rr_2);

    //     const CCTK_REAL deltakerr2_2 = bh_mass*bh_mass - bh_spin2 ;
    //     const CCTK_REAL deltakerr  = sqrt(deltakerr2_2) ;

    //     const CCTK_REAL costh_2  = y1_2/rr_2 ;
    //     const CCTK_REAL costh2_2 = costh_2*costh_2 ;
    //     const CCTK_REAL sinth2_2 = 1. - costh2_2 ;
    //     const CCTK_REAL sinth_2  = sqrt(sinth2_2) ;

    //     // const CCTK_REAL R_x    = gamma*x1_2/rr_2 ;
    //     const CCTK_REAL R_x    = x1_2/rr_2 ;
    //     const CCTK_REAL R_y    = y1_2/rr_2 ;
    //     const CCTK_REAL R_z    = z1_2/rr_2 ;

    //     // const CCTK_REAL x_R    = gamma*x1_2/rr_2 ;
    //     // const CCTK_REAL x_R    = x1_2/rr_2 ;
    //     // const CCTK_REAL y_R    = y1_2/rr_2 ;
    //     // const CCTK_REAL z_R    = z1_2/rr_2 ;

    //     const CCTK_REAL sinth2ph_x = z1_2/rr2_2 ;
    //     const CCTK_REAL sinth2ph_y = 0;
    //     const CCTK_REAL sinth2ph_z = -x1_2/rr2_2 ;


    //     // const CCTK_REAL sinthth_x  = z1_2*gamma*x1_2/(rr_2*rr2_2) ;
    //     const CCTK_REAL sinthth_x  = x1_2*y1_2/(rr_2*rr2_2) ; 
    //     const CCTK_REAL sinthth_y  = -rho2_2/(rr_2*rr2_2) ; 
    //     const CCTK_REAL sinthth_z  = y1_2*z1_2/(rr_2*rr2_2) ; 

  
    //     const CCTK_REAL rBL    = rr_2 + bh_mass + 0.25*deltakerr2_2 / rr_2 ;   // Boyer-Lindquist coordinate r

    //     const CCTK_REAL RRrBL  = rr2_2 + rr_2*bh_mass + 0.25*deltakerr2_2 ;

    //     const CCTK_REAL rho2kerr   = rBL*rBL + bh_spin2 * costh2_2 ;
    //     const CCTK_REAL rhokerr    = sqrt(rho2kerr) ;

    //     const CCTK_REAL sigma  = (2.*bh_mass*rBL)/rho2kerr;
    //     const CCTK_REAL hh     = (1 + sigma) / (RRrBL*RRrBL + rr2_2*bh_spin*bh_spin * costh2_2) ;

    //     const CCTK_REAL psi4_2 = rho2kerr / rr2_2 ;
    //     const CCTK_REAL psi2_2 = sqrt(psi4_2) ;
    //     psi1_2 = sqrt(psi2_2) ;
    //     // const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
    //     // const CCTK_REAL psi2_1 = sqrt(psi4_1);
    //     // psi1_1 = sqrt(psi2_1);

    //     // non-axisymmetric perturbation.
    //     /* pert = 1. + AA * (x1_2*x1_2 - y1_2*y1_2)/(bh_mass*bh_mass) * exp( -2.*rr2_2/deltakerr2_2 ) ; */
        
    //     alpha0  = (rr_2 + 0.5*deltakerr)*(rr_2 - 0.5*deltakerr) / rr_2 * \
    //              1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2_2)) ;
    //     const CCTK_REAL alpha02 = alpha0*alpha0 ;

       

    //     // add non-axisymmetric perturbation on conformal factor
    //     // NOTE: the perturbation is only taken into account for the 3-metric grid functions (not extrinsic curvature, lapse, ...)
    //     const CCTK_REAL argpert_cf = (rr_1 - R0pert_conf_fac)/Sigmapert_conf_fac;
    //     const CCTK_REAL pert_cf = 1. + Apert_conf_fac * (x1_1*x1_1*gamma2 - y1_1*y1_1)*mu*mu * exp( -0.5*argpert_cf*argpert_cf );

    //     const CCTK_REAL conf_fac = psi4_1 * pert_cf;

    //     // 3-metric
    //     // gxx[ind] = psi4_2*(1+bh_spin2*hh*y1_2*y1_2) + conf_fac * (1. + h_rho2_1 * sinph * sinph) - 1;
    //     // gxy[ind] = -psi4_2*bh_spin2*hh*y1_2*x1_2 - conf_fac * h_rho2_1 * sinph * cosph;
    //     // gxz[ind] = 0;
    //     // gyy[ind] = psi4_2 * ( 1. + bh_spin2 * hh * x1_2*x1_2) + conf_fac * (1. + h_rho2_1 * cosph * cosph) - 1;
    //     // gyz[ind] = 0;
    //     // gzz[ind] = psi4_2 + conf_fac - 1;


    //     gxx[ind] = psi4_2*(1+bh_spin2*hh*z1_2*z1_2) + Gb[1][1] - 1;
    //     gxy[ind] = 0 + Gb[1][2];
    //     gxz[ind] = - psi4_2*bh_spin2*hh*z1_2*x1_2 + Gb[1][3];
    //     gyy[ind] = psi4_2 + Gb[2][2] - 1;
    //     gyz[ind] = 0 + Gb[2][3];
    //     gzz[ind] = psi4_2 * ( 1. + bh_spin2 * hh * x1_2*x1_2) +  Gb[3][3] - 1;


    //     CCTK_REAL gammaB[4][4];
    //     for (int a = 0; a < 4; ++a) {
    //       for (int b = 0; b < 4; ++b) {
    //         gammaB[a][b] = 0.0;
    //       }
    //     }
    //     gammaB[1][1] = psi4_2*(1+bh_spin2*hh*z1_2*z1_2);
    //     gammaB[1][2] = 0;
    //     gammaB[1][3] = - psi4_2*bh_spin2*hh*z1_2*x1_2;
    //     gammaB[2][1] = gammaB[1][2];
    //     gammaB[2][2] = psi4_2;
    //     gammaB[2][3] = 0;
    //     gammaB[3][1] = gammaB[1][3];
    //     gammaB[3][2] = gammaB[2][3];
    //     gammaB[3][3] = psi4_2 * ( 1. + bh_spin2 * hh * x1_2*x1_2); 

    //     CCTK_REAL gammaB_inv[4][4];
    //     for (int a = 0; a < 4; ++a) {
    //       for (int b = 0; b < 4; ++b) {
    //         gammaB_inv[a][b] = 0.0;
    //       }
    //     }
    //     CCTK_REAL det_gammaB_inv =
    //         gammaB[1][1]*(gammaB[2][2]*gammaB[3][3] - gammaB[2][3]*gammaB[3][2])
    //       - gammaB[1][2]*(gammaB[2][1]*gammaB[3][3] - gammaB[2][3]*gammaB[3][1])
    //       + gammaB[1][3]*(gammaB[2][1]*gammaB[3][2] - gammaB[2][2]*gammaB[3][1]);

    //     gammaB_inv[1][1] =  (gammaB[2][2]*gammaB[3][3] - gammaB[2][3]*gammaB[3][2]) / det_gammaB_inv;
    //     gammaB_inv[1][2] = -(gammaB[1][2]*gammaB[3][3] - gammaB[1][3]*gammaB[3][2]) / det_gammaB_inv;
    //     gammaB_inv[1][3] =  (gammaB[1][2]*gammaB[2][3] - gammaB[1][3]*gammaB[2][2]) / det_gammaB_inv;
    //     gammaB_inv[2][1] = -(gammaB[2][1]*gammaB[3][3] - gammaB[2][3]*gammaB[3][1]) / det_gammaB_inv;
    //     gammaB_inv[2][2] =  (gammaB[1][1]*gammaB[3][3] - gammaB[1][3]*gammaB[3][1]) / det_gammaB_inv;
    //     gammaB_inv[2][3] = -(gammaB[1][1]*gammaB[2][3] - gammaB[1][3]*gammaB[2][1]) / det_gammaB_inv;
    //     gammaB_inv[3][1] =  (gammaB[2][1]*gammaB[3][2] - gammaB[2][2]*gammaB[3][1]) / det_gammaB_inv;
    //     gammaB_inv[3][2] = -(gammaB[1][1]*gammaB[3][2] - gammaB[1][2]*gammaB[3][1]) / det_gammaB_inv;
    //     gammaB_inv[3][3] =  (gammaB[1][1]*gammaB[2][2] - gammaB[1][2]*gammaB[2][1]) / det_gammaB_inv;

    //     /*
    //       d/drho = rho_1/r * d/dr  +    z1_1/r^2 * d/dth
    //       d/dz   =   z1_1/r * d/dr  -  rho_1/r^2 * d/dth

    //       Kxx = 0.5 * 2xy/rho_1        * exp(2F2-F0_1) * dW/drho   = 0.5 * rho_1 * sin(2phi) * exp(2F2-F0_1) * dW/drho
    //       Kyy = - Kxx
    //       Kzz = 0
    //       Kxy =-0.5 * (x1_1*gamma^2-y1_1^2)/rho_1  * exp(2F2-F0_1) * dW/drho   = 0.5 * rho_1 * cos(2phi) * exp(2F2-F0_1) * dW/drho
    //       Kxz = 0.5 * y1_1 * exp(2F2-F0_1) * dW/dz
    //       Kyz =-0.5 * x1_1*gamma * exp(2F2-F0_1) * dW/dz
    //     */

    //     /*
    //       Close to the axis and the origin, Kij = 0.
    //       The "coordinate" part of the expressions above behave like rho_1 (or r).
    //       Let's first consider a threshold of rho_1 < 1e-8. The sphere r < 1e-8 is included in this cylinder.
    //       In this case, we just set d/drho and d/dz = 0 as proxies.
    //     */


    //     check_nan_or_inf("gxx",gxx[ind]);
    //     check_nan_or_inf("gxy",gxy[ind]);
    //     check_nan_or_inf("gxz",gxz[ind]);
    //     check_nan_or_inf("gyy",gyy[ind]);
    //     check_nan_or_inf("gyz",gyz[ind]);
    //     check_nan_or_inf("gzz",gzz[ind]);


    //     const CCTK_REAL HF     = - bh_spin2*bh_spin * alpha0 * sigma/rhokerr * costh_2  ;  // we are dividing by sinth2_2
    //     const CCTK_REAL Athph  = HF / rr_2 ;                                        // we are dividing by sinth_2

    //     const CCTK_REAL aux    =  rho2kerr * (rBL*rBL - bh_spin2) + 2.*rBL*rBL * (rBL*rBL + bh_spin2);

    //     const CCTK_REAL HE     = bh_spin*bh_mass * aux / (rhokerr*rhokerr*rhokerr) * 
    //              1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2_2)) ;

    //     const CCTK_REAL ARph   = HE / rr2_2 ;                                       // we are dividing by sinth2_2


    //     const CCTK_REAL Axx = 2.*ARph *  R_x * sinth2ph_x                     +  2.*Athph *  sinthth_x * sinth2ph_x ;
    //     const CCTK_REAL Axy =    ARph * (R_x * sinth2ph_y + R_y * sinth2ph_x) +     Athph * (sinthth_x * sinth2ph_y + sinthth_y * sinth2ph_x) ;
    //     const CCTK_REAL Axz =    ARph * (R_x * sinth2ph_z + R_z * sinth2ph_x) +     Athph * (sinthth_x * sinth2ph_z + sinthth_z * sinth2ph_x) ; 
    //     const CCTK_REAL Ayy = 2.*ARph *  R_y * sinth2ph_y                     +  2.*Athph *  sinthth_y * sinth2ph_y ;
    //     const CCTK_REAL Ayz =    ARph * (R_y * sinth2ph_z + R_z * sinth2ph_y) +     Athph * (sinthth_y * sinth2ph_z + sinthth_z * sinth2ph_y) ;
    //     const CCTK_REAL Azz = 2.*ARph *  R_z * sinth2ph_z                     +  2.*Athph *  sinthth_z * sinth2ph_z ;

    //     CCTK_REAL dW_drho, dW_dz;
    //     const CCTK_REAL exp_auxi = exp(2. * F2_1[ind] - F0_1[ind]);

    //     if (rho_1 < 1e-8) {
    //       dW_drho = 0.;
    //       dW_dz   = 0.;
    //     }
    //     else {
    //       dW_drho = rho_1/rr_1 * dW_dr_1[ind]  +   z1_1/rr2_1 * dW_dth_1[ind];
    //       dW_dz   =  z1_1/rr_1 * dW_dr_1[ind]  -  rho_1/rr2_1 * dW_dth_1[ind];
    //     }

    //     CCTK_REAL gamma_final[4][4];
    //     for (int a = 0; a < 4; ++a) {
    //       for (int b = 0; b < 4; ++b) {
    //         gamma_final[a][b] = 0.0;
    //       }
    //     }
    //     gamma_final[1][1] = gxx[ind];
    //     gamma_final[1][2] = gxy[ind];
    //     gamma_final[1][3] = gxz[ind];
    //     gamma_final[2][1] = gxy[ind];
    //     gamma_final[2][2] = gyy[ind];
    //     gamma_final[2][3] = gyz[ind];
    //     gamma_final[3][1] = gxz[ind];
    //     gamma_final[3][2] = gyz[ind];
    //     gamma_final[3][3] = gzz[ind];


    //     CCTK_REAL gamma_final_inv[4][4];
    //     for (int a = 0; a < 4; ++a) {
    //       for (int b = 0; b < 4; ++b) {
    //         gamma_final_inv[a][b] = 0.0;
    //       }
    //     }
    //     CCTK_REAL det_gamma_final =
    //         gamma_final[1][1]*(gamma_final[2][2]*gamma_final[3][3] - gamma_final[2][3]*gamma_final[3][2])
    //       - gamma_final[1][2]*(gamma_final[2][1]*gamma_final[3][3] - gamma_final[2][3]*gamma_final[3][1])
    //       + gamma_final[1][3]*(gamma_final[2][1]*gamma_final[3][2] - gamma_final[2][2]*gamma_final[3][1]);

    //     gamma_final_inv[1][1] =  (gamma_final[2][2]*gamma_final[3][3] - gamma_final[2][3]*gamma_final[3][2]) / det_gamma_final;
    //     gamma_final_inv[1][2] = -(gamma_final[1][2]*gamma_final[3][3] - gamma_final[1][3]*gamma_final[3][2]) / det_gamma_final;
    //     gamma_final_inv[1][3] =  (gamma_final[1][2]*gamma_final[2][3] - gamma_final[1][3]*gamma_final[2][2]) / det_gamma_final;
    //     gamma_final_inv[2][1] = -(gamma_final[2][1]*gamma_final[3][3] - gamma_final[2][3]*gamma_final[3][1]) / det_gamma_final;
    //     gamma_final_inv[2][2] =  (gamma_final[1][1]*gamma_final[3][3] - gamma_final[1][3]*gamma_final[3][1]) / det_gamma_final;
    //     gamma_final_inv[2][3] = -(gamma_final[1][1]*gamma_final[2][3] - gamma_final[1][3]*gamma_final[2][1]) / det_gamma_final;
    //     gamma_final_inv[3][1] =  (gamma_final[2][1]*gamma_final[3][2] - gamma_final[2][2]*gamma_final[3][1]) / det_gamma_final;
    //     gamma_final_inv[3][2] = -(gamma_final[1][1]*gamma_final[3][2] - gamma_final[1][2]*gamma_final[3][1]) / det_gamma_final;
    //     gamma_final_inv[3][3] =  (gamma_final[1][1]*gamma_final[2][2] - gamma_final[1][2]*gamma_final[2][1]) / det_gamma_final;


    //     CCTK_REAL K_B[4][4]; // extrinsic curvature
    //     for (int a = 0; a < 4; ++a) {
    //       for (int b = 0; b < 4; ++b) {
    //         K_B[a][b] = 0.0;
    //       }
    //     } //K_0\mu might not be zero but irrelevant for what i want to compute

    //     K_B[1][1] =  Axx / psi2_2;
    //     K_B[1][2] =  Axy / psi2_2;
    //     K_B[1][3] =  Axz / psi2_2;
    //     K_B[2][1] =  Axy / psi2_2;
    //     K_B[2][2] =  Ayy / psi2_2;
    //     K_B[2][3] =  Ayz / psi2_2;
    //     K_B[3][1] =  Axz / psi2_2;
    //     K_B[3][2] =  Ayz / psi2_2;
    //     K_B[3][3] =  Azz / psi2_2;


    //     CCTK_REAL Kfinal[4][4]; // extrinsic curvature
    //     for (int a = 0; a < 4; ++a) {
    //       for (int b = 0; b < 4; ++b) {
    //         Kfinal[a][b] = 0.0;
    //       }
    //     }

    //     for (int i = 1; i < 4; ++i) {
    //       for (int j = 1; j < 4; ++j) {
    //         CCTK_REAL sum1 = 0.0;
    //         CCTK_REAL sum2 = 0.0;
    //         for (int m = 1; m < 4; ++m) {
    //           for (int n = 1; n < 4; ++n) {
    //             sum1 += gamma_final[m][i] * (K_A[j][n] * gammaA_inv[n][m] + K_B[j][n] * gammaB_inv[n][m]);
    //             sum2 += gamma_final[m][j] * (K_A[i][n] * gammaA_inv[n][m] + K_B[i][n] * gammaB_inv[n][m]);
    //           }
    //         }
    //       Kfinal[i][j] = 0.5 * (sum1 + sum2);
    //       }
    //     }



    //     kxx[ind] = Kfinal[1][1];
    //     kxy[ind] = Kfinal[1][2];
    //     kxz[ind] = Kfinal[1][3];
    //     kyy[ind] = Kfinal[2][2];
    //     kyz[ind] = Kfinal[2][3];
    //     kzz[ind] = Kfinal[3][3];

    //     check_nan_or_inf("kxx",kxx[ind]);
    //     check_nan_or_inf("kxy",kxy[ind]);
    //     check_nan_or_inf("kxz",kxz[ind]);
    //     check_nan_or_inf("kyy",kyy[ind]);
    //     check_nan_or_inf("kyz",kyz[ind]);
    //     check_nan_or_inf("kzz",kzz[ind]);  

          
    // } // end if bh_spin_direction == "y1_1"


        // lapse value (field initialization below)
        // No lapse regularization needed for the BS, the lapse is non-zero
        //pre boost since i want to compute unboosted quantities
        const CCTK_REAL alph = exp(F0_1[ind]);// + alpha0 - 1;

    


        // let's add a perturbation to the Proca field as well
        // NOTE: the perturbation is added directed to every instance of e^{i m \varphi}, hence its derivatives are not taken into account
        // TODO (?): Design perturbation more generically as ~ cos((m+1)\varphi)
        const CCTK_REAL argpert_Proca = (rr_1 - R0pert_Proca)/Sigmapert_Proca;
        const CCTK_REAL pert_Proca = 1. + Apert_Proca * (x1_1*x1_1*gamma2 - y1_1*y1_1)*mu*mu * exp( -0.5*argpert_Proca*argpert_Proca ); //ignorar por agora


        // ----- Proca fields -----

        // TODO: check what happens with divisions by rr_1 and sinth_1, can we work around them?

        // Real and imaginay part of the harmonic dependence: exp[i(m\varphi - \omega t)] rotation not implemented as of yet
        const CCTK_REAL harm_re = (coswt * cosmph + sinwt * sinmph) * pert_Proca;
        const CCTK_REAL harm_im = (coswt * sinmph - sinwt * cosmph) * pert_Proca;

        // No need to change the radial component, R and r coincide

        CCTK_REAL A1_unboosted[4]; //A_\mu real part
        CCTK_REAL A2_unboosted[4]; //A_\mu imag part

        // A_t
        A1_unboosted[0] = V_1[ind] * sinwt; 
        A2_unboosted[0] = V_1[ind] * coswt;

        // A_x
        A1_unboosted[1] = x1_1*gamma/rr_1 * H1r_1[ind] * harm_re + costh_1*cosph/rr_1 * H2_1[ind] * harm_re + sinph/rr_1 * H3_1[ind] * harm_im;
        A2_unboosted[1] = x1_1*gamma/rr_1 * H1r_1[ind] * harm_im + costh_1*cosph/rr_1 * H2_1[ind] * harm_im - sinph/rr_1 * H3_1[ind] * harm_re;
        
        // A_y
        A1_unboosted[2] = y1_1/rr_1 * H1r_1[ind] * harm_re + costh_1*sinph/rr_1 * H2_1[ind] * harm_re - cosph/rr_1 * H3_1[ind] * harm_im;
        A2_unboosted[2] = y1_1/rr_1 * H1r_1[ind] * harm_im + costh_1*sinph/rr_1 * H2_1[ind] * harm_im + cosph/rr_1 * H3_1[ind] * harm_re;
        
        // A_z
        A1_unboosted[3] = (z1_1/rr_1 * H1r_1[ind] - sinth_1/rr_1 * H2_1[ind]) * harm_re;
        A2_unboosted[3] = (z1_1/rr_1 * H1r_1[ind] - sinth_1/rr_1 * H2_1[ind]) * harm_im;

        const CCTK_REAL dH1r_dr_1 = dH1_dr_1[ind]/rr_1 - H1r_1[ind]/rr_1;

        // Build unboosted field-strength tensor F_{mu nu} (only 0i components from time/spatial derivatives)
        CCTK_REAL F1_unb[4][4], F2_unb[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            F1_unb[a][b] = 0.0;
            F2_unb[a][b] = 0.0;
          }
        }

        const CCTK_REAL dV_dx = dV_dr_1[ind]*R_x_1 + dV_dth_1[ind]*th_x_1;
        const CCTK_REAL dV_dy = dV_dr_1[ind]*R_y_1 + dV_dth_1[ind]*th_y_1;
        const CCTK_REAL dV_dz = dV_dr_1[ind]*R_z_1 + dV_dth_1[ind]*th_z_1;

        const CCTK_REAL dH1r_dx = dH1r_dr_1 * R_x_1 + dH1_dr_1[ind]*th_x_1;
        const CCTK_REAL dH1r_dy = dH1r_dr_1 * R_y_1 + dH1_dr_1[ind]*th_y_1;
        const CCTK_REAL dH1r_dz = dH1r_dr_1 * R_z_1 + dH1_dr_1[ind]*th_z_1;

        const CCTK_REAL dH2_dx = dH2_dr_1[ind] * R_x_1 + dH2_dr_1[ind]*th_x_1;
        const CCTK_REAL dH2_dy = dH2_dr_1[ind] * R_y_1 + dH2_dr_1[ind]*th_y_1;
        const CCTK_REAL dH2_dz = dH2_dr_1[ind] * R_z_1 + dH2_dr_1[ind]*th_z_1;

        const CCTK_REAL dH3_dx = dH3_dr_1[ind] * R_x_1 + dH3_dr_1[ind]*th_x_1;
        const CCTK_REAL dH3_dy = dH3_dr_1[ind] * R_y_1 + dH3_dr_1[ind]*th_y_1;
        const CCTK_REAL dH3_dz = dH3_dr_1[ind] * R_z_1 + dH3_dr_1[ind]*th_z_1;

        const CCTK_REAL dA1x_dt = omega_BS * ( x1_1*gamma/rr_1 * H1r_1[ind] * harm_im + costh_1*cosph/rr_1 * H2_1[ind] * harm_im - sinph/rr_1 * H3_1[ind] * harm_re );
        const CCTK_REAL dA1y_dt = omega_BS * ( y1_1/rr_1 * H1r_1[ind] * harm_im + costh_1*sinph/rr_1 * H2_1[ind] * harm_im + cosph/rr_1 * H3_1[ind] * harm_re );
        const CCTK_REAL dA1z_dt = omega_BS * ( (z1_1/rr_1 * H1r_1[ind] - sinth_1/rr_1 * H2_1[ind]) * harm_im );

        const CCTK_REAL dA2x_dt = -omega_BS * ( x1_1*gamma/rr_1 * H1r_1[ind] * harm_re + costh_1*cosph/rr_1 * H2_1[ind] * harm_re + sinph/rr_1 * H3_1[ind] * harm_im );
        const CCTK_REAL dA2y_dt = -omega_BS * ( y1_1/rr_1 * H1r_1[ind] * harm_re + costh_1*sinph/rr_1 * H2_1[ind] * harm_re - cosph/rr_1 * H3_1[ind] * harm_im );
        const CCTK_REAL dA2z_dt = -omega_BS * ( (z1_1/rr_1 * H1r_1[ind] - sinth_1/rr_1 * H2_1[ind]) * harm_re );
        

        const CCTK_REAL dA1t_dx = dV_dx * sinwt;
        const CCTK_REAL dA1t_dy = dV_dy * sinwt;
        const CCTK_REAL dA1t_dz = dV_dz * sinwt;

        const CCTK_REAL dA2t_dx = dV_dx * coswt;
        const CCTK_REAL dA2t_dy = dV_dy * coswt;
        const CCTK_REAL dA2t_dz = dV_dz * coswt;

        const CCTK_REAL dA1x_dx = (-(sinwt*(H3_1[ind]*(d_sinph_dx*rr_1 - \
                                  R_x_1*sinph) + rr_1*sinph*(dH3_dr_1[ind]*R_x_1 + \
                                  th_x_1*dH3_dth_1[ind]))) + \
                                  coswt*(H1r_1[ind]*(rr_1 - R_x_1*x1_1*gamma) + \
                                  (H2_1[ind]*(d_cosph_dx*rr_1 - cosph*R_x_1) + \
                                  cosph*rr_1*(dH2_dr_1[ind]*R_x_1 + \
                                  dH2_dth_1[ind]*th_x_1))*costh_1 + \
                                  rr_1*(dH1_dr_1[ind]*R_x_1*x1_1*gamma + dH1_dth_1[ind]*th_x_1*x1_1*gamma + \
                                  cosph*H2_1[ind]*(-sinth_1*th_x_1))))/rr2_1;
        const CCTK_REAL dA1x_dy = (-(sinwt*(H3_1[ind]*(d_sinph_dy*rr_1 - \
                                  R_y_1*sinph) + rr_1*sinph*(dH3_dr_1[ind]*R_y_1 + \
                                  th_y_1*dH3_dth_1[ind]))) + \
                                  coswt*((-H1r_1[ind] + \
                                  dH1_dr_1[ind]*rr_1)*R_y_1*x1_1*gamma + \
                                  dH1_dth_1[ind]*rr_1*th_y_1*x1_1*gamma + \
                                  (H2_1[ind]*(d_cosph_dy*rr_1 - cosph*R_y_1) + \
                                  cosph*rr_1*(dH2_dr_1[ind]*R_y_1 + \
                                  dH2_dth_1[ind]*th_y_1))*costh_1 + \
                                  cosph*H2_1[ind]*rr_1*(-sinth_1*th_y_1)))/rr2_1;
        const CCTK_REAL dA1x_dz = (R_z_1*((-H1r_1[ind] + \
                                  dH1_dr_1[ind]*rr_1)*x1_1*gamma*coswt + (H3_1[ind] - \
                                  dH3_dr_1[ind]*rr_1)*sinph*sinwt) + \
                                  rr_1*th_z_1*(dH1_dth_1[ind]*x1_1*gamma*coswt - \
                                  sinph*sinwt*dH3_dth_1[ind]) + \
                                  cosph*coswt*(rr_1*(dH2_dr_1[ind]*R_z_1 + \
                                  dH2_dth_1[ind]*th_z_1)*costh_1 + \
                                  H2_1[ind]*(-(R_z_1*costh_1) + \
                                  rr_1*(-sinth_1*th_z_1))))/rr2_1;

        const CCTK_REAL dA1y_dx = ((H3_1[ind]*(d_cosph_dx*rr_1 - cosph*R_x_1) + \
                                  cosph*rr_1*(dH3_dr_1[ind]*R_x_1 + dH3_dth_1[ind]*\
                                  th_x_1))*sinwt + \
                                  coswt*((-H1r_1[ind] + dH1_dr_1[ind]*rr_1\
                                  )*R_x_1*y1_1 + dH1_dth_1[ind]*rr_1*th_x_1*y1_1 + \
                                  (H2_1[ind]*(d_sinph_dx*rr_1 - R_x_1*sinph) + \
                                  rr_1*sinph*(dH2_dr_1[ind]*R_x_1 + dH2_dth_1[ind]*\
                                  th_x_1))*costh_1 + \
                                  H2_1[ind]*rr_1*sinph*(-sinth_1*th_x_1)))/rr2_1;
        const CCTK_REAL dA1y_dy = ((H3_1[ind]*(d_cosph_dy*rr_1 - cosph*R_y_1) + \
                                  cosph*rr_1*(dH3_dr_1[ind]*R_y_1 + dH3_dth_1[ind]*\
                                  th_y_1))*sinwt + \
                                  coswt*(H1r_1[ind]*(rr_1 - R_y_1*y1_1) + (\
                                  H2_1[ind]*(d_sinph_dy*rr_1 - R_y_1*sinph) + \
                                  rr_1*sinph*(dH2_dr_1[ind]*R_y_1 + dH2_dth_1[ind]*\
                                  th_y_1))*costh_1 + rr_1*(dH1_dr_1[ind]*R_y_1*y1_1 + \
                                  dH1_dth_1[ind]*th_y_1*y1_1 + \
                                  H2_1[ind]*sinph*(-sinth_1*th_y_1))))/rr2_1;
        const CCTK_REAL dA1y_dz = (rr_1*th_z_1*(coswt*(dH1_dth_1[ind]*y1_1 + \
                                  dH2_dth_1[ind]*sinph*costh_1) + \
                                  cosph*dH3_dth_1[ind]*sinwt) + \
                                  R_z_1*(coswt*(-(H1r_1[ind]*y1_1) + \
                                  rr_1*(dH1_dr_1[ind]*y1_1 + dH2_dr_1[ind]*sinph*costh_1)) + \
                                  cosph*(-H3_1[ind] + dH3_dr_1[ind]*rr_1)*sinwt) \
                                  + H2_1[ind]*sinph*coswt*(-(R_z_1*costh_1) + \
                                  rr_1*(-sinth_1*th_z_1)))/rr2_1;

        const CCTK_REAL dA1z_dx = (coswt*(-(H1r_1[ind]*R_x_1*z1_1) + \
                                  rr_1*(R_x_1*(dH1_dr_1[ind]*z1_1 - \
                                  dH2_dr_1[ind]*sinth_1) + th_x_1*(dH1_dth_1[ind]*z1_1 - \
                                  dH2_dth_1[ind]*sinth_1)) + \
                                  H2_1[ind]*(R_x_1*sinth_1 - \
                                  rr_1*(costh_1*th_x_1))))/rr2_1;
        const CCTK_REAL dA1z_dy = (coswt*(-(H1r_1[ind]*R_y_1*z1_1) + \
                                  rr_1*(R_y_1*(dH1_dr_1[ind]*z1_1 - \
                                  dH2_dr_1[ind]*sinth_1) + th_y_1*(dH1_dth_1[ind]*z1_1 - \
                                  dH2_dth_1[ind]*sinth_1)) + \
                                  H2_1[ind]*(R_y_1*sinth_1 - \
                                  rr_1*(costh_1*th_y_1))))/rr2_1;
        const CCTK_REAL dA1z_dz = (coswt*(H1r_1[ind]*(rr_1 - R_z_1*z1_1) + \
                                  rr_1*(R_z_1*(dH1_dr_1[ind]*z1_1 - \
                                  dH2_dr_1[ind]*sinth_1) + th_z_1*(dH1_dth_1[ind]*z1_1 - \
                                  dH2_dth_1[ind]*sinth_1)) + \
                                  H2_1[ind]*(R_z_1*sinth_1 - \
                                  rr_1*(costh_1*th_z_1))))/rr2_1;

        const CCTK_REAL dA2x_dx = (-((H3_1[ind]*(d_sinph_dx*rr_1 - R_x_1*sinph) + \
                                  rr_1*sinph*(dH3_dr_1[ind]*R_x_1 + \
                                  dH3_dth_1[ind]*th_x_1))*coswt) + \
                                  sinth_1*(H1r_1[ind]*(-rr_1 + R_x_1*x1_1*gamma) - \
                                  (H2_1[ind]*(d_cosph_dx*rr_1 - cosph*R_x_1) + \
                                  cosph*rr_1*(dH2_dr_1[ind]*R_x_1 + \
                                  dH2_dth_1[ind]*th_x_1))*costh_1 - \
                                  rr_1*(dH1_dr_1[ind]*R_x_1*x1_1*gamma + dH1_dth_1[ind]*th_x_1*x1_1*gamma + \
                                  cosph*H2_1[ind]*(-sinth_1*th_x_1))))/rr2_1;                  
        const CCTK_REAL dA2x_dy = (-((H3_1[ind]*(d_sinph_dy*rr_1 - R_y_1*sinph) + \
                                  rr_1*sinph*(dH3_dr_1[ind]*R_y_1 + dH3_dth_1[ind]*\
                                  th_y_1))*coswt) + \
                                  sinth_1*(H1r_1[ind]*R_y_1*x1_1*gamma - \
                                  (H2_1[ind]*(d_cosph_dy*rr_1 - cosph*R_y_1) + \
                                  cosph*rr_1*(dH2_dr_1[ind]*R_y_1 + dH2_dth_1[ind]*\
                                  th_y_1))*costh_1 - rr_1*(dH1_dr_1[ind]*R_y_1*x1_1*gamma + \
                                  dH1_dth_1[ind]*th_y_1*x1_1*gamma + \
                                  cosph*H2_1[ind]*(-sinth_1*th_y_1))))/rr2_1;
        const CCTK_REAL dA2x_dz = (-(rr_1*th_z_1*(dH3_dth_1[ind]*sinph*coswt + dH1_dth_1[ind]*x1_1*gamma*sinth_1)) + \
                                  R_z_1*((H3_1[ind] - \
                                  dH3_dr_1[ind]*rr_1)*sinph*coswt + \
                                  (H1r_1[ind] - dH1_dr_1[ind]*rr_1)*x1_1*gamma*sinwt) - cosph*sinth_1*(rr_1*(dH2_dr_1[ind]\
                                  *R_z_1 + dH2_dth_1[ind]*th_z_1)*costh_1 + \
                                  H2_1[ind]*(-(R_z_1*costh_1) + \
                                  rr_1*(-sinth_1*th_z_1))))/rr2_1;
        
        const CCTK_REAL dA2y_dx = ((H3_1[ind]*(d_cosph_dx*rr_1 - cosph*R_x_1) + \
                                  cosph*rr_1*(dH3_dr_1[ind]*R_x_1 + \
                                  dH3_dth_1[ind]*th_x_1))*coswt + \
                                  sinth_1*(H1r_1[ind]*R_x_1*y1_1 - \
                                  (H2_1[ind]*(d_sinph_dx*rr_1 - R_x_1*sinph) + \
                                  rr_1*sinph*(dH2_dr_1[ind]*R_x_1 + \
                                  dH2_dth_1[ind]*th_x_1))*costh_1 - \
                                  rr_1*(dH1_dr_1[ind]*R_x_1*y1_1 + dH1_dth_1[ind]*th_x_1*y1_1 + \
                                  H2_1[ind]*sinph*(-sinth_1*th_x_1))))/rr2_1;
        const CCTK_REAL dA2y_dy = ((H3_1[ind]*(d_cosph_dy*rr_1 - cosph*R_y_1) + \
                                  cosph*rr_1*(dH3_dr_1[ind]*R_y_1 + dH3_dth_1[ind]*\
                                  th_y_1))*coswt + \
                                  sinth_1*(H1r_1[ind]*(-rr_1 + R_y_1*y1_1) - \
                                  (H2_1[ind]*(d_sinph_dy*rr_1 - R_y_1*sinph) + \
                                  rr_1*sinph*(dH2_dr_1[ind]*R_y_1 + dH2_dth_1[ind]*\
                                  th_y_1))*costh_1 - rr_1*(dH1_dr_1[ind]*R_y_1*y1_1 + \
                                  dH1_dth_1[ind]*th_y_1*y1_1 + \
                                  H2_1[ind]*sinph*(-sinth_1*th_y_1))))/rr2_1;
        const CCTK_REAL dA2y_dz = (-(rr_1*th_z_1*(-(cosph*dH3_dth_1[ind]*coswt) + (dH1_dth_1[ind]*y1_1 + \
                                  dH2_dth_1[ind]*sinph*costh_1)*sinth_1)) \
                                  + R_z_1*(cosph*(-H3_1[ind] + \
                                  dH3_dr_1[ind]*rr_1)*coswt + \
                                  (H1r_1[ind]*y1_1 - rr_1*(dH1_dr_1[ind]*y1_1 + \
                                  dH2_dr_1[ind]*sinph*costh_1))*sinth_1) \
                                  + H2_1[ind]*sinph*sinth_1*(R_z_1*costh_1 - rr_1*(-sinth_1*th_z_1)))/rr2_1;
                                                              
        const CCTK_REAL dA2z_dx = (sinth_1*(H1r_1[ind]*R_x_1*z1_1 + \
                                  rr_1*(R_x_1*(-(dH1_dr_1[ind]*z1_1) + \
                                  dH2_dr_1[ind]*sinth_1) + th_x_1*(-(dH1_dth_1[ind]*z1_1) \
                                  + dH2_dth_1[ind]*sinth_1)) + \
                                  H2_1[ind]*(-(R_x_1*sinth_1) + \
                                  rr_1*(costh_1*th_x_1))))/rr2_1;
        const CCTK_REAL dA2z_dy = (sinth_1*(H1r_1[ind]*R_y_1*z1_1 + \
                                  rr_1*(R_y_1*(-(dH1_dr_1[ind]*z1_1) + \
                                  dH2_dr_1[ind]*sinth_1) + th_y_1*(-(dH1_dth_1[ind]*z1_1) \
                                  + dH2_dth_1[ind]*sinth_1)) + \
                                  H2_1[ind]*(-(R_y_1*sinth_1) + \
                                  rr_1*(costh_1*th_y_1))))/rr2_1;
        const CCTK_REAL dA2z_dz = (sinth_1*(H1r_1[ind]*(-rr_1 + R_z_1*z1_1) \
                                  + rr_1*(R_z_1*(-(dH1_dr_1[ind]*z1_1) + \
                                  dH2_dr_1[ind]*sinth_1) + th_z_1*(-(dH1_dth_1[ind]*z1_1) \
                                  + dH2_dth_1[ind]*sinth_1)) + \
                                  H2_1[ind]*(-(R_z_1*sinth_1) + \
                                  rr_1*(costh_1*th_z_1))))/rr2_1;



        // // Spatial derivatives of A_0 = V_1 * {sinwt, coswt}
        // const CCTK_REAL dV_dx = dV_dr_1[ind]*R_x_1 + dV_dth_1[ind]*th_x_1;
        // const CCTK_REAL dV_dy = dV_dr_1[ind]*R_y_1 + dV_dth_1[ind]*th_y_1;
        // const CCTK_REAL dV_dz = dV_dr_1[ind]*R_z_1 + dV_dth_1[ind]*th_z_1;

        // F_{0i} = d_t A_i - d_i A_0  (covariant indices)
        F1_unb[0][1] = dA1x_dt - dA1t_dx; F1_unb[1][0] = -F1_unb[0][1];
        F1_unb[0][2] = dA1y_dt - dA1t_dy; F1_unb[2][0] = -F1_unb[0][2];
        F1_unb[0][3] = dA1z_dt - dA1t_dz; F1_unb[3][0] = -F1_unb[0][3];
        F1_unb[1][2] = dA1y_dx - dA1x_dy; F1_unb[2][1] = -F1_unb[1][2];
        F1_unb[1][3] = dA1z_dx - dA1x_dz; F1_unb[3][1] = -F1_unb[1][3];
        F1_unb[2][3] = dA1z_dy - dA1y_dz; F1_unb[3][2] = -F1_unb[2][3];


        F2_unb[0][1] = dA2x_dt - dA2t_dx; F2_unb[1][0] = -F2_unb[0][1];
        F2_unb[0][2] = dA2y_dt - dA2t_dy; F2_unb[2][0] = -F2_unb[0][2];
        F2_unb[0][3] = dA2z_dt - dA2t_dz; F2_unb[3][0] = -F2_unb[0][3];
        F2_unb[1][2] = dA2y_dx - dA2x_dy; F2_unb[2][1] = -F2_unb[1][2];
        F2_unb[1][3] = dA2z_dx - dA2x_dz; F2_unb[3][1] = -F2_unb[1][3];
        F2_unb[2][3] = dA2z_dy - dA2y_dz; F2_unb[3][2] = -F2_unb[2][3];


        CCTK_REAL A1_boosted[4]; //A_\mu real part
        CCTK_REAL A2_boosted[4]; //A_\mu imag part
        // Boosted components
        for (int a = 0; a < 4; ++a) {
          A1_boosted[a] = 0.0;
          A2_boosted[a] = 0.0;
          for (int mu = 0; mu < 4; ++mu) {
            A1_boosted[a] += invLambda[mu][a] * A1_unboosted[mu];
            A2_boosted[a] += invLambda[mu][a] * A2_unboosted[mu];
          }
        }

        /* store spatial components */
        A1x[ind] = A1_boosted[1];
        A1y[ind] = A1_boosted[2];
        A1z[ind] = A1_boosted[3];

        A2x[ind] = A2_boosted[1];
        A2y[ind] = A2_boosted[2];
        A2z[ind] = A2_boosted[3];
        

        /*
          A_\phi = -n^\mu A_\mu = - (A_t + W_1*A_ph)/alpha
                = -i * e^{i (m ph_1 - w t)} * (V_1 + W_1 H3_1 sinth_1) / alpha
        */
        // Aphi1[ind] = (V_1[ind] + W_1[ind] * sinth_1 * H3_1[ind]) / alph * harm_im;
        // Aphi2[ind] =-(V_1[ind] + W_1[ind] * sinth_1 * H3_1[ind]) / alph * harm_re;

        Aphi1[ind] = 1 / alpha1 * ( - A1_boosted[0] + betaup1[1] * A1_boosted[1] + betaup1[2] * A1_boosted[2] + betaup1[3] * A1_boosted[3]);
        Aphi2[ind] = 1 / alpha1 * ( - A2_boosted[0] + betaup1[1] * A2_boosted[1] + betaup1[2] * A2_boosted[2] + betaup1[3] * A2_boosted[3]);


        //------ Boosted field-strength tensor F_{mu nu} ------
        CCTK_REAL F1_boosted[4][4], F2_boosted[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            F1_boosted[a][b] = 0.0;
            F2_boosted[a][b] = 0.0;
            for (int mu = 0; mu < 4; ++mu) {
              for (int nu = 0; nu < 4; ++nu) {
                F1_boosted[a][b] += invLambda[mu][a] * invLambda[nu][b] * F1_unb[mu][nu];
                F2_boosted[a][b] += invLambda[mu][a] * invLambda[nu][b] * F2_unb[mu][nu];
              }
            }
          }
        }


        // ----- Electric fields -----
        
        // First, E_i in (R, th, ph_1) coordinates
        CCTK_REAL E1d_r, E2d_r, E1d_th, E2d_th, E1d_ph_o_sinth, E2d_ph_o_sinth;

        // E_r
        /*
          E_r = i * e^{i(m phi - w t)} / alpha * [- (m*W_1 - w) H1r_1 + dV/dr + W_1 sinth_1 dH3/dr]
        */

        E1d_r = -(- (mm * W_1[ind] - omega_BS) * H1r_1[ind] + dV_dr_1[ind] + W_1[ind] * sinth_1 * dH3_dr_1[ind]) / alph * harm_im;
        E2d_r =  (- (mm * W_1[ind] - omega_BS) * H1r_1[ind] + dV_dr_1[ind] + W_1[ind] * sinth_1 * dH3_dr_1[ind]) / alph * harm_re;

        // E_th
        /*
          E_th = i * e^{i(m phi - w t)} / alpha * [- (m*W_1 - w) * H2_1 + dV/dth + W_1 * d(H3_1 * sinth_1)/dth]
        */
        E1d_th = -(- (mm * W_1[ind] - omega_BS) * H2_1[ind] + dV_dth_1[ind] + W_1[ind] * (sinth_1 * dH3_dth_1[ind] + costh_1 * H3_1[ind])) / alph * harm_im;
        E2d_th =  (- (mm * W_1[ind] - omega_BS) * H2_1[ind] + dV_dth_1[ind] + W_1[ind] * (sinth_1 * dH3_dth_1[ind] + costh_1 * H3_1[ind])) / alph * harm_re;

        // E_ph / sinth_1
        /*
          E_ph = - (m * V_1 + w * H3_1 * sinth_1) / alpha * e^{i(m phi - w t)}

          We include here the division by sinth_1 which arises when computing E^ph_1.
          
          on the axis sinth_1=0, we need to regularize the division by sinth_1
          We have V_1(theta=0,pi) = 0, so with l'Hôpital's rule
          V_1 / sinth_1 ~ \pm dV/dth    for theta = 0, pi resp.
        */
        if (fabs(sinth_1) < 1e-8) {
          // TODO: see how to deal with this if now we allow rr_1==0
          const CCTK_INT zsign = (costh_1>=0) ? 1 : -1; // costh_1==0 shouldn't happen on the axis for a grid point, this would mean rr_1==0 too...

          E1d_ph_o_sinth = - (mm * zsign * dV_dth_1[ind] + omega_BS * H3_1[ind]) / alph * harm_re;
          E2d_ph_o_sinth = - (mm * zsign * dV_dth_1[ind] + omega_BS * H3_1[ind]) / alph * harm_im;
        
        } else {
          E1d_ph_o_sinth = - (mm * V_1[ind] / sinth_1 + omega_BS * H3_1[ind]) / alph * harm_re;
          E2d_ph_o_sinth = - (mm * V_1[ind] / sinth_1 + omega_BS * H3_1[ind]) / alph * harm_im;
        }


        // E^i components
        // Spherical auxiliaries
        
        // E^r/r = e^{-2*F1_1} * E_r / r
        const CCTK_REAL E1u_r_o_r = exp(-2*F1_1[ind]) * E1d_r / rr_1;
        const CCTK_REAL E2u_r_o_r = exp(-2*F1_1[ind]) * E2d_r / rr_1;

        // E^th = e^{-2*F1_1} / r^2 * E_th
        const CCTK_REAL E1u_th = exp(-2*F1_1[ind]) / rr2_1 * E1d_th;
        const CCTK_REAL E2u_th = exp(-2*F1_1[ind]) / rr2_1 * E2d_th;

        // E^ph_1 = e^{-2*F2_1} / r^2 / sinth_1^2 * E_ph
        // We compute r * sinth_1 * E^ph_1. The other division by sinth_1 is managed with E_ph above 
        const CCTK_REAL rsinthE1u_ph = exp(-2*F2_1[ind]) / rr_1 * E1d_ph_o_sinth;
        const CCTK_REAL rsinthE2u_ph = exp(-2*F2_1[ind]) / rr_1 * E2d_ph_o_sinth;


        // Finally Cartesian components
        CCTK_REAL E1_unboosted[4]; //E^\mu real part
        CCTK_REAL E2_unboosted[4]; //E^\mu imag part

        // E^t
        E1_unboosted[0] = 0.;
        E2_unboosted[0] = 0.;
        // E^x1_1*gamma
        E1_unboosted[1] = x1_1 * gamma * E1u_r_o_r + z1_1 * cosph * E1u_th - sinph * rsinthE1u_ph;
        E2_unboosted[1] = x1_1 * gamma * E2u_r_o_r + z1_1 * cosph * E2u_th - sinph * rsinthE2u_ph;
        // E^y1_1
        E1_unboosted[2] = y1_1 * E1u_r_o_r + z1_1 * sinph * E1u_th + cosph * rsinthE1u_ph;
        E2_unboosted[2] = y1_1 * E2u_r_o_r + z1_1 * sinph * E2u_th + cosph * rsinthE2u_ph;
        // E^z1_1
        E1_unboosted[3] = z1_1 * E1u_r_o_r - rho_1 * E1u_th;
        E2_unboosted[3] = z1_1 * E2u_r_o_r - rho_1 * E2u_th;


        //Boosted components

        //cannot be boosted since its foliation dependent
        //furthermore it was also wrong for a generic 4-vector boost.
        //must compute the magnetic field?
        CCTK_REAL E1_boosted[4]; //E_\mu real part
        CCTK_REAL E2_boosted[4]; //E_\mu imag part
        for (int a = 0; a < 4; ++a) {
          E1_boosted[a] = 0.0;
          E2_boosted[a] = 0.0;
          // for (int mu = 0; mu < 4; ++mu) {
          //   E1_boosted[a] += Lambda[a][mu] * E1_unboosted[mu];
          //   E2_boosted[a] += Lambda[a][mu] * E2_unboosted[mu];
          // }
        }
        E1_boosted[1] = 1 / alpha1 * (F1_boosted[1][0] - betaup1[1] * F1_boosted[1][1] - betaup1[2] * F1_boosted[1][2] - betaup1[3] * F1_boosted[1][3]);
        E1_boosted[2] = 1 / alpha1 * (F1_boosted[2][0] - betaup1[1] * F1_boosted[2][1] - betaup1[2] * F1_boosted[2][2] - betaup1[3] * F1_boosted[2][3]);
        E1_boosted[3] = 1 / alpha1 * (F1_boosted[3][0] - betaup1[1] * F1_boosted[3][1] - betaup1[2] * F1_boosted[3][2] - betaup1[3] * F1_boosted[3][3]);

        E2_boosted[1] = 1 / alpha1 * (F2_boosted[1][0] - betaup1[1] * F2_boosted[1][1] - betaup1[2] * F2_boosted[1][2] - betaup1[3] * F2_boosted[1][3]);
        E2_boosted[2] = 1 / alpha1 * (F2_boosted[2][0] - betaup1[1] * F2_boosted[2][1] - betaup1[2] * F2_boosted[2][2] - betaup1[3] * F2_boosted[2][3]);
        E2_boosted[3] = 1 / alpha1 * (F2_boosted[3][0] - betaup1[1] * F2_boosted[3][1] - betaup1[2] * F2_boosted[3][2] - betaup1[3] * F2_boosted[3][3]);

        /* store spatial components */ 
        E1x[ind] = E1_boosted[1];
        E1y[ind] = E1_boosted[2];
        E1z[ind] = E1_boosted[3];

        E2x[ind] = E2_boosted[1];
        E2y[ind] = E2_boosted[2];
        E2z[ind] = E2_boosted[3];



        // zero-initialize constraint damping variable Zeta
        Zeta1[ind] = 0;
        Zeta2[ind] = 0;


        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1_1 + psi1_2 - 1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "ProcaBS")) {
          alp[ind] = alph;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // shift
        if (CCTK_EQUALS(initial_shift, "ProcaBS")) {
          betax[ind] =  W_1[ind] * y1_1;
          betay[ind] = -W_1[ind] * x1_1*gamma;
          betaz[ind] =  0.;
        }


      } /* for i */
    }   /* for j */
  }     /* for k */

  free(F1_1); free(F2_1); free(F0_1); free(W_1);
  free(H1r_1); free(H2_1); free(H3_1); free(V_1);
  free(dW_dr_1); free(dW_dth_1);
  free(dH3_dr_1); free(dH3_dth_1);
  free(dV_dr_1); free(dV_dth_1);

  return;
}
