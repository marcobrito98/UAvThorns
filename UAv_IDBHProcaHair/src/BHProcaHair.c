
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


void UAv_IDBHProcaHair(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

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

  CCTK_REAL *F1_extd, *F2_extd, *F0_extd, *Wbar_extd, *H2_extd, *H3_extd, *V_extd;
  F1_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  F2_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  F0_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  Wbar_extd  = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  H2_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  H3_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  V_extd     = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // We'll use A_r rather than the input A_x
  // Notation here: A ~ H1 dx + ... = H1r dr + ...    (input file has the first one: H1)
  CCTK_REAL *H1r_extd;
  H1r_extd          = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // now we need to take the derivatives of the Wbar function and store their values

  CCTK_REAL *dWbar_dr_extd, *dWbar_dth_extd, *d2Wbar_drth_extd;
  dWbar_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dWbar_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  d2Wbar_drth_extd = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  
  // Same for H3 and V
  
  CCTK_REAL *dH3_dr_extd, *dH3_dth_extd, *d2H3_drth_extd;
  dH3_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dH3_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  d2H3_drth_extd = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  
  CCTK_REAL *dV_dr_extd, *dV_dth_extd, *d2V_drth_extd;
  dV_dr_extd    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dV_dth_extd   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  d2V_drth_extd = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));


  // // Some auxi file for debug
  // FILE* debugfile = fopen ("testdebug.txt", "w");
  // if (debugfile == NULL) {
  //   CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
  //   "Unable to open file %s\n", "testdebug.txt");
  // } else {
  //   CCTK_VInfo(CCTK_THORNSTRING, "Write test file %s", "testdebug.txt");
  // }

  const CCTK_REAL oodX       = 1. / dX;
  // const CCTK_REAL oodXsq     = oodX * oodX;
  const CCTK_REAL oodX12     = 1. / (12. * dX);
  const CCTK_REAL oodth12    = 1. / (12. * dtheta);
  const CCTK_REAL oodXsq12   = oodX * oodX12;
  const CCTK_REAL oodXdth4   = 1. / (4.  * dX * dtheta);
  const CCTK_REAL oodXdth144 = 1. / (144. * dX * dtheta);
  // const CCTK_REAL oodXsqdth2 = 1. / (2.  * dX * dX * dtheta);
  const CCTK_REAL oodXsqdth144 = oodX * oodXdth144;


  // First loop on z>=0 half-space (i.e. input values of 0 <= theta <= pi/2)

  for (int jj = 0; jj < Ntheta; jj++) {
    for (int i = 0; i < NX; i++) {

      CCTK_INT j, jm1, jm2, jp1, jp2, jp3, jp4;
      /* let's use the fact that the solution is axi-symmetric (and that
         theta[0] = 0) for the boundary points in j */
      if (jj == 0) {
        j   = jj;
        jp1 = jj+1;
        jp2 = jj+2;
        jp3 = jj+3;
        jp4 = jj+4;
        jm1 = jj+1;
        jm2 = jj+2;
      } else if (jj == 1) {
        j   = jj;
        jp1 = jj+1;
        jp2 = jj+2;
        jp3 = jj+3;
        jp4 = jj+4;
        jm1 = jj-1;
        jm2 = jj;
      } else if (jj == Ntheta - 4) { // Shouldn't be needed here, but just in case
        j   = jj;
        jm1 = jj-1;
        jm2 = jj-2;
        jp1 = jj+1;
        jp2 = jj+2;
        jp3 = jj+3;
        jp4 = jj+2; 
      } else if (jj == Ntheta - 3) { // Shouldn't be needed here, but just in case
        j   = jj;
        jm1 = jj-1;
        jm2 = jj-2;
        jp1 = jj+1;
        jp2 = jj+2;
        jp3 = jj+1;
        jp4 = jj; 
      } else if (jj == Ntheta - 2) {
        j   = jj;
        jm1 = jj-1;
        jm2 = jj-2;
        jp1 = jj+1;
        jp2 = jj;
        jp3 = jj-1; // Shouldn't be needed here, but just in case
        jp4 = jj-2; // Shouldn't be needed here, but just in case
      } else if (jj == Ntheta - 1) {
        j   = jj;
        jm1 = jj-1;
        jm2 = jj-2;
        jp1 = jj-1;
        jp2 = jj-2;
        jp3 = jj-3; // Shouldn't be needed here, but just in case
        jp4 = jj-4; // Shouldn't be needed here, but just in case
      } else {
        j   = jj;
        jp1 = jj+1;
        jp2 = jj+2;
        jp3 = jj+3; // Shouldn't be needed here, but just in case
        jp4 = jj+4; // Shouldn't be needed here, but just in case
        jm1 = jj-1;
        jm2 = jj-2;
      }

      const CCTK_INT ind    = i + j*NX;

      const CCTK_INT indim1 = i-1 + j*NX;
      const CCTK_INT indip1 = i+1 + j*NX;
      const CCTK_INT indim2 = i-2 + j*NX;
      const CCTK_INT indip2 = i+2 + j*NX;
      const CCTK_INT indip3 = i+3 + j*NX;
      const CCTK_INT indip4 = i+4 + j*NX;
      const CCTK_INT indip5 = i+5 + j*NX;

      const CCTK_INT indim2jm1 = i-2 + jm1*NX;
      const CCTK_INT indim2jm2 = i-2 + jm2*NX;
      const CCTK_INT indim2jp1 = i-2 + jp1*NX;
      const CCTK_INT indim2jp2 = i-2 + jp2*NX;
      const CCTK_INT indim2jp3 = i-2 + jp3*NX;
      const CCTK_INT indim2jp4 = i-2 + jp4*NX;

      const CCTK_INT indim1jm1 = i-1 + jm1*NX;
      const CCTK_INT indim1jm2 = i-1 + jm2*NX;
      const CCTK_INT indim1jp1 = i-1 + jp1*NX;
      const CCTK_INT indim1jp2 = i-1 + jp2*NX;
      const CCTK_INT indim1jp3 = i-1 + jp3*NX;
      const CCTK_INT indim1jp4 = i-1 + jp4*NX;

      const CCTK_INT indjm1 = i + jm1*NX;
      const CCTK_INT indjm2 = i + jm2*NX;
      const CCTK_INT indjp1 = i + jp1*NX;
      const CCTK_INT indjp2 = i + jp2*NX;
      const CCTK_INT indjp3 = i + jp3*NX;
      const CCTK_INT indjp4 = i + jp4*NX;

      const CCTK_INT indip1jm1 = i+1 + jm1*NX;
      const CCTK_INT indip1jm2 = i+1 + jm2*NX;
      const CCTK_INT indip1jp1 = i+1 + jp1*NX;
      const CCTK_INT indip1jp2 = i+1 + jp2*NX;
      const CCTK_INT indip1jp3 = i+1 + jp3*NX;
      const CCTK_INT indip1jp4 = i+1 + jp4*NX;

      const CCTK_INT indip2jm1 = i+2 + jm1*NX;
      const CCTK_INT indip2jm2 = i+2 + jm2*NX;
      const CCTK_INT indip2jp1 = i+2 + jp1*NX;
      const CCTK_INT indip2jp2 = i+2 + jp2*NX;
      const CCTK_INT indip2jp3 = i+2 + jp3*NX;
      const CCTK_INT indip2jp4 = i+2 + jp4*NX;
      
      const CCTK_INT indip3jm1 = i+3 + jm1*NX;
      const CCTK_INT indip3jm2 = i+3 + jm2*NX;
      const CCTK_INT indip3jp1 = i+3 + jp1*NX;
      const CCTK_INT indip3jp2 = i+3 + jp2*NX;
      const CCTK_INT indip3jp3 = i+3 + jp3*NX;
      const CCTK_INT indip3jp4 = i+3 + jp4*NX;
      
      const CCTK_INT indip4jm1 = i+4 + jm1*NX;
      const CCTK_INT indip4jm2 = i+4 + jm2*NX;
      const CCTK_INT indip4jp1 = i+4 + jp1*NX;
      const CCTK_INT indip4jp2 = i+4 + jp2*NX;
      const CCTK_INT indip4jp3 = i+4 + jp3*NX;
      const CCTK_INT indip4jp4 = i+4 + jp4*NX;

      const CCTK_INT indip5jm1 = i+5 + jm1*NX;
      const CCTK_INT indip5jm2 = i+5 + jm2*NX;
      const CCTK_INT indip5jp1 = i+5 + jp1*NX;
      const CCTK_INT indip5jp2 = i+5 + jp2*NX;
      const CCTK_INT indip5jp3 = i+5 + jp3*NX;
      const CCTK_INT indip5jp4 = i+5 + jp4*NX;
      

      // Just copy input values of ansatz functions
      F1_extd[ind]   = F1_in[ind];
      F2_extd[ind]   = F2_in[ind];
      F0_extd[ind]   = F0_in[ind];
      Wbar_extd[ind] = Wbar_in[ind];
      H2_extd[ind]   = H2_in[ind];
      H3_extd[ind]   = H3_in[ind];
      V_extd[ind]    = V_in[ind];


      const CCTK_REAL lX = X[i];
      /* const CCTK_REAL lth = theta[j]; */
      /* printf("X[%3d] = %lf\n", i, lX); */


      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL Wbar_th = (-Wbar_in[indjp2] + 8 * Wbar_in[indjp1] - 8 * Wbar_in[indjm1] + Wbar_in[indjm2]) *
        oodth12;

      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL H3_th = (-H3_in[indjp2] + 8 * H3_in[indjp1] - 8 * H3_in[indjm1] + H3_in[indjm2]) *
        oodth12;


      // Apparently dV/dth (th=0) != 0, which can't be captured by centered finite differences and th=0 symmetry
      CCTK_REAL V_th;
      if (jj==0) {
        // 1st derivative with 4th order accuracy (forward stencils)
        V_th = (- 25 * V_in[ind] + 48 * V_in[indjp1] - 36 * V_in[indjp2] + 16 * V_in[indjp3] - 3 * V_in[indjp4]) *
          oodth12;
      } else if (jj==1) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (- 3 * V_in[indjm1] - 10 * V_in[ind] + 18 * V_in[indjp1] - 6 * V_in[indjp2] + V_in[indjp3]) * 
          oodth12;
      } else {
        // 1st derivative with 4th order accuracy (centered stencils)
        V_th = (-V_in[indjp2] + 8 * V_in[indjp1] - 8 * V_in[indjm1] + V_in[indjm2]) *
          oodth12;
      }

      CCTK_REAL Wbar_X, Wbar_Xth;
      CCTK_REAL H3_X, H3_Xth;
      CCTK_REAL V_X, V_Xth;
      CCTK_REAL H1_X = 0.;

      /*
      Regarding finite differencing orders: plotting dWbar_dr, d2Wbar_drth, there were small discontinuities near r=r_H
      during tests with the previous 2nd order accuracy for i==0 and i==1.
      Those alleviate when moving to 4th order accuracy.

      For i==NX-1 and i==NX-2, we keep 2nd order for now. The issue is not appearing as clearly,
      and they represent points which are physically far, so maybe better to keep the computation more local.
      */

      if (i == 0) {
        /* this point is X == 0, r == rH, R == rH/4. dWbar_dX goes to zero here. but
           since we're interested in dWbar_dr, and since drxdr diverges (here), we
           will use L'Hopital's rule. for that, we will write instead the 2nd
           derivative */

        // 2nd derivative with 4th order accuracy (forward stencils)
        Wbar_X = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] + 214 * Wbar_in[indip2] 
                  - 156 * Wbar_in[indip3] + 61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) * oodXsq12;
        
        // mixed derivatives with 4th order accuracy (central stencils in j (1st der) and forward in i (2nd der))
        Wbar_Xth = ( 
              -  45 * Wbar_in[indjp2] +  154 * Wbar_in[indip1jp2] -  214 * Wbar_in[indip2jp2]
                  +  156 * Wbar_in[indip3jp2] -  61 * Wbar_in[indip4jp2] + 10 * Wbar_in[indip5jp2]
              + 360 * Wbar_in[indjp1] - 1232 * Wbar_in[indip1jp1] + 1712 * Wbar_in[indip2jp1]
                  - 1248 * Wbar_in[indip3jp1] + 488 * Wbar_in[indip4jp1] - 80 * Wbar_in[indip5jp1]
              - 360 * Wbar_in[indjm1] + 1232 * Wbar_in[indip1jm1] - 1712 * Wbar_in[indip2jm1]
                  + 1248 * Wbar_in[indip3jm1] - 488 * Wbar_in[indip4jm1] + 80 * Wbar_in[indip5jm1] 
              +  45 * Wbar_in[indjm2] -  154 * Wbar_in[indip1jm2] +  214 * Wbar_in[indip2jm2]
                  -  156 * Wbar_in[indip3jm2] +  61 * Wbar_in[indip4jm2] - 10 * Wbar_in[indip5jm2]
        ) * oodXsqdth144;
        
        // Same for H3

        // 2nd derivative with 4th order accuracy (forward stencils)
        H3_X = (45 * H3_in[ind] - 154 * H3_in[indip1] + 214 * H3_in[indip2] 
                  - 156 * H3_in[indip3] + 61 * H3_in[indip4] - 10 * H3_in[indip5]) * oodXsq12;
        
        // mixed derivatives with 4th order accuracy (central stencils in j (1st der) and forward in i (2nd der))
        H3_Xth = ( 
              -  45 * H3_in[indjp2] +  154 * H3_in[indip1jp2] -  214 * H3_in[indip2jp2]
                  +  156 * H3_in[indip3jp2] -  61 * H3_in[indip4jp2] + 10 * H3_in[indip5jp2]
              + 360 * H3_in[indjp1] - 1232 * H3_in[indip1jp1] + 1712 * H3_in[indip2jp1]
                  - 1248 * H3_in[indip3jp1] + 488 * H3_in[indip4jp1] - 80 * H3_in[indip5jp1]
              - 360 * H3_in[indjm1] + 1232 * H3_in[indip1jm1] - 1712 * H3_in[indip2jm1]
                  + 1248 * H3_in[indip3jm1] - 488 * H3_in[indip4jm1] + 80 * H3_in[indip5jm1] 
              +  45 * H3_in[indjm2] -  154 * H3_in[indip1jm2] +  214 * H3_in[indip2jm2]
                  -  156 * H3_in[indip3jm2] +  61 * H3_in[indip4jm2] - 10 * H3_in[indip5jm2]
        ) * oodXsqdth144;

        // V

        // 2nd derivative with 4th order accuracy (forward stencils)
        V_X = (45 * V_in[ind] - 154 * V_in[indip1] + 214 * V_in[indip2] 
                  - 156 * V_in[indip3] + 61 * V_in[indip4] - 10 * V_in[indip5]) * oodXsq12;
        
        // For V, near the axis, we need non-centered stencils because dV/dth != 0
        // mixed derivatives with 4th order accuracy (forward stencils in i (2nd der))
        if (jj==0) {
          // 1st derivative with 4th order forward stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth = ( 
                - 1125 * V_in[ind]    + 3850 * V_in[indip1]    -  5350 * V_in[indip2]
                    + 3900 * V_in[indip3]    - 1525 * V_in[indip4]    + 250 * V_in[indip5]
                + 2160 * V_in[indjp1] - 7392 * V_in[indip1jp1] + 10272 * V_in[indip2jp1]
                    - 7488 * V_in[indip3jp1] + 2928 * V_in[indip4jp1] - 480 * V_in[indip5jp1]
                - 1620 * V_in[indjp2] + 5544 * V_in[indip1jp2] -  7704 * V_in[indip2jp2]
                    + 5616 * V_in[indip3jp2] - 2196 * V_in[indip4jp2] + 360 * V_in[indip5jp2]
                +  720 * V_in[indjp3] - 2464 * V_in[indip1jp3] +  3424 * V_in[indip2jp3]
                    - 2496 * V_in[indip3jp3] +  976 * V_in[indip4jp3] - 160 * V_in[indip5jp3]
                -  135 * V_in[indjp4] +  462 * V_in[indip1jp4] -   642 * V_in[indip2jp4]
                    +  468 * V_in[indip3jp4] -  183 * V_in[indip4jp4] +  30 * V_in[indip5jp4]
          ) * oodXsqdth144;
        } else if (jj==1) {
          // 1st derivative with 4th order stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth = ( 
                - 135 * V_in[indjm1] +  462 * V_in[indip1jm1] -  642 * V_in[indip2jm1]
                    +  468 * V_in[indip3jm1] -  183 * V_in[indip4jm1] +  30 * V_in[indip5jm1]
                - 450 * V_in[ind]    + 1540 * V_in[indip1]    - 2140 * V_in[indip2]
                    + 1560 * V_in[indip3]    -  610 * V_in[indip4]    + 100 * V_in[indip5]
                + 810 * V_in[indjp1] - 2772 * V_in[indip1jp1] + 3852 * V_in[indip2jp1]
                    - 2808 * V_in[indip3jp1] + 1098 * V_in[indip4jp1] - 180 * V_in[indip5jp1]
                - 270 * V_in[indjp2] +  924 * V_in[indip1jp2] - 1284 * V_in[indip2jp2]
                    +  936 * V_in[indip3jp2] -  366 * V_in[indip4jp2] +  60 * V_in[indip5jp2]
                +  45 * V_in[indjp3] -  154 * V_in[indip1jp3] +  214 * V_in[indip2jp3]
                    -  156 * V_in[indip3jp3] +   61 * V_in[indip4jp3] -  10 * V_in[indip5jp3]
          ) * oodXsqdth144;
        } else {
          V_Xth = ( 
            // 1st derivative with 4th order central stencils for j
                -  45 * V_in[indjp2] +  154 * V_in[indip1jp2] -  214 * V_in[indip2jp2]
                    +  156 * V_in[indip3jp2] -  61 * V_in[indip4jp2] + 10 * V_in[indip5jp2]
                + 360 * V_in[indjp1] - 1232 * V_in[indip1jp1] + 1712 * V_in[indip2jp1]
                    - 1248 * V_in[indip3jp1] + 488 * V_in[indip4jp1] - 80 * V_in[indip5jp1]
                - 360 * V_in[indjm1] + 1232 * V_in[indip1jm1] - 1712 * V_in[indip2jm1]
                    + 1248 * V_in[indip3jm1] - 488 * V_in[indip4jm1] + 80 * V_in[indip5jm1] 
                +  45 * V_in[indjm2] -  154 * V_in[indip1jm2] +  214 * V_in[indip2jm2]
                    -  156 * V_in[indip3jm2] +  61 * V_in[indip4jm2] - 10 * V_in[indip5jm2]
          ) * oodXsqdth144;
        }


        // A_r

        // For X=0 <=> x=0 <=> r=r_h, we need to compute dA_x/dX (cf below)
        // 1st derivative with 4th order accuracy (forward stencils)
        H1_X = (- 25 * H1_in[ind] + 48 * H1_in[indip1] - 36 * H1_in[indip2] + 16 * H1_in[indip3] - 3 * H1_in[indip4]) * oodX12;
        
      
      
      } else if (i == 1) {
        // 1st derivative, 4th order accuracy
        Wbar_X = (- 3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] - 6 * Wbar_in[indip2] + Wbar_in[indip3]) * oodX12;

        // 1st derivative, 4th order accuracy (central stencils in j)
        Wbar_Xth = (
             3 * Wbar_in[indim1jp2] + 10 * Wbar_in[indjp2] -  18 * Wbar_in[indip1jp2] +  6 * Wbar_in[indip2jp2] -     Wbar_in[indip3jp2]
          - 24 * Wbar_in[indim1jp1] - 80 * Wbar_in[indjp1] + 144 * Wbar_in[indip1jp1] - 48 * Wbar_in[indip2jp1] + 8 * Wbar_in[indip3jp1]
          + 24 * Wbar_in[indim1jm1] + 80 * Wbar_in[indjm1] - 144 * Wbar_in[indip1jm1] + 48 * Wbar_in[indip2jm1] - 8 * Wbar_in[indip3jm1]
          -  3 * Wbar_in[indim1jm2] - 10 * Wbar_in[indjm2] +  18 * Wbar_in[indip1jm2] -  6 * Wbar_in[indip2jm2] +     Wbar_in[indip3jm2]
        ) * oodXdth144;
        
        // 1st derivative, 4th order accuracy
        H3_X = (- 3 * H3_in[indim1] - 10 * H3_in[ind] + 18 * H3_in[indip1] - 6 * H3_in[indip2] + H3_in[indip3]) * oodX12;

        // 1st derivative, 4th order accuracy (central stencils in j)
        H3_Xth = (
             3 * H3_in[indim1jp2] + 10 * H3_in[indjp2] -  18 * H3_in[indip1jp2] +  6 * H3_in[indip2jp2] -     H3_in[indip3jp2]
          - 24 * H3_in[indim1jp1] - 80 * H3_in[indjp1] + 144 * H3_in[indip1jp1] - 48 * H3_in[indip2jp1] + 8 * H3_in[indip3jp1]
          + 24 * H3_in[indim1jm1] + 80 * H3_in[indjm1] - 144 * H3_in[indip1jm1] + 48 * H3_in[indip2jm1] - 8 * H3_in[indip3jm1]
          -  3 * H3_in[indim1jm2] - 10 * H3_in[indjm2] +  18 * H3_in[indip1jm2] -  6 * H3_in[indip2jm2] +     H3_in[indip3jm2]
        ) * oodXdth144;
        
        // 1st derivative, 4th order accuracy
        V_X = (- 3 * V_in[indim1] - 10 * V_in[ind] + 18 * V_in[indip1] - 6 * V_in[indip2] + V_in[indip3]) * oodX12;

        // For V, near the axis, we need non-centered stencils because dV/dth != 0
        // mixed 1st derivatives with 4th order accuracy
        if (jj==0) {
          // 1st derivative with 4th order forward stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth = (
               75 * V_in[indim1]    + 250 * V_in[ind]    - 450 * V_in[indip1]    + 150 * V_in[indip2]    - 25 * V_in[indip3]
            - 144 * V_in[indim1jp1] - 480 * V_in[indjp1] + 864 * V_in[indip1jp1] - 288 * V_in[indip2jp1] + 48 * V_in[indip3jp1]
            + 108 * V_in[indim1jp2] + 360 * V_in[indjp2] - 648 * V_in[indip1jp2] + 216 * V_in[indip2jp2] - 36 * V_in[indip3jp2]
            -  48 * V_in[indim1jp3] - 160 * V_in[indjp3] + 288 * V_in[indip1jp3] -  96 * V_in[indip2jp3] + 16 * V_in[indip3jp3]  
            +   9 * V_in[indim1jp4] +  30 * V_in[indjp4] -  54 * V_in[indip1jp4] +  18 * V_in[indip2jp4] -  3 * V_in[indip3jp4]
          ) * oodXdth144;
        } else if (jj==1) {
          // 1st derivative with 4th order stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth = (
               9 * V_in[indim1jm1] +  30 * V_in[indjm1] -  54 * V_in[indip1jm1] +  18 * V_in[indip2jm1] -  3 * V_in[indip3jm1]
            + 30 * V_in[indim1]    + 100 * V_in[ind]    - 180 * V_in[indip1]    +  60 * V_in[indip2]    - 10 * V_in[indip3]
            - 54 * V_in[indim1jp1] - 180 * V_in[indjp1] + 324 * V_in[indip1jp1] - 108 * V_in[indip2jp1] + 18 * V_in[indip3jp1]
            + 18 * V_in[indim1jp2] +  60 * V_in[indjp2] - 108 * V_in[indip1jp2] +  36 * V_in[indip2jp2] -  6 * V_in[indip3jp2]
            -  3 * V_in[indim1jp3] -  10 * V_in[indjp3] +  18 * V_in[indip1jp3] -   6 * V_in[indip2jp3] +      V_in[indip3jp3]
          ) * oodXdth144;
        } else {
          // 1st derivative with 4th order central stencils for j
          V_Xth = (
               3 * V_in[indim1jp2] + 10 * V_in[indjp2] -  18 * V_in[indip1jp2] +  6 * V_in[indip2jp2] -     V_in[indip3jp2]
            - 24 * V_in[indim1jp1] - 80 * V_in[indjp1] + 144 * V_in[indip1jp1] - 48 * V_in[indip2jp1] + 8 * V_in[indip3jp1]
            + 24 * V_in[indim1jm1] + 80 * V_in[indjm1] - 144 * V_in[indip1jm1] + 48 * V_in[indip2jm1] - 8 * V_in[indip3jm1]
            -  3 * V_in[indim1jm2] - 10 * V_in[indjm2] +  18 * V_in[indip1jm2] -  6 * V_in[indip2jm2] +     V_in[indip3jm2]
          ) * oodXdth144;
        }



      } else if (i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;

        Wbar_Xth  = ( Wbar_in[indip1jp1] - Wbar_in[indip1jm1] - Wbar_in[indim1jp1] + Wbar_in[indim1jm1] ) * oodXdth4;
        
        // 1st derivative with 2nd order accuracy (central stencils)
        H3_X = (-H3_in[indim1] + H3_in[indip1]) * 0.5 * oodX;

        H3_Xth  = ( H3_in[indip1jp1] - H3_in[indip1jm1] - H3_in[indim1jp1] + H3_in[indim1jm1] ) * oodXdth4;

        // 1st derivative with 2nd order accuracy (central stencils)
        V_X = (-V_in[indim1] + V_in[indip1]) * 0.5 * oodX;

        // For V, near the axis, we need non-centered stencils because dV/dth != 0
        // mixed 1st derivatives with 2nd order accuracy
        if (jj==0) {
          // 1st derivative with 2nd order forward stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth  = ( - 3*V_in[indip1] + 4*V_in[indip1jp1] - V_in[indip1jp2] 
                     + 3*V_in[indim1] - 4*V_in[indim1jp1] + V_in[indim1jp2] ) * oodXdth4;
        } else {
          // 1st derivative with 2nd order central stencils for j
          V_Xth  = ( V_in[indip1jp1] - V_in[indip1jm1] - V_in[indim1jp1] + V_in[indim1jm1] ) * oodXdth4;
        }


      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4*Wbar_in[indim1] + 3*Wbar_in[ind]) * 0.5 * oodX;
        Wbar_Xth = 0.; // we don't actually use this variable at large r, so just
                       // set it to zero
        
        // 1st derivative with 2nd order accuracy (backward stencils)
        H3_X = (H3_in[indim2] - 4*H3_in[indim1] + 3*H3_in[ind]) * 0.5 * oodX;
        H3_Xth = 0.; // we don't actually use this variable at large r, so just
                     // set it to zero

        // 1st derivative with 2nd order accuracy (backward stencils)
        V_X = (V_in[indim2] - 4*V_in[indim1] + 3*V_in[ind]) * 0.5 * oodX;
        V_Xth = 0.; // we don't actually use this variable at large r, so just
                    // set it to zero

      } else {
        // 4th order accurate stencils
        Wbar_X    = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
        Wbar_Xth  = (
            -Wbar_in[indim2jp2] +  8*Wbar_in[indim1jp2] -  8*Wbar_in[indip1jp2] +   Wbar_in[indip2jp2]
         + 8*Wbar_in[indim2jp1] - 64*Wbar_in[indim1jp1] + 64*Wbar_in[indip1jp1] - 8*Wbar_in[indip2jp1]
         - 8*Wbar_in[indim2jm1] + 64*Wbar_in[indim1jm1] - 64*Wbar_in[indip1jm1] + 8*Wbar_in[indip2jm1]
         +   Wbar_in[indim2jm2] -  8*Wbar_in[indim1jm2] +  8*Wbar_in[indip1jm2] -   Wbar_in[indip2jm2] ) * oodXdth144;

        // 4th order accurate stencils
        H3_X    = (-H3_in[indip2] + 8 * H3_in[indip1] - 8 * H3_in[indim1] + H3_in[indim2]) * oodX12;
        H3_Xth  = (
            -H3_in[indim2jp2] +  8*H3_in[indim1jp2] -  8*H3_in[indip1jp2] +   H3_in[indip2jp2]
         + 8*H3_in[indim2jp1] - 64*H3_in[indim1jp1] + 64*H3_in[indip1jp1] - 8*H3_in[indip2jp1]
         - 8*H3_in[indim2jm1] + 64*H3_in[indim1jm1] - 64*H3_in[indip1jm1] + 8*H3_in[indip2jm1]
         +   H3_in[indim2jm2] -  8*H3_in[indim1jm2] +  8*H3_in[indip1jm2] -   H3_in[indip2jm2] ) * oodXdth144;

        // 4th order accurate stencils
        V_X    = (-V_in[indip2] + 8 * V_in[indip1] - 8 * V_in[indim1] + V_in[indim2]) * oodX12;
        // For V, near the axis, we need non-centered stencils because dV/dth != 0
        // mixed 1st derivatives with 4th order accuracy
        if (jj==0) {
          // 1st derivative with 4th order forward stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth  = (
            - 25 * V_in[indim2]    + 200 * V_in[indim1]    - 200 * V_in[indip1]    + 25 * V_in[indip2]
            + 48 * V_in[indim2jp1] - 384 * V_in[indim1jp1] + 384 * V_in[indip1jp1] - 48 * V_in[indip2jp1]
            - 36 * V_in[indim2jp2] + 288 * V_in[indim1jp2] - 288 * V_in[indip1jp2] + 36 * V_in[indip2jp2]
            + 16 * V_in[indim2jp3] - 128 * V_in[indim1jp3] + 128 * V_in[indip1jp3] - 16 * V_in[indip2jp3]
            -  3 * V_in[indim2jp4] +  24 * V_in[indim1jp4] -  24 * V_in[indip1jp4] +  3 * V_in[indip2jp4]
          ) * oodXdth144;
        } else if (jj==1) {
          // 1st derivative with 4th order stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth  = (
            -  3 * V_in[indim2jm1] +  24 * V_in[indim1jm1] -  24 * V_in[indip1jm1] +  3 * V_in[indip2jm1]
            - 10 * V_in[indim2]    +  80 * V_in[indim1]    -  80 * V_in[indip1]    + 10 * V_in[indip2]
            + 18 * V_in[indim2jp1] - 144 * V_in[indim1jp1] + 144 * V_in[indip1jp1] - 18 * V_in[indip2jp1]
            -  6 * V_in[indim2jp2] +  48 * V_in[indim1jp2] -  48 * V_in[indip1jp2] +  6 * V_in[indip2jp2]
            +      V_in[indim2jp3] -   8 * V_in[indim1jp3] +   8 * V_in[indip1jp3] -      V_in[indip2jp3]
          ) * oodXdth144;
        } else {
          // 1st derivative with 4th order central stencils for j
          // /!\ Reverse order of j points compared to the rest
          V_Xth  = (
              -V_in[indim2jp2] +  8*V_in[indim1jp2] -  8*V_in[indip1jp2] +   V_in[indip2jp2]
           + 8*V_in[indim2jp1] - 64*V_in[indim1jp1] + 64*V_in[indip1jp1] - 8*V_in[indip2jp1]
           - 8*V_in[indim2jm1] + 64*V_in[indim1jm1] - 64*V_in[indip1jm1] + 8*V_in[indip2jm1]
           +   V_in[indim2jm2] -  8*V_in[indim1jm2] +  8*V_in[indip1jm2] -   V_in[indip2jm2] ) * oodXdth144;
        }
      }

      // from the X coordinate used in the input files to the x coordinate
      // We need to be careful at X == 1 for radial derivatives (coordinate change is singular)
      if (i == NX - 1) {
        dWbar_dr_extd[ind]    = 0.; // Sensibly, dWbar_dr_in should vanish (dXdr == 0, Wbar_X bounded)
        d2Wbar_drth_extd[ind] = 0.; // Wbar_Xth is set to 0 above anyway

        // Same for H3, V
        dH3_dr_extd[ind]    = 0.;
        d2H3_drth_extd[ind] = 0.;

        dV_dr_extd[ind]    = 0.;
        d2V_drth_extd[ind] = 0.;


        H1r_extd[ind] = 0.; // A_r = 0 at infinity

      } else {
        const CCTK_REAL rx = C0*lX/(1. - lX);
        // from the x coordinate to the metric coordinate r
        // const CCTK_REAL rr = sqrt(rH*rH + rx*rx);

        // corresponding derivatives
        const CCTK_REAL dXdrx = 1./(C0 + rx) - rx/((C0 + rx)*(C0 + rx));

        CCTK_REAL drxdr;
        if (i == 0) { // rx == 0 (X == 0)
          // Including here additional C0 contribution from L'Hôpital's rule: d(dr/dx)/dX = dx/dX * d2r/dx2 = C0 / rH at i==0
          drxdr = rH / C0;
        } else {
          drxdr = sqrt(rH*rH + rx*rx)/rx;
        }
        const CCTK_REAL dXdr = dXdrx * drxdr;

        dWbar_dr_extd[ind]    = dXdr * Wbar_X;
        d2Wbar_drth_extd[ind] = dXdr * Wbar_Xth;
        
        dH3_dr_extd[ind]      = dXdr * H3_X;
        d2H3_drth_extd[ind]   = dXdr * H3_Xth;
        
        dV_dr_extd[ind]       = dXdr * V_X;
        d2V_drth_extd[ind]    = dXdr * V_Xth;

        // A_r
        /*
          A_r = dx/dr * A_x = r/x * A_x
          
          At x=0, we need to regularize
            A_r ~ (dA_x/dX) / d(dr/dx)/dX    
          the latter term being "drxdr" computed above to include L'Hôpital's rule in that case
        */
        if (i==0) { // rx == 0 (X == 0)
          H1r_extd[ind] = drxdr * H1_X;
        } else {
          H1r_extd[ind] = drxdr * H1_in[ind];
        }
      } // if/else i == NX - 1

      dWbar_dth_extd[ind]   = Wbar_th;
      dH3_dth_extd[ind]     = H3_th;
      dV_dth_extd[ind]      = V_th;

      // fprintf (debugfile, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f ", 
      //             Wbar_in[ind], dWbar_dr_in[ind], dWbar_dth_in[ind], d2Wbar_drth_in[ind],
      //             H3_in[ind], dH3_dr_in[ind], dH3_dth_in[ind], d2H3_drth_in[ind],
      //             V_in[ind], dV_dr_in[ind], dV_dth_in[ind], d2V_drth_in[ind]);
    } // for i
    // fprintf (debugfile, "\n");
  } // for jj

  // fclose(debugfile);


  // Second loop on z<0 half-space (completion by symmetry)

  // Even parity: F1, F2, F0, Wbar, H1, H3, V and their r derivatives
  // Odd parity:  H2, and theta derivatives of even functions

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
      H1r_extd[ind]      = H1r_extd[indsym];
      
      Wbar_extd[ind]     = Wbar_extd[indsym];
      H3_extd[ind]       = H3_extd[indsym];
      V_extd[ind]        = V_extd[indsym];

      dWbar_dr_extd[ind] = dWbar_dr_extd[indsym];
      dH3_dr_extd[ind]   = dH3_dr_extd[indsym];
      dV_dr_extd[ind]    = dV_dr_extd[indsym];

      // Odd
      H2_extd[ind]          = - H2_extd[indsym];

      dWbar_dth_extd[ind]   = - dWbar_dth_extd[indsym];
      dH3_dth_extd[ind]     = - dH3_dth_extd[indsym];
      dV_dth_extd[ind]      = - dV_dth_extd[indsym];
      
      d2Wbar_drth_extd[ind] = - d2Wbar_drth_extd[indsym];
      d2H3_drth_extd[ind]   = - d2H3_drth_extd[indsym];
      d2V_drth_extd[ind]    = - d2V_drth_extd[indsym];

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

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;

        CCTK_REAL RR  = sqrt(RR2);
        /* note that there are divisions by RR in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        // from (quasi-)isotropic coordinate R to the metric coordinate r
        const CCTK_REAL rr = RR * (1. + 0.25 * rH / RR) * (1. + 0.25 * rH / RR);

        // from the metric coordinate r to the x coordinate
        const CCTK_REAL rx = sqrt(rr*rr - rH*rH);

        // and finally to the X radial coordinate (used in input files)
        const CCTK_REAL lX = rx / (C0 + rx);

        const CCTK_REAL ltheta = acos( z1/RR );

        X_g[ind]     = lX;
        theta_g[ind] = ltheta;
      }
    }
  }

  /* now for the interpolation */

  const CCTK_INT N_dims  = 2;   // 2-D interpolation

  const CCTK_INT N_input_arrays  = 17;
  const CCTK_INT N_output_arrays = 17;

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

  /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
  input_arrays[0] = (const void *) F1_extd;
  input_arrays[1] = (const void *) F2_extd;
  input_arrays[2] = (const void *) F0_extd;
  input_arrays[3] = (const void *) Wbar_extd;
  input_arrays[4] = (const void *) dWbar_dr_extd;
  input_arrays[5] = (const void *) dWbar_dth_extd;
  input_arrays[6] = (const void *) d2Wbar_drth_extd;
  input_arrays[7] = (const void *) H1r_extd;
  input_arrays[8] = (const void *) H2_extd;
  input_arrays[9] = (const void *) H3_extd;
  input_arrays[10]= (const void *) V_extd;
  input_arrays[11]= (const void *) dH3_dr_extd;
  input_arrays[12]= (const void *) dH3_dth_extd;
  input_arrays[13]= (const void *) d2H3_drth_extd;
  input_arrays[14]= (const void *) dV_dr_extd;
  input_arrays[15]= (const void *) dV_dth_extd;
  input_arrays[16]= (const void *) d2V_drth_extd;

  /* output arrays */
  void *output_arrays[N_output_arrays];
  CCTK_INT output_array_type_codes[N_output_arrays];
  CCTK_REAL *F1, *F2, *F0, *Wbar, *H1r, *H2, *H3, *V;
  CCTK_REAL *dWbar_dr, *dWbar_dth, *d2Wbar_drth;
  CCTK_REAL *dH3_dr, *dH3_dth, *d2H3_drth;
  CCTK_REAL *dV_dr, *dV_dth, *d2V_drth;

  F1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F0          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  Wbar        = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dWbar_dr    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dWbar_dth   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  d2Wbar_drth = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H1r         = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  H3          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  V           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dr      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dth     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  d2H3_drth   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dr       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dth      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  d2V_drth    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

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

  output_arrays[0] = (void *) F1;
  output_arrays[1] = (void *) F2;
  output_arrays[2] = (void *) F0;
  output_arrays[3] = (void *) Wbar;
  output_arrays[4] = (void *) dWbar_dr;
  output_arrays[5] = (void *) dWbar_dth;
  output_arrays[6] = (void *) d2Wbar_drth;
  output_arrays[7] = (void *) H1r;
  output_arrays[8] = (void *) H2;
  output_arrays[9] = (void *) H3;
  output_arrays[10]= (void *) V;
  output_arrays[11]= (void *) dH3_dr;
  output_arrays[12]= (void *) dH3_dth;
  output_arrays[13]= (void *) d2H3_drth;
  output_arrays[14]= (void *) dV_dr;
  output_arrays[15]= (void *) dV_dth;
  output_arrays[16]= (void *) d2V_drth;


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
  
  free(F1_extd); free(F2_extd); free(F0_extd); free(Wbar_extd);
  free(H1r_extd); free(H2_extd); free(H3_extd); free(V_extd);
  free(dWbar_dr_extd); free(dWbar_dth_extd); free(d2Wbar_drth_extd);
  free(dH3_dr_extd); free(dH3_dth_extd); free(d2H3_drth_extd);
  free(dV_dr_extd); free(dV_dth_extd); free(d2V_drth_extd);


  /* printf("F1 = %g\n", F1[0]); */
  /* printf("F2 = %g\n", F2[0]); */
  /* printf("F0 = %g\n", F0[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W = %g\n", W[0]); */


  /* now we finally write the metric and all 3+1 quantities. first we write the
     3-metric and extrinsic curvature, then Proca fields, then lapse and shift */

  const CCTK_REAL tt = cctk_time;
  const CCTK_REAL omega = mm * OmegaH;

  const CCTK_REAL coswt = cos(omega * tt);
  const CCTK_REAL sinwt = sin(omega * tt);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        /* note that there are divisions by RR in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        const CCTK_REAL RR  = sqrt(RR2);

        // from (quasi-)isotropic coordinate R to the metric coordinate r
        const CCTK_REAL rr = RR * (1. + 0.25 * rH / RR) * (1. + 0.25 * rH / RR);
        const CCTK_REAL rrP = pow(rr, Wbar_r_power);
        const CCTK_REAL rrPp1 = rr*rrP;

        /*
        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        const CCTK_REAL rho  = sqrt(rho2);
        */

        const CCTK_REAL costh  = z1/RR;
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

        // from Wbar to W function, Wbar = r^p * W
        const CCTK_REAL W        = Wbar[ind] / rrP;
        const CCTK_REAL dW_dth   = dWbar_dth[ind] / rrP;
        const CCTK_REAL dW_dr    = dWbar_dr[ind] / rrP - Wbar_r_power * Wbar[ind] / rrPp1;
        const CCTK_REAL d2W_drth = d2Wbar_drth[ind] / rrP - Wbar_r_power * dWbar_dth[ind] / rrPp1;

        // add non-axisymmetric perturbation on conformal factor
        // NOTE: the perturbation is only taken into account for the 3-metric grid functions (not extrinsic curvature, lapse, ...)
        const CCTK_REAL argpert_cf = (RR - R0pert_conf_fac)/Sigmapert_conf_fac;
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
          KRph/(R sin(th)^2)  = - 1/2 exp(2F2-F0) (1 + rH/(4R))^6 R dW/dr
          Kthph/(R sin(th))^3 = - 1/2 exp(2F2-F0) (1 + rH/(4R))^5 dW/dth / sin(th) 1/(R - rH/4)
        */

        // KRph/(R sin(th)^2)
        const CCTK_REAL KRph_o_Rsinth2 = -0.5 * exp(2. * F2[ind] - F0[ind]) * aux6 * RR * dW_dr;

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
          dWdth_o_sinth_den = (eps_o_1meps - eps_o_1meps*eps_o_1meps) * d2W_drth / sinth;
        else
          dWdth_o_sinth_den = dW_dth / (den * sinth);
        
        // Kthph/(R sin(th))^3
        const CCTK_REAL Kthph_o_R3sinth3 = -0.5 * exp(2. * F2[ind] - F0[ind]) * aux5 * dWdth_o_sinth_den;

        // extrinsic curvature
        kxx[ind] = 2.*KRph_o_Rsinth2 *  x1 * sinth2ph_x                     +  2.*Kthph_o_R3sinth3 *  Rsinthth_x * R2sinth2ph_x;
        kxy[ind] =    KRph_o_Rsinth2 * (x1 * sinth2ph_y + y1 * sinth2ph_x)  +     Kthph_o_R3sinth3 * (Rsinthth_x * R2sinth2ph_y + Rsinthth_y * R2sinth2ph_x);
        kxz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_x   +     Kthph_o_R3sinth3 *                              Rsinthth_z * R2sinth2ph_x;
        kyy[ind] = 2.*KRph_o_Rsinth2 *  y1 * sinth2ph_y                     +  2.*Kthph_o_R3sinth3 *  Rsinthth_y * R2sinth2ph_y;
        kyz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_y   +     Kthph_o_R3sinth3 *                              Rsinthth_z * R2sinth2ph_y;
        kzz[ind] = 0.;

        
        // lapse value (field initialization below)
        const CCTK_REAL alph = exp(F0[ind]) * den / (RR + 0.25*rH);


        // let's add a perturbation to the Proca field as well
        // NOTE: the perturbation is added directed to every instance of e^{i m \varphi}, hence its derivatives are not taken into account
        // TODO (?): Design perturbation more generically as ~ cos((m+1)\varphi)
        const CCTK_REAL argpert_Proca = (RR - R0pert_Proca)/Sigmapert_Proca;
        const CCTK_REAL pert_Proca = 1. + Apert_Proca * (x1*x1 - y1*y1)*mu*mu * exp( -0.5*argpert_Proca*argpert_Proca );


        // ----- Proca fields -----

        // Real and imaginay part of the harmonic dependence: exp[i(m\varphi - \omega t)]
        const CCTK_REAL harm_re = (coswt * cosmph + sinwt * sinmph) * pert_Proca;
        const CCTK_REAL harm_im = (coswt * sinmph - sinwt * cosmph) * pert_Proca;

        // Radial component in R coordinate
        // A_R = dr/dR * A_r 
        // dr/dR = (1-rH/(4R))*(1+rH/(4R)) = eps * aux
        const CCTK_REAL H1R = eps * aux * H1r[ind];


        // A_x
        A1x[ind] = x1/RR * H1R * harm_re + costh*cosph/RR * H2[ind] * harm_re + sinph/RR * H3[ind] * harm_im;
        A2x[ind] = x1/RR * H1R * harm_im + costh*cosph/RR * H2[ind] * harm_im - sinph/RR * H3[ind] * harm_re;
        
        // A_y
        A1y[ind] = y1/RR * H1R * harm_re + costh*sinph/RR * H2[ind] * harm_re - cosph/RR * H3[ind] * harm_im;
        A2y[ind] = y1/RR * H1R * harm_im + costh*sinph/RR * H2[ind] * harm_im + cosph/RR * H3[ind] * harm_re;
        
        // A_z
        A1z[ind] = (z1/RR * H1R - sinth/RR * H2[ind]) * harm_re;
        A2z[ind] = (z1/RR * H1R - sinth/RR * H2[ind]) * harm_im;

        // A_\phi
        /*
          A_\phi = -n^\mu A_\mu = - (A_t + W*A_ph)/alpha
                = -i * e^{i (m ph - w t)} * (V + W H3 sinth) / alpha
          
          if at R ~ rH/4 we need to regularize the division by R - rH/4
          That's the same as for the extrinsic curvature above, with f(R) = V + W H3 sinth  (= 0 on horizon)
        */
        if (fabs(eps) < 1e-4) {
          // f(R) = V + W H3 sinth
          const CCTK_REAL df_dr = dV_dr[ind] + dW_dr * H3[ind] * sinth + W * dH3_dr[ind] * sinth;
          const CCTK_REAL reg  = exp(-F0[ind]) * (RR + 0.25*rH) * (eps_o_1meps - eps_o_1meps*eps_o_1meps) * df_dr;

          Aphi1[ind] =  reg * harm_im;
          Aphi2[ind] = -reg * harm_re;
        } else {
          Aphi1[ind] = (V[ind] + W * sinth * H3[ind]) / alph * harm_im;
          Aphi2[ind] =-(V[ind] + W * sinth * H3[ind]) / alph * harm_re;
        }

        // ----- Electric fields -----
        
        // First, E_i in (R, th, ph) coordinates
        CCTK_REAL E1d_R, E2d_R, E1d_th, E2d_th, E1d_ph_o_sinth, E2d_ph_o_sinth;

        // E_R
        /*
          E_R = i * e^{i(m phi - w t)} * aux^2 * e^{-F_0} * [-m * (W - OmegaH) H1r + dV/dr + W sinth dH3/dr]

          aux contribution is canceled with inverse metric below.
        */

        E1d_R = -exp(-F0[ind]) * (- mm * (W - OmegaH) * H1r[ind] + dV_dr[ind] + W * sinth * dH3_dr[ind]) * harm_im;
        E2d_R =  exp(-F0[ind]) * (- mm * (W - OmegaH) * H1r[ind] + dV_dr[ind] + W * sinth * dH3_dr[ind]) * harm_re;

        // E_th
        /*
          E_th = i * e^{i(m phi - w t)} / alpha * [- m * (W - OmegaH) * H2 + dV/dth + W * d(H3 * sinth)/dth]
          if at R ~ rH/4 we need to regularize the division by R - rH/4
          That's the same as for A_\phi, with 
            1) f(R) = W - OmegaH
            2) f(R) = dV/dth + W d(H3 sinth)/dth   (= 0 on horizon, it coincides with d(V + OmegaH H3 sinth)/dth there)
        */
        if (fabs(eps) < 1e-4) {
          // f(R) = dV/dth + W d(H3 sinth)/dth
          const CCTK_REAL df_dr = d2V_drth[ind] + dW_dr * (sinth * dH3_dth[ind] + costh * H3[ind]) + W * (sinth * d2H3_drth[ind] + costh * dH3_dr[ind]);
          const CCTK_REAL reg  = exp(-F0[ind]) * (RR + 0.25*rH) * (eps_o_1meps - eps_o_1meps*eps_o_1meps) * (- mm * H2[ind] * dW_dr + df_dr );

          E1d_th = -reg * harm_im;
          E2d_th =  reg * harm_re;
        } else {
          E1d_th = -(- mm * (W - OmegaH) * H2[ind] + dV_dth[ind] + W * (sinth * dH3_dth[ind] + costh * H3[ind])) / alph * harm_im;
          E2d_th =  (- mm * (W - OmegaH) * H2[ind] + dV_dth[ind] + W * (sinth * dH3_dth[ind] + costh * H3[ind])) / alph * harm_re;
        }

        // E_ph / sinth
        /*
          E_ph = - m * (V + OmegaH * H3 * sinth) / alpha * e^{i(m phi - w t)}

          We include here the division by sinth which arises when computing E^ph.
          
          if at R ~ rH/4 we need to regularize the division by R - rH/4
          That's the same as for A_\phi, with f(R) = V + OmegaH H3 sinth   (= 0 on horizon)

          on the axis sinth=0, we need to regularize the division by sinth
          We have V(theta=0,pi) = 0, so with l'Hôpital's rule
          V / sinth ~ \pm dV/dth    for theta = 0, pi resp.

          On the horizon we need to combine both. In that case, f(R) = dV/dth + OmegaH H3   (corresponds to d(V + OmegaH H3 sinth)/dth = 0)

          Using abusive shorthands for Ahor = \xi^\mu A_\mu    \propto    V + OmegaH H3 sinth
        */
        if (fabs(eps) < 1e-4 && fabs(sinth) < 1e-8) {
          const CCTK_INT zsign = (costh>=0) ? 1 : -1; // costh==0 shouldn't happen on the axis for a grid point, this would mean RR==0 too...

          const CCTK_REAL d2Ahor_drth = zsign * d2V_drth[ind] + OmegaH * dH3_dr[ind];
          const CCTK_REAL reg  = exp(-F0[ind]) * (RR + 0.25*rH) * (eps_o_1meps - eps_o_1meps*eps_o_1meps) * d2Ahor_drth;

          E1d_ph_o_sinth = - mm * reg * harm_re;
          E2d_ph_o_sinth = - mm * reg * harm_im;

        } else if (fabs(eps) < 1e-4) {
          const CCTK_REAL dAhor_dr_o_sinth = dV_dr[ind] / sinth + OmegaH * dH3_dr[ind];
          const CCTK_REAL reg  = exp(-F0[ind]) * (RR + 0.25*rH) * (eps_o_1meps - eps_o_1meps*eps_o_1meps) * dAhor_dr_o_sinth;

          E1d_ph_o_sinth = - mm * reg * harm_re;
          E2d_ph_o_sinth = - mm * reg * harm_im;
        
        } else if (fabs(sinth) < 1e-8) {
          const CCTK_INT zsign = (costh>=0) ? 1 : -1; // costh==0 shouldn't happen on the axis for a grid point, this would mean RR==0 too...

          E1d_ph_o_sinth = - mm * (zsign * dV_dth[ind] + OmegaH * H3[ind]) / alph * harm_re;
          E2d_ph_o_sinth = - mm * (zsign * dV_dth[ind] + OmegaH * H3[ind]) / alph * harm_im;
        
        } else {
          E1d_ph_o_sinth = - mm * (V[ind] / sinth + OmegaH * H3[ind]) / alph * harm_re;
          E2d_ph_o_sinth = - mm * (V[ind] / sinth + OmegaH * H3[ind]) / alph * harm_im;
        }


        // E^i components
        // Spherical auxiliaries
        
        // E^R/R = e^{-2*F1} / aux^4 * E_R / R
        // aux contribution not included in E_R above
        const CCTK_REAL E1u_R_o_R = exp(-2*F1[ind]) / aux2 * E1d_R / RR;
        const CCTK_REAL E2u_R_o_R = exp(-2*F1[ind]) / aux2 * E2d_R / RR;

        // E^th = e^{-2*F1} / aux^4 / R^2 * E_th
        const CCTK_REAL E1u_th = exp(-2*F1[ind]) / RR2 / aux4 * E1d_th;
        const CCTK_REAL E2u_th = exp(-2*F1[ind]) / RR2 / aux4 * E2d_th;

        // E^ph = e^{-2*F2} / aux^4 / R^2 / sinth^2 * E_ph
        // We compute R * sinth * E^ph. The other division by sinth is managed with E_ph above 
        const CCTK_REAL RsinthE1u_ph = exp(-2*F2[ind]) / aux4 / RR * E1d_ph_o_sinth;
        const CCTK_REAL RsinthE2u_ph = exp(-2*F2[ind]) / aux4 / RR * E2d_ph_o_sinth;


        // Finally Cartesian components
        // E^x
        E1x[ind] = x1 * E1u_R_o_R + z1 * cosph * E1u_th - sinph * RsinthE1u_ph;
        E2x[ind] = x1 * E2u_R_o_R + z1 * cosph * E2u_th - sinph * RsinthE2u_ph;

        // E^y
        E1y[ind] = y1 * E1u_R_o_R + z1 * sinph * E1u_th + cosph * RsinthE1u_ph;
        E2y[ind] = y1 * E2u_R_o_R + z1 * sinph * E2u_th + cosph * RsinthE2u_ph;

        // E^z
        E1z[ind] = z1 * E1u_R_o_R - RR * sinth * E1u_th;
        E2z[ind] = z1 * E2u_R_o_R - RR * sinth * E2u_th;




        // zero-initialize constraint damping variable Zeta
        Zeta1[ind] = 0;
        Zeta2[ind] = 0;


        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "ProcaHairyBH")) {
          alp[ind] = alph;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // shift
        if (CCTK_EQUALS(initial_shift, "ProcaHairyBH")) {
          betax[ind] =  W * y1;
          betay[ind] = -W * x1;
          betaz[ind] =  0.;
        }

      } /* for i */
    }   /* for j */
  }     /* for k */

  free(F1); free(F2); free(F0); free(Wbar);
  free(H1r); free(H2); free(H3); free(V);
  free(dWbar_dr); free(dWbar_dth); free(d2Wbar_drth);
  free(dH3_dr); free(dH3_dth); free(d2H3_drth);
  free(dV_dr); free(dV_dth); free(d2V_drth);

  return;
}
