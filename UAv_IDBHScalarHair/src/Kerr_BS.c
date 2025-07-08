//este e o codigo copiado do kerrnewman.
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

void UAv_ID_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL []);


void check_nan_or_inf(const char* var_name, double value) {
    if (isnan(value)) {
        fprintf(stderr, "Error: %s is NaN\n", var_name);
        abort(); // Break execution
    } else if (isinf(value)) {
        fprintf(stderr, "Error: %s is Inf\n", var_name);
        abort(); // Break execution
    }
}


void UAv_ID_Kerr_BS(CCTK_ARGUMENTS)
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

  CCTK_REAL *Xtmp, *thtmp, *F1_in, *F2_in, *F0_in, *phi0_in, *Wbar_in;
  Xtmp     = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  thtmp    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F1_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F2_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F0_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  phi0_in  = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  Wbar_in  = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));

  // we get the data from the input file
  UAv_ID_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, phi0_in, Wbar_in);

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


  // now we need to take the derivatives of the Wbar function
  // Then we convert to W and store the values

  CCTK_REAL *W_in, *dW_dr_in, *dW_dth_in;
  W_in        = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dr_in    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dth_in   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  const CCTK_REAL oodX       = 1. / dX;
  // const CCTK_REAL oodXsq     = oodX * oodX;
  const CCTK_REAL oodX12     = 1. / (12. * dX);
  const CCTK_REAL oodXsq12   = oodX * oodX12;
  const CCTK_REAL oodth12    = 1. / (12. * dtheta);

  for (int jj = 0; jj < Ntheta; jj++) {
    for (int i = 0; i < NX; i++) {

      CCTK_INT j, jm1, jm2, jp1, jp2;
      /* let's use the fact that the solution is axi-symmetric (and that
         theta[0] = 0) for the boundary points in j */
      if (jj == 0) {
        j   = jj;
        jp1 = jj+1;
        jp2 = jj+2;
        jm1 = jj+1;
        jm2 = jj+2;
      } else if (jj == 1) {
        j   = jj;
        jp1 = jj+1;
        jp2 = jj+2;
        jm1 = jj-1;
        jm2 = jj;
      } else if (jj == Ntheta - 2) {
        j   = jj;
        jm1 = jj-1;
        jm2 = jj-2;
        jp1 = jj+1;
        jp2 = jj;
      } else if (jj == Ntheta - 1) {
        j   = jj;
        jm1 = jj-1;
        jm2 = jj-2;
        jp1 = jj-1;
        jp2 = jj-2;
      } else {
        j   = jj;
        jp1 = jj+1;
        jp2 = jj+2;
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

      const CCTK_INT indjm1 = i + jm1*NX;
      const CCTK_INT indjm2 = i + jm2*NX;
      const CCTK_INT indjp1 = i + jp1*NX;
      const CCTK_INT indjp2 = i + jp2*NX;


      const CCTK_REAL lX = X[i];
      /* const CCTK_REAL lth = theta[j]; */
      /* printf("X[%3d] = %lf\n", i, lX); */


      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL Wbar_th = (-Wbar_in[indjp2] + 8 * Wbar_in[indjp1] - 8 * Wbar_in[indjm1] + Wbar_in[indjm2]) *
        oodth12;

      CCTK_REAL Wbar_X;
      CCTK_REAL Wbar_XX = 0.; // Used for r=0 (i==0), if Wbar_r_power == 2.

      /*
      Regarding finite differencing orders: plotting W and dW_dr, there were small discontinuities near r=0
      during tests with the previous 2nd order accuracy for i==0 and i==1.
      Those vanish when moving to 4th order accuracy.

      For i==NX-1 and i==NX-2, we keep 2nd order for now. The issue is not appearing as clearly,
      and they represent points which are physically far, so maybe better to keep the computation more local.
      */

      if (i == 0) {
        /* For the Boson Star, there's no issue, dWbar/dX != 0 at X==0, and x and r coordinates coincide. */

        // 1st derivative with 4th order accuracy (forward stencils)
        Wbar_X =(- 25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] + 16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) * oodX12;

        if (Wbar_r_power == 2) {
          // If Wbar = r^2 * W, to compute W(r=0), we need to compute Wbar_XX.
          // 2nd derivative with 4th order accuracy (forward stencils)
          Wbar_XX = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] + 214 * Wbar_in[indip2] 
                    - 156 * Wbar_in[indip3] + 61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) * oodXsq12;
        }

      } else if (i == 1 ) {
        // 1st derivative, 4th order accuracy
        Wbar_X = (- 3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] - 6 * Wbar_in[indip2] + Wbar_in[indip3]) * oodX12;

      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4*Wbar_in[indim1] + 3*Wbar_in[ind]) * 0.5 * oodX;

      } else if (i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;

      } else {
        // 4th order accurate stencils
        Wbar_X    = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
      
      }

      // From the X coordinate used in the input files to the r coordinate (coincides with x for the Boson Star, rH=0).
      // We also do the conversion from Wbar to W here, to tackle r = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  r == 0
      if (i == 0) {
        // At r=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_in[ind]    = 0.; 
        dW_dth_in[ind]   = 0.; 

        // For W we need more care depending on the power
        switch (Wbar_r_power)
        {
        case 0:   // Wbar = W
          W_in[ind]        = Wbar_in[ind];
          break;
        
        case 1:   // Wbar = r * W
          /*
          dWbar/dr = W + r * dW/dr
                   = W + 0 * 0     at r=0
          
          dWbar/dr = dWbar/dX * dX/dr
          dX/dr = C/(C+r)^2 = 1/C  at r=0
          */
          W_in[ind]        = Wbar_X / C0;
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
          W_in[ind]        = (0.5 * Wbar_XX - Wbar_X)/(C0*C0);
          break;
        
        default:  // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }
      }

      // We need to be careful at X == 1 (r == infty) for radial derivatives (coordinate change is singular)
      else if (i == NX - 1) {
        // W -> 0 for r -> infty
        W_in[ind]        = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction paper) also gives:
        dW_dr_in[ind]    = 0.; 
        dW_dth_in[ind]   = 0.; 

      } else {
        const CCTK_REAL rr = C0*lX/(1. - lX);

        // corresponding derivatives
        // const CCTK_REAL dXdr = 1./(C0 + rr) - rr/((C0 + rr)*(C0 + rr));
        const CCTK_REAL dXdr = C0/((C0 + rr)*(C0 + rr));

        const CCTK_REAL Wbar_r = dXdr * Wbar_X;
        
        // Now translate from Wbar to W
        switch (Wbar_r_power) // We could put a generic power for the computation here I guess...
        {
        case 0:   // Wbar = W
          W_in[ind]        = Wbar_in[ind];
          dW_dr_in[ind]    = Wbar_r;
          dW_dth_in[ind]   = Wbar_th;
          break;
        
        case 1:   // Wbar = r * W
          W_in[ind]        = Wbar_in[ind] / rr;
          dW_dr_in[ind]    = (Wbar_r - W_in[ind]) / rr; // dW/dr  =  1/r * dWbar/dr - Wbar / r^2  =  (dWbar/dr - W) / r
          dW_dth_in[ind]   = Wbar_th / rr;
          break;
        
        case 2: ; // Wbar = r^2 * W
          // empty statement after case to prevent compilation error on some gcc versions...
          const CCTK_REAL rr2_2 = rr*rr;
          W_in[ind]        = Wbar_in[ind] / rr2_2;
          dW_dr_in[ind]    = Wbar_r / rr2_2 - 2 * W_in[ind] / rr; // dW/dr  =  1/r^2 * dWbar/dr - 2 * Wbar / r^3  =  1/r^2 * dWbar/dr - 2 * W / r
          dW_dth_in[ind]   = Wbar_th / rr2_2;
          break;
        
        default:  // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }
      } // if/else i==...
    
    } // for i
  } // for jj


  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *X_g_1, *theta_g_1;
  X_g_1     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  theta_g_1 = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  // CCTK_REAL *X_g_2, *theta_g_2;
  // X_g_2     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // theta_g_2 = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1_1  = x[ind] - x0;
        const CCTK_REAL y1_1  = y[ind] - y0;
        const CCTK_REAL z1_1  = z[ind] - z0;

        const CCTK_REAL rr2_1 = x1_1*x1_1 + y1_1*y1_1 + z1_1*z1_1;

        CCTK_REAL rr_1  = sqrt(rr2_1);
        /* For the Boson Star, x, r and R coordinates coincide (rH=0). */
        
        // From r to the X radial coordinate (used in input files)
        const CCTK_REAL lX_1 = rr_1 / (C0 + rr_1);

        CCTK_REAL ltheta_1 = rr_1 < 1e-16 ? 0 : acos( z1_1/rr_1 );    // There should be at most one point in the grid with rr~0. Not sure about the threshold.
        if (ltheta_1 > 0.5*M_PI)    // symmetry along the equatorial plane
          ltheta_1 = M_PI - ltheta_1;

        X_g_1[ind]     = lX_1;
        theta_g_1[ind] = ltheta_1;
        // const CCTK_REAL x1_2  = x[ind] - x0_2;
        // const CCTK_REAL y1_2  = y[ind] - y0_2;
        // const CCTK_REAL z1_2  = z[ind] - z0_2;

        // const CCTK_REAL rr2_2 = x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;

        // CCTK_REAL rr_2  = sqrt(rr2_2);
        // /* For the Boson Star, x, r and R coordinates coincide (rH=0). */
        
        // // From r to the X radial coordinate (used in input files)
        // const CCTK_REAL lX_2 = rr_2 / (C0 + rr_2);

        // CCTK_REAL ltheta_2 = rr_2 < 1e-16 ? 0 : acos( z1_2/rr_2);    // There should be at most one point in the grid with rr~0. Not sure about the threshold.
        // if (ltheta_2 > 0.5*M_PI)    // symmetry along the equatorial plane
        //   ltheta_2 = M_PI - ltheta_2;

        // X_g_2[ind]     = lX_2;
        // theta_g_2[ind] = ltheta_2;
      }
    }
  }

  /* now for the interpolation */

  const CCTK_INT N_dims  = 2;   // 2-D interpolation

  const CCTK_INT N_input_arrays  = 7;
  const CCTK_INT N_output_arrays = 7;

  /* origin and stride of the input coordinates. with this Cactus reconstructs
     the whole X and theta array. */
  CCTK_REAL origin[N_dims];
  CCTK_REAL delta [N_dims];
  origin[0] = X[0];  origin[1] = theta[0];
  delta[0]  = dX;    delta[1]  = dtheta;

  /* points onto which we want to interpolate, ie, the grid points themselves in
     (X, theta) coordinates (computed above) */
  const void *interp_coords_1[N_dims];
  interp_coords_1[0] = (const void *) X_g_1;
  interp_coords_1[1] = (const void *) theta_g_1;

  // const void *interp_coords_2[N_dims];
  // interp_coords_2[0] = (const void *) X_g_2;
  // interp_coords_2[1] = (const void *) theta_g_2;


  /* input arrays */
  const void *input_arrays[N_input_arrays];
  CCTK_INT input_array_type_codes[N_input_arrays];
  CCTK_INT input_array_dims[N_dims];
  input_array_dims[0] = NX;
  input_array_dims[1] = Ntheta;

  input_array_type_codes[0] = CCTK_VARIABLE_REAL;
  input_array_type_codes[1] = CCTK_VARIABLE_REAL;
  input_array_type_codes[2] = CCTK_VARIABLE_REAL;
  input_array_type_codes[3] = CCTK_VARIABLE_REAL;
  input_array_type_codes[4] = CCTK_VARIABLE_REAL;
  input_array_type_codes[5] = CCTK_VARIABLE_REAL;
  input_array_type_codes[6] = CCTK_VARIABLE_REAL;

  /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
  input_arrays[0] = (const void *) F1_in;
  input_arrays[1] = (const void *) F2_in;
  input_arrays[2] = (const void *) F0_in;
  input_arrays[3] = (const void *) phi0_in;
  input_arrays[4] = (const void *) W_in;
  input_arrays[5] = (const void *) dW_dr_in;
  input_arrays[6] = (const void *) dW_dth_in;

  /* output arrays */
  void *output_arrays_1[N_output_arrays];
  CCTK_INT output_array_type_codes_1[N_output_arrays];
  void *output_arrays_2[N_output_arrays];
  CCTK_INT output_array_type_codes_2[N_output_arrays];
  CCTK_REAL *F1_1, *F2_1, *F0_1, *phi0_1, *W_1;
  CCTK_REAL *dW_dr_1, *dW_dth_1;

  F1_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F2_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F0_1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  phi0_1        = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  W_1           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dr_1       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dth_1      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  
  // F1_2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // F2_2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // F0_2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // phi0_2        = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // W_2           = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // dW_dr_2       = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // dW_dth_2      = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  output_array_type_codes_1[0] = CCTK_VARIABLE_REAL;
  output_array_type_codes_1[1] = CCTK_VARIABLE_REAL;
  output_array_type_codes_1[2] = CCTK_VARIABLE_REAL;
  output_array_type_codes_1[3] = CCTK_VARIABLE_REAL;
  output_array_type_codes_1[4] = CCTK_VARIABLE_REAL;
  output_array_type_codes_1[5] = CCTK_VARIABLE_REAL;
  output_array_type_codes_1[6] = CCTK_VARIABLE_REAL;

  // output_array_type_codes_2[0] = CCTK_VARIABLE_REAL;
  // output_array_type_codes_2[1] = CCTK_VARIABLE_REAL;   // nao posso tirar este senao da interpolation error.
  // output_array_type_codes_2[2] = CCTK_VARIABLE_REAL;
  // output_array_type_codes_2[3] = CCTK_VARIABLE_REAL;
  // output_array_type_codes_2[4] = CCTK_VARIABLE_REAL;
  // output_array_type_codes_2[5] = CCTK_VARIABLE_REAL;
  // output_array_type_codes_2[6] = CCTK_VARIABLE_REAL;

  output_arrays_1[0] = (void *) F1_1;
  output_arrays_1[1] = (void *) F2_1;
  output_arrays_1[2] = (void *) F0_1;
  output_arrays_1[3] = (void *) phi0_1;
  output_arrays_1[4] = (void *) W_1;
  output_arrays_1[5] = (void *) dW_dr_1;
  output_arrays_1[6] = (void *) dW_dth_1;

  // output_arrays_2[0] = (void *) F1_2;
  // output_arrays_2[1] = (void *) F2_2;
  // output_arrays_2[2] = (void *) F0_2;
  // output_arrays_2[3] = (void *) phi0_2;   // sem isto da segmentation fault, mas assim ele põe os dados do ficheiro.
  // output_arrays_2[4] = (void *) W_2;
  // output_arrays_2[5] = (void *) dW_dr_2;
  // output_arrays_2[6] = (void *) dW_dth_2;







  /* handle and settings for the interpolation routine */
  int operator_handle, param_table_handle;
  operator_handle    = CCTK_InterpHandle("Lagrange polynomial interpolation");
  param_table_handle = Util_TableCreateFromString("order=4 boundary_extrapolation_tolerance={0.1 1.0 0.05 0.05}");

  CCTK_INFO("Interpolating result...");

  /* do the actual interpolation, and check for error returns */
  int status_1 = CCTK_InterpLocalUniform(N_dims, operator_handle,
                                       param_table_handle,
                                       origin, delta,
                                       N_interp_points,
                                       CCTK_VARIABLE_REAL,
                                       interp_coords_1,
                                       N_input_arrays, input_array_dims,
                                       input_array_type_codes,
                                       input_arrays,
                                       N_output_arrays, output_array_type_codes_1,
                                       output_arrays_1);
  if (status_1 < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
    "interpolation screwed up!");
  }

  //   /* do the actual interpolation, and check for error returns */
  // int status_2 = CCTK_InterpLocalUniform(N_dims, operator_handle,
  //                                      param_table_handle,
  //                                      origin, delta,
  //                                      N_interp_points,
  //                                      CCTK_VARIABLE_REAL,
  //                                      interp_coords_2,
  //                                      N_input_arrays, input_array_dims,
  //                                      input_array_type_codes,
  //                                      input_arrays,
  //                                      N_output_arrays, output_array_type_codes_2,
  //                                      output_arrays_2);
  // if (status_2 < 0) {
  //   CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
  //   "interpolation screwed up!");
  // }

  free(X_g_1); free(theta_g_1);
  free(Xtmp); free(thtmp);
  free(F1_in); free(F2_in); free(F0_in); free(phi0_in); free(Wbar_in);
  free(W_in); free(dW_dr_in); free(dW_dth_in);


  /* printf("F1 = %g\n", F1[0]); */
  /* printf("F2 = %g\n", F2[0]); */
  /* printf("F0 = %g\n", F0[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W = %g\n", W[0]); */


  /* now we finally write the metric and all 3+1 quantities. first we write the
     3-metric, lapse and scalar fields */
  /* For the Boson Star, in order to avoid unneeded regularizations and divisions,
      we express K_ij in terms of dW/drho and dW/dz.
      For points close to the axis and the origin, K_ij = 0.
  */

  const CCTK_REAL tt = cctk_time;

  const CCTK_REAL coswt = cos(omega_BS * tt);
  const CCTK_REAL sinwt = sin(omega_BS * tt);

  const CCTK_REAL bh_spin2= bh_spin*bh_spin;
  const CCTK_REAL bh_mass2= bh_mass*bh_mass;

  const CCTK_REAL rBLp  = bh_mass + sqrt( bh_mass2 - bh_spin2 );
  const CCTK_REAL rBLm  = bh_mass - sqrt( bh_mass2 - bh_spin2 );

  const CCTK_REAL horizon_radius = 0.5*sqrt(bh_mass2-bh_spin2);

  printf("cctk_lsh[0] = %d\n",cctk_lsh[0]);
  printf("cctk_lsh[1] = %d\n",cctk_lsh[1]);
  printf("cctk_lsh[2] = %d\n",cctk_lsh[2]);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        //Boson Star A

        const CCTK_REAL x1_1  = x[ind] - x0;
        const CCTK_REAL y1_1  = y[ind] - y0;
        const CCTK_REAL z1_1  = z[ind] - z0;

        // For the Boson Star, r = R, no coordinate change needed.
        const CCTK_REAL rr2_1 = x1_1*x1_1 + y1_1*y1_1 + z1_1*z1_1;
        const CCTK_REAL rr_1  = sqrt(rr2_1);

        const CCTK_REAL rho2_1 = x1_1*x1_1 + y1_1*y1_1;
        const CCTK_REAL rho_1  = sqrt(rho2_1);
        

        
        const CCTK_REAL ph_1 = atan2(y1_1, x1_1);
        // If x1_2=y1_2=0, should return 0? The other metric functions should vanish anyway to make sure that this doesn't matter,
        // but can this lead to nan depending on the C implementation?

        const CCTK_REAL cosph_1  = cos(ph_1);
        const CCTK_REAL sinph_1  = sin(ph_1);

        const CCTK_REAL cosmph_1 = cos(mm*ph_1);
        const CCTK_REAL sinmph_1 = sin(mm*ph_1);

        const CCTK_REAL h_rho2_1 = exp(2. * (F2_1[ind] - F1_1[ind])) - 1.;

        //Black Hole B

        const CCTK_REAL x1_2  = x[ind] - x0_2;
        const CCTK_REAL y1_2  = y[ind] - y0_2;
        const CCTK_REAL z1_2  = z[ind] - z0_2;

        const CCTK_REAL bh_v2 = bh_v * bh_v;
        const CCTK_REAL gamma2 = 1. / (1. - bh_v2);
        const CCTK_REAL gamma = sqrt(gamma2);
        // printf("gamma = %e\n", gamma);
 
        CCTK_REAL rr2_2 = x1_2*x1_2*gamma2 + y1_2*y1_2 + z1_2*z1_2;
        if( rr2_2 < pow( eps_r, 2 ) ) {
        rr2_2 = pow( eps_r, 2 );
        }
        const CCTK_REAL rr_2  = sqrt(rr2_2);


        CCTK_REAL rho2_2 = x1_2*x1_2*gamma2 + y1_2*y1_2;
        if( rho2_2 < pow( eps_r, 2 ) ){
        rho2_2 = pow( eps_r, 2 );
        }
        const CCTK_REAL rho_2  = sqrt(rho2_2);
        

        const CCTK_REAL theta_2 = acos(z1_2/rr_2);

        const CCTK_REAL deltakerr2_2 = bh_mass*bh_mass - bh_spin2 ;
        const CCTK_REAL deltakerr  = sqrt(deltakerr2_2) ;

        const CCTK_REAL costh  = z1_2/rr_2 ;
        const CCTK_REAL costh2 = costh*costh ;
        const CCTK_REAL sinth2 = 1. - costh2 ;
        const CCTK_REAL sinth  = sqrt(sinth2) ;

        const CCTK_REAL R_x    = x1_2*gamma/rr_2 ;
        const CCTK_REAL R_y    = y1_2/rr_2 ;
        const CCTK_REAL R_z    = z1_2/rr_2 ;

        // // const CCTK_REAL x_R    = gamma*x1_2/rr_2 ;
        // const CCTK_REAL x_R    = x1_2*gamma/rr_2 ;
        // const CCTK_REAL y_R    = y1_2/rr_2 ;
        // const CCTK_REAL z_R    = z1_2/rr_2 ;

        // const CCTK_REAL sinth2ph_x = -y1_2/rr2_2 ;
        // // const CCTK_REAL sinth2ph_y =  gamma*x1_2/rr2_2 ;
        // const CCTK_REAL sinth2ph_y =  x1_2*gamma/rr2_2 ;


        // // const CCTK_REAL sinthth_x  = z1_2*gamma*x1_2/(rr_2*rr2_2) ;
        // const CCTK_REAL sinthth_x  = z1_2*x1_2/(rr_2*rr2_2) ; 
        // const CCTK_REAL sinthth_y  = z1_2*y1_2/(rr_2*rr2_2) ; 
        // const CCTK_REAL sinthth_z  = -sinth2/rr_2 ; 

        // // const CCTK_REAL sinthx_th  = gamma*x1_2 * costh ;
        // const CCTK_REAL sinthx_th  = x1_2*gamma * costh ;
        // const CCTK_REAL sinthy_th  = y1_2 * costh ;
        // const CCTK_REAL sinthz_th  = -rr_2 * sinth2 ;


        const CCTK_REAL costh2_x = -2*x1_2*gamma*pow(z1_2,2)/pow(rr2_2,2);
        const CCTK_REAL costh2_y = -2*y1_2*pow(z1_2,2)/pow(rr2_2,2);
        const CCTK_REAL costh2_z = 2*rho2_2*z1_2/pow(rr2_2,2);


        const CCTK_REAL rBL    = rr_2 + bh_mass + 0.25*(bh_mass2-bh_spin2) / rr_2 ;   // Boyer-Lindquist coordinate r

        // const CCTK_REAL RRrBL  = rr2_2 + rr_2*bh_mass + 0.25*(bh_mass2-bh_spin2);

        CCTK_REAL rho2kerr   = rBL*rBL + bh_spin2 * costh2 ;
        // if( rho2kerr < pow( eps_r, 2 ) ) {
        //   rho2kerr = pow( eps_r, 2 );
        // }
        // const CCTK_REAL rhokerr    = sqrt(rho2kerr) ;

        const CCTK_REAL sigma  = (2.*bh_mass*rBL)/rho2kerr;
        const CCTK_REAL hh     = (1 + sigma) / (rho2kerr * rr2_2);

        // NaN/Inf checks for sigma and hh
        check_nan_or_inf("sigma", sigma);
        check_nan_or_inf("hh", hh);

        const CCTK_REAL psi4_2 = rho2kerr / rr2_2 ;
        const CCTK_REAL psi2_2 = sqrt(psi4_2) ;
        const CCTK_REAL psi1_2 = sqrt(psi2_2) ;
        const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        const CCTK_REAL psi2_1 = sqrt(psi4_1);
        const CCTK_REAL psi1_1 = sqrt(psi2_1);

        // non-axisymmetric perturbation.
        /* pert = 1. + AA * (x1_2*x1_2 - y1_2*y1_2)/(bh_mass*bh_mass) * exp( -2.*rr2_2/deltakerr2_2 ) ; */
        
        const CCTK_REAL alpha0  = (rr_2 + horizon_radius)*(rr_2 - horizon_radius) / rr_2 * \
                 1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2)) ;
        const CCTK_REAL alpha02 = alpha0*alpha0;


        //fazer as derivadas das funcoes auxiliares

        const CCTK_REAL dr_dR = 1 + (bh_spin2 - bh_mass*bh_mass)/(4*rr2_2);
        
        const CCTK_REAL drho2kerr_dx = 2*rBL*dr_dR*R_x + bh_spin2*costh2_x;
        const CCTK_REAL drho2kerr_dy = 2*rBL*dr_dR*R_y + bh_spin2*costh2_y;
        const CCTK_REAL drho2kerr_dz = 2*rBL*dr_dR*R_z + bh_spin2*costh2_z;

        const CCTK_REAL dsigma_dx = (2*bh_mass*(rho2kerr*R_x*dr_dR - \
                                 rBL*drho2kerr_dx))/pow(rho2kerr,2);
        const CCTK_REAL dsigma_dy = (2*bh_mass*(rho2kerr*R_y*dr_dR - \
                                 rBL*drho2kerr_dy))/pow(rho2kerr,2);
        const CCTK_REAL dsigma_dz = (2*bh_mass*(rho2kerr*R_z*dr_dR - \
                                 rBL*drho2kerr_dz))/pow(rho2kerr,2);


        const CCTK_REAL dpsi4_2_dx = (-2*rho2kerr*R_x + rr_2*drho2kerr_dx)/pow(rr_2,3);
        const CCTK_REAL dpsi4_2_dy = (-2*rho2kerr*R_y + rr_2*drho2kerr_dy)/pow(rr_2,3);
        const CCTK_REAL dpsi4_2_dz = (-2*rho2kerr*R_z + rr_2*drho2kerr_dz)/pow(rr_2,3);



        const CCTK_REAL bphi = -bh_spin * sigma * sinth2;
        // const CCTK_REAL gammaphiphi = psi4_2 * rr2_2 * sinth2 * (1 + bh_spin2 * hh * rr2_2 * sinth2);
        const CCTK_REAL gammaphiphi = rho2kerr*rho2_2/rr2_2 *(1 + bh_spin2 * hh * rho2_2);
        const CCTK_REAL bphiup = bphi / gammaphiphi;

        // Check for NaN or Inf in these quantities
        check_nan_or_inf("bphi", bphi);
        check_nan_or_inf("gammaphiphi", gammaphiphi);
        check_nan_or_inf("bphiup", bphiup);



        CCTK_REAL betad[4] = {0.0, 0.0, 0.0, 0.0};
        betad[1] = -y1_2/rho2_2 * bphi;
        betad[2] =  x1_2*gamma/rho2_2 * bphi;
        betad[3] = 0.;

        // Check for NaN or Inf in betad components
        check_nan_or_inf("betad[1]", betad[1]);
        check_nan_or_inf("betad[2]", betad[2]);
        check_nan_or_inf("betad[3]", betad[3]);

        //To change the spin direction, change the indices accordingly.
        const CCTK_REAL dbetadphi_dx = (-2*bh_spin*pow(z1_2,2)*sigma*R_x)/pow(rr_2,3) - bh_spin*(1 - pow(z1_2,2)/pow(rr_2,2))*dsigma_dx;
        const CCTK_REAL dbetadphi_dy = (-2*bh_spin*pow(z1_2,2)*sigma*R_y)/pow(rr_2,3) - bh_spin*(1 - pow(z1_2,2)/pow(rr_2,2))*dsigma_dy;
        const CCTK_REAL dbetadphi_dz = (bh_spin*(2*z1_2*sigma*(rr_2 - z1_2*R_z) + rr_2*(-rho2_2)*dsigma_dz))/pow(rr_2,3);


        const CCTK_REAL dhh_dx = (-((1 + sigma)*(2*rho2kerr*R_x + \
                                 rr_2*drho2kerr_dx))+rr_2*rho2kerr*dsigma_dx)/(pow(rr_2,3)*pow(rho2kerr,2));
        const CCTK_REAL dhh_dy = (-((1 + sigma)*(2*rho2kerr*R_y + \
                                 rr_2*drho2kerr_dy))+rr_2*rho2kerr*dsigma_dy)/(pow(rr_2,3)*pow(rho2kerr,2));
        const CCTK_REAL dhh_dz = (-((1 + sigma)*(2*rho2kerr*R_z + \
                                 rr_2*drho2kerr_dz))+rr_2*rho2kerr*dsigma_dz)/(pow(rr_2,3)*pow(rho2kerr,2));

        
        const CCTK_REAL dalpha_dx = (2*pow(rr_2,2)*(pow(rBL,2) + \
                                    pow(bh_spin,2)*(1 - (-1 + costh2)*sigma))*R_x + \
                                    2*pow(horizon_radius,2)*(pow(rBL,2) + \
                                    pow(bh_spin,2)*(1 - (-1 + costh2)*sigma))*R_x + \
                                    rr_2*pow(horizon_radius,2)*(2*rBL*R_x*dr_dR + \
                                    pow(bh_spin,2)*(-(sigma*costh2_x) - (-1 + costh2)*dsigma_dx)) + \
                                    pow(rr_2,3)*(-2*rBL*R_x*dr_dR + pow(bh_spin,2)*(sigma*costh2_x + (-1 \
                                    + costh2)*dsigma_dx)))/(2.*pow(rr_2,2)*pow(pow(rBL,2) + pow(bh_spin,2)*(1 - (-1 +costh2)*sigma),1.5));
        const CCTK_REAL dalpha_dy = (2*pow(rr_2,2)*(pow(rBL,2) + pow(bh_spin,2)*(1 - (-1 + \
                                    costh2)*sigma))*R_y + 2*pow(horizon_radius,2)*(pow(rBL,2) + \
                                    pow(bh_spin,2)*(1 - (-1 + costh2)*sigma))*R_y + \
                                    rr_2*pow(horizon_radius,2)*(2*rBL*R_y*dr_dR + \
                                    pow(bh_spin,2)*(-(sigma*costh2_y) - (-1 + costh2)*dsigma_dy)) + \
                                    pow(rr_2,3)*(-2*rBL*R_y*dr_dR+ pow(bh_spin,2)*(sigma*costh2_y + (-1 \
                                    + costh2)*dsigma_dy)))/(2.*pow(rr_2,2)*pow(pow(rBL,2) + \
                                    pow(bh_spin,2)*(1 - (-1 + costh2)*sigma),1.5));
        const CCTK_REAL dalpha_dz = (2*pow(rr_2,2)*(pow(rBL,2) + pow(bh_spin,2)*(1 - (-1 + \
                                    costh2)*sigma))*R_z + 2*pow(horizon_radius,2)*(pow(rBL,2) + \
                                    pow(bh_spin,2)*(1 - (-1 + costh2)*sigma))*R_z + \
                                    rr_2*pow(horizon_radius,2)*(2*rBL*R_z*dr_dR + \
                                    pow(bh_spin,2)*(-(sigma*costh2_z) - (-1 + costh2)*dsigma_dz)) + \
                                    pow(rr_2,3)*(-2*rBL*R_z*dr_dR+ pow(bh_spin,2)*(sigma*costh2_z + (-1 \
                                    + costh2)*dsigma_dz)))/(2.*pow(rr_2,2)*pow(pow(rBL,2) + \
                                    pow(bh_spin,2)*(1 - (-1 + costh2)*sigma),1.5));

                                    
        // Check for NaN or Inf in all these quantities
        check_nan_or_inf("dr_dR", dr_dR);
        check_nan_or_inf("costh2_x", costh2_x);
        check_nan_or_inf("costh2_y", costh2_y);
        check_nan_or_inf("costh2_z", costh2_z);

        check_nan_or_inf("dsigma_dx", dsigma_dx);
        check_nan_or_inf("dsigma_dy", dsigma_dy);
        check_nan_or_inf("dsigma_dz", dsigma_dz);

        check_nan_or_inf("dpsi4_2_dx", dpsi4_2_dx);
        check_nan_or_inf("dpsi4_2_dy", dpsi4_2_dy);
        check_nan_or_inf("dpsi4_2_dz", dpsi4_2_dz);

        check_nan_or_inf("dhh_dx", dhh_dx);
        check_nan_or_inf("dhh_dy", dhh_dy);
        check_nan_or_inf("dhh_dz", dhh_dz);

        check_nan_or_inf("dalpha_dx", dalpha_dx);
        check_nan_or_inf("dalpha_dy", dalpha_dy);
        check_nan_or_inf("dalpha_dz", dalpha_dz);


        check_nan_or_inf("drho2kerr_dx", drho2kerr_dx);
        check_nan_or_inf("drho2kerr_dy", drho2kerr_dy);
        check_nan_or_inf("drho2kerr_dz", drho2kerr_dz);

        check_nan_or_inf("dbetadphi_dx", drho2kerr_dx);
        check_nan_or_inf("dbetadphi_dy", drho2kerr_dy);
        check_nan_or_inf("dbetadphi_dz", drho2kerr_dz);


        CCTK_REAL dbetad[4][4];
        // Initialize g to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            dbetad[i][j] = 0.0; 
        // Compute derivatives of the beta vector. To change the spin direction, change the indices accordingly.
        dbetad[1][1] = (2*x1_2*gamma*y1_2*bphi - y1_2*rho2_2*dbetadphi_dx)/pow(rho2_2,2);

        dbetad[1][2] = ((-pow(x1_2,2)*gamma2 + pow(y1_2,2))*bphi - y1_2*(rho2_2)*dbetadphi_dy)/pow(rho2_2,2);

        dbetad[1][3] = -((y1_2*dbetadphi_dz)/(rho2_2));

        dbetad[2][1] = ((-pow(x1_2,2)*gamma2 + pow(y1_2,2))*bphi + x1_2*rho2_2*dbetadphi_dx)/pow(rho2_2,2);

        dbetad[2][2] = (x1_2*gamma*(-2*y1_2*bphi + rho2_2*dbetadphi_dy))/pow(rho2_2,2);
        dbetad[2][3] = (x1_2*gamma*dbetadphi_dz)/(rho2_2);
        dbetad[3][1] = 0;
        dbetad[3][2] = 0;
        dbetad[3][3] = 0;
        
        // Check for NaN or Inf in all dbetad[i][j] components
        for (int ii = 1; ii <= 3; ++ii) {
          for (int jj = 1; jj <= 3; ++jj) {
            char dbetad_name[32];
            snprintf(dbetad_name, sizeof(dbetad_name), "dbetad[%d][%d]", ii, jj);
            check_nan_or_inf(dbetad_name, dbetad[ii][jj]);
          }
        }


        //capital Gs refer to the unboosted frame.

        const CCTK_REAL Gtt = -alpha02 + bphi*bphiup;
        const CCTK_REAL Gxt = betad[1];
        const CCTK_REAL Gxx = psi4_2*(1+bh_spin2*hh*y1_2*y1_2);
        const CCTK_REAL Gxy = -psi4_2*bh_spin2*hh*y1_2*gamma*x1_2;
        const CCTK_REAL Gty = betad[2]; 
        

       

        gxx[ind] = gamma2*Gxx + 2*gamma2*bh_v*Gxt + gamma2*bh_v2*Gtt;
        gxy[ind] = gamma*Gxy+gamma*bh_v*Gty;
        gxz[ind] = 0;
        // gyy[ind] = (1 + bh_spin2*pow(x1_2,2)*gamma2*hh)*psi4_2; //+ psi4_1* (1. + h_rho2_1 * cosph_1 * cosph_1) - 1;
        gyy[ind] = psi4_2 * ( 1. + bh_spin2 * hh * gamma2*x1_2*x1_2 );
        gyz[ind] = 0;
        gzz[ind] = psi4_2; // + psi4_1 - 1;


        check_nan_or_inf("gxx",gxx[ind]);
        check_nan_or_inf("gxy",gxy[ind]);
        check_nan_or_inf("gxz",gxz[ind]);
        check_nan_or_inf("gyy",gyy[ind]);
        check_nan_or_inf("gyz",gyz[ind]);
        check_nan_or_inf("gzz",gzz[ind]);




         // Create an array g to store the metric components at each grid point
        CCTK_REAL g[4][4];

        // Initialize g to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            g[i][j] = 0.0; 

        g[0][0] = pow(gamma,2)*(-pow(alpha0,2) - 2*bh_v*betad[1] + (pow(betad[1],2) + \
                  pow(betad[2],2) + pow(bh_spin,2)*pow(y1_2*betad[2] + \
                  x1_2*betad[1]*gamma,2)*hh)/((1 + pow(bh_spin,2)*(rho2_2)*hh)*psi4_2) \
                  + pow(bh_v,2)*(1 + pow(bh_spin,2)*pow(y1_2,2)*hh)*psi4_2);
        g[1][1] = gxx[ind];
        g[1][2] = gxy[ind];
        g[1][3] = gxz[ind];
        g[2][1] = g[1][2];
        g[2][2] = gyy[ind];
        g[2][3] = gyz[ind];
        g[3][1] = gxz[ind];
        g[3][2] = gyz[ind];
        g[3][3] = gzz[ind];

        
        // Check for NaN or Inf in all g[i][j] components
        for (int i = 0; i <= 3; ++i) {
          for (int j = 0; j <= 3; ++j) {
            char g_name[32];
            snprintf(g_name, sizeof(g_name), "g[%d][%d]", i, j);
            check_nan_or_inf(g_name, g[i][j]);
          }
        }
        

        //from here on the derivatives already take into account the gammas correctly.
        
        CCTK_REAL dg[4][4][4]; // dg[i][j][k] = \partial_k g_{ij}
        // Initialize dg to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            for (int kk = 0; kk < 4; ++kk)
              dg[ii][jj][kk] = 0.0;

        // Example: dg[1][1][1] = dgxx_dx, dg[1][1][2] = dgxx_dy, etc.
        dg[1][1][1] = (pow(gamma,3)*(2*pow(bh_v,2)*betad[2]*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*psi4_2*(pow(bh_spin,2)*x1_2*y1_2*gamma*hh*\
                      dbetad[1][1] + (1 + pow(bh_spin,2)*pow(y1_2,2)*hh)*dbetad[2][1]) + \
                      pow(1 + pow(bh_spin,2)*(rho2_2)*hh,2)*pow(psi4_2,2)*(-2*bh_v*(bh_v*\
                      alpha0*dalpha_dx + dbetad[1][1]) + \
                      pow(bh_spin,2)*pow(y1_2,2)*psi4_2*dhh_dx + (1 + \
                      pow(bh_spin,2)*pow(y1_2,2)*hh)*dpsi4_2_dx) + \
                      pow(bh_v,2)*pow(betad[2],2)*(-(pow(bh_spin,2)*x1_2*gamma*psi4_2*(2*hh*\
                      (1 + pow(bh_spin,2)*pow(y1_2,2)*hh) + x1_2*gamma*dhh_dx)) - (1 + \
                      pow(bh_spin,2)*pow(y1_2,2)*hh)*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dx) + \
                      pow(bh_v,2)*pow(betad[1],2)*(pow(bh_spin,2)*pow(y1_2,2)*psi4_2*(2*pow(\
                      bh_spin,2)*x1_2*gamma*pow(hh,2) - dhh_dx) - (1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dx) + 2*pow(bh_v,2)*betad[1]*((1 \
                      + pow(bh_spin,2)*(rho2_2)*hh)*psi4_2*((1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*dbetad[1][1] + \
                      pow(bh_spin,2)*x1_2*y1_2*gamma*hh*dbetad[2][1]) + \
                      pow(bh_spin,2)*y1_2*betad[2]*(psi4_2*(hh + \
                      pow(bh_spin,2)*(pow(y1_2,2) - pow(x1_2,2)*pow(gamma,2))*pow(hh,2) + \
                      x1_2*gamma*dhh_dx) - x1_2*gamma*hh*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dx))))/(pow(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh,2)*pow(psi4_2,2)); // dgxx_dx
        dg[1][1][2] = (pow(gamma,2)*(2*pow(bh_v,2)*betad[2]*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*psi4_2*(pow(bh_spin,2)*x1_2*y1_2*gamma*hh*\
                      dbetad[1][2] + (1 + pow(bh_spin,2)*pow(y1_2,2)*hh)*dbetad[2][2]) + \
                      pow(bh_v,2)*pow(betad[2],2)*(pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*\
                      psi4_2*(2*pow(bh_spin,2)*y1_2*pow(hh,2) - dhh_dy) - (1 + \
                      pow(bh_spin,2)*pow(y1_2,2)*hh)*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dy) + \
                      pow(bh_v,2)*pow(betad[1],2)*(-(pow(bh_spin,2)*y1_2*psi4_2*(2*hh*(1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh) + y1_2*dhh_dy)) - (1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dy) + 2*pow(bh_v,2)*betad[1]*((1 \
                      + pow(bh_spin,2)*(rho2_2)*hh)*psi4_2*((1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*dbetad[1][2] + \
                      pow(bh_spin,2)*x1_2*y1_2*gamma*hh*dbetad[2][2]) + \
                      pow(bh_spin,2)*x1_2*betad[2]*gamma*(psi4_2*(hh + \
                      pow(bh_spin,2)*(-rho2_2)*pow(hh,2) + y1_2*dhh_dy) - y1_2*hh*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dy)) + pow(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh,2)*pow(psi4_2,2)*(-2*bh_v*(bh_v*alpha0*\
                      dalpha_dy + dbetad[1][2]) + dpsi4_2_dy + \
                      pow(bh_spin,2)*y1_2*(y1_2*psi4_2*dhh_dy + hh*(2*psi4_2 + \
                      y1_2*dpsi4_2_dy)))))/(pow(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh,2)*pow(psi4_2,2)); // dgxx_dy
        dg[1][1][3] = (pow(gamma,2)*(-2*pow(bh_v,2)*alpha0*pow(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh,2)*pow(psi4_2,2)*dalpha_dz + \
                      2*pow(bh_spin,2)*pow(bh_v,2)*x1_2*y1_2*betad[2]*gamma*hh*psi4_2*\
                      dbetad[1][3] + \
                      2*pow(bh_spin,4)*pow(bh_v,2)*x1_2*pow(y1_2,3)*betad[2]*gamma*pow(hh,2)\
                      *psi4_2*dbetad[1][3] + \
                      2*pow(bh_spin,4)*pow(bh_v,2)*pow(x1_2,3)*y1_2*betad[2]*pow(gamma,3)*\
                      pow(hh,2)*psi4_2*dbetad[1][3] - 2*bh_v*pow(psi4_2,2)*dbetad[1][3] - \
                      4*pow(bh_spin,2)*bh_v*pow(y1_2,2)*hh*pow(psi4_2,2)*dbetad[1][3] - \
                      4*pow(bh_spin,2)*bh_v*pow(x1_2,2)*pow(gamma,2)*hh*pow(psi4_2,2)*\
                      dbetad[1][3] - \
                      2*pow(bh_spin,4)*bh_v*pow(y1_2,4)*pow(hh,2)*pow(psi4_2,2)*dbetad[1][3]\
                       - 4*pow(bh_spin,4)*bh_v*pow(x1_2,2)*pow(y1_2,2)*pow(gamma,2)*pow(hh,\
                      2)*pow(psi4_2,2)*dbetad[1][3] - \
                      2*pow(bh_spin,4)*bh_v*pow(x1_2,4)*pow(gamma,4)*pow(hh,2)*pow(psi4_2,2)\
                      *dbetad[1][3] + 2*pow(bh_v,2)*betad[2]*psi4_2*dbetad[2][3] + \
                      4*pow(bh_spin,2)*pow(bh_v,2)*pow(y1_2,2)*betad[2]*hh*psi4_2*dbetad[2][\
                      3] + 2*pow(bh_spin,2)*pow(bh_v,2)*pow(x1_2,2)*betad[2]*pow(gamma,2)*\
                      hh*psi4_2*dbetad[2][3] + \
                      2*pow(bh_spin,4)*pow(bh_v,2)*pow(y1_2,4)*betad[2]*pow(hh,2)*psi4_2*\
                      dbetad[2][3] + \
                      2*pow(bh_spin,4)*pow(bh_v,2)*pow(x1_2,2)*pow(y1_2,2)*betad[2]*pow(\
                      gamma,2)*pow(hh,2)*psi4_2*dbetad[2][3] - \
                      pow(bh_spin,2)*pow(bh_v,2)*pow(x1_2,2)*pow(betad[2],2)*pow(gamma,2)*\
                      psi4_2*dhh_dz + pow(bh_spin,2)*pow(y1_2,2)*pow(psi4_2,3)*dhh_dz + \
                      2*pow(bh_spin,4)*pow(y1_2,4)*hh*pow(psi4_2,3)*dhh_dz + \
                      2*pow(bh_spin,4)*pow(x1_2,2)*pow(y1_2,2)*pow(gamma,2)*hh*pow(psi4_2,3)\
                      *dhh_dz + pow(bh_spin,6)*pow(y1_2,6)*pow(hh,2)*pow(psi4_2,3)*dhh_dz + \
                      2*pow(bh_spin,6)*pow(x1_2,2)*pow(y1_2,4)*pow(gamma,2)*pow(hh,2)*pow(\
                      psi4_2,3)*dhh_dz + \
                      pow(bh_spin,6)*pow(x1_2,4)*pow(y1_2,2)*pow(gamma,4)*pow(hh,2)*pow(\
                      psi4_2,3)*dhh_dz + (1 + pow(bh_spin,2)*pow(y1_2,2)*hh)*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*(-(pow(bh_v,2)*pow(betad[2],2)) + (1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*pow(psi4_2,2))*dpsi4_2_dz + \
                      2*pow(bh_v,2)*betad[1]*(psi4_2*((1 + pow(bh_spin,2)*(rho2_2)*hh)*((1 \
                      + pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*dbetad[1][3] + \
                      pow(bh_spin,2)*x1_2*y1_2*gamma*hh*dbetad[2][3]) + \
                      pow(bh_spin,2)*x1_2*y1_2*betad[2]*gamma*dhh_dz) - \
                      pow(bh_spin,2)*x1_2*y1_2*betad[2]*gamma*hh*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dz) + \
                      pow(bh_v,2)*pow(betad[1],2)*(-(pow(bh_spin,2)*pow(y1_2,2)*psi4_2*dhh_dz) - (1 + pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh)*dpsi4_2_dz)))/(pow(1 + \
                      pow(bh_spin,2)*(rho2_2)*hh,2)*pow(psi4_2,2)); // dgxx_dz
        dg[1][2][1] = -(pow(gamma,2)*(bh_v*dbetad[2][1] + \
                      pow(bh_spin,2)*y1_2*(x1_2*gamma*psi4_2*dhh_dx + hh*(psi4_2 + x1_2*gamma*dpsi4_2_dx)))); // dgxy_dx
        dg[1][2][2] = -(gamma*(bh_v*dbetad[2][2] + \
                      pow(bh_spin,2)*x1_2*gamma*(y1_2*psi4_2*dhh_dy + hh*(psi4_2 + y1_2*dpsi4_2_dy)))); // dgxy_dy
        dg[1][2][3] = -(gamma*(bh_v*dbetad[2][3] + \
                      pow(bh_spin,2)*x1_2*y1_2*gamma*(psi4_2*dhh_dz + hh*dpsi4_2_dz))); // dgxy_dz
        dg[1][3][1] = 0.0; // dgxz_dx
        dg[1][3][2] = 0.0; // dgxz_dy
        dg[1][3][3] = 0.0; // dgxz_dz
        dg[2][2][1] = gamma*(pow(bh_spin,2)*x1_2*gamma*psi4_2*(2*hh + x1_2*gamma*dhh_dx) + \
                      (1 + pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*dpsi4_2_dx); // dgyy_dx
        dg[2][2][2] = pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*psi4_2*dhh_dy + (1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*dpsi4_2_dy; // dgyy_dy
        dg[2][2][3] = pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*psi4_2*dhh_dz + (1 + \
                      pow(bh_spin,2)*pow(x1_2,2)*pow(gamma,2)*hh)*dpsi4_2_dz; // dgyy_dz
        dg[2][3][1] = 0.0; // dgyz_dx
        dg[2][3][2] = 0.0; // dgyz_dy
        dg[2][3][3] = 0.0; // dgyz_dz
        dg[3][3][1] = gamma*dpsi4_2_dx; // dgzz_dx
        dg[3][3][2] = dpsi4_2_dy; // dgzz_dy
        dg[3][3][3] = dpsi4_2_dz; // dgzz_dz
        // Set symmetric components
        dg[2][1][1] = dg[1][2][1]; // dgyx_dx = dgxy_dx
        dg[2][1][2] = dg[1][2][2]; // dgyx_dy = dgxy_dy
        dg[2][1][3] = dg[1][2][3]; // dgyx_dz = dgxy_dz
        dg[3][1][1] = dg[1][3][1]; // dgzx_dx = dgxz_dx
        dg[3][1][2] = dg[1][3][2]; // dgzx_dy = dgxz_dy
        dg[3][1][3] = dg[1][3][3]; // dgzx_dz = dgxz_dz
        dg[3][2][1] = dg[2][3][1]; // dgzy_dx = dgyz_dx
        dg[3][2][2] = dg[2][3][2]; // dgzy_dy = dgyz_dy
        dg[3][2][3] = dg[2][3][3]; // dgzy_dz = dgyz_dz



        
        //time derivatives. nas coordenadas adptadas à foliação do boost a métrica é estacionária.
        dg[1][1][0] = dg[1][1][1]*bh_v*gamma; // ∂g_xx/∂t
        dg[1][2][0] = dg[1][2][1]*bh_v*gamma; // ∂g_xy/∂t
        dg[1][3][0] = dg[1][3][1]*bh_v*gamma; // ∂g_xz/∂t
        dg[2][1][0] = dg[1][2][0]; // ∂g_yx/∂t
        dg[2][2][0] = dg[2][2][1]*bh_v*gamma; // ∂g_yy/∂t
        dg[2][3][0] = dg[2][3][1]*bh_v*gamma; // ∂g_yz/∂t
        dg[3][1][0] = dg[1][3][0]; // ∂g_zx/∂t
        dg[3][2][0] = dg[2][3][0]; // ∂g_zy/∂t
        dg[3][3][0] = dg[3][3][1]*bh_v*gamma; // ∂g_zz/∂t




        // Check for NaN or Inf in all metric derivatives
        for (int ii = 1; ii <= 3; ++ii) {
          for (int jj = 1; jj <= 3; ++jj) {
            for (int kk = 0; kk <= 3; ++kk) {
              char dg_name[32];
              snprintf(dg_name, sizeof(dg_name), "dg[%d][%d][%d]", ii, jj, kk);
              check_nan_or_inf(dg_name, dg[ii][jj][kk]);
            }
          }
        }



        // Compute inverse metric g^{ij} (spatial part only)
        CCTK_REAL det_g =
            g[1][1]*(g[2][2]*g[3][3] - g[2][3]*g[3][2])
          - g[1][2]*(g[2][1]*g[3][3] - g[2][3]*g[3][1])
          + g[1][3]*(g[2][1]*g[3][2] - g[2][2]*g[3][1]);

        CCTK_REAL g_inv[4][4]; // Inverse metric
        // Initialize g_inv to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            g_inv[i][j] = 0.0;

        if (fabs(det_g) < 1e-12) {
            // Abort execution due to singular metric
            CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Singular spatial metric (det_g = %e) at grid point (%d,%d,%d). Aborting.", det_g, i, j, k);
            abort();
        } else {
            g_inv[1][1] =  (g[2][2]*g[3][3] - g[2][3]*g[3][2]) / det_g;
            g_inv[1][2] = -(g[1][2]*g[3][3] - g[1][3]*g[3][2]) / det_g;
            g_inv[1][3] =  (g[1][2]*g[2][3] - g[1][3]*g[2][2]) / det_g;
            g_inv[2][1] = -(g[2][1]*g[3][3] - g[2][3]*g[3][1]) / det_g;
            g_inv[2][2] =  (g[1][1]*g[3][3] - g[1][3]*g[3][1]) / det_g;
            g_inv[2][3] = -(g[1][1]*g[2][3] - g[1][3]*g[2][1]) / det_g;
            g_inv[3][1] =  (g[2][1]*g[3][2] - g[2][2]*g[3][1]) / det_g;
            g_inv[3][2] = -(g[1][1]*g[3][2] - g[1][2]*g[3][1]) / det_g;
            g_inv[3][3] =  (g[1][1]*g[2][2] - g[1][2]*g[2][1]) / det_g;
        }

        // CCTK_REAL betaup[4];
        //           betaup[1] = g_inv[1][1] * betad[1] + g_inv[1][2] * betad[2] + g_inv[1][3] * betad[3];
        //           betaup[2] = g_inv[2][1] * betad[1] + g_inv[2][2] * betad[2] + g_inv[2][3] * betad[3];
        //           betaup[3] = g_inv[3][1] * betad[1] + g_inv[3][2] * betad[2] + g_inv[3][3] * betad[3];

        // // Check for NaN or Inf in betaup components
        // check_nan_or_inf("betaup[1]", betaup[1]);
        // check_nan_or_inf("betaup[2]", betaup[2]);
        // check_nan_or_inf("betaup[3]", betaup[3]);
        //estava mal porque os shifts agora sao diferentes depois do boost.

        
        // Christoffel symbols of the spatial metric (only spatial indices 1..3)
        CCTK_REAL Gamma[4][4][4]; // Gamma^i_{jk}
        // Initialize to zero
        for (int i = 1; i <= 3; ++i)
          for (int j = 1; j <= 3; ++j)
            for (int k = 1; k <= 3; ++k)
              Gamma[i][j][k] = 0.0;


        // Compute Christoffel symbols (spatial part)
        // Gamma^i_{jk} = 0.5 * g^{il} (dg_{lj}/dx^k + dg_{lk}/dx^j - dg_{jk}/dx^l)
        for (int i = 1; i <= 3; ++i) {
          for (int j = 1; j <= 3; ++j) {
            for (int k = 1; k <= 3; ++k) {
              for (int l = 1; l <= 3; ++l) {
                Gamma[i][j][k] += 0.5 * g_inv[i][l] * (dg[l][j][k] + dg[l][k][j] - dg[j][k][l]);
              }
            }
          }
        }
        // Check for NaN or Inf in all Christoffel symbols
        for (int ii = 1; ii <= 3; ++ii) {
          for (int jj = 1; jj <= 3; ++jj) {
            for (int kk = 1; kk <= 3; ++kk) {
              char gamma_name[32];
              snprintf(gamma_name, sizeof(gamma_name), "Gamma[%d][%d][%d]", ii, jj, kk);
              check_nan_or_inf(gamma_name, Gamma[ii][jj][kk]);
            }
          }
        }

        CCTK_REAL new_betad[4];
        // Initialize g to zero
        for (int i = 0; i < 4; ++i)
            new_betad[i] = 0.0; 

        new_betad[1] = pow(gamma,2)*(betad[1] + bh_v*(pow(alpha0,2) + bh_v*betad[1] - \
                       (pow(betad[1],2) + pow(betad[2],2) + pow(bh_spin,2)*pow(y1_2*betad[2] \
                       + x1_2*betad[1]*gamma,2)*hh)/((1 + \
                       pow(bh_spin,2)*(rho2_2)*hh)*psi4_2) - (1 + pow(bh_spin,2)*pow(y1_2,2)*hh)*psi4_2));
        new_betad[2] =  gamma*(betad[2] + pow(bh_spin,2)*bh_v*x1_2*y1_2*gamma*hh*psi4_2);
        new_betad[3] = 0.;
        check_nan_or_inf("new_betad[1]", new_betad[1]);
        check_nan_or_inf("new_betad[2]", new_betad[2]);
        check_nan_or_inf("new_betad[3]", new_betad[3]);


        CCTK_REAL new_betaup[4];
        // Initialize new_betaup to zero
        for (int i = 0; i < 4; ++i)
            new_betaup[i] = 0.0;

        new_betaup[1] = g_inv[1][1] * new_betad[1] + g_inv[1][2] * new_betad[2] + g_inv[1][3] * new_betad[3];
        new_betaup[2] = g_inv[2][1] * new_betad[1] + g_inv[2][2] * new_betad[2] + g_inv[2][3] * new_betad[3];
        new_betaup[3] = g_inv[3][1] * new_betad[1] + g_inv[3][2] * new_betad[2] + g_inv[3][3] * new_betad[3];

        // Check for NaN or Inf in betaup components
        check_nan_or_inf("betaup[1]", new_betaup[1]);
        check_nan_or_inf("betaup[2]", new_betaup[2]);
        check_nan_or_inf("betaup[3]", new_betaup[3]);





        // CCTK_REAL new_lapse = sqrt(-g[0][0] + betad[1]*betaup[1] + betad[2]*betaup[2] + betad[3]*betaup[3]);
        double lapse_arg = -g[0][0] + new_betad[1]*new_betaup[1] + new_betad[2]*new_betaup[2] + new_betad[3]*new_betaup[3];
        if (lapse_arg < 0) {
            fprintf(stderr, "Negative argument in sqrt for new_lapse: %.9e\n", lapse_arg);
            // print more context here
        }
        CCTK_REAL new_lapse = sqrt(lapse_arg);
        if (new_lapse < SMALL){
            new_lapse = SMALL;
        }
        if (isnan(new_lapse)) {
        fprintf(stderr, "Error: %s is NaN\n", "new_lapse");
        fprintf(stderr, "g00 = %.9e \n", g[0][0]);
        fprintf(stderr, "beta2 = %.9e \n", new_betad[1]*new_betaup[1] + new_betad[2]*new_betaup[2] + new_betad[3]*new_betaup[3]);
        fprintf(stderr, "Error: new_lapse is nan at grid point (%d,%d,%d)\n", x[CCTK_GFINDEX3D (cctkGH, i, j, k)], y[CCTK_GFINDEX3D (cctkGH, i, j, k)], z[CCTK_GFINDEX3D (cctkGH, i, j, k)]);
        abort(); // Break execution
        }


    
      

        CCTK_REAL dW_drho_1, dW_dz_1;
        const CCTK_REAL exp_auxi_1 = exp(2. * F2_1[ind] - F0_1[ind]);

        if (rho_1 < 1e-8) {
          dW_drho_1 = 0.;
          dW_dz_1   = 0.;
        }
        else {
          dW_drho_1 = rho_1/rr_1 * dW_dr_1[ind]  +   z1_1/rr2_1 * dW_dth_1[ind];
          dW_dz_1   =  z1_1/rr_1 * dW_dr_1[ind]  -  rho_1/rr2_1 * dW_dth_1[ind];
        }


        // Compute covariant derivatives of the shift
        CCTK_REAL Dbetad[4][4];
        // Initialize Dbetad to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            Dbetad[i][j] = 0.0;
          
        for (int i = 1; i <= 3; ++i) {
          for (int j = 1; j <= 3; ++j) {
            Dbetad[i][j] = dbetad[i][j]; //aqui deve ser a derivada dos shifts depois do boost.
            for (int k = 1; k <= 3; ++k)
              Dbetad[i][j] -= Gamma[k][i][j] * betad[k]; // e aqui também.
            // Check for NaN or Inf in Dbetad
            char Dbetad_name[32];
            snprintf(Dbetad_name, sizeof(Dbetad_name), "Dbetad[%d][%d]", i, j);
            check_nan_or_inf(Dbetad_name, Dbetad[i][j]);
          }
        }




      // Compute extrinsic curvature K_{ij}
      kxx[ind] = 0.5 / new_lapse * (Dbetad[1][1] + Dbetad[1][1] - dg[1][1][0]);
      kxy[ind] = 0.5 / new_lapse * (Dbetad[1][2] + Dbetad[2][1] - dg[1][2][0]);
      kxz[ind] = 0.5 / new_lapse * (Dbetad[1][3] + Dbetad[3][1] - dg[1][3][0]);
      kyy[ind] = 0.5 / new_lapse * (Dbetad[2][2] + Dbetad[2][2] - dg[2][2][0]);
      kyz[ind] = 0.5 / new_lapse * (Dbetad[2][3] + Dbetad[3][2] - dg[2][3][0]);
      kzz[ind] = 0.5 / new_lapse * (Dbetad[3][3] + Dbetad[3][3] - dg[3][3][0]);


     

        
        check_nan_or_inf("kxx",kxx[ind]);
        check_nan_or_inf("kxy",kxy[ind]);
        check_nan_or_inf("kxz",kxz[ind]);
        check_nan_or_inf("kyy",kyy[ind]);
        check_nan_or_inf("kyz",kyz[ind]);
        check_nan_or_inf("kzz",kzz[ind]);  

        // // let's add a perturbation to the scalar field as well
        // const CCTK_REAL argpert_phi_1 = (rr_1 - R0pert_phi)/Sigmapert_phi;
        // const CCTK_REAL pert_phi_1 = 1. + Apert_phi * (x1_1*x1_1 - y1_1*y1_1)*mu*mu * exp( -0.5*argpert_phi_1*argpert_phi_1);
        
        // const CCTK_REAL argpert_phi_2 = (rr_2 - R0pert_phi)/Sigmapert_phi;
        // const CCTK_REAL pert_phi_2 = 1. + Apert_phi * (x1_2*x1_2 - y1_2*y1_2)*mu*mu * exp( -0.5*argpert_phi_2*argpert_phi_2);

        const CCTK_REAL phi0_l_1 = phi0_1[ind];// * pert_phi_1;
        const CCTK_REAL phi0_l_2 = 0.;

        // scalar fields

              //star 1
        CCTK_REAL phi1_1  = phi0_l_1 * (coswt * cosmph_1 + sinwt * sinmph_1);
        CCTK_REAL phi2_1  = phi0_l_1 * (coswt * sinmph_1 - sinwt * cosmph_1);
        

              //BH 2
        CCTK_REAL phi1_2  = 0;//phi0_l_2 * (coswt * cosmph_2 + sinwt * sinmph_2);
        CCTK_REAL phi2_2  = 0;//phi0_l_2 * (coswt * sinmph_2 - sinwt * cosmph_2);
        /////////////////////////////

        phi1[ind]  = phi1_1 + phi1_2;
        phi2[ind]  = phi2_1 + phi2_2;

        

        const CCTK_REAL alph_1 = exp(F0_1[ind]);
        const CCTK_REAL alph_2 = new_lapse;

        // No regularization needed for the BS, the lapse is non-zero

        CCTK_REAL Kphi1_1 = 0.5 * (mm * W_1[ind] - omega_BS) / alph_1 * phi2_1;
        CCTK_REAL Kphi2_1 = 0.5 * (omega_BS - mm * W_1[ind]) / alph_1 * phi1_1;

        CCTK_REAL Kphi1_2 = 0.;
        CCTK_REAL Kphi2_2 = 0.;


        Kphi1[ind] = Kphi1_1 + Kphi1_2;
        Kphi2[ind] = Kphi2_1 + Kphi2_2;


        

        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1_1 + psi1_2 - 1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "Kerr_BS")) {
          alp[ind] = alph_1 + alph_2 - 1;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // CCTK_REAL Delt  = rBL*rBL + bh_spin2 - 2 * bh_mass * rBL;
        // CCTK_REAL fctFF = ( rBL*rBL + bh_spin2 ) * ( rBL*rBL + bh_spin2 ) - Delt * bh_spin2 * sinth2;
        // bphi = 2.0 * bh_spin * bh_mass * rBL / fctFF;

        // shift
        if (CCTK_EQUALS(initial_shift, "Kerr_BS")) {
          betax[ind] =  betad[1] ;//por enquato o shift da bs é zero pois é estática.
          betay[ind] =  betad[2];
          betaz[ind] =   0.;
        }

      } /* for i */
    }   /* for j */
  }     /* for k */


  free(F1_1); free(F2_1); free(F0_1); free(phi0_1); free(W_1);
  free(dW_dr_1); free(dW_dth_1);

  // free(F1_2); free(F2_2); free(F0_2); free(phi0_2); free(W_2);
  // free(dW_dr_2); free(dW_dth_2);

  return;
}
