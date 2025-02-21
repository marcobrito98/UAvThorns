
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

        const CCTK_REAL rr2_2_1 = x1_1*x1_1 + y1_1*y1_1 + z1_1*z1_1;

        CCTK_REAL rr_1  = sqrt(rr2_2_1);
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



  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        //Boson Star A

        const CCTK_REAL x1_1  = x[ind] - x0;
        const CCTK_REAL y1_1  = y[ind] - y0;
        const CCTK_REAL z1_1  = z[ind] - z0;

        // For the Boson Star, r = R, no coordinate change needed.
        const CCTK_REAL rr2_2_1 = x1_1*x1_1 + y1_1*y1_1 + z1_1*z1_1;
        const CCTK_REAL rr_1  = sqrt(rr2_2_1);

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
        const CCTK_REAL bh_spin2 = bh_spin*bh_spin;
        const CCTK_REAL gamma2 = 1. / (1. - bh_v2);
        const CCTK_REAL gamma = sqrt(gamma2);
        // const CCTK_REAL rr2_2 = gamma2*x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;
        // const CCTK_REAL rr_2  = sqrt(rr2_2);
        const CCTK_REAL rr2_2 = x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;
        const CCTK_REAL rr_2  = sqrt(rr2_2);


        // const CCTK_REAL rho2_2 = gamma2*x1_2*x1_2 + y1_2*y1_2;
        // const CCTK_REAL rho_2  = sqrt(rho2_2);
        const CCTK_REAL rho2_2 = x1_2*x1_2 + y1_2*y1_2;
        const CCTK_REAL rho_2  = sqrt(rho2_2);
        

        const CCTK_REAL theta_2 = acos(z1_2/rr_2);

        const CCTK_REAL deltakerr2_2 = bh_mass*bh_mass - bh_spin2 ;
        const CCTK_REAL deltakerr  = sqrt(deltakerr2_2) ;

        const CCTK_REAL costh  = z1_2/rr_2 ;
        const CCTK_REAL costh2 = costh*costh ;
        const CCTK_REAL sinth2 = 1. - costh2 ;
        const CCTK_REAL sinth  = sqrt(sinth2) ;

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
        const CCTK_REAL sinthth_z  = -sinth2/rr_2 ; 

        // const CCTK_REAL sinthx_th  = gamma*x1_2 * costh ;
        const CCTK_REAL sinthx_th  = x1_2 * costh ;
        const CCTK_REAL sinthy_th  = y1_2 * costh ;
        const CCTK_REAL sinthz_th  = -rr_2 * sinth2 ;


        const CCTK_REAL rBL    = rr_2 + bh_mass + 0.25*deltakerr2_2 / rr_2 ;   // Boyer-Lindquist coordinate r

        const CCTK_REAL RRrBL  = rr2_2 + rr_2*bh_mass + 0.25*deltakerr2_2 ;

        const CCTK_REAL rho2kerr   = rBL*rBL + bh_spin2 * costh2 ;
        const CCTK_REAL rhokerr    = sqrt(rho2kerr) ;

        const CCTK_REAL sigma  = (2.*bh_mass*rBL)/rho2kerr;
        const CCTK_REAL hh     = (1 + sigma) / (RRrBL*RRrBL + rr2_2*bh_spin*bh_spin * costh2) ;

        const CCTK_REAL psi4_2 = rho2kerr / rr2_2 ;
        const CCTK_REAL psi2_2 = sqrt(psi4_2) ;
        const CCTK_REAL psi1_2 = sqrt(psi2_2) ;
        const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        const CCTK_REAL psi2_1 = sqrt(psi4_1);
        const CCTK_REAL psi1_1 = sqrt(psi2_1);

        // non-axisymmetric perturbation.
        /* pert = 1. + AA * (x1_2*x1_2 - y1_2*y1_2)/(bh_mass*bh_mass) * exp( -2.*rr2_2/deltakerr2_2 ) ; */
        
        const CCTK_REAL alpha0  = (rr_2 + 0.5*deltakerr)*(rr_2 - 0.5*deltakerr) / rr_2 * \
                 1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2)) ;
        const CCTK_REAL alpha02 = alpha0*alpha0 ;

        const CCTK_REAL dr_dR = 1 + (bh_spin2 - bh_mass*bh_mass)/(4*rr2_2);
        const CCTK_REAL delta_metric = rBL*rBL-2*bh_mass*rBL+bh_spin2;
        const CCTK_REAL betadphi = -bh_spin*sigma*sinth2;
        const CCTK_REAL dbetadphi_dth = -(4*bh_spin*bh_mass*rBL*(bh_spin2+rBL*rBL)*sinth*costh)/pow(rho2kerr,2);
        const CCTK_REAL dbetadphi_dR = dr_dR*2*bh_spin*bh_mass*(rBL*rBL-bh_spin2*costh2)*sinth2/pow(rho2kerr,2);

        const CCTK_REAL gammaphiphi= psi4_2*rr2_2*sinth2*(1 + bh_spin2*hh*rr2_2*sinth2);
        //const CCTK_REAL dgammaphiphi_dth= -2*bh_spin2*delta_metric*sinth*costh/rr2_2;
        const CCTK_REAL dgammaphiphi_dth= (delta_metric+8*bh_mass*rBL*pow(bh_spin2+rBL*rBL,2)/pow(bh_spin2+2*rBL*rBL+bh_spin2*(costh2-sinth2),2))*2*costh*sinth;
        //const CCTK_REAL dgammaphiphi_dth= 4*bh_spin2*bh_mass*rBL*(bh_spin2+rBL*rBL)*costh*sinth/pow(rho2kerr,2);
        // const CCTK_REAL dgammaphiphi_dR= dr_dR*2*(rr_2*(2*rBL*(bh_spin2+rBL*rBL)+bh_spin2*(bh_mass-rBL)*sinth2))/(rr2_2*rr_2) + \
        //                                 2*(-pow(bh_spin2+rBL*rBL,2)+bh_spin2*delta_metric*sinth2)/(rr2_2*rr_2);

        const CCTK_REAL dgammaphiphi_dR= dr_dR*(2*rBL*(bh_spin2+rBL*rBL)*(rBL*rBL+bh_spin2*(costh2-sinth2))*sinth2 + \
                                      2*bh_spin2*(rBL*(bh_spin2-bh_mass*rBL)+bh_spin2*(bh_mass-rBL)*costh2)*sinth2*sinth2)/pow(rho2kerr,2);
        const CCTK_REAL betauphi = betadphi/gammaphiphi;
        const CCTK_REAL dbetauphi_dth = (gammaphiphi*dbetadphi_dth - betadphi*dgammaphiphi_dth)/pow(gammaphiphi,2);
        const CCTK_REAL dbetauphi_dR = (gammaphiphi*dbetadphi_dR - betadphi*dgammaphiphi_dR)/pow(gammaphiphi,2);
       


        //capital Gs refer to the unboosted frame.

        // const CCTK_REAL Gtt = -alpha02 + betadphi*betauphi;
        // const CCTK_REAL Gxt = bh_spin*sigma*y1_2/rr2_2;
        // const CCTK_REAL Gxx = psi4_2*(1+bh_spin2*hh*y1_2*y1_2);
        // // const CCTK_REAL Gxy = -psi4_2*bh_spin2*hh*y1_2*gamma*x1_2; 
        // // const CCTK_REAL Gty = -bh_spin*sigma*gamma*x1_2/rr2_2;
        // const CCTK_REAL Gxy = -psi4_2*bh_spin2*hh*y1_2*x1_2; 
        // const CCTK_REAL Gty = -bh_spin*sigma*x1_2/rr2_2;
        // // const CCTK_REAL fff = bh_mass/(bh_spin-bh_spin);

        // printf("%.6f",fff);



        // check_nan_or_inf("betauphi",betauphi);
        // check_nan_or_inf("betadphi",betadphi);
        // check_nan_or_inf("gammaphiphi",gammaphiphi);
        // check_nan_or_inf("dbetadphi_dth",dbetadphi_dth);
        // check_nan_or_inf("dbetadphi_dR",dbetadphi_dR);
        // check_nan_or_inf("dbetauphi_dth",dbetauphi_dth);
        // check_nan_or_inf("dbetauphi_dR",dbetauphi_dR);
       
        // check_nan_or_inf("Gtt",Gtt);
        // check_nan_or_inf("Gxt",Gxt);
        // check_nan_or_inf("Gxx",Gxx);
        // check_nan_or_inf("Gty",Gty);
        // check_nan_or_inf("Gxy",Gxy);
        // check_nan_or_inf("delta_metric",delta_metric);
        // check_nan_or_inf("alpha0",alpha0);
        // check_nan_or_inf("hh",hh);
        // check_nan_or_inf("sigma",sigma);
        // check_nan_or_inf("psi4_2",psi4_2);

        // check_nan_or_inf("fff",fff);


        // 3-metric
        // gxx[ind] = gamma2*Gxx + 2*gamma2*bh_v*Gxt + gamma2*bh_v2*Gtt;
        gxx[ind] =pow(psi1_1+psi1_2-1,4)*(1. + bh_spin2*hh*y1_2*y1_2);
        // gxy[ind] = gamma*Gxy+gamma*bh_v*Gty;
        gxy[ind] = -pow(psi1_2,4)*bh_spin2*hh*y1_2*x1_2;// deve ser o termo problematico. talvez testar a outra sobreposição (que nao contruibui para esta componente)
        gxz[ind] = 0;
        // gyy[ind] = psi4_2 * ( 1. + bh_spin2 * hh * gamma2*x1_2*x1_2 );
        gyy[ind] = pow(psi1_1+psi1_2-1,4)* (1. + bh_spin2 * hh * x1_2*x1_2) ;
        gyz[ind] = 0;
        gzz[ind] = pow(psi1_1+psi1_2-1,4);

        check_nan_or_inf("gxx",gxx[ind]);
        check_nan_or_inf("gxy",gxy[ind]);
        check_nan_or_inf("gxz",gxz[ind]);
        check_nan_or_inf("gyy",gyy[ind]);
        check_nan_or_inf("gyz",gyz[ind]);
        check_nan_or_inf("gzz",gzz[ind]);


        const CCTK_REAL HF     = - bh_spin2*bh_spin * alpha0 * sigma/rhokerr * costh  ;  // we are dividing by sinth2
        const CCTK_REAL Athph  = HF / rr_2 ;                                        // we are dividing by sinth

        const CCTK_REAL aux    =  rho2kerr * (rBL*rBL - bh_spin2) + 2.*rBL*rBL * (rBL*rBL + bh_spin2);

        const CCTK_REAL HE     = bh_spin*bh_mass * aux / (rhokerr*rhokerr*rhokerr) * 
                 1. / sqrt(rBL*rBL + bh_spin2 * ( 1. + sigma*sinth2)) ;

        const CCTK_REAL ARph   = HE / rr2_2 ;                                       // we are dividing by sinth2



        // //capital Ks refer to the unboosted frame.
        // const CCTK_REAL Ktht = betadphi*dbetauphi_dth/(-2*alpha0);
        // const CCTK_REAL KRt = betadphi*dbetauphi_dR/(-2*alpha0);

        // // const CCTK_REAL Kxt = R_x*KRt + gamma*x1_2*z1_2/(rho_2*rr2_2) * Ktht;
        // const CCTK_REAL Kxt = R_x*KRt + x1_2*z1_2/(rho_2*rr2_2) * Ktht;
        // const CCTK_REAL Kyt = R_y*KRt + y1_2*z1_2/(rho_2*rr2_2) * Ktht;
        // const CCTK_REAL Kzt = R_z*KRt + rho_2/rr2_2 * Ktht;


        const CCTK_REAL Axx = 2.*ARph *  R_x * sinth2ph_x                     +  2.*Athph *  sinthth_x * sinth2ph_x ;
        const CCTK_REAL Axy =    ARph * (R_x * sinth2ph_y + R_y * sinth2ph_x) +     Athph * (sinthth_x * sinth2ph_y + sinthth_y * sinth2ph_x) ;
        const CCTK_REAL Axz =    ARph *                     R_z * sinth2ph_x  +     Athph *                           sinthth_z * sinth2ph_x  ; 
        const CCTK_REAL Ayy = 2.*ARph *  R_y * sinth2ph_y                     +  2.*Athph *  sinthth_y * sinth2ph_y ;
        const CCTK_REAL Ayz =    ARph *                     R_z * sinth2ph_y  +     Athph *                           sinthth_z * sinth2ph_y  ;


        //K esta errado. tinha usado foliacao errada.

        // extrinsic curvature (this will be zero due to W=0 for the boson star. only BH matters) No caso em repouso, estara correto? verificar.
        kxx[ind] = Axx / psi2_2;
        kxy[ind] = Axy / psi2_2;
        kxz[ind] = Axz / psi2_2;
        kyy[ind] = Ayy / psi2_2;
        kyz[ind] = Ayz / psi2_2;
        kzz[ind] = 0.;

        
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
        const CCTK_REAL alph_2 = alpha0;

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
        else if (CCTK_EQUALS(initial_lapse, "BH_BS")) {
          alp[ind] = alph_1 + alph_2 - 1;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // shift
        if (CCTK_EQUALS(initial_shift, "BH_BS")) {
          betax[ind] =  0.;
          betay[ind] = 0.;
          betaz[ind] =  0.;
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
