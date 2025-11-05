// este e o codigo copiado do kerrnewman.
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

void UAv_ID_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL[], CCTK_REAL[],
                      CCTK_REAL[], CCTK_REAL[], CCTK_REAL[], CCTK_REAL[], CCTK_REAL[]);

void check_nan_or_inf(const char *var_name, double value)
{
  if (isnan(value))
  {
    fprintf(stderr, "Error: %s is NaN\n", var_name);
    abort(); // Break execution
  }
  else if (isinf(value))
  {
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

  CCTK_INT NF;     // NF will be the actual size of the arrays
  CCTK_INT NX;     // NX will be the number of X points
  CCTK_INT Ntheta; // Ntheta will be the number of theta points

  CCTK_REAL *Xtmp, *thtmp, *F1_in, *F2_in, *F0_in, *phi0_in, *Wbar_in;
  Xtmp = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  thtmp = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  F1_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  F2_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  F0_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  phi0_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  Wbar_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));

  // we get the data from the input file
  UAv_ID_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, phi0_in, Wbar_in);

  Ntheta = NF / NX;

  CCTK_VInfo(CCTK_THORNSTRING, "NX     = %d", NX);
  CCTK_VInfo(CCTK_THORNSTRING, "Ntheta = %d", Ntheta);
  CCTK_VInfo(CCTK_THORNSTRING, "NF     = %d", NF);

  // now we create arrays with the X and theta coordinates
  CCTK_REAL X[NX], theta[Ntheta];
  for (int i = 0; i < NX; i++)
  {
    X[i] = Xtmp[i];
    /* printf("X[%3d] = %lf\n", i, X[i]); */
  }
  for (int i = 0; i < Ntheta; i++)
  {
    theta[i] = thtmp[i * NX];
    /* printf("theta[%3d] = %lf\n", i, theta[i]); */
  }

  // the spacing in each coordinate is
  const CCTK_REAL dX = (X[NX - 1] - X[0]) / (NX - 1);
  const CCTK_REAL dtheta = (theta[Ntheta - 1] - theta[0]) / (Ntheta - 1);

  /* printf("dX     = %e\n", dX); */
  /* printf("dtheta = %e\n", dtheta); */

  // make sure spacing is uniform in the provided grid
  for (int i = 1; i < NX; i++)
  {
    if (fabs(X[i] - X[i - 1] - dX) > SMALL)
    {
      printf("i = %d\n", i);
      CCTK_WARN(0, "X grid is not uniformly spaced. Aborting.");
    }
  }
  for (int j = 1; j < Ntheta; j++)
  {
    if (fabs(theta[j] - theta[j - 1] - dtheta) > SMALL)
      CCTK_WARN(0, "theta grid is not uniformly spaced. Aborting.");
  }

  // now we need to take the derivatives of the Wbar function
  // Then we convert to W and store the values

  CCTK_REAL *W_in, *dW_dr_in, *dW_dth_in;
  W_in = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dW_dr_in = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dW_dth_in = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  const CCTK_REAL oodX = 1. / dX;
  // const CCTK_REAL oodXsq     = oodX * oodX;
  const CCTK_REAL oodX12 = 1. / (12. * dX);
  const CCTK_REAL oodXsq12 = oodX * oodX12;
  const CCTK_REAL oodth12 = 1. / (12. * dtheta);

  for (int jj = 0; jj < Ntheta; jj++)
  {
    for (int i = 0; i < NX; i++)
    {

      CCTK_INT j, jm1, jm2, jp1, jp2;
      /* let's use the fact that the solution is axi-symmetric (and that
         theta[0] = 0) for the boundary points in j */
      if (jj == 0)
      {
        j = jj;
        jp1 = jj + 1;
        jp2 = jj + 2;
        jm1 = jj + 1;
        jm2 = jj + 2;
      }
      else if (jj == 1)
      {
        j = jj;
        jp1 = jj + 1;
        jp2 = jj + 2;
        jm1 = jj - 1;
        jm2 = jj;
      }
      else if (jj == Ntheta - 2)
      {
        j = jj;
        jm1 = jj - 1;
        jm2 = jj - 2;
        jp1 = jj + 1;
        jp2 = jj;
      }
      else if (jj == Ntheta - 1)
      {
        j = jj;
        jm1 = jj - 1;
        jm2 = jj - 2;
        jp1 = jj - 1;
        jp2 = jj - 2;
      }
      else
      {
        j = jj;
        jp1 = jj + 1;
        jp2 = jj + 2;
        jm1 = jj - 1;
        jm2 = jj - 2;
      }

      const CCTK_INT ind = i + j * NX;

      const CCTK_INT indim1 = i - 1 + j * NX;
      const CCTK_INT indip1 = i + 1 + j * NX;
      const CCTK_INT indim2 = i - 2 + j * NX;
      const CCTK_INT indip2 = i + 2 + j * NX;
      const CCTK_INT indip3 = i + 3 + j * NX;
      const CCTK_INT indip4 = i + 4 + j * NX;
      const CCTK_INT indip5 = i + 5 + j * NX;

      const CCTK_INT indjm1 = i + jm1 * NX;
      const CCTK_INT indjm2 = i + jm2 * NX;
      const CCTK_INT indjp1 = i + jp1 * NX;
      const CCTK_INT indjp2 = i + jp2 * NX;

      const CCTK_REAL lX = X[i];
      /* const CCTK_REAL lth = theta[j]; */
      /* printf("X[%3d] = %lf\n", i, lX); */

      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL Wbar_th = (-Wbar_in[indjp2] + 8 * Wbar_in[indjp1] - 8 * Wbar_in[indjm1] + Wbar_in[indjm2]) *
                                oodth12;

      CCTK_REAL Wbar_X;
      CCTK_REAL Wbar_XX = 0.; // Used for rhor=0 (i==0), if Wbar_r_power == 2.

      /*
      Regarding finite differencing orders: plotting W and dW_dr, there were small discontinuities near rhor=0
      during tests with the previous 2nd order accuracy for i==0 and i==1.
      Those vanish when moving to 4th order accuracy.

      For i==NX-1 and i==NX-2, we keep 2nd order for now. The issue is not appearing as clearly,
      and they represent points which are physically far, so maybe better to keep the computation more local.
      */

      if (i == 0)
      {
        /* For the Boson Star, there's no issue, dWbar/dX != 0 at X==0, and x and rhor coordinates coincide. */

        // 1st derivative with 4th order accuracy (forward stencils)
        Wbar_X = (-25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] + 16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) * oodX12;

        if (Wbar_r_power == 2)
        {
          // If Wbar = rhor^2 * W, to compute W(rhor=0), we need to compute Wbar_XX.
          // 2nd derivative with 4th order accuracy (forward stencils)
          Wbar_XX = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] + 214 * Wbar_in[indip2] - 156 * Wbar_in[indip3] + 61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) * oodXsq12;
        }
      }
      else if (i == 1)
      {
        // 1st derivative, 4th order accuracy
        Wbar_X = (-3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] - 6 * Wbar_in[indip2] + Wbar_in[indip3]) * oodX12;
      }
      else if (i == NX - 1)
      {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4 * Wbar_in[indim1] + 3 * Wbar_in[ind]) * 0.5 * oodX;
      }
      else if (i == NX - 2)
      {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;
      }
      else
      {
        // 4th order accurate stencils
        Wbar_X = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
      }

      // From the X coordinate used in the input files to the rhor coordinate (coincides with x for the Boson Star, rH=0).
      // We also do the conversion from Wbar to W here, to tackle rhor = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  rhor == 0
      if (i == 0)
      {
        // At rhor=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_in[ind] = 0.;
        dW_dth_in[ind] = 0.;

        // For W we need more care depending on the power
        switch (Wbar_r_power)
        {
        case 0: // Wbar = W
          W_in[ind] = Wbar_in[ind];
          break;

        case 1: // Wbar = rhor * W
          /*
          dWbar/dr = W + rhor * dW/dr
                   = W + 0 * 0     at rhor=0

          dWbar/dr = dWbar/dX * dX/dr
          dX/dr = C/(C+rhor)^2 = 1/C  at rhor=0
          */
          W_in[ind] = Wbar_X / C0;
          break;

        case 2: // Wbar = rhor^2 * W
          /*
          dWbar/dr   = 2r * W + rhor^2 * dW/dr
          d2Wbar/dr2 = 2  * W + 4r  * dW/dr + rhor^2 * d2W/dr2
                     = 2  * W + 0 + 0        at rhor=0

          d2Wbar/dr2 = d2Wbar/dX2 * (dX/dr)^2 + dWbar/dX * d2X/dr2
          dX/dr   =   C/(C+rhor)^2 =  1/C      at rhor=0
          d2X/dr2 = -2C/(C+rhor)^3 = -2/C^2    at rhor=0

          W (rhor=0) = 1/C^2 * [1/2 * d2Wbar/dX^2 (X=0)  -  dWbar/dX (X=0)]
          */
          W_in[ind] = (0.5 * Wbar_XX - Wbar_X) / (C0 * C0);
          break;

        default: // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }
      }

      // We need to be careful at X == 1 (rhor == infty) for radial derivatives (coordinate change is singular)
      else if (i == NX - 1)
      {
        // W -> 0 for rhor -> infty
        W_in[ind] = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction paper) also gives:
        dW_dr_in[ind] = 0.;
        dW_dth_in[ind] = 0.;
      }
      else
      {
        const CCTK_REAL rr = C0 * lX / (1. - lX);

        // corresponding derivatives
        // const CCTK_REAL dXdr = 1./(C0 + rr) - rr/((C0 + rr)*(C0 + rr));
        const CCTK_REAL dXdr = C0 / ((C0 + rr) * (C0 + rr));

        const CCTK_REAL Wbar_r = dXdr * Wbar_X;

        // Now translate from Wbar to W
        switch (Wbar_r_power) // We could put a generic power for the computation here I guess...
        {
        case 0: // Wbar = W
          W_in[ind] = Wbar_in[ind];
          dW_dr_in[ind] = Wbar_r;
          dW_dth_in[ind] = Wbar_th;
          break;

        case 1: // Wbar = rhor * W
          W_in[ind] = Wbar_in[ind] / rr;
          dW_dr_in[ind] = (Wbar_r - W_in[ind]) / rr; // dW/dr  =  1/rhor * dWbar/dr - Wbar / rhor^2  =  (dWbar/dr - W) / rhor
          dW_dth_in[ind] = Wbar_th / rr;
          break;

        case 2:; // Wbar = rhor^2 * W
          // empty statement after case to prevent compilation error on some gcc versions...
          const CCTK_REAL rr2_2 = rr * rr;
          W_in[ind] = Wbar_in[ind] / rr2_2;
          dW_dr_in[ind] = Wbar_r / rr2_2 - 2 * W_in[ind] / rr; // dW/dr  =  1/rhor^2 * dWbar/dr - 2 * Wbar / rhor^3  =  1/rhor^2 * dWbar/dr - 2 * W / rhor
          dW_dth_in[ind] = Wbar_th / rr2_2;
          break;

        default: // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }
      } // if/else i==...

    } // for i
  } // for jj

  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]; // total points

  CCTK_REAL *X_g_1, *theta_g_1;
  X_g_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  theta_g_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));

  // CCTK_REAL *X_g_2, *theta_g_2;
  // X_g_2     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  // theta_g_2 = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  for (int k = 0; k < cctk_lsh[2]; ++k)
  {
    for (int j = 0; j < cctk_lsh[1]; ++j)
    {
      for (int i = 0; i < cctk_lsh[0]; ++i)
      {

        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_REAL x1_1 = x[ind] - x0;
        const CCTK_REAL y1_1 = y[ind] - y0;
        const CCTK_REAL z1_1 = z[ind] - z0;

        const CCTK_REAL rr2_1 = x1_1 * x1_1 + y1_1 * y1_1 + z1_1 * z1_1;

        CCTK_REAL rr_1 = sqrt(rr2_1);
        /* For the Boson Star, x, rhor and R coordinates coincide (rH=0). */

        // From rhor to the X radial coordinate (used in input files)
        const CCTK_REAL lX_1 = rr_1 / (C0 + rr_1);

        CCTK_REAL ltheta_1 = rr_1 < 1e-16 ? 0 : acos(z1_1 / rr_1); // There should be at most one point in the grid with rr~0. Not sure about the threshold.
        if (ltheta_1 > 0.5 * M_PI)                                 // symmetry along the equatorial plane
          ltheta_1 = M_PI - ltheta_1;

        X_g_1[ind] = lX_1;
        theta_g_1[ind] = ltheta_1;
        // const CCTK_REAL x1_2  = x[ind] - x0_2;
        // const CCTK_REAL y1_2  = y[ind] - y0_2;
        // const CCTK_REAL z1_2  = z[ind] - z0_2;

        // const CCTK_REAL rr2_2 = x1_2*x1_2 + y1_2*y1_2 + z1_2*z1_2;

        // CCTK_REAL rr_2  = sqrt(rr2_2);
        // /* For the Boson Star, x, rhor and R coordinates coincide (rH=0). */

        // // From rhor to the X radial coordinate (used in input files)
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

  const CCTK_INT N_dims = 2; // 2-D interpolation

  const CCTK_INT N_input_arrays = 7;
  const CCTK_INT N_output_arrays = 7;

  /* origin and stride of the input coordinates. with this Cactus reconstructs
     the whole X and theta array. */
  CCTK_REAL origin[N_dims];
  CCTK_REAL delta[N_dims];
  origin[0] = X[0];
  origin[1] = theta[0];
  delta[0] = dX;
  delta[1] = dtheta;

  /* points onto which we want to interpolate, ie, the grid points themselves in
     (X, theta) coordinates (computed above) */
  const void *interp_coords_1[N_dims];
  interp_coords_1[0] = (const void *)X_g_1;
  interp_coords_1[1] = (const void *)theta_g_1;

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
  input_arrays[0] = (const void *)F1_in;
  input_arrays[1] = (const void *)F2_in;
  input_arrays[2] = (const void *)F0_in;
  input_arrays[3] = (const void *)phi0_in;
  input_arrays[4] = (const void *)W_in;
  input_arrays[5] = (const void *)dW_dr_in;
  input_arrays[6] = (const void *)dW_dth_in;

  /* output arrays */
  void *output_arrays_1[N_output_arrays];
  CCTK_INT output_array_type_codes_1[N_output_arrays];
  void *output_arrays_2[N_output_arrays];
  CCTK_INT output_array_type_codes_2[N_output_arrays];
  CCTK_REAL *F1_1, *F2_1, *F0_1, *phi0_1, *W_1;
  CCTK_REAL *dW_dr_1, *dW_dth_1;

  F1_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  F2_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  F0_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  phi0_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  W_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dr_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dth_1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));

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

  output_arrays_1[0] = (void *)F1_1;
  output_arrays_1[1] = (void *)F2_1;
  output_arrays_1[2] = (void *)F0_1;
  output_arrays_1[3] = (void *)phi0_1;
  output_arrays_1[4] = (void *)W_1;
  output_arrays_1[5] = (void *)dW_dr_1;
  output_arrays_1[6] = (void *)dW_dth_1;

  /* handle and settings for the interpolation routine */
  int operator_handle, param_table_handle;
  operator_handle = CCTK_InterpHandle("Lagrange polynomial interpolation");
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
  if (status_1 < 0)
  {
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

  free(X_g_1);
  free(theta_g_1);
  free(Xtmp);
  free(thtmp);
  free(F1_in);
  free(F2_in);
  free(F0_in);
  free(phi0_in);
  free(Wbar_in);
  free(W_in);
  free(dW_dr_in);
  free(dW_dth_in);

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

  const CCTK_REAL bh_spin2 = bh_spin * bh_spin;
  const CCTK_REAL bh_mass2 = bh_mass * bh_mass;

  const CCTK_REAL rBLp = bh_mass + sqrt(bh_mass2 - bh_spin2);
  const CCTK_REAL rBLm = bh_mass - sqrt(bh_mass2 - bh_spin2);

  const CCTK_REAL horizon_radius = 0.5 * sqrt(bh_mass2 - bh_spin2);

  // printf("cctk_lsh[0] = %d\n",cctk_lsh[0]);
  // printf("cctk_lsh[1] = %d\n",cctk_lsh[1]);
  // printf("cctk_lsh[2] = %d\n",cctk_lsh[2]);

  for (int k = 0; k < cctk_lsh[2]; ++k)
  {
    for (int j = 0; j < cctk_lsh[1]; ++j)
    {
      for (int i = 0; i < cctk_lsh[0]; ++i)
      {

        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // Boson Star A

        const CCTK_REAL x1_1 = x[ind] - x0;
        const CCTK_REAL y1_1 = y[ind] - y0;
        const CCTK_REAL z1_1 = z[ind] - z0;

        // For the Boson Star, rhor = R, no coordinate change needed.
        const CCTK_REAL rr2_1 = x1_1 * x1_1 + y1_1 * y1_1 + z1_1 * z1_1;
        const CCTK_REAL rr_1 = sqrt(rr2_1);

        const CCTK_REAL rho2_1 = x1_1 * x1_1 + y1_1 * y1_1;
        const CCTK_REAL rho_1 = sqrt(rho2_1);

        const CCTK_REAL ph_1 = atan2(y1_1, x1_1);
        // If x1_2=y1_2=0, should return 0? The other metric functions should vanish anyway to make sure that this doesn't matter,
        // but can this lead to nan depending on the C implementation?

        const CCTK_REAL cosph_1 = cos(ph_1);
        const CCTK_REAL sinph_1 = sin(ph_1);

        const CCTK_REAL cosmph_1 = cos(mm * ph_1);
        const CCTK_REAL sinmph_1 = sin(mm * ph_1);

        const CCTK_REAL h_rho2_1 = exp(2. * (F2_1[ind] - F1_1[ind])) - 1.;

        // Black Hole B

        const CCTK_REAL x1_2 = x[ind] - x0_2;
        const CCTK_REAL y1_2 = y[ind] - y0_2;
        const CCTK_REAL z1_2 = z[ind] - z0_2;

        const CCTK_REAL bh_v2 = bh_v * bh_v;
        const CCTK_REAL gamma2 = 1. / (1. - bh_v2);
        const CCTK_REAL gamma = sqrt(gamma2);

        CCTK_REAL rr2_2 = x1_2 * x1_2 * gamma2 + y1_2 * y1_2 + z1_2 * z1_2;
        if (rr2_2 < pow(eps_r, 2))
        {
          rr2_2 = pow(eps_r, 2);
        }
        const CCTK_REAL rr_2 = sqrt(rr2_2);

        CCTK_REAL rho2_2 = x1_2 * x1_2 * gamma2 + y1_2 * y1_2;
        if (rho2_2 < pow(eps_r, 2))
        {
          rho2_2 = pow(eps_r, 2);
        }
        const CCTK_REAL rho_2 = sqrt(rho2_2);

        const CCTK_REAL rho3_2 = rho2_2 * rho_2;

        // R0pert2 = (rr_2 - R0pert)*(rr_2 - R0pert) ;

        const CCTK_REAL costh = z1_2 / rr_2;
        const CCTK_REAL costh2 = costh * costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth = sqrt(sinth2);

        const CCTK_REAL R_x = x1_2 / rr_2;
        const CCTK_REAL R_y = y1_2 / rr_2;
        const CCTK_REAL R_z = z1_2 / rr_2;

        const CCTK_REAL x_R = x1_2 / rr_2;
        const CCTK_REAL y_R = y1_2 / rr_2;
        const CCTK_REAL z_R = z1_2 / rr_2;

        const CCTK_REAL sinth2ph_x = -y1_2 / rr2_2;
        const CCTK_REAL sinth2ph_y = x1_2 / rr2_2;

        const CCTK_REAL sinthth_x = z1_2 * x1_2 / (rr_2 * rr2_2);
        const CCTK_REAL sinthth_y = z1_2 * y1_2 / (rr_2 * rr2_2);
        const CCTK_REAL sinthth_z = -sinth2 / rr_2;

        const CCTK_REAL sinthx_th = x1_2 * costh;
        const CCTK_REAL sinthy_th = y1_2 * costh;
        const CCTK_REAL sinthz_th = -rr_2 * sinth2;

        CCTK_REAL rBL, rBL2;
        rBL = rr_2 * (1.0 + 0.25 * rBLp / rr_2) * (1.0 + 0.25 * rBLp / rr_2);
        rBL2 = rBL * rBL;

        CCTK_REAL drBLdR;
        drBLdR = 1.0 - rBLp * rBLp / (16.0 * rr2_2);

        CCTK_REAL Delt, Sigm, Sigm2, fctFF;
        Delt = rBL2 + bh_spin2 - 2 * bh_mass * rBL;
        Sigm = rBL2 + bh_spin2 * costh2;
        Sigm2 = Sigm * Sigm;
        fctFF = (rBL2 + bh_spin2) * (rBL2 + bh_spin2) - Delt * bh_spin2 * sinth2;

        const CCTK_REAL psi4_2 = Sigm / rr2_2;
        const CCTK_REAL psi2_2 = sqrt(psi4_2);
        const CCTK_REAL psi1_2 = sqrt(psi2_2);
        const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        const CCTK_REAL psi2_1 = sqrt(psi4_1);
        const CCTK_REAL psi1_1 = sqrt(psi2_1);

        CCTK_REAL fctGG, fctHH, dfctGG_dR, dfctHH_dR, dfctHH_dth;
        fctGG = rBLm / (rr2_2 * (rBL - rBLm));
        dfctGG_dR = -rBLm * (-rBLm + rBL + rr_2 * drBLdR) / (rr2_2 * pow(rBLm - rBL, 2));
        fctHH = (2.0 * bh_mass * rBL + Sigm) / (rr2_2 * Sigm2);
        dfctHH_dR = (-2 * (bh_spin2 * costh2 + rBL2) * (bh_spin2 * costh2 + 2 * bh_mass * rBL + rBL2) - 2 * rr_2 * (-(bh_mass * bh_spin2 * costh2) + bh_spin2 * costh2 * rBL + 3 * bh_mass * rBL2 + pow(rBL, 3)) * drBLdR) / (pow(rr_2, 3) * pow(bh_spin2 * costh2 + rBL2, 3));
        dfctHH_dth = (2 * bh_spin2 * costh * (bh_spin2 * costh2 + 4 * bhmass * rBL + pow(rBL, 2)) * sinth) / (rr2_2 * pow(bh_spin2 * costh2 + rBL2, 3));

        CCTK_REAL detgij;
        detgij = pow(psi4_2, 3) * (1.0 + rr2_2 * fctGG) * (1.0 + bh_spin2 * rho2_2 * fctHH);

        /*----------------------------------*/

        /*=== initialize gauge functions ===*/
        /*----------------------------------*/

        const CCTK_REAL alpha0 = (4.0 * rr_2 - rBLp) * sqrt(rBL - rBLm) / sqrt(16.0 * rr_2 * (rBL2 + bh_spin2 * (1.0 + 2.0 * bh_mass * rBL * sinth2 / Sigm)));
        const CCTK_REAL alpha02 = alpha0 * alpha0;
        const CCTK_REAL bphiup = -2.0 * bh_mass * bh_spin * rBL / fctFF;
        const CCTK_REAL bphi = bphiup * psi4_2 * sinth2;

        CCTK_REAL G[4][4];
        // Initialize G to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            G[i][j] = 0.0;

        G[0][0] = -alpha02 + bphi * bphiup;
        G[0][1] = -y1_2/rho2_2 * bphi;
        G[0][2] = x1_2/rho2_2 * bphi;
        G[0][3] = 0;
        G[1][0] = G[0][1];
        G[2][0] = G[0][2];
        G[3][0] = G[0][3];
        G[1][1] = psi4_2 * (1.0 + x1_2 * x1_2 * fctGG + bh_spin2 * y1_2 * y1_2 * fctHH);
        G[1][2] = psi4_2 * (x1_2 * y1_2 * fctGG - bh_spin2 * x1_2 * y1_2 * fctHH);
        G[1][3] = psi4_2 * (x1_2 * z1_2 * fctGG);
        G[2][1] = G[1][2];
        G[2][2] = psi4_2 * (1.0 + y1_2 * y1_2 * fctGG + bh_spin2 * x1_2 * x1_2 * fctHH);
        G[2][3] = psi4_2 * (y1_2 * z1_2 * fctGG);
        G[3][1] = G[1][3];
        G[3][2] = G[2][3];
        G[3][3] = psi4_2 * (1.0 + z1_2 * z1_2 * fctGG);

        // const CCTK_REAL Gtt = G[0][0];
        // const CCTK_REAL Gxt = G[0][1];
        // const CCTK_REAL Gxx = G[1][1];
        // const CCTK_REAL Gxy = G[1][2];
        // const CCTK_REAL Gty = G[0][2];
        // const CCTK_REAL Gyy = G[2][2];

        CCTK_REAL G_inv[4][4];
        // Initialize G_inv to zero. Inverse of 3 metric
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            G_inv[i][j] = 0.0;

        /* Invert the spatial 3x3 block of G into G_inv */
        {
          const CCTK_REAL G11 = G[1][1];
          const CCTK_REAL G12 = G[1][2];
          const CCTK_REAL G13 = G[1][3];
          const CCTK_REAL G22 = G[2][2];
          const CCTK_REAL G23 = G[2][3];
          const CCTK_REAL G33 = G[3][3];

          const CCTK_REAL det =
              G11 * (G22 * G33 - G23 * G23) - G12 * (G12 * G33 - G13 * G23) + G13 * (G12 * G23 - G13 * G22);

          const CCTK_REAL inv_det = 1.0 / det;

          G_inv[1][1] = (G22 * G33 - G23 * G23) * inv_det;
          G_inv[1][2] = (G13 * G23 - G12 * G33) * inv_det;
          G_inv[1][3] = (G12 * G23 - G13 * G22) * inv_det;

          G_inv[2][1] = G_inv[1][2];
          G_inv[2][2] = (G11 * G33 - G13 * G13) * inv_det;
          G_inv[2][3] = (G13 * G12 - G11 * G23) * inv_det;

          G_inv[3][1] = G_inv[1][3];
          G_inv[3][2] = G_inv[2][3];
          G_inv[3][3] = (G11 * G22 - G12 * G12) * inv_det;
        }

        // CCTK_REAL betaup[4];
        // betaup[0] = 0.0;
        // betaup[1] = G_inv[1][1] * betad[1] + G_inv[1][2] * betad[2] + G_inv[1][3] * betad[3];
        // betaup[2] = G_inv[2][1] * betad[1] + G_inv[2][2] * betad[2] + G_inv[2][3] * betad[3];
        // betaup[3] = G_inv[3][1] * betad[1] + G_inv[3][2] * betad[2] + G_inv[3][3] * betad[3];

        // Derivatives of the metric functions
        CCTK_REAL dG[4][4][4];
        for (int a = 0; a < 4; ++a)
        {
          for (int b = 0; b < 4; ++b)
          {
            for (int c = 0; c < 4; ++c)
            {
              dG[a][b][c] = 0.0;
            }
          }
        }
        // dG[a][b][c] = dG_ab/dx^c

        dG[1][1][1] = pow(bh_spin, 2) * pow(y1_2, 2) * psi4_2 * dhh_dx + (1 + pow(bh_spin, 2) * pow(y1_2, 2) * hh) * dpsi4_2_dx;
        dG[1][1][2] = pow(bh_spin, 2) * y1_2 * psi4_2 * (2 * hh + y1_2 * dhh_dy) + (1 + pow(bh_spin, 2) * pow(y1_2, 2) * hh) * dpsi4_2_dy;
        dG[1][1][3] = pow(bh_spin, 2) * pow(y1_2, 2) * psi4_2 * dhh_dz + (1 + pow(bh_spin, 2) * pow(y1_2, 2) * hh) * dpsi4_2_dz;
        dG[1][2][1] = -(pow(bh_spin, 2) * y1_2 * (x1_2 * gamma * psi4_2 * dhh_dx + hh * (psi4_2 + x1_2 * gamma * dpsi4_2_dx)));
        dG[1][2][2] = -(pow(bh_spin, 2) * x1_2 * gamma * (y1_2 * psi4_2 * dhh_dy + hh * (psi4_2 + y1_2 * dpsi4_2_dy)));
        dG[1][2][3] = -(pow(bh_spin, 2) * x1_2 * gamma * y1_2 * (psi4_2 * dhh_dz + hh * dpsi4_2_dz));

        // dG[1][3][i] = 0

        dG[2][2][1] = pow(bh_spin, 2) * x1_2 * gamma * psi4_2 * (2 * hh + x1_2 * gamma * dhh_dx) + (1 +
                                                                                                    pow(bh_spin, 2) * pow(x1_2 * gamma, 2) * hh) *
                                                                                                       dpsi4_2_dx;
        dG[2][2][2] = pow(bh_spin, 2) * pow(x1_2 * gamma, 2) * psi4_2 * dhh_dy + (1 +
                                                                                  pow(bh_spin, 2) * pow(x1_2 * gamma, 2) * hh) *
                                                                                     dpsi4_2_dy;
        dG[2][2][3] = pow(bh_spin, 2) * pow(x1_2 * gamma, 2) * psi4_2 * dhh_dz + (1 +
                                                                                  pow(bh_spin, 2) * pow(x1_2 * gamma, 2) * hh) *
                                                                                     dpsi4_2_dz;

        // dG23_dx^i = 0

        dG[3][3][1] = dpsi4_2_dx;
        dG[3][3][2] = dpsi4_2_dy;
        dG[3][3][3] = dpsi4_2_dz;

        // dG[mu][nu][0] = 0
        dG[0][0][1] = -2 * alpha0 * dalpha_dx - (pow(bh_spin, 2) * (pow(betad[1], 2) + pow(betad[2], 2) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * hh) * (2 * x1_2 * gamma * hh + (rho2_2)*dhh_dx)) / (pow(1 + pow(bh_spin, 2) * (rho2_2)*hh, 2) * psi4_2) + (2 * betad[1] * dbetad[1][1] + 2 * betad[2] * dbetad[2][1] + 2 * pow(bh_spin, 2) * (x1_2 * gamma * betad[1] + y1_2 * betad[2]) * hh * (betad[1] + x1_2 * gamma * dbetad[1][1] + y1_2 * dbetad[2][1]) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * dhh_dx) / ((1 + pow(bh_spin, 2) * (rho2_2)*hh) * psi4_2) - ((pow(betad[1], 2) + pow(betad[2], 2) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * hh) * dpsi4_2_dx) / ((1 + pow(bh_spin, 2) * (rho2_2)*hh) * pow(psi4_2, 2));
        dG[0][0][2] = -2 * alpha0 * dalpha_dy - (pow(bh_spin, 2) * (pow(betad[1], 2) + pow(betad[2], 2) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * hh) * (2 * y1_2 * hh + (rho2_2)*dhh_dy)) / (pow(1 + pow(bh_spin, 2) * (rho2_2)*hh, 2) * psi4_2) + (2 * betad[1] * dbetad[1][2] + 2 * betad[2] * dbetad[2][2] + 2 * pow(bh_spin, 2) * (x1_2 * gamma * betad[1] + y1_2 * betad[2]) * hh * (betad[2] + x1_2 * gamma * dbetad[1][2] + y1_2 * dbetad[2][2]) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * dhh_dy) / ((1 + pow(bh_spin, 2) * (rho2_2)*hh) * psi4_2) - ((pow(betad[1], 2) + pow(betad[2], 2) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * hh) * dpsi4_2_dy) / ((1 + pow(bh_spin, 2) * (rho2_2)*hh) * pow(psi4_2, 2));
        dG[0][0][3] = -2 * alpha0 * dalpha_dz - (pow(bh_spin, 2) * (rho2_2) * (pow(betad[1], 2) + pow(betad[2], 2) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * hh) * dhh_dz) / (pow(1 + pow(bh_spin, 2) * (rho2_2)*hh, 2) * psi4_2) + (2 * betad[1] * dbetad[1][3] + 2 * betad[2] * dbetad[2][3] + 2 * pow(bh_spin, 2) * (x1_2 * gamma * betad[1] + y1_2 * betad[2]) * hh * (x1_2 * gamma * dbetad[1][3] + y1_2 * dbetad[2][3]) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * dhh_dz) / ((1 + pow(bh_spin, 2) * (rho2_2)*hh) * psi4_2) - ((pow(betad[1], 2) + pow(betad[2], 2) + pow(bh_spin, 2) * pow(x1_2 * gamma * betad[1] + y1_2 * betad[2], 2) * hh) * dpsi4_2_dz) / ((1 + pow(bh_spin, 2) * (rho2_2)*hh) * pow(psi4_2, 2));
        dG[0][1][1] = dbetad[1][1];
        dG[0][1][2] = dbetad[1][2];
        dG[0][1][3] = dbetad[1][3];
        dG[0][2][1] = dbetad[2][1];
        dG[0][2][2] = dbetad[2][2];
        dG[0][2][3] = dbetad[2][3];
        dG[0][3][1] = dbetad[3][1];
        dG[0][3][2] = dbetad[3][2];
        dG[0][3][3] = dbetad[3][3];
        // symmetries
        dG[2][1][1] = dG[1][2][1];
        dG[2][1][2] = dG[1][2][2];
        dG[2][1][3] = dG[1][2][3];
        dG[3][1][1] = dG[1][3][1];
        dG[3][1][2] = dG[1][3][2];
        dG[3][1][3] = dG[1][3][3];
        dG[3][2][1] = dG[2][3][1];
        dG[3][2][2] = dG[2][3][2];
        dG[3][2][3] = dG[2][3][3];
        dG[1][0][1] = dG[0][1][1];
        dG[1][0][2] = dG[0][1][2];
        dG[1][0][3] = dG[0][1][3];
        dG[2][0][1] = dG[0][2][1];
        dG[2][0][2] = dG[0][2][2];
        dG[2][0][3] = dG[0][2][3];
        dG[3][0][1] = dG[0][3][1];
        dG[3][0][2] = dG[0][3][2];
        dG[3][0][3] = dG[0][3][3];

        CCTK_REAL invLambda[4][4];
        for (int a = 0; a < 4; ++a)
        {
          for (int b = 0; b < 4; ++b)
          {
            invLambda[a][b] = 0.0;
          }
        }
        invLambda[0][0] = gamma;
        invLambda[0][1] = gamma * bh_v;
        invLambda[1][0] = gamma * bh_v;
        invLambda[1][1] = gamma;
        invLambda[2][2] = 1.;
        invLambda[3][3] = 1.;

        CCTK_REAL Lambda[4][4];
        for (int a = 0; a < 4; ++a)
        {
          for (int b = 0; b < 4; ++b)
          {
            Lambda[a][b] = 0.0;
          }
        }
        Lambda[0][0] = gamma;
        Lambda[0][1] = -gamma * bh_v;
        Lambda[1][0] = -gamma * bh_v;
        Lambda[1][1] = gamma;
        Lambda[2][2] = 1.;
        Lambda[3][3] = 1.;

        CCTK_REAL Gb[4][4];
        for (int a = 0; a < 4; ++a)
        {
          for (int b = 0; b < 4; ++b)
          {
            Gb[a][b] = 0.0;
          }
        }

        for (int a = 0; a < 4; ++a)
        {
          for (int b = 0; b < 4; ++b)
          {
            CCTK_REAL sum = 0.0;
            for (int mu = 0; mu < 4; ++mu)
              for (int nu = 0; nu < 4; ++nu)
                sum += invLambda[mu][a] * invLambda[nu][b] * G[mu][nu];
            Gb[a][b] = sum;
          }
        }

        // depois de testado fazer aqui a sobreposição
        gxx[ind] = Gb[1][1];
        gxy[ind] = Gb[1][2];
        gxz[ind] = Gb[1][3];
        gyy[ind] = Gb[2][2];
        gyz[ind] = Gb[2][3];
        gzz[ind] = Gb[3][3];

        check_nan_or_inf("gxx", gxx[ind]);
        check_nan_or_inf("gxy", gxy[ind]);
        check_nan_or_inf("gxz", gxz[ind]);
        check_nan_or_inf("gyy", gyy[ind]);
        check_nan_or_inf("gyz", gyz[ind]);
        check_nan_or_inf("gzz", gzz[ind]);

        CCTK_REAL dg[4][4][4]; // dg[i][j][k] = \partial_k g_{ij} (boosted metric)
        // Initialize dg to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            for (int kk = 0; kk < 4; ++kk)
              dg[ii][jj][kk] = 0.0;

        for (int a = 0; a < 4; ++a)
        {
          for (int b = 0; b < 4; ++b)
          {
            for (int c = 0; c < 4; ++c)
            {
              CCTK_REAL sum = 0.0;
              for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                  for (int lam = 0; lam < 4; ++lam)
                    sum += invLambda[mu][a] * invLambda[nu][b] * invLambda[lam][c] * dG[mu][nu][lam];
              dg[a][b][c] = sum;
            }
          }
        }

        // Check for NaN or Inf in all metric derivatives
        for (int ii = 1; ii <= 3; ++ii)
        {
          for (int jj = 1; jj <= 3; ++jj)
          {
            for (int kk = 0; kk <= 3; ++kk)
            {
              char dg_name[32];
              snprintf(dg_name, sizeof(dg_name), "dg[%d][%d][%d]", ii, jj, kk);
              check_nan_or_inf(dg_name, dg[ii][jj][kk]);
            }
          }
        }
        //////////////////////////////////////////////////////////////////////////////////////////

        // // Compute inverse metric g^{ij} (spatial part only)
        // CCTK_REAL det_g =
        //     g[1][1]*(g[2][2]*g[3][3] - g[2][3]*g[3][2])
        //   - g[1][2]*(g[2][1]*g[3][3] - g[2][3]*g[3][1])
        //   + g[1][3]*(g[2][1]*g[3][2] - g[2][2]*g[3][1]);

        // CCTK_REAL g_inv[4][4]; // Inverse metric
        // // Initialize g_inv to zero
        // for (int i = 0; i < 4; ++i)
        //   for (int j = 0; j < 4; ++j)
        //     g_inv[i][j] = 0.0;

        // if (fabs(det_g) < 1e-12) {
        //     // Abort execution due to singular metric
        //     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        //            "Singular spatial metric (det_g = %e) at grid point (%d,%d,%d). Aborting.", det_g, i, j, k);
        //     abort();
        // } else {
        //     g_inv[1][1] =  (g[2][2]*g[3][3] - g[2][3]*g[3][2]) / det_g;
        //     g_inv[1][2] = -(g[1][2]*g[3][3] - g[1][3]*g[3][2]) / det_g;
        //     g_inv[1][3] =  (g[1][2]*g[2][3] - g[1][3]*g[2][2]) / det_g;
        //     g_inv[2][1] = -(g[2][1]*g[3][3] - g[2][3]*g[3][1]) / det_g;
        //     g_inv[2][2] =  (g[1][1]*g[3][3] - g[1][3]*g[3][1]) / det_g;
        //     g_inv[2][3] = -(g[1][1]*g[2][3] - g[1][3]*g[2][1]) / det_g;
        //     g_inv[3][1] =  (g[2][1]*g[3][2] - g[2][2]*g[3][1]) / det_g;
        //     g_inv[3][2] = -(g[1][1]*g[3][2] - g[1][2]*g[3][1]) / det_g;
        //     g_inv[3][3] =  (g[1][1]*g[2][2] - g[1][2]*g[2][1]) / det_g;
        // }

        // // Christoffel symbols of the spatial metric (only spatial indices 1..3)
        // CCTK_REAL Gamma[4][4][4]; // Gamma^i_{jk}
        // // Initialize to zero
        // for (int i = 1; i <= 3; ++i)
        //   for (int j = 1; j <= 3; ++j)
        //     for (int k = 1; k <= 3; ++k)
        //       Gamma[i][j][k] = 0.0;

        // // Compute Christoffel symbols (spatial part)
        // // Gamma^i_{jk} = 0.5 * g^{il} (dg_{lj}/dx^k + dg_{lk}/dx^j - dg_{jk}/dx^l)
        // for (int i = 1; i <= 3; ++i) {
        //   for (int j = 1; j <= 3; ++j) {
        //     for (int k = 1; k <= 3; ++k) {
        //       for (int l = 1; l <= 3; ++l) {
        //         Gamma[i][j][k] += 0.5 * g_inv[i][l] * (dg[l][j][k] + dg[l][k][j] - dg[j][k][l]);
        //       }
        //     }
        //   }
        // }
        // // Check for NaN or Inf in all Christoffel symbols
        // for (int ii = 1; ii <= 3; ++ii) {
        //   for (int jj = 1; jj <= 3; ++jj) {
        //     for (int kk = 1; kk <= 3; ++kk) {
        //       char gamma_name[32];
        //       snprintf(gamma_name, sizeof(gamma_name), "Gamma[%d][%d][%d]", ii, jj, kk);
        //       check_nan_or_inf(gamma_name, Gamma[ii][jj][kk]);
        //     }
        //   }
        // }

        // CCTK_REAL new_betad[4];
        // // Initialize g to zero
        // for (int i = 0; i < 4; ++i)
        //     new_betad[i] = 0.0;

        // new_betad[1] = pow(gamma,2)*(betad[1] + bh_v*(pow(alpha0,2) + bh_v*betad[1] - \
        //                (pow(betad[1],2) + pow(betad[2],2) + pow(bh_spin,2)*pow(y1_2*betad[2] \
        //                + x1_2*betad[1]*gamma,2)*hh)/((1 + \
        //                pow(bh_spin,2)*(rho2_2)*hh)*psi4_2) - (1 + pow(bh_spin,2)*pow(y1_2,2)*hh)*psi4_2));
        // new_betad[2] =  gamma*(betad[2] + pow(bh_spin,2)*bh_v*x1_2*y1_2*gamma*hh*psi4_2);
        // new_betad[3] = 0.;
        // check_nan_or_inf("new_betad[1]", new_betad[1]);
        // check_nan_or_inf("new_betad[2]", new_betad[2]);
        // check_nan_or_inf("new_betad[3]", new_betad[3]);

        // CCTK_REAL new_betaup[4];
        // // Initialize new_betaup to zero
        // for (int i = 0; i < 4; ++i)
        //     new_betaup[i] = 0.0;

        // new_betaup[1] = g_inv[1][1] * new_betad[1] + g_inv[1][2] * new_betad[2] + g_inv[1][3] * new_betad[3];
        // new_betaup[2] = g_inv[2][1] * new_betad[1] + g_inv[2][2] * new_betad[2] + g_inv[2][3] * new_betad[3];
        // new_betaup[3] = g_inv[3][1] * new_betad[1] + g_inv[3][2] * new_betad[2] + g_inv[3][3] * new_betad[3];

        // // Check for NaN or Inf in betaup components
        // check_nan_or_inf("new_betaup[1]", new_betaup[1]);
        // check_nan_or_inf("new_betaup[2]", new_betaup[2]);
        // check_nan_or_inf("new_betaup[3]", new_betaup[3]);

        // // CCTK_REAL new_lapse = sqrt(-g[0][0] + betad[1]*betaup[1] + betad[2]*betaup[2] + betad[3]*betaup[3]);
        // CCTK_REAL lapse_arg = -g[0][0] + new_betad[1]*new_betaup[1] + new_betad[2]*new_betaup[2] + new_betad[3]*new_betaup[3];
        // if (lapse_arg < pow(SMALL, 4)) {
        //     // fprintf(stderr, "Negative argument in sqrt for new_lapse: %.9e\n", lapse_arg);
        //     lapse_arg = 0; // so that is kij become nan my new approximation is lacking something.
        // }
        // CCTK_REAL new_lapse = sqrt(lapse_arg);

        // if (isnan(new_lapse)) {
        // fprintf(stderr, "Error: %s is NaN\n", "new_lapse");
        // fprintf(stderr, "g00 = %.9e \n", g[0][0]);
        // fprintf(stderr, "beta2 = %.9e \n", new_betad[1]*new_betaup[1] + new_betad[2]*new_betaup[2] + new_betad[3]*new_betaup[3]);
        // // fprintf(stderr, "Error: new_lapse is nan at grid point (%d,%d,%d)\n", x[CCTK_GFINDEX3D (cctkGH, i, j, k)], y[CCTK_GFINDEX3D (cctkGH, i, j, k)], z[CCTK_GFINDEX3D (cctkGH, i, j, k)]);
        // fprintf(stderr, "Error: new_lapse is nan at grid point (%lf,%lf,%lf)\n", x1_2, y1_2, z1_2);

        // abort(); // Break execution
        // }
        // // if ((rr_2 < horizon_radius + 1e-5) && (rr_2 > horizon_radius - 1e-5 )) {
        // //   fprintf(stderr, "Warning: new_lapse at grid point (%lf,%lf,%lf)\n", x1_2, y1_2, z1_2);
        // //   fprintf(stderr, "new_lapse = %e\n", new_lapse);
        // //   // abort(); // Break execution
        // // }

        CCTK_REAL dW_drho_1, dW_dz_1;
        const CCTK_REAL exp_auxi_1 = exp(2. * F2_1[ind] - F0_1[ind]);

        if (rho_1 < 1e-8)
        {
          dW_drho_1 = 0.;
          dW_dz_1 = 0.;
        }
        else
        {
          dW_drho_1 = rho_1 / rr_1 * dW_dr_1[ind] + z1_1 / rr2_1 * dW_dth_1[ind];
          dW_dz_1 = z1_1 / rr_1 * dW_dr_1[ind] - rho_1 / rr2_1 * dW_dth_1[ind];
        }

        // // Compute covariant derivatives of the shift
        // CCTK_REAL Dbetad[4][4];
        // // Initialize Dbetad to zero
        // for (int i = 0; i < 4; ++i)
        //   for (int j = 0; j < 4; ++j)
        //     Dbetad[i][j] = 0.0;

        // for (int i = 1; i <= 3; ++i) {
        //   for (int j = 1; j <= 3; ++j) {
        //     Dbetad[i][j] = dbetad[i][j]; //aqui deve ser a derivada dos shifts depois do boost.
        //     for (int k = 1; k <= 3; ++k)
        //       Dbetad[i][j] -= Gamma[k][i][j] * betad[k]; // e aqui também.
        //     // Check for NaN or Inf in Dbetad
        //     char Dbetad_name[32];
        //     snprintf(Dbetad_name, sizeof(Dbetad_name), "Dbetad[%d][%d]", i, j);
        //     check_nan_or_inf(Dbetad_name, Dbetad[i][j]);
        //   }
        // }

        // Compute extrinsic curvature K_{ij}
        //  if (rr_2 < horizon_radius + eps_r || rr_2 > horizon_radius - eps_r) {
        //    // Example: recompute with new coordinates (replace new_x1_2, etc. with your values)
        //    KerrVars hor = compute_kerr_vars(x1_2, y1_2, z1_2, bh_v, gamma, eps_r, bh_mass, bh_spin, bh_mass2, bh_spin2, horizon_radius);
        //    // Use new_vars.rr_2, new_vars.rho_2, etc. as needed

        //   CCTK_REAL dBB_dr[4][4]; // dBB[i][j] = \partial_r B_{ij} indices refer to spherical coordinates
        //   // Initialize dBB to zero
        //   for (int ii = 0; ii < 4; ++ii)
        //     for (int jj = 0; jj < 4; ++jj)
        //       dBB_dr[ii][jj] = 0.0;

        //   dBB_dr[1][1] = ;
        //   dBB_dr[1][2] = ;
        //   dBB_dr[1][3] = ;
        //   dBB_dr[2][2] = ;
        //   dBB_dr[2][3] = ;
        //   dBB_dr[3][3] = ;

        //   const CCTK_REAL dnew_lapse_dr = r*sqrt((pow(-1 + pow(bh_v,2),2)*pow(gamma,2)*pow(hor.psi4_2,2)*pow(hor.sinth,2)*pow(1 \
        //                                   + pow(bh_spin,2)*pow(r,2)*hor.hh*pow(hor.sinth,2),2))/pow(hor.psi4_2*(\
        //                                   r*hor.sinth + pow(bh_spin,2)*pow(r,3)*hor.hh*pow(hor.sinth,3)) + \
        //                                   bh_v*hor.bphi*hor.sinph,2))*(hor.costh*hor.dalpha_dz + \
        //                                   hor.sinth*(hor.sinph*hor.dalpha_dy + hor.cosph*hor.dalpha_dx));

        //   const CCTK_REAL bxx[ind] = ;
        //   const CCTK_REAL bxy[ind] = ;
        //   const CCTK_REAL bxz[ind] = ;
        //   const CCTK_REAL byy[ind] = ;
        //   const CCTK_REAL byz[ind] = ;
        //   const CCTK_REAL bzz[ind] = ;

        //   kxx[ind] = ;
        //   kxy[ind] = ;
        //   kxz[ind] = ;
        //   kyy[ind] = ;
        //   kyz[ind] = ;
        //   kzz[ind] = ;

        // } else{

        const CCTK_REAL HF = -bh_spin2 * bh_spin * alpha0 * sigma / rhokerr * costh; // we are dividing by sinth2
        const CCTK_REAL Athph = HF / rr_2;                                           // we are dividing by sinth
        const CCTK_REAL aux = rho2kerr * (rBL * rBL - bh_spin2) + 2. * rBL * rBL * (rBL * rBL + bh_spin2);
        const CCTK_REAL HE = bh_spin * bh_mass * aux / (rhokerr * rhokerr * rhokerr) *
                             1. / sqrt(rBL * rBL + bh_spin2 * (1. + sigma * sinth2));
        const CCTK_REAL ARph = HE / rr2_2; // we are dividing by sinth2

        // capital Ks refer to the unboosted frame. Original quantities.

        const CCTK_REAL dbeta2_dx = -((pow(bh_spin, 2) * sigma * (rr_2 * (rho2_2) * (1 + pow(bh_spin, 2) * hh * (rho2_2)) * sigma * dpsi4_2_dx + psi4_2 * (sigma * (2 * (-2 * pow(z1_2, 2) + pow(rr_2, 2) + 2 * pow(bh_spin, 2) * hh * pow(-rho2_2, 2)) * R_x + pow(bh_spin, 2) * rr_2 * pow(-rho2_2, 2) * dhh_dx) - 2 * rr_2 * (rho2_2) * (1 + pow(bh_spin, 2) * hh * (-pow(z1_2, 2) + pow(rr_2, 2))) * dsigma_dx))) / (pow(psi4_2, 2) * pow(rr_2, 5) * pow(1 + pow(bh_spin, 2) * hh * (rho2_2), 2)));
        const CCTK_REAL dbeta2_dy = -((pow(bh_spin, 2) * sigma * (rr_2 * (rho2_2) * (1 + pow(bh_spin, 2) * hh * (rho2_2)) * sigma * dpsi4_2_dy + psi4_2 * (sigma * (2 * (-2 * pow(z1_2, 2) + pow(rr_2, 2) + 2 * pow(bh_spin, 2) * hh * pow(-rho2_2, 2)) * R_y + pow(bh_spin, 2) * rr_2 * pow(-rho2_2, 2) * dhh_dy) - 2 * rr_2 * (rho2_2) * (1 + pow(bh_spin, 2) * hh * (-pow(z1_2, 2) + pow(rr_2, 2))) * dsigma_dy))) / (pow(psi4_2, 2) * pow(rr_2, 5) * pow(1 + pow(bh_spin, 2) * hh * (rho2_2), 2)));
        const CCTK_REAL dbeta2_dz = (pow(bh_spin, 2) * sigma * (sigma * (psi4_2 * (-2 * z1_2 * rr_2 - 2 * (-2 * pow(z1_2, 2) + pow(rr_2, 2) + 2 * pow(bh_spin, 2) * hh * pow(-rho2_2, 2)) * R_z - pow(bh_spin, 2) * rr_2 * pow(-rho2_2, 2) * dhh_dz) - rr_2 * (rho2_2) * (1 + pow(bh_spin, 2) * hh * (rho2_2)) * dpsi4_2_dz) + 2 * psi4_2 * rr_2 * (rho2_2) * (1 + pow(bh_spin, 2) * hh * (rho2_2)) * dsigma_dz)) / (pow(psi4_2, 2) * pow(rr_2, 5) * pow(1 + pow(bh_spin, 2) * hh * (rho2_2), 2));
        check_nan_or_inf("dbeta2_dx", dbeta2_dx);
        check_nan_or_inf("dbeta2_dy", dbeta2_dy);
        check_nan_or_inf("dbeta2_dz", dbeta2_dz);

        const CCTK_REAL dbetadphi_dth = -(4 * bh_spin * bh_mass * rBL * (bh_spin2 + rBL * rBL) * sinth * costh) / pow(rho2kerr, 2);
        const CCTK_REAL dbetadphi_dR = dr_dR * 2 * bh_spin * bh_mass * (rBL * rBL - bh_spin2 * costh2) * sinth2 / pow(rho2kerr, 2);
        const CCTK_REAL delta_metric = rBL * rBL - 2 * bh_mass * rBL + bh_spin2;
        // const CCTK_REAL gammaphiphi= psi4_2*rr2_2*sinth2*(1 + bh_spin2*hh*rr2_2*sinth2);
        // const CCTK_REAL dgammaphiphi_dth= -2*bh_spin2*delta_metric*sinth*costh/rr2_2;
        const CCTK_REAL dgammaphiphi_dth = (delta_metric + 8 * bh_mass * rBL * pow(bh_spin2 + rBL * rBL, 2) / pow(bh_spin2 + 2 * rBL * rBL + bh_spin2 * (costh2 - sinth2), 2)) * 2 * costh * sinth;
        const CCTK_REAL dgammaphiphi_dR = dr_dR * 2 * (rr_2 * (2 * rBL * (bh_spin2 + rBL * rBL) + bh_spin2 * (bh_mass - rBL) * sinth2)) / (rr2_2 * rr_2) +
                                          2 * (-pow(bh_spin2 + rBL * rBL, 2) + bh_spin2 * delta_metric * sinth2) / (rr2_2 * rr_2);

        const CCTK_REAL dbetauphi_dth = (gammaphiphi * dbetadphi_dth - bphi * dgammaphiphi_dth) / pow(gammaphiphi, 2);
        const CCTK_REAL dbetauphi_dR = (gammaphiphi * dbetadphi_dR - bphi * dgammaphiphi_dR) / pow(gammaphiphi, 2);

        const CCTK_REAL Ktt = 0.5 * (dbeta2_dx * betaup[1] + dbeta2_dy * betaup[2] + dbeta2_dz * betaup[3]) / alpha0; // time derivatives here are zero since we are computing the original Ktt. nevertheless its argument now for the computation of the new quatities is gamma x.
        const CCTK_REAL Ktht = bphi * dbetauphi_dth / (-2 * alpha0);
        const CCTK_REAL KRt = bphi * dbetauphi_dR / (-2 * alpha0);

        const CCTK_REAL Kxt = R_x * KRt + x1_2 * gamma * z1_2 / (rho_2 * rr2_2) * Ktht;
        const CCTK_REAL Kyt = R_y * KRt + y1_2 * z1_2 / (rho_2 * rr2_2) * Ktht;
        const CCTK_REAL Kzt = R_z * KRt - rho_2 / rr2_2 * Ktht;

        const CCTK_REAL Axx = 2. * ARph * R_x * sinth2ph_x + 2. * Athph * sinthth_x * sinth2ph_x;
        const CCTK_REAL Axy = ARph * (R_x * sinth2ph_y + R_y * sinth2ph_x) + Athph * (sinthth_x * sinth2ph_y + sinthth_y * sinth2ph_x);
        const CCTK_REAL Axz = ARph * R_z * sinth2ph_x + Athph * sinthth_z * sinth2ph_x;
        const CCTK_REAL Ayy = 2. * ARph * R_y * sinth2ph_y + 2. * Athph * sinthth_y * sinth2ph_y;
        const CCTK_REAL Ayz = ARph * R_z * sinth2ph_y + Athph * sinthth_z * sinth2ph_y;
        // aparentemente nao preciso de mexer aqui com gammas no K de kerr original. só depois é que o argumento de altera.

        CCTK_REAL first_term[4][4];
        // Initialize first_term to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            first_term[ii][jj] = 0.0;

        first_term[1][1] = gamma2 * Axx / psi2_2 + bh_v2 * gamma2 * Ktt + 2 * Kxt * bh_v * gamma2;
        first_term[1][2] = gamma * Axy / psi2_2 + gamma * bh_v * Kyt;
        first_term[1][3] = gamma * Axz / psi2_2 + gamma * bh_v * Kzt;
        first_term[2][1] = first_term[1][2]; // symmetric component;
        first_term[2][2] = Ayy / psi2_2;
        first_term[2][3] = Ayz / psi2_2;
        first_term[3][1] = first_term[1][3]; // symmetric component;
        first_term[3][2] = first_term[2][3]; // symmetric component;
        first_term[3][3] = 0.0;

        CCTK_REAL second_term[4][4];
        // Initialize second_term to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            second_term[ii][jj] = 0.0;
        second_term[1][1] = -0.5 * bh_v * betaup[1] / alpha0 * dg[1][1][0];
        second_term[1][2] = -0.5 * bh_v * betaup[1] / alpha0 * dg[1][2][0];
        second_term[1][3] = -0.5 * bh_v * betaup[1] / alpha0 * dg[1][3][0];
        second_term[2][1] = second_term[1][2]; // symmetric component;
        second_term[2][2] = -0.5 * bh_v * betaup[1] / alpha0 * dg[2][2][0];
        second_term[2][3] = -0.5 * bh_v * betaup[1] / alpha0 * dg[2][3][0];
        second_term[3][1] = second_term[1][3]; // symmetric component;
        second_term[3][2] = second_term[2][3]; // symmetric component;
        second_term[3][3] = -0.5 * bh_v * betaup[1] / alpha0 * dg[3][3][0];

        CCTK_REAL third_term[4][4];
        // Initialize third_term to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            third_term[ii][jj] = 0.0;

        third_term[1][1] = 0.5 * bh_v / alpha0 * dg[1][1][1];
        third_term[1][2] = 0.5 * bh_v / alpha0 * dg[1][2][1];
        third_term[1][3] = 0.5 * bh_v / alpha0 * dg[1][3][1];
        third_term[2][1] = third_term[1][2]; // symmetric component
        third_term[2][2] = 0.5 * bh_v / alpha0 * dg[2][2][1];
        third_term[2][3] = 0.5 * bh_v / alpha0 * dg[2][3][1];
        third_term[3][1] = third_term[1][3]; // symmetric component
        third_term[3][2] = third_term[2][3]; // symmetric component
        third_term[3][3] = 0.5 * bh_v / alpha0 * dg[3][3][1];

        CCTK_REAL forth_term[4][4];
        // Initialize third_term to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            forth_term[ii][jj] = 0.0;

        forth_term[1][1] = -0.5 * betaup[2] / alpha0 * (gamma - 1) * dg[1][1][2];
        forth_term[1][2] = -0.5 * betaup[2] / alpha0 * (gamma - 1) * dg[1][2][2];
        forth_term[1][3] = -0.5 * betaup[2] / alpha0 * (gamma - 1) * dg[1][3][2];
        forth_term[2][1] = forth_term[1][2]; // symmetric component
        forth_term[2][2] = -0.5 * betaup[2] / alpha0 * (gamma - 1) * dg[2][2][2];
        forth_term[2][3] = -0.5 * betaup[2] / alpha0 * (gamma - 1) * dg[2][3][2];
        forth_term[3][1] = forth_term[1][3]; // symmetric component
        forth_term[3][2] = forth_term[2][3]; // symmetric component
        forth_term[3][3] = -0.5 * betaup[2] / alpha0 * (gamma - 1) * dg[3][3][2];

        CCTK_REAL fifth_term[4][4];
        // Initialize third_term to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            fifth_term[ii][jj] = 0.0;

        fifth_term[1][1] = -0.5 * betaup[3] / alpha0 * (gamma - 1) * dg[1][1][3];
        fifth_term[1][2] = -0.5 * betaup[3] / alpha0 * (gamma - 1) * dg[1][2][3];
        fifth_term[1][3] = -0.5 * betaup[3] / alpha0 * (gamma - 1) * dg[1][3][3];
        fifth_term[2][1] = fifth_term[1][2]; // symmetric component
        fifth_term[2][2] = -0.5 * betaup[3] / alpha0 * (gamma - 1) * dg[2][2][3];
        fifth_term[2][3] = -0.5 * betaup[3] / alpha0 * (gamma - 1) * dg[2][3][3];
        fifth_term[3][1] = fifth_term[1][3]; // symmetric component
        fifth_term[3][2] = fifth_term[2][3]; // symmetric component
        fifth_term[3][3] = -0.5 * betaup[3] / alpha0 * (gamma - 1) * dg[3][3][3];

        kxx[ind] = gamma * (first_term[1][1] + second_term[1][1] + third_term[1][1]) + forth_term[1][1] + fifth_term[1][1];
        kxy[ind] = gamma * (first_term[1][2] + second_term[1][2] + third_term[1][2]) + forth_term[1][2] + fifth_term[1][2];
        kxz[ind] = gamma * (first_term[1][3] + second_term[1][3] + third_term[1][3]) + forth_term[1][3] + fifth_term[1][3];
        kyy[ind] = gamma * (first_term[2][2] + second_term[2][2] + third_term[2][2]) + forth_term[2][2] + fifth_term[2][2];
        kyz[ind] = gamma * (first_term[2][3] + second_term[2][3] + third_term[2][3]) + forth_term[2][3] + fifth_term[2][3];
        kzz[ind] = gamma * (first_term[3][3] + second_term[3][3] + third_term[3][3]) + forth_term[3][3] + fifth_term[3][3];

        // }

        check_nan_or_inf("kxx", kxx[ind]);
        check_nan_or_inf("kxy", kxy[ind]);
        check_nan_or_inf("kxz", kxz[ind]);
        check_nan_or_inf("kyy", kyy[ind]);
        check_nan_or_inf("kyz", kyz[ind]);
        check_nan_or_inf("kzz", kzz[ind]);

        const CCTK_REAL phi0_l_1 = phi0_1[ind]; // * pert_phi_1;
        const CCTK_REAL phi0_l_2 = 0.;

        // scalar fields

        // star 1
        CCTK_REAL phi1_1 = phi0_l_1 * (coswt * cosmph_1 + sinwt * sinmph_1);
        CCTK_REAL phi2_1 = phi0_l_1 * (coswt * sinmph_1 - sinwt * cosmph_1);

        // BH 2
        CCTK_REAL phi1_2 = 0; // phi0_l_2 * (coswt * cosmph_2 + sinwt * sinmph_2);
        CCTK_REAL phi2_2 = 0; // phi0_l_2 * (coswt * sinmph_2 - sinwt * cosmph_2);
        /////////////////////////////

        phi1[ind] = phi1_1 + phi1_2;
        phi2[ind] = phi2_1 + phi2_2;

        const CCTK_REAL alph_1 = exp(F0_1[ind]);
        const CCTK_REAL alph_2 = 1.0;

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
        else if (CCTK_EQUALS(initial_lapse, "Kerr_BS"))
        {
          alp[ind] = alph_1 + alph_2 - 1;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // CCTK_REAL Delt  = rBL*rBL + bh_spin2 - 2 * bh_mass * rBL;
        // CCTK_REAL fctFF = ( rBL*rBL + bh_spin2 ) * ( rBL*rBL + bh_spin2 ) - Delt * bh_spin2 * sinth2;
        // bphi = 2.0 * bh_spin * bh_mass * rBL / fctFF;

        // shift
        if (CCTK_EQUALS(initial_shift, "Kerr_BS"))
        {
          betax[ind] = betad[1]; // por enquato o shift da bs é zero pois é estática.
          betay[ind] = betad[2];
          betaz[ind] = 0.;
        }

      } /* for i */
    } /* for j */
  } /* for k */

  free(F1_1);
  free(F2_1);
  free(F0_1);
  free(phi0_1);
  free(W_1);
  free(dW_dr_1);
  free(dW_dth_1);

  // free(F1_2); free(F2_2); free(F0_2); free(phi0_2); free(W_2);
  // free(dW_dr_2); free(dW_dth_2);

  return;
}
