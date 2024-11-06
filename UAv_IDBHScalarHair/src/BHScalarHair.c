
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

void UAv_IDBHScalarHair_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL []);


void UAv_IDBHScalarHair(CCTK_ARGUMENTS)
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
  UAv_IDBHScalarHair_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, phi0_in, Wbar_in);

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


  // now we need to take the derivatives of the Wbar function and store their values

  CCTK_REAL *dWbar_dr_in, *dWbar_dth_in, *d2Wbar_drth_in;
  dWbar_dr_in    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dWbar_dth_in   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  d2Wbar_drth_in = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  const CCTK_REAL oodX       = 1. / dX;
  const CCTK_REAL oodXsq     = oodX * oodX;
  const CCTK_REAL oodX12     = 1. / (12. * dX);
  const CCTK_REAL oodth12    = 1. / (12. * dtheta);
  const CCTK_REAL oodXdth4   = 1. / (4.  * dX * dtheta);
  const CCTK_REAL oodXdth144 = 1. / (144. * dX * dtheta);
  const CCTK_REAL oodXsqdth2 = 1. / (2.  * dX * dX * dtheta);

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

      const CCTK_INT indim2jm1 = i-2 + jm1*NX;
      const CCTK_INT indim2jm2 = i-2 + jm2*NX;
      const CCTK_INT indim2jp1 = i-2 + jp1*NX;
      const CCTK_INT indim2jp2 = i-2 + jp2*NX;

      const CCTK_INT indim1jm1 = i-1 + jm1*NX;
      const CCTK_INT indim1jm2 = i-1 + jm2*NX;
      const CCTK_INT indim1jp1 = i-1 + jp1*NX;
      const CCTK_INT indim1jp2 = i-1 + jp2*NX;

      const CCTK_INT indjm1 = i + jm1*NX;
      const CCTK_INT indjm2 = i + jm2*NX;
      const CCTK_INT indjp1 = i + jp1*NX;
      const CCTK_INT indjp2 = i + jp2*NX;

      const CCTK_INT indip1jm1 = i+1 + jm1*NX;
      const CCTK_INT indip1jm2 = i+1 + jm2*NX;
      const CCTK_INT indip1jp1 = i+1 + jp1*NX;
      const CCTK_INT indip1jp2 = i+1 + jp2*NX;

      const CCTK_INT indip2jm1 = i+2 + jm1*NX;
      const CCTK_INT indip2jm2 = i+2 + jm2*NX;
      const CCTK_INT indip2jp1 = i+2 + jp1*NX;
      const CCTK_INT indip2jp2 = i+2 + jp2*NX;
      
      const CCTK_INT indip3jm1 = i+3 + jm1*NX;
      //const CCTK_INT indip3jm2 = i+3 + jm2*NX;
      const CCTK_INT indip3jp1 = i+3 + jp1*NX;
      //const CCTK_INT indip3jp2 = i+3 + jp2*NX;


      const CCTK_REAL lX = X[i];
      /* const CCTK_REAL lth = theta[j]; */
      /* printf("X[%3d] = %lf\n", i, lX); */


      // 1st derivative with 4th order accuracy (central stencils)
      const CCTK_REAL Wbar_th = (-Wbar_in[indjp2] + 8 * Wbar_in[indjp1] - 8 * Wbar_in[indjm1] + Wbar_in[indjm2]) *
        oodth12;

      CCTK_REAL Wbar_X, Wbar_Xth;

      if (i == 1 || i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;

        Wbar_Xth  = ( Wbar_in[indip1jp1] - Wbar_in[indip1jm1] - Wbar_in[indim1jp1] + Wbar_in[indim1jm1] ) * oodXdth4;

      } else if (i == 0) {
        /* this point is X == 0, r == rH, R == rH/4. dWbar_dX goes to zero here. but
           since we're interested in dWbar_dr, and since drxdr diverges (here), we
           will use L'Hopital's rule. for that, we will write instead the 2nd
           derivative */

        // 2nd derivative with 2nd order accuracy (forward stencils)
        Wbar_X = (2*Wbar_in[ind] - 5*Wbar_in[indip1] + 4*Wbar_in[indip2] - Wbar_in[indip3]) * oodXsq;

        // mixed derivatives with 2nd order accuracy (central stencils in j (1st der) and forward in i (2nd der))
        Wbar_Xth = ( 2*Wbar_in[indjp1] - 2*Wbar_in[indjm1] - 5*Wbar_in[indip1jp1] + 5*Wbar_in[indip1jm1]
                    + 4*Wbar_in[indip2jp1] - 4*Wbar_in[indip2jm1] - Wbar_in[indip3jp1] + Wbar_in[indip3jm1]) * oodXsqdth2;
        
      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4*Wbar_in[indim1] + 3*Wbar_in[ind]) * 0.5 * oodX;
        Wbar_Xth = 0.; // we don't actually use this variable at large r, so just
                    // set it to zero

      } else {
        // 4th order accurate stencils
        Wbar_X    = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
        Wbar_Xth  = (
            -Wbar_in[indim2jp2] +  8*Wbar_in[indim1jp2] -  8*Wbar_in[indip1jp2] +   Wbar_in[indip2jp2]
         + 8*Wbar_in[indim2jp1] - 64*Wbar_in[indim1jp1] + 64*Wbar_in[indip1jp1] - 8*Wbar_in[indip2jp1]
         - 8*Wbar_in[indim2jm1] + 64*Wbar_in[indim1jm1] - 64*Wbar_in[indip1jm1] + 8*Wbar_in[indip2jm1]
         +   Wbar_in[indim2jm2] -  8*Wbar_in[indim1jm2] +  8*Wbar_in[indip1jm2] -   Wbar_in[indip2jm2] ) * oodXdth144;

      }

      // from the X coordinate used in the input files to the x coordinate
      // We need to be careful at X == 1 for radial derivatives (coordinate change is singular)
      if (i == NX - 1) {
        dWbar_dr_in[ind]    = 0.; // Sensibly, dWbar_dr_in should vanish (dXdr == 0, Wbar_X bounded)
        d2Wbar_drth_in[ind] = 0.; // Wbar_Xth is set to 0 above anyway

      } else {
        const CCTK_REAL rx = C0*lX/(1. - lX);
        // from the x coordinate to the metric coordinate r
        // const CCTK_REAL rr = sqrt(rH*rH + rx*rx);

        // corresponding derivatives
        const CCTK_REAL dXdrx = 1./(C0 + rx) - rx/((C0 + rx)*(C0 + rx));

        CCTK_REAL drxdr;
        if (i == 0) { // rx == 0 (X == 0)
          drxdr = sqrt(rH*rH + rx*rx);
        } else {
          drxdr = sqrt(rH*rH + rx*rx)/rx;
        }
        const CCTK_REAL dXdr = dXdrx * drxdr;

        dWbar_dr_in[ind]    = dXdr * Wbar_X;
        d2Wbar_drth_in[ind] = dXdr * Wbar_Xth;
      }

      dWbar_dth_in[ind]   = Wbar_th;
    }
  }


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

        CCTK_REAL ltheta = acos( z1/RR );
        if (ltheta > 0.5*M_PI)    // symmetry along the equatorial plane
          ltheta = M_PI - ltheta;

        X_g[ind]     = lX;
        theta_g[ind] = ltheta;
      }
    }
  }

  /* now for the interpolation */

  const CCTK_INT N_dims  = 2;   // 2-D interpolation

  const CCTK_INT N_input_arrays  = 8;
  const CCTK_INT N_output_arrays = 8;

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
  input_array_dims[1] = Ntheta;

  input_array_type_codes[0] = CCTK_VARIABLE_REAL;
  input_array_type_codes[1] = CCTK_VARIABLE_REAL;
  input_array_type_codes[2] = CCTK_VARIABLE_REAL;
  input_array_type_codes[3] = CCTK_VARIABLE_REAL;
  input_array_type_codes[4] = CCTK_VARIABLE_REAL;
  input_array_type_codes[5] = CCTK_VARIABLE_REAL;
  input_array_type_codes[6] = CCTK_VARIABLE_REAL;
  input_array_type_codes[7] = CCTK_VARIABLE_REAL;

  /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
  input_arrays[0] = (const void *) F1_in;
  input_arrays[1] = (const void *) F2_in;
  input_arrays[2] = (const void *) F0_in;
  input_arrays[3] = (const void *) phi0_in;
  input_arrays[4] = (const void *) Wbar_in;
  input_arrays[5] = (const void *) dWbar_dr_in;
  input_arrays[6] = (const void *) dWbar_dth_in;
  input_arrays[7] = (const void *) d2Wbar_drth_in;

  /* output arrays */
  void *output_arrays[N_output_arrays];
  CCTK_INT output_array_type_codes[N_output_arrays];
  CCTK_REAL *F1, *F2, *F0, *phi0, *Wbar;
  CCTK_REAL *dWbar_dr, *dWbar_dth, *d2Wbar_drth;

  F1          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F2          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F0          = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  phi0        = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  Wbar        = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dWbar_dr    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dWbar_dth   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  d2Wbar_drth = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  output_array_type_codes[0] = CCTK_VARIABLE_REAL;
  output_array_type_codes[1] = CCTK_VARIABLE_REAL;
  output_array_type_codes[2] = CCTK_VARIABLE_REAL;
  output_array_type_codes[3] = CCTK_VARIABLE_REAL;
  output_array_type_codes[4] = CCTK_VARIABLE_REAL;
  output_array_type_codes[5] = CCTK_VARIABLE_REAL;
  output_array_type_codes[6] = CCTK_VARIABLE_REAL;
  output_array_type_codes[7] = CCTK_VARIABLE_REAL;

  output_arrays[0] = (void *) F1;
  output_arrays[1] = (void *) F2;
  output_arrays[2] = (void *) F0;
  output_arrays[3] = (void *) phi0;
  output_arrays[4] = (void *) Wbar;
  output_arrays[5] = (void *) dWbar_dr;
  output_arrays[6] = (void *) dWbar_dth;
  output_arrays[7] = (void *) d2Wbar_drth;


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
  free(F1_in); free(F2_in); free(F0_in); free(phi0_in); free(Wbar_in);
  free(dWbar_dr_in); free(dWbar_dth_in); free(d2Wbar_drth_in);


  /* printf("F1 = %g\n", F1[0]); */
  /* printf("F2 = %g\n", F2[0]); */
  /* printf("F0 = %g\n", F0[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W = %g\n", W[0]); */


  /* now we finally write the metric and all 3+1 quantities. first we write the
     3-metric, lapse and scalar fields */

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

        // from Wbar to W function, Wbar = r^p * W
        const CCTK_REAL W        = Wbar[ind] / rrP;
        const CCTK_REAL dW_dth   = dWbar_dth[ind] / rrP;
        const CCTK_REAL dW_dr    = dWbar_dr[ind] / rrP - Wbar_r_power * Wbar[ind] / rrPp1;
        const CCTK_REAL d2W_drth = d2Wbar_drth[ind] / rrP - Wbar_r_power * dWbar_dth[ind] / rrPp1;

        // add non-axisymmetric perturbation on conformal factor
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

          

        // let's add a perturbation to the scalar field as well
        const CCTK_REAL argpert_phi = (RR - R0pert_phi)/Sigmapert_phi;
        const CCTK_REAL pert_phi = 1. + Apert_phi * (x1*x1 - y1*y1)*mu*mu * exp( -0.5*argpert_phi*argpert_phi );

        const CCTK_REAL phi0_l = phi0[ind] * pert_phi;

        // scalar fields
        phi1[ind]  = phi0_l * (coswt * cosmph + sinwt * sinmph);
        phi2[ind]  = phi0_l * (coswt * sinmph - sinwt * cosmph);

        const CCTK_REAL alph = exp(F0[ind]) * (RR - 0.25*rH) / (RR + 0.25*rH);

        // if at R ~ rH/4 we need to regularize the division by R - rH/4
        // That's the same as for the extrinsic curvature above, with f(R) = W - OmegaH
        if (fabs(eps) < 1e-4) {
          const CCTK_REAL reg  = exp(-F0[ind]) * (RR + 0.25*rH) * (eps_o_1meps - eps_o_1meps*eps_o_1meps) * dW_dr;

          Kphi1[ind] =  0.5 * mm * reg * phi2[ind];
          Kphi2[ind] = -0.5 * mm * reg * phi1[ind];
        } else {
          Kphi1[ind] = 0.5 * mm * (W - OmegaH) / alph * phi2[ind];
          Kphi2[ind] = 0.5 * mm * (OmegaH - W) / alph * phi1[ind];
        }

        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "HairyBH")) {
          alp[ind] = alph;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // shift
        if (CCTK_EQUALS(initial_shift, "HairyBH")) {
          betax[ind] =  W * y1;
          betay[ind] = -W * x1;
          betaz[ind] =  0.;
        }

      } /* for i */
    }   /* for j */
  }     /* for k */

  free(F1); free(F2); free(F0); free(phi0); free(Wbar);
  free(dWbar_dr); free(dWbar_dth); free(d2Wbar_drth);

  return;
}
