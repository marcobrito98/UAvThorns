// este e o codigo copiado do kerrnewman.
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include <math.h>
#include <mpfr.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define SMALL (1.e-9)
#define MPFR_PREC 1024               // ~308 decimal digits
#define HORIZON_PROXIMITY_FACTOR 0.2 // compute with high precision within 20% of horizon

void UAv_ID_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL[], CCTK_REAL[],
                      CCTK_REAL[], CCTK_REAL[], CCTK_REAL[], CCTK_REAL[],
                      CCTK_REAL[]);

void check_nan_or_inf(const char *var_name, double value) {
  if (isnan(value)) {
    fprintf(stderr, "Error: %s is NaN\n", var_name);
    abort(); // Break execution
  } else if (isinf(value)) {
    fprintf(stderr, "Error: %s is Inf\n", var_name);
    abort(); // Break execution
  }
}

// Numerically stable inversion of a symmetric positive-definite 3x3 matrix
// (indices 1..3) Uses diagonal scaling and long-double accumulators, then
// Cholesky: M = S^{-1} A S^{-1}, A = L L^T Returns 1 on success, 0 on failure
// (non-SPD or singular).
static inline int invert_spd3x3(const CCTK_REAL M[4][4], CCTK_REAL Minv[4][4]) {
  // Enforce symmetry (guards against tiny asymmetries)
  long double m11 = (long double)M[1][1];
  long double m22 = (long double)M[2][2];
  long double m33 = (long double)M[3][3];
  long double m12 = 0.5L * ((long double)M[1][2] + (long double)M[2][1]);
  long double m13 = 0.5L * ((long double)M[1][3] + (long double)M[3][1]);
  long double m23 = 0.5L * ((long double)M[2][3] + (long double)M[3][2]);

  // Diagonal scaling to improve conditioning: A = S M S, S = diag(1/sqrt(mii))
  const long double eps = 1e-300L;
  if (!(m11 > 0 && m22 > 0 && m33 > 0))
    return 0;
  long double s1 = 1.0L / sqrtl(fmaxl(m11, eps));
  long double s2 = 1.0L / sqrtl(fmaxl(m22, eps));
  long double s3 = 1.0L / sqrtl(fmaxl(m33, eps));

  long double a11 = s1 * s1 * m11;
  long double a22 = s2 * s2 * m22;
  long double a33 = s3 * s3 * m33;
  long double a12 = s1 * s2 * m12;
  long double a13 = s1 * s3 * m13;
  long double a23 = s2 * s3 * m23;

  // Cholesky A = L L^T (lower triangular L)
  if (!(a11 > 0))
    return 0;
  long double L11 = sqrtl(a11);

  long double L21 = a12 / L11;
  long double L31 = a13 / L11;

  long double t22 = a22 - L21 * L21;
  if (!(t22 > 0))
    return 0;
  long double L22 = sqrtl(t22);

  long double L32 = (a23 - L31 * L21) / L22;

  long double t33 = a33 - L31 * L31 - L32 * L32;
  if (!(t33 > 0))
    return 0;
  long double L33 = sqrtl(t33);

  // Invert L (lower triangular): compute Linv so that Linv * L = I
  long double Linv11 = 1.0L / L11;
  long double Linv21 = -L21 * (Linv11 / L22);
  long double Linv22 = 1.0L / L22;
  long double Linv31 = -(L31 * Linv11 + L32 * Linv21) / L33;
  long double Linv32 = -L32 * (Linv22 / L33);
  long double Linv33 = 1.0L / L33;

  // A^{-1} = (L^{-T} L^{-1}) = (Linv^T * Linv)
  long double i11 = Linv11 * Linv11 + Linv21 * Linv21 + Linv31 * Linv31;
  long double i12 =
      Linv21 * Linv22 + Linv31 * Linv32 + Linv11 * 0.0L; // explicit for clarity
  long double i13 = Linv31 * Linv33 + 0.0L;              // since Linv is lower
  long double i22 = Linv22 * Linv22 + Linv32 * Linv32;
  long double i23 = Linv32 * Linv33;
  long double i33 = Linv33 * Linv33;

  // Undo scaling: M^{-1} = S * A^{-1} * S
  long double S1 = s1, S2 = s2, S3 = s3;
  long double mInv11 = S1 * S1 * i11;
  long double mInv12 = S1 * S2 * i12;
  long double mInv13 = S1 * S3 * i13;
  long double mInv22 = S2 * S2 * i22;
  long double mInv23 = S2 * S3 * i23;
  long double mInv33 = S3 * S3 * i33;

  Minv[1][1] = (CCTK_REAL)mInv11;
  Minv[1][2] = (CCTK_REAL)mInv12;
  Minv[1][3] = (CCTK_REAL)mInv13;
  Minv[2][1] = Minv[1][2];
  Minv[2][2] = (CCTK_REAL)mInv22;
  Minv[2][3] = (CCTK_REAL)mInv23;
  Minv[3][1] = Minv[1][3];
  Minv[3][2] = Minv[2][3];
  Minv[3][3] = (CCTK_REAL)mInv33;

  return 1;
}

void UAv_ID_Kerr_BS(CCTK_ARGUMENTS) {
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
  for (int i = 0; i < NX; i++) {
    X[i] = Xtmp[i];
    /* printf("X[%3d] = %lf\n", i, X[i]); */
  }
  for (int i = 0; i < Ntheta; i++) {
    theta[i] = thtmp[i * NX];
    /* printf("theta[%3d] = %lf\n", i, theta[i]); */
  }

  // the spacing in each coordinate is
  const CCTK_REAL dX = (X[NX - 1] - X[0]) / (NX - 1);
  const CCTK_REAL dtheta = (theta[Ntheta - 1] - theta[0]) / (Ntheta - 1);

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
  W_in = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dW_dr_in = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dW_dth_in = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  const CCTK_REAL oodX = 1. / dX;
  // const CCTK_REAL oodXsq     = oodX * oodX;
  const CCTK_REAL oodX12 = 1. / (12. * dX);
  const CCTK_REAL oodXsq12 = oodX * oodX12;
  const CCTK_REAL oodth12 = 1. / (12. * dtheta);

  for (int jj = 0; jj < Ntheta; jj++) {
    for (int i = 0; i < NX; i++) {

      CCTK_INT j, jm1, jm2, jp1, jp2;
      /* let's use the fact that the solution is axi-symmetric (and that
              theta[0] = 0) for the boundary points in j */
      if (jj == 0) {
        j = jj;
        jp1 = jj + 1;
        jp2 = jj + 2;
        jm1 = jj + 1;
        jm2 = jj + 2;
      } else if (jj == 1) {
        j = jj;
        jp1 = jj + 1;
        jp2 = jj + 2;
        jm1 = jj - 1;
        jm2 = jj;
      } else if (jj == Ntheta - 2) {
        j = jj;
        jm1 = jj - 1;
        jm2 = jj - 2;
        jp1 = jj + 1;
        jp2 = jj;
      } else if (jj == Ntheta - 1) {
        j = jj;
        jm1 = jj - 1;
        jm2 = jj - 2;
        jp1 = jj - 1;
        jp2 = jj - 2;
      } else {
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

      if (i == 0) {
        /* For the Boson Star, there's no issue, dWbar/dX != 0 at X==0, and x and rhor coordinates coincide. */

        // 1st derivative with 4th order accuracy (forward stencils)
        Wbar_X = (-25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] + 16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) * oodX12;

        if (Wbar_r_power == 2) {
          // If Wbar = rhor^2 * W, to compute W(rhor=0), we need to compute Wbar_XX.
          // 2nd derivative with 4th order accuracy (forward stencils)
          Wbar_XX = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] + 214 * Wbar_in[indip2] - 156 * Wbar_in[indip3] + 61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) * oodXsq12;
        }

      } else if (i == 1) {
        // 1st derivative, 4th order accuracy
        Wbar_X = (-3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] - 6 * Wbar_in[indip2] + Wbar_in[indip3]) * oodX12;

      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4 * Wbar_in[indim1] + 3 * Wbar_in[ind]) * 0.5 * oodX;

      } else if (i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;

      } else {
        // 4th order accurate stencils
        Wbar_X = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;
      }

      // From the X coordinate used in the input files to the rhor coordinate (coincides with x for the Boson Star, rH=0).
      // We also do the conversion from Wbar to W here, to tackle rhor = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  rhor == 0
      if (i == 0) {
        // At rhor=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_in[ind] = 0.;
        dW_dth_in[ind] = 0.;

        // For W we need more care depending on the power
        switch (Wbar_r_power) {
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
      else if (i == NX - 1) {
        // W -> 0 for rhor -> infty
        W_in[ind] = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction paper) also gives:
        dW_dr_in[ind] = 0.;
        dW_dth_in[ind] = 0.;

      } else {
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

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

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

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

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
        if (rr2_2 < pow(eps_r, 2)) {
          rr2_2 = pow(eps_r, 2);
        }
        const CCTK_REAL rr_2 = sqrt(rr2_2);

        CCTK_REAL rho2_2 = x1_2 * x1_2 * gamma2 + y1_2 * y1_2;
        if (rho2_2 < pow(eps_r, 2)) {
          rho2_2 = pow(eps_r, 2);
        }
        const CCTK_REAL rho_2 = sqrt(rho2_2);

        const CCTK_REAL theta_2 = acos(z1_2 / rr_2);

        const CCTK_REAL deltakerr2_2 = bh_mass * bh_mass - bh_spin2;
        const CCTK_REAL deltakerr = sqrt(deltakerr2_2);

        const CCTK_REAL costh = z1_2 / rr_2;
        const CCTK_REAL costh2 = costh * costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth = sqrt(sinth2);

        const CCTK_REAL R_x = x1_2 * gamma / rr_2;
        const CCTK_REAL R_y = y1_2 / rr_2;
        const CCTK_REAL R_z = z1_2 / rr_2;

        const CCTK_REAL x_R = x1_2 * gamma / rr_2;
        const CCTK_REAL y_R = y1_2 / rr_2;
        const CCTK_REAL z_R = z1_2 / rr_2;

        const CCTK_REAL sinth2ph_x = -y1_2 / rr2_2;
        const CCTK_REAL sinth2ph_y = gamma * x1_2 / rr2_2;

        const CCTK_REAL sinthth_x = z1_2 * gamma * x1_2 / (rr_2 * rr2_2);
        const CCTK_REAL sinthth_y = z1_2 * y1_2 / (rr_2 * rr2_2);
        const CCTK_REAL sinthth_z = -sinth2 / rr_2;

        const CCTK_REAL sinthx_th = gamma * x1_2 * costh;
        const CCTK_REAL sinthy_th = y1_2 * costh;
        const CCTK_REAL sinthz_th = -rr_2 * sinth2;

        const CCTK_REAL costh2_x = -2 * x1_2 * gamma * pow(z1_2, 2) / pow(rr2_2, 2);
        const CCTK_REAL costh2_y = -2 * y1_2 * pow(z1_2, 2) / pow(rr2_2, 2);
        const CCTK_REAL costh2_z = 2 * rho2_2 * z1_2 / pow(rr2_2, 2);

        const CCTK_REAL rBL = rr_2 + bh_mass + 0.25 * (bh_mass2 - bh_spin2) / rr_2; // Boyer-Lindquist coordinate rhor

        // const CCTK_REAL RRrBL  = rr2_2 + rr_2*bh_mass + 0.25*(bh_mass2-bh_spin2);

        CCTK_REAL rho2kerr = rBL * rBL + bh_spin2 * costh2;
        // if( rho2kerr < pow( eps_r, 2 ) ) {
        //   rho2kerr = pow( eps_r, 2 );
        // }
        const CCTK_REAL rhokerr = sqrt(rho2kerr);

        const CCTK_REAL sigma = (2. * bh_mass * rBL) / rho2kerr;
        const CCTK_REAL hh = (1 + sigma) / (rho2kerr * rr2_2);

        // NaN/Inf checks for sigma and hh
        check_nan_or_inf("sigma", sigma);
        check_nan_or_inf("hh", hh);

        const CCTK_REAL psi4_2 = rho2kerr / rr2_2;
        const CCTK_REAL psi2_2 = sqrt(psi4_2);
        const CCTK_REAL psi1_2 = sqrt(psi2_2);
        const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        const CCTK_REAL psi2_1 = sqrt(psi4_1);
        const CCTK_REAL psi1_1 = sqrt(psi2_1);

        // non-axisymmetric perturbation.
        /* pert = 1. + AA * (x1_2*x1_2 - y1_2*y1_2)/(bh_mass*bh_mass) * exp( -2.*rr2_2/deltakerr2_2 ) ; */

        CCTK_REAL alpha0 = (rr_2 + horizon_radius) * (rr_2 - horizon_radius) / rr_2 *
                           1. / sqrt(rBL * rBL + bh_spin2 * (1. + sigma * sinth2));
        // check_nan_or_inf("alpha0", alpha0);
        // if (alpha0 < SMALL) {
        //   // If alpha0 is too small, we set it to zero to avoid division by zero.
        //   // This is a safeguard against numerical issues, especially near the horizon.
        //   // It should not affect the physics significantly, as alpha0 is a small correction.
        //   // In the context of the Kerr metric, this corresponds to the case where we are very
        //   // close to the horizon, where the metric becomes singular.
        //   alpha0 = SMALL;
        // }
        // if (alpha0 < SMALL) {
        //   // If alpha0 is too small, we set it to zero to avoid division by zero.
        //   // This is a safeguard against numerical issues, especially near the horizon.
        //   // It should not affect the physics significantly, as alpha0 is a small correction.
        //   // In the context of the Kerr metric, this corresponds to the case where we are very
        //   // close to the horizon, where the metric becomes singular.
        //   alpha0 = SMALL;
        // }
        const CCTK_REAL alpha02 = alpha0 * alpha0;

        // fazer as derivadas das funcoes auxiliares

        const CCTK_REAL dr_dR = 1 + (bh_spin2 - bh_mass * bh_mass) / (4 * rr2_2);

        const CCTK_REAL drho2kerr_dx = 2 * rBL * dr_dR * R_x + bh_spin2 * costh2_x;
        const CCTK_REAL drho2kerr_dy = 2 * rBL * dr_dR * R_y + bh_spin2 * costh2_y;
        const CCTK_REAL drho2kerr_dz = 2 * rBL * dr_dR * R_z + bh_spin2 * costh2_z;

        const CCTK_REAL dsigma_dx = (2 * bh_mass * (rho2kerr * R_x * dr_dR - rBL * drho2kerr_dx)) / pow(rho2kerr, 2);
        const CCTK_REAL dsigma_dy = (2 * bh_mass * (rho2kerr * R_y * dr_dR - rBL * drho2kerr_dy)) / pow(rho2kerr, 2);
        const CCTK_REAL dsigma_dz = (2 * bh_mass * (rho2kerr * R_z * dr_dR - rBL * drho2kerr_dz)) / pow(rho2kerr, 2);

        const CCTK_REAL dpsi4_2_dx = (-2 * rho2kerr * R_x + rr_2 * drho2kerr_dx) / pow(rr_2, 3);
        const CCTK_REAL dpsi4_2_dy = (-2 * rho2kerr * R_y + rr_2 * drho2kerr_dy) / pow(rr_2, 3);
        const CCTK_REAL dpsi4_2_dz = (-2 * rho2kerr * R_z + rr_2 * drho2kerr_dz) / pow(rr_2, 3);

        const CCTK_REAL bphi = -bh_spin * sigma * sinth2;
        // const CCTK_REAL gammaphiphi = psi4_2 * rr2_2 * sinth2 * (1 + bh_spin2 * hh * rr2_2 * sinth2);
        const CCTK_REAL gammaphiphi = rho2kerr * rho2_2 / rr2_2 * (1 + bh_spin2 * hh * rho2_2);
        const CCTK_REAL bphiup = bphi / gammaphiphi;

        // Check for NaN or Inf in these quantities
        check_nan_or_inf("bphi", bphi);
        check_nan_or_inf("gammaphiphi", gammaphiphi);
        check_nan_or_inf("bphiup", bphiup);

        CCTK_REAL betadrest[4] = {0.0, 0.0, 0.0, 0.0};
        betadrest[1] = -y1_2 / rho2_2 * bphi;
        betadrest[2] = x1_2 * gamma / rho2_2 * bphi;
        betadrest[3] = 0.;

        // Check for NaN or Inf in betadrest components
        check_nan_or_inf("betadrest[1]", betadrest[1]);
        check_nan_or_inf("betadrest[2]", betadrest[2]);
        check_nan_or_inf("betadrest[3]", betadrest[3]);

        // To change the spin direction, change the indices accordingly.
        const CCTK_REAL dbetadrestphi_dx = (-2 * bh_spin * pow(z1_2, 2) * sigma * R_x) / pow(rr_2, 3) - bh_spin * (1 - pow(z1_2, 2) / pow(rr_2, 2)) * dsigma_dx;
        const CCTK_REAL dbetadrestphi_dy = (-2 * bh_spin * pow(z1_2, 2) * sigma * R_y) / pow(rr_2, 3) - bh_spin * (1 - pow(z1_2, 2) / pow(rr_2, 2)) * dsigma_dy;
        const CCTK_REAL dbetadrestphi_dz = (bh_spin * (2 * z1_2 * sigma * (rr_2 - z1_2 * R_z) + rr_2 * (-rho2_2) * dsigma_dz)) / pow(rr_2, 3);

        CCTK_REAL dbetadrest[4][4];
        // Initialize g to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            dbetadrest[i][j] = 0.0;
        // Compute derivatives of the betadrest vector. To change the spin direction, change the indices accordingly.
        dbetadrest[1][1] = (2 * x1_2 * gamma * y1_2 * bphi - y1_2 * rho2_2 * dbetadrestphi_dx) / pow(rho2_2, 2);
        dbetadrest[1][2] = ((-pow(x1_2, 2) * gamma2 + pow(y1_2, 2)) * bphi - y1_2 * (rho2_2)*dbetadrestphi_dy) / pow(rho2_2, 2);
        dbetadrest[1][3] = -((y1_2 * dbetadrestphi_dz) / (rho2_2));
        dbetadrest[2][1] = ((-pow(x1_2, 2) * gamma2 + pow(y1_2, 2)) * bphi + x1_2 * gamma * rho2_2 * dbetadrestphi_dx) / pow(rho2_2, 2);
        dbetadrest[2][2] = (x1_2 * gamma * (-2 * y1_2 * bphi + rho2_2 * dbetadrestphi_dy)) / pow(rho2_2, 2);
        dbetadrest[2][3] = (x1_2 * gamma * dbetadrestphi_dz) / (rho2_2);
        dbetadrest[3][1] = 0;
        dbetadrest[3][2] = 0;
        dbetadrest[3][3] = 0;

        // Check for NaN or Inf in all dbetadrest[i][j] components
        for (int ii = 1; ii <= 3; ++ii) {
          for (int jj = 1; jj <= 3; ++jj) {
            char dbetadrest_name[32];
            snprintf(dbetadrest_name, sizeof(dbetadrest_name), "dbetadrest[%d][%d]", ii, jj);
            check_nan_or_inf(dbetadrest_name, dbetadrest[ii][jj]);
          }
        }

        const CCTK_REAL dhh_dx = (-((1 + sigma) * (2 * rho2kerr * R_x +
                                                   rr_2 * drho2kerr_dx)) +
                                  rr_2 * rho2kerr * dsigma_dx) /
                                 (pow(rr_2, 3) * pow(rho2kerr, 2));
        const CCTK_REAL dhh_dy = (-((1 + sigma) * (2 * rho2kerr * R_y +
                                                   rr_2 * drho2kerr_dy)) +
                                  rr_2 * rho2kerr * dsigma_dy) /
                                 (pow(rr_2, 3) * pow(rho2kerr, 2));
        const CCTK_REAL dhh_dz = (-((1 + sigma) * (2 * rho2kerr * R_z +
                                                   rr_2 * drho2kerr_dz)) +
                                  rr_2 * rho2kerr * dsigma_dz) /
                                 (pow(rr_2, 3) * pow(rho2kerr, 2));

        const CCTK_REAL dalpha_dx = (2 * pow(rr_2, 2) * (pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma)) * R_x +
                                     2 * pow(horizon_radius, 2) * (pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma)) * R_x +
                                     rr_2 * pow(horizon_radius, 2) * (2 * rBL * R_x * dr_dR + bh_spin2 * (-(sigma * costh2_x) - (-1 + costh2) * dsigma_dx)) +
                                     pow(rr_2, 3) * (-2 * rBL * R_x * dr_dR + bh_spin2 * (sigma * costh2_x + (-1 + costh2) * dsigma_dx))) /
                                    (2. * pow(rr_2, 2) * pow(pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma), 1.5));
        const CCTK_REAL dalpha_dy = (2 * pow(rr_2, 2) * (pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma)) * R_y + 2 * pow(horizon_radius, 2) * (pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma)) * R_y +
                                     rr_2 * pow(horizon_radius, 2) * (2 * rBL * R_y * dr_dR + bh_spin2 * (-(sigma * costh2_y) - (-1 + costh2) * dsigma_dy)) +
                                     pow(rr_2, 3) * (-2 * rBL * R_y * dr_dR + bh_spin2 * (sigma * costh2_y + (-1 + costh2) * dsigma_dy))) /
                                    (2. * pow(rr_2, 2) * pow(pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma), 1.5));
        const CCTK_REAL dalpha_dz = (2 * pow(rr_2, 2) * (pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma)) * R_z + 2 * pow(horizon_radius, 2) * (pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma)) * R_z +
                                     rr_2 * pow(horizon_radius, 2) * (2 * rBL * R_z * dr_dR + bh_spin2 * (-(sigma * costh2_z) - (-1 + costh2) * dsigma_dz)) +
                                     pow(rr_2, 3) * (-2 * rBL * R_z * dr_dR + bh_spin2 * (sigma * costh2_z + (-1 + costh2) * dsigma_dz))) /
                                    (2. * pow(rr_2, 2) * pow(pow(rBL, 2) + bh_spin2 * (1 - (-1 + costh2) * sigma), 1.5));

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

        check_nan_or_inf("dbetadrestphi_dx", drho2kerr_dx);
        check_nan_or_inf("dbetadrestphi_dy", drho2kerr_dy);
        check_nan_or_inf("dbetadrestphi_dz", drho2kerr_dz);

        // capital Gs refer to the unboosted frame.

        CCTK_REAL G[4][4];
        // Initialize G to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            G[i][j] = 0.0;

        G[0][0] = -alpha02 + bphi * bphiup;
        G[0][1] = betadrest[1];
        G[0][2] = betadrest[2];
        G[0][3] = betadrest[3];
        G[1][0] = G[0][1];
        G[2][0] = G[0][2];
        G[3][0] = G[0][3];
        G[1][1] = psi4_2 * (1 + bh_spin2 * hh * y1_2 * y1_2);
        G[1][2] = -psi4_2 * bh_spin2 * hh * y1_2 * x1_2 * gamma;
        G[1][3] = 0.0;
        G[2][1] = G[1][2];
        G[2][2] = psi4_2 * (1 + bh_spin2 * hh * x1_2 * x1_2 * gamma2);
        G[2][3] = 0.0;
        G[3][1] = G[1][3];
        G[3][2] = G[2][3];
        G[3][3] = psi4_2;

        // CCTK_REAL G3[4][4];
        // G3[0][0] = bphi * bphiup;
        // G3[0][1] = G[0][1];
        // G3[0][2] = G[0][2];
        // G3[0][3] = G[0][3];
        // G3[1][0] = G3[0][1];
        // G3[2][0] = G3[0][2];
        // G3[3][0] = G3[0][3];
        // G3[1][1] = G[1][1];
        // G3[1][2] = G[1][2];
        // G3[1][3] = G[1][3];
        // G3[2][1] = G[2][1];
        // G3[2][2] = G[2][2];
        // G3[2][3] = G[2][3];
        // G3[3][1] = G[3][1];
        // G3[3][2] = G[3][2];
        // G3[3][3] = G[3][3];

        // // const CCTK_REAL Gtt = G[0][0];
        // // const CCTK_REAL Gxt = G[0][1];
        // // const CCTK_REAL Gxx = G[1][1];
        // // const CCTK_REAL Gxy = G[1][2];
        // // const CCTK_REAL Gty = G[0][2];
        // // const CCTK_REAL Gyy = G[2][2];

        // CCTK_REAL G3_inv[4][4];
        // // Initialize G3_inv to zero. Inverse of 3 metric
        // for (int i = 0; i < 4; ++i)
        //   for (int j = 0; j < 4; ++j)
        //     G3_inv[i][j] = 0.0;

        // G3_inv[1][1] = (1 + bh_spin2 * pow(x1_2 * gamma, 2) * hh) / ((1 + bh_spin2 * (rho2_2)*hh) * psi4_2);
        // G3_inv[1][2] = (bh_spin2 * x1_2 * gamma * y1_2 * hh) / ((1 + bh_spin2 * (rho2_2)*hh) * psi4_2);
        // G3_inv[1][3] = 0;
        // G3_inv[2][1] = G3_inv[1][2];
        // G3_inv[2][2] = (1 + bh_spin2 * pow(y1_2, 2) * hh) / ((1 + bh_spin2 * (rho2_2)*hh) * psi4_2);
        // G3_inv[2][3] = 0;
        // G3_inv[3][1] = G3_inv[1][3];
        // G3_inv[3][2] = G3_inv[2][3];
        // G3_inv[3][3] = 1 / psi4_2;

        // CCTK_REAL dG3_inv[4][4][4];
        // // Initialize dG3_inv to zero
        // for (int a = 0; a < 4; ++a) {
        //   for (int b = 0; b < 4; ++b) {
        //     for (int c = 0; c < 4; ++c) {
        //       dG3_inv[a][b][c] = 0.0;
        //     }
        //   }
        // }
        // // dG3_inv[a][b][c] = dG^ab/dx^c
        // dG3_inv[1][1][1] = (bh_spin2 * pow(y1_2, 2) * psi4_2 * (2 * bh_spin2 * x1_2 * gamma * pow(hh, 2) - dhh_dx) - (1 + bh_spin2 * pow(x1_2 * gamma, 2) * hh) * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dx) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2));
        // dG3_inv[1][1][2] = -((bh_spin2 * y1_2 * psi4_2 * (2 * hh * (1 + bh_spin2 * pow(x1_2 * gamma, 2) * hh) + y1_2 * dhh_dy) + (1 + bh_spin2 * pow(x1_2 * gamma, 2) * hh) * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dy) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2)));
        // dG3_inv[1][1][3] = -((bh_spin2 * pow(y1_2, 2) * psi4_2 * dhh_dz + (1 + bh_spin2 * pow(x1_2 * gamma, 2) * hh) * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dz) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2)));
        // dG3_inv[1][2][1] = (bh_spin2 * y1_2 * (psi4_2 * (hh + bh_spin2 * (-pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * pow(hh, 2) + x1_2 * gamma * dhh_dx) - x1_2 * gamma * hh * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dx)) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2));
        // dG3_inv[1][2][2] = (bh_spin2 * x1_2 * gamma * (psi4_2 * (hh + bh_spin2 * (x1_2 * gamma - y1_2) * (x1_2 * gamma + y1_2) * pow(hh, 2) + y1_2 * dhh_dy) - y1_2 * hh * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dy)) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2));
        // dG3_inv[1][2][3] = (bh_spin2 * x1_2 * gamma * y1_2 * (psi4_2 * dhh_dz - hh * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dz)) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2));
        // dG3_inv[2][2][1] = -((bh_spin2 * x1_2 * gamma * psi4_2 * (2 * hh * (1 + bh_spin2 * pow(y1_2, 2) * hh) + x1_2 * gamma * dhh_dx) + (1 + bh_spin2 * pow(y1_2, 2) * hh) * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dx) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2)));
        // dG3_inv[2][2][2] = (bh_spin2 * pow(x1_2 * gamma, 2) * psi4_2 * (2 * bh_spin2 * y1_2 * pow(hh, 2) - dhh_dy) - (1 + bh_spin2 * pow(y1_2, 2) * hh) * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dy) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2));
        // dG3_inv[2][2][3] = -((bh_spin2 * pow(x1_2 * gamma, 2) * psi4_2 * dhh_dz + (1 + bh_spin2 * pow(y1_2, 2) * hh) * (1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh) * dpsi4_2_dz) / (pow(1 + bh_spin2 * (pow(x1_2 * gamma, 2) + pow(y1_2, 2)) * hh, 2) * pow(psi4_2, 2)));
        // dG3_inv[3][3][1] = -dpsi4_2_dx / (psi4_2 * psi4_2);
        // dG3_inv[3][3][2] = -dpsi4_2_dy / (psi4_2 * psi4_2);
        // dG3_inv[3][3][3] = -dpsi4_2_dz / (psi4_2 * psi4_2);
        // // dG3_inv[mu][nu][0] = 0
        // // Symmetries
        // dG3_inv[2][1][1] = dG3_inv[1][2][1];
        // dG3_inv[2][1][2] = dG3_inv[1][2][2];
        // dG3_inv[2][1][3] = dG3_inv[1][2][3];
        // // actually identically equal to zero anyway
        // dG3_inv[3][1][1] = dG3_inv[1][3][1];
        // dG3_inv[3][1][2] = dG3_inv[1][3][2];
        // dG3_inv[3][1][3] = dG3_inv[1][3][3];
        // dG3_inv[3][2][1] = dG3_inv[2][3][1];
        // dG3_inv[3][2][2] = dG3_inv[2][3][2];
        // dG3_inv[3][2][3] = dG3_inv[2][3][3];

        // CCTK_REAL betaup[4];
        // betaup[0] = 0.0;
        // betaup[1] = G3_inv[1][1] * betadrest[1] + G3_inv[1][2] * betadrest[2] + G3_inv[1][3] * betadrest[3];
        // betaup[2] = G3_inv[2][1] * betadrest[1] + G3_inv[2][2] * betadrest[2] + G3_inv[2][3] * betadrest[3];
        // betaup[3] = G3_inv[3][1] * betadrest[1] + G3_inv[3][2] * betadrest[2] + G3_inv[3][3] * betadrest[3];

        // CCTK_REAL dbetaup[4][4];
        // // Initialize dbetaup to zero
        // for (int i = 0; i < 4; ++i)
        //   for (int j = 0; j < 4; ++j)
        //     dbetaup[i][j] = 0.0;

        // dbetaup[1][1] = dG3_inv[1][1][1] * betadrest[1] + G3_inv[1][1] * dbetadrest[1][1] +
        //                 dG3_inv[1][2][1] * betadrest[2] + G3_inv[1][2] * dbetadrest[2][1] +
        //                 dG3_inv[1][3][1] * betadrest[3] + G3_inv[1][3] * dbetadrest[3][1];
        // dbetaup[1][2] = dG3_inv[1][1][2] * betadrest[1] + G3_inv[1][1] * dbetadrest[1][2] +
        //                 dG3_inv[1][2][2] * betadrest[2] + G3_inv[1][2] * dbetadrest[2][2] +
        //                 dG3_inv[1][3][2] * betadrest[3] + G3_inv[1][3] * dbetadrest[3][2];
        // dbetaup[1][3] = dG3_inv[1][1][3] * betadrest[1] + G3_inv[1][1] * dbetadrest[1][3] +
        //                 dG3_inv[1][2][3] * betadrest[2] + G3_inv[1][2] * dbetadrest[2][3] +
        //                 dG3_inv[1][3][3] * betadrest[3] + G3_inv[1][3] * dbetadrest[3][3];
        // dbetaup[2][1] = dG3_inv[2][1][1] * betadrest[1] + G3_inv[2][1] * dbetadrest[1][1] +
        //                 dG3_inv[2][2][1] * betadrest[2] + G3_inv[2][2] * dbetadrest[2][1] +
        //                 dG3_inv[2][3][1] * betadrest[3] + G3_inv[2][3] * dbetadrest[3][1];
        // dbetaup[2][2] = dG3_inv[2][1][2] * betadrest[1] + G3_inv[2][1] * dbetadrest[1][2] +
        //                 dG3_inv[2][2][2] * betadrest[2] + G3_inv[2][2] * dbetadrest[2][2] +
        //                 dG3_inv[2][3][2] * betadrest[3] + G3_inv[2][3] * dbetadrest[3][2];
        // dbetaup[2][3] = dG3_inv[2][1][3] * betadrest[1] + G3_inv[2][1] * dbetadrest[1][3] +
        //                 dG3_inv[2][2][3] * betadrest[2] + G3_inv[2][2] * dbetadrest[2][3] +
        //                 dG3_inv[2][3][3] * betadrest[3] + G3_inv[2][3] * dbetadrest[3][3];
        // dbetaup[3][1] = dG3_inv[3][1][1] * betadrest[1] + G3_inv[3][1] * dbetadrest[1][1] +
        //                 dG3_inv[3][2][1] * betadrest[2] + G3_inv[3][2] * dbetadrest[2][1] +
        //                 dG3_inv[3][3][1] * betadrest[3] + G3_inv[3][3] * dbetadrest[3][1];
        // dbetaup[3][2] = dG3_inv[3][1][2] * betadrest[1] + G3_inv[3][1] * dbetadrest[1][2] +
        //                 dG3_inv[3][2][2] * betadrest[2] + G3_inv[3][2] * dbetadrest[2][2] +
        //                 dG3_inv[3][3][2] * betadrest[3] + G3_inv[3][3] * dbetadrest[3][2];
        // dbetaup[3][3] = dG3_inv[3][1][3] * betadrest[1] + G3_inv[3][1] * dbetadrest[1][3] +
        //                 dG3_inv[3][2][3] * betadrest[2] + G3_inv[3][2] * dbetadrest[2][3] +
        //                 dG3_inv[3][3][3] * betadrest[3] + G3_inv[3][3] * dbetadrest[3][3];

        // const CCTK_REAL dbeta2_dx = -((bh_spin2 * sigma * (rr_2 * (rho2_2) * (1 + bh_spin2 * hh * (rho2_2)) * sigma * dpsi4_2_dx + psi4_2 * (sigma * (2 * (-2 * pow(z1_2, 2) + pow(rr_2, 2) + 2 * bh_spin2 * hh * pow(-rho2_2, 2)) * R_x + bh_spin2 * rr_2 * pow(-rho2_2, 2) * dhh_dx) - 2 * rr_2 * (rho2_2) * (1 + bh_spin2 * hh * (-pow(z1_2, 2) + pow(rr_2, 2))) * dsigma_dx))) / (pow(psi4_2, 2) * pow(rr_2, 5) * pow(1 + bh_spin2 * hh * (rho2_2), 2)));
        // const CCTK_REAL dbeta2_dy = -((bh_spin2 * sigma * (rr_2 * (rho2_2) * (1 + bh_spin2 * hh * (rho2_2)) * sigma * dpsi4_2_dy + psi4_2 * (sigma * (2 * (-2 * pow(z1_2, 2) + pow(rr_2, 2) + 2 * bh_spin2 * hh * pow(-rho2_2, 2)) * R_y + bh_spin2 * rr_2 * pow(-rho2_2, 2) * dhh_dy) - 2 * rr_2 * (rho2_2) * (1 + bh_spin2 * hh * (-pow(z1_2, 2) + pow(rr_2, 2))) * dsigma_dy))) / (pow(psi4_2, 2) * pow(rr_2, 5) * pow(1 + bh_spin2 * hh * (rho2_2), 2)));
        // const CCTK_REAL dbeta2_dz = (bh_spin2 * sigma * (sigma * (psi4_2 * (-2 * z1_2 * rr_2 - 2 * (-2 * pow(z1_2, 2) + pow(rr_2, 2) + 2 * bh_spin2 * hh * pow(-rho2_2, 2)) * R_z - bh_spin2 * rr_2 * pow(-rho2_2, 2) * dhh_dz) - rr_2 * (rho2_2) * (1 + bh_spin2 * hh * (rho2_2)) * dpsi4_2_dz) + 2 * psi4_2 * rr_2 * (rho2_2) * (1 + bh_spin2 * hh * (rho2_2)) * dsigma_dz)) / (pow(psi4_2, 2) * pow(rr_2, 5) * pow(1 + bh_spin2 * hh * (rho2_2), 2));
        // check_nan_or_inf("dbeta2_dx", dbeta2_dx);
        // check_nan_or_inf("dbeta2_dy", dbeta2_dy);
        // check_nan_or_inf("dbeta2_dz", dbeta2_dz);

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

        dG[1][1][1] = bh_spin2 * pow(y1_2, 2) * psi4_2 * dhh_dx + (1 + bh_spin2 * pow(y1_2, 2) * hh) * dpsi4_2_dx;

        dG[1][1][2] = bh_spin2 * y1_2 * psi4_2 * (2 * hh + y1_2 * dhh_dy) + (1 + bh_spin2 * pow(y1_2, 2) * hh) * dpsi4_2_dy;

        dG[1][1][3] = bh_spin2 * pow(y1_2, 2) * psi4_2 * dhh_dz + (1 + bh_spin2 * pow(y1_2, 2) * hh) * dpsi4_2_dz;

        dG[1][2][1] = -(bh_spin2 * y1_2 * (x1_2 * gamma * psi4_2 * dhh_dx + hh * (psi4_2 + x1_2 * gamma * dpsi4_2_dx)));

        dG[1][2][2] = -(bh_spin2 * x1_2 * gamma * (y1_2 * psi4_2 * dhh_dy + hh * (psi4_2 + y1_2 * dpsi4_2_dy)));

        dG[1][2][3] = -(bh_spin2 * x1_2 * gamma * y1_2 * (psi4_2 * dhh_dz + hh * dpsi4_2_dz));

        // dG[1][3][i] = 0

        dG[2][2][1] = bh_spin2 * x1_2 * gamma * psi4_2 * (2 * hh + x1_2 * gamma * dhh_dx) + (1 +
                                                                                             bh_spin2 * pow(x1_2 * gamma, 2) * hh) *
                                                                                                dpsi4_2_dx;
        dG[2][2][2] = bh_spin2 * pow(x1_2 * gamma, 2) * psi4_2 * dhh_dy + (1 +
                                                                           bh_spin2 * pow(x1_2 * gamma, 2) * hh) *
                                                                              dpsi4_2_dy;
        dG[2][2][3] = bh_spin2 * pow(x1_2 * gamma, 2) * psi4_2 * dhh_dz + (1 +
                                                                           bh_spin2 * pow(x1_2 * gamma, 2) * hh) *
                                                                              dpsi4_2_dz;

        // dG23_dx^i = 0

        dG[3][3][1] = dpsi4_2_dx;
        dG[3][3][2] = dpsi4_2_dy;
        dG[3][3][3] = dpsi4_2_dz;

        // dG[mu][nu][0] = 0
        dG[0][0][1] = -2 * alpha0 * dalpha_dx - (bh_spin2 * (pow(betadrest[1], 2) + pow(betadrest[2], 2) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * hh) * (2 * x1_2 * gamma * hh + (rho2_2)*dhh_dx)) / (pow(1 + bh_spin2 * (rho2_2)*hh, 2) * psi4_2) + (2 * betadrest[1] * dbetadrest[1][1] + 2 * betadrest[2] * dbetadrest[2][1] + 2 * bh_spin2 * (x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2]) * hh * (betadrest[1] + x1_2 * gamma * dbetadrest[1][1] + y1_2 * dbetadrest[2][1]) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * dhh_dx) / ((1 + bh_spin2 * (rho2_2)*hh) * psi4_2) - ((pow(betadrest[1], 2) + pow(betadrest[2], 2) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * hh) * dpsi4_2_dx) / ((1 + bh_spin2 * (rho2_2)*hh) * pow(psi4_2, 2));
        dG[0][0][2] = -2 * alpha0 * dalpha_dy - (bh_spin2 * (pow(betadrest[1], 2) + pow(betadrest[2], 2) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * hh) * (2 * y1_2 * hh + (rho2_2)*dhh_dy)) / (pow(1 + bh_spin2 * (rho2_2)*hh, 2) * psi4_2) + (2 * betadrest[1] * dbetadrest[1][2] + 2 * betadrest[2] * dbetadrest[2][2] + 2 * bh_spin2 * (x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2]) * hh * (betadrest[2] + x1_2 * gamma * dbetadrest[1][2] + y1_2 * dbetadrest[2][2]) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * dhh_dy) / ((1 + bh_spin2 * (rho2_2)*hh) * psi4_2) - ((pow(betadrest[1], 2) + pow(betadrest[2], 2) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * hh) * dpsi4_2_dy) / ((1 + bh_spin2 * (rho2_2)*hh) * pow(psi4_2, 2));
        dG[0][0][3] = -2 * alpha0 * dalpha_dz - (bh_spin2 * (rho2_2) * (pow(betadrest[1], 2) + pow(betadrest[2], 2) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * hh) * dhh_dz) / (pow(1 + bh_spin2 * (rho2_2)*hh, 2) * psi4_2) + (2 * betadrest[1] * dbetadrest[1][3] + 2 * betadrest[2] * dbetadrest[2][3] + 2 * bh_spin2 * (x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2]) * hh * (x1_2 * gamma * dbetadrest[1][3] + y1_2 * dbetadrest[2][3]) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * dhh_dz) / ((1 + bh_spin2 * (rho2_2)*hh) * psi4_2) - ((pow(betadrest[1], 2) + pow(betadrest[2], 2) + bh_spin2 * pow(x1_2 * gamma * betadrest[1] + y1_2 * betadrest[2], 2) * hh) * dpsi4_2_dz) / ((1 + bh_spin2 * (rho2_2)*hh) * pow(psi4_2, 2));
        dG[0][1][1] = dbetadrest[1][1];
        dG[0][1][2] = dbetadrest[1][2];
        dG[0][1][3] = dbetadrest[1][3];
        dG[0][2][1] = dbetadrest[2][1];
        dG[0][2][2] = dbetadrest[2][2];
        dG[0][2][3] = dbetadrest[2][3];
        dG[0][3][1] = dbetadrest[3][1];
        dG[0][3][2] = dbetadrest[3][2];
        dG[0][3][3] = dbetadrest[3][3];
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

        // Boosting the metric

        CCTK_REAL invLambda[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
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
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
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
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            Gb[a][b] = 0.0;
          }
        }

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
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

        // Invert the spatial 3x3 block of the boosted metric Gb into Gb3_inv
        CCTK_REAL Gb3_inv[4][4];
        for (int a = 0; a < 4; ++a)
          for (int b = 0; b < 4; ++b)
            Gb3_inv[a][b] = 0.0;

        {

          // Build symmetric 3x3 block M from Gb
          CCTK_REAL M[4][4] = {{0}};
          M[1][1] = Gb[1][1];
          M[1][2] = 0.5 * (Gb[1][2] + Gb[2][1]);
          M[1][3] = 0.5 * (Gb[1][3] + Gb[3][1]);
          M[2][1] = M[1][2];
          M[2][2] = Gb[2][2];
          M[2][3] = 0.5 * (Gb[2][3] + Gb[3][2]);
          M[3][1] = M[1][3];
          M[3][2] = M[2][3];
          M[3][3] = Gb[3][3];

          if (!invert_spd3x3(M, Gb3_inv)) {
            CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Gb_spatial not SPD or ill-conditioned at (%d,%d,%d).",
                       i, j, k);
          }

          // Optional checks
          check_nan_or_inf("Gb3_inv[1][1]", Gb3_inv[1][1]);
          check_nan_or_inf("Gb3_inv[1][2]", Gb3_inv[1][2]);
          check_nan_or_inf("Gb3_inv[1][3]", Gb3_inv[1][3]);
          check_nan_or_inf("Gb3_inv[2][2]", Gb3_inv[2][2]);
          check_nan_or_inf("Gb3_inv[2][3]", Gb3_inv[2][3]);
          check_nan_or_inf("Gb3_inv[3][3]", Gb3_inv[3][3]);

          // Verify Gb3_inv * Gb_spatial ≈ I (indices 1..3)
          {
            const CCTK_REAL tol = SMALL;
            for (int i3 = 1; i3 <= 3; ++i3) {
              for (int j3 = 1; j3 <= 3; ++j3) {
                CCTK_REAL s = 0.0;
                for (int k3 = 1; k3 <= 3; ++k3)
                  s += Gb3_inv[i3][k3] * Gb[k3][j3];
                const CCTK_REAL delta = (i3 == j3) ? 1.0 : 0.0;
                if (fabs(s - delta) > tol) {
                  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                             "Gb3_inv check failed at (%d,%d): %g (tol=%g)", i3,
                             j3, (double)s, (double)tol);
                  break;
                }
              }
            }
          }
        }

        CCTK_REAL beta[4];
        beta[0] = NAN;
        beta[1] = Gb[1][0];
        beta[2] = Gb[2][0];
        beta[3] = Gb[3][0];

        CCTK_REAL betaup[4];
        betaup[0] = 0.0;
        betaup[1] = Gb3_inv[1][1] * beta[1] + Gb3_inv[1][2] * beta[2] +
                    Gb3_inv[1][3] * beta[3];
        betaup[2] = Gb3_inv[2][1] * beta[1] + Gb3_inv[2][2] * beta[2] +
                    Gb3_inv[2][3] * beta[3];
        betaup[3] = Gb3_inv[3][1] * beta[1] + Gb3_inv[3][2] * beta[2] +
                    Gb3_inv[3][3] * beta[3];

        // Derivatives of the boosted metric
        CCTK_REAL dg[4][4][4]; // dg[i][j][k] = \partial_k g_{ij} (boosted metric)
        // Initialize dg to zero
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 4; ++jj)
            for (int kk = 0; kk < 4; ++kk)
              dg[ii][jj][kk] = 0.0;

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              CCTK_REAL sum = 0.0;
              for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                  for (int lam = 0; lam < 4; ++lam)
                    sum += invLambda[mu][a] * invLambda[nu][b] * invLambda[lam][c] * dG[mu][nu][lam];
              dg[a][b][c] = sum;
            }
          }
        }

        CCTK_REAL new_alpha2 = -Gb[0][0] + betaup[1] * beta[1] +
                               betaup[2] * beta[2] + betaup[3] * beta[3];
        CCTK_REAL new_alpha = sqrt(new_alpha2);

        CCTK_REAL dW_drho_1, dW_dz_1;
        const CCTK_REAL exp_auxi_1 = exp(2. * F2_1[ind] - F0_1[ind]);

        if (rho_1 < 1e-8) {
          dW_drho_1 = 0.;
          dW_dz_1 = 0.;
        } else {
          dW_drho_1 = rho_1 / rr_1 * dW_dr_1[ind] + z1_1 / rr2_1 * dW_dth_1[ind];
          dW_dz_1 = z1_1 / rr_1 * dW_dr_1[ind] - rho_1 / rr2_1 * dW_dth_1[ind];
        }

        CCTK_REAL K_B[4][4]; // extrinsic curvature
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            K_B[a][b] = 0.0;
          }
        } // K_0\mu might not be zero but irrelevant for what i want to compute

        for (int a = 1; a < 4; ++a) {
          for (int b = 1; b < 4; ++b) {
            CCTK_REAL sum1 = 0.0;
            CCTK_REAL sum2 = 0.0;
            CCTK_REAL sum3 = 0.0;
            for (int c = 1; c < 4; ++c) {
              sum1 += betaup[c] * dg[a][b][c];
              sum2 += betaup[c] * dg[b][c][a];
              sum3 += betaup[c] * dg[a][c][b];
            }
            K_B[a][b] = -0.5 / new_alpha * (dg[a][b][0] - sum1 - (dg[0][b][a] - sum2) - (dg[0][a][b] - sum3));
          }
        }

        kxx[ind] = K_B[1][1];
        kxy[ind] = K_B[1][2];
        kxz[ind] = K_B[1][3];
        kyy[ind] = K_B[2][2];
        kyz[ind] = K_B[2][3];
        kzz[ind] = K_B[3][3];

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
          betax[ind] = betadrest[1]; // por enquato o shift da bs é zero pois é estática.
          betay[ind] = betadrest[2];
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