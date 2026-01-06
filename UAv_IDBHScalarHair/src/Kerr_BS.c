// este e o codigo copiado do kerrnewman.
#include "UAv_Derivatives.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define SMALL (1.e-9)

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
  UAv_ID_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, phi0_in,
                   Wbar_in);

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
      const CCTK_REAL Wbar_th = (-Wbar_in[indjp2] + 8 * Wbar_in[indjp1] -
                                 8 * Wbar_in[indjm1] + Wbar_in[indjm2]) *
                                oodth12;

      CCTK_REAL Wbar_X;
      CCTK_REAL Wbar_XX = 0.; // Used for rhor=0 (i==0), if Wbar_r_power == 2.

      /*
      Regarding finite differencing orders: plotting W and dW_dr, there were
      small discontinuities near rhor=0 during tests with the previous 2nd order
      accuracy for i==0 and i==1. Those vanish when moving to 4th order
      accuracy.

      For i==NX-1 and i==NX-2, we keep 2nd order for now. The issue is not
      appearing as clearly, and they represent points which are physically far,
      so maybe better to keep the computation more local.
      */

      if (i == 0) {
        /* For the Boson Star, there's no issue, dWbar/dX != 0 at X==0, and x
         * and rhor coordinates coincide. */

        // 1st derivative with 4th order accuracy (forward stencils)
        Wbar_X =
            (-25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] +
             16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) *
            oodX12;

        if (Wbar_r_power == 2) {
          // If Wbar = rhor^2 * W, to compute W(rhor=0), we need to compute
          // Wbar_XX. 2nd derivative with 4th order accuracy (forward stencils)
          Wbar_XX = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] +
                     214 * Wbar_in[indip2] - 156 * Wbar_in[indip3] +
                     61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) *
                    oodXsq12;
        }
      } else if (i == 1) {
        // 1st derivative, 4th order accuracy
        Wbar_X =
            (-3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] -
             6 * Wbar_in[indip2] + Wbar_in[indip3]) *
            oodX12;
      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4 * Wbar_in[indim1] + 3 * Wbar_in[ind]) *
                 0.5 * oodX;
      } else if (i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        Wbar_X = (-Wbar_in[indim1] + Wbar_in[indip1]) * 0.5 * oodX;
      } else {
        // 4th order accurate stencils
        Wbar_X = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] +
                  Wbar_in[indim2]) *
                 oodX12;
      }

      // From the X coordinate used in the input files to the rhor coordinate
      // (coincides with x for the Boson Star, rH=0). We also do the conversion
      // from Wbar to W here, to tackle rhor = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  rhor == 0
      if (i == 0) {
        // At rhor=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_in[ind] = 0.;
        dW_dth_in[ind] = 0.;

        // For W we need more care depending on the pow
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

        default: // As of writing, this should be prevented by the scope of the
                 // parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Unknown value of Wbar_r_power: %d. Aborting.",
                     Wbar_r_power);
          break;
        }
      }

      // We need to be careful at X == 1 (rhor == infty) for radial derivatives
      // (coordinate change is singular)
      else if (i == NX - 1) {
        // W -> 0 for rhor -> infty
        W_in[ind] = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction
        // paper) also gives:
        dW_dr_in[ind] = 0.;
        dW_dth_in[ind] = 0.;
      } else {
        const CCTK_REAL rr = C0 * lX / (1. - lX);

        // corresponding derivatives
        // const CCTK_REAL dXdr = 1./(C0 + rr) - rr/((C0 + rr)*(C0 + rr));
        const CCTK_REAL dXdr = C0 / ((C0 + rr) * (C0 + rr));

        const CCTK_REAL Wbar_r = dXdr * Wbar_X;

        // Now translate from Wbar to W
        switch (Wbar_r_power) // We could put a generic pow for the computation
                              // here I guess...
        {
        case 0: // Wbar = W
          W_in[ind] = Wbar_in[ind];
          dW_dr_in[ind] = Wbar_r;
          dW_dth_in[ind] = Wbar_th;
          break;

        case 1: // Wbar = rhor * W
          W_in[ind] = Wbar_in[ind] / rr;
          dW_dr_in[ind] =
              (Wbar_r - W_in[ind]) / rr; // dW/dr  =  1/rhor * dWbar/dr - Wbar /
                                         // rhor^2  =  (dWbar/dr - W) / rhor
          dW_dth_in[ind] = Wbar_th / rr;
          break;

        case 2:; // Wbar = rhor^2 * W
          // empty statement after case to prevent compilation error on some gcc
          // versions...
          const CCTK_REAL rr2_2 = rr * rr;
          W_in[ind] = Wbar_in[ind] / rr2_2;
          dW_dr_in[ind] =
              Wbar_r / rr2_2 -
              2 * W_in[ind] /
                  rr; // dW/dr  =  1/rhor^2 * dWbar/dr - 2 * Wbar / rhor^3  =
                      // 1/rhor^2 * dWbar/dr - 2 * W / rhor
          dW_dth_in[ind] = Wbar_th / rr2_2;
          break;

        default: // As of writing, this should be prevented by the scope of the
                 // parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Unknown value of Wbar_r_power: %d. Aborting.",
                     Wbar_r_power);
          break;
        }
      } // if/else i==...

    } // for i
  } // for jj

  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points =
      cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]; // total points

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

        CCTK_REAL ltheta_1 =
            rr_1 < 1e-16
                ? 0
                : acos(z1_1 /
                       rr_1);      // There should be at most one point in the grid
                                   // with rr~0. Not sure about the threshold.
        if (ltheta_1 > 0.5 * M_PI) // symmetry along the equatorial plane
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

        // CCTK_REAL ltheta_2 = rr_2 < 1e-16 ? 0 : acos( z1_2/rr_2);    // There
        // should be at most one point in the grid with rr~0. Not sure about the
        // threshold. if (ltheta_2 > 0.5*M_PI)    // symmetry along the
        // equatorial plane
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
  param_table_handle = Util_TableCreateFromString(
      "order=4 boundary_extrapolation_tolerance={0.1 1.0 0.05 0.05}");

  CCTK_INFO("Interpolating result...");

  /* do the actual interpolation, and check for error returns */
  int status_1 = CCTK_InterpLocalUniform(
      N_dims, operator_handle, param_table_handle, origin, delta,
      N_interp_points, CCTK_VARIABLE_REAL, interp_coords_1, N_input_arrays,
      input_array_dims, input_array_type_codes, input_arrays, N_output_arrays,
      output_array_type_codes_1, output_arrays_1);
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
  //                                      N_output_arrays,
  //                                      output_array_type_codes_2,
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
  /* For the Boson Star, in order to avoid unneeded regularizations and
     divisions, we express K_ij in terms of dW/drho and dW/dz. For points close
     to the axis and the origin, K_ij = 0.
  */

  const CCTK_REAL tt = cctk_time;

  const CCTK_REAL coswt = cos(omega_BS * tt);
  const CCTK_REAL sinwt = sin(omega_BS * tt);

  const CCTK_REAL bh_spin2 = bh_spin * bh_spin;
  const CCTK_REAL bh_mass2 = bh_mass * bh_mass;

  const CCTK_REAL rBLp = bh_mass + sqrt(bh_mass2 - bh_spin2);
  const CCTK_REAL rBLm = bh_mass - sqrt(bh_mass2 - bh_spin2);

  // const CCTK_REAL horizon_radius = 0.5 * sqrt(bh_mass2 - bh_spin2);

  // printf("cctk_lsh[0] = %d\n",cctk_lsh[0]);
  // printf("cctk_lsh[1] = %d\n",cctk_lsh[1]);
  // printf("cctk_lsh[2] = %d\n",cctk_lsh[2]);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // Boson Star A ignorar por agora

        const CCTK_REAL x1_1 = x[ind] - x0;
        const CCTK_REAL y1_1 = y[ind] - y0;
        const CCTK_REAL z1_1 = z[ind] - z0;

        // For the Boson Star, rhor = R, no coordinate change needed.
        const CCTK_REAL rr2_1 = x1_1 * x1_1 + y1_1 * y1_1 + z1_1 * z1_1;
        const CCTK_REAL rr_1 = sqrt(rr2_1);

        const CCTK_REAL rho2_1 = x1_1 * x1_1 + y1_1 * y1_1;
        const CCTK_REAL rho_1 = sqrt(rho2_1);

        const CCTK_REAL ph_1 = atan2(y1_1, x1_1);
        // If x1_2=y1_2=0, should return 0? The other metric functions should
        // vanish anyway to make sure that this doesn't matter, but can this
        // lead to nan depending on the C implementation?

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

        // All quantities evaluated at point (x*gamma,y,z)

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

        const CCTK_REAL rho3_2 = rho2_2 * rho_2;

        // R0pert2 = (rr_2 - R0pert)*(rr_2 - R0pert) ;

        const CCTK_REAL costh = z1_2 / rr_2;
        const CCTK_REAL costh2 = costh * costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth = sqrt(sinth2);

        const CCTK_REAL ph_2 = atan2(y1_2, x1_2);
        // If x1_2=y1_2=0, should return 0? The other metric functions should
        // vanish anyway to make sure that this doesn't matter, but can this
        // lead to nan depending on the C implementation?

        const CCTK_REAL cosph_2 = cos(ph_2);
        const CCTK_REAL sinph_2 = sin(ph_2);

        const CCTK_REAL R_x = x1_2 * gamma / rr_2;
        const CCTK_REAL R_y = y1_2 / rr_2;
        const CCTK_REAL R_z = z1_2 / rr_2;

        const CCTK_REAL x_R = x1_2 * gamma / rr_2;
        const CCTK_REAL y_R = y1_2 / rr_2;
        const CCTK_REAL z_R = z1_2 / rr_2;

        const CCTK_REAL th_x = costh * R_x / rho_2;
        const CCTK_REAL th_y = costh * R_y / rho_2;
        const CCTK_REAL th_z = -rho_2 / rr2_2;
        // auxilliary quantities
        CCTK_REAL rBL, rBL2;
        rBL = rr_2 * (1.0 + 0.25 * rBLp / rr_2) * (1.0 + 0.25 * rBLp / rr_2);
        rBL2 = rBL * rBL;
        CCTK_REAL Delt, Sigm, Sigm2, fctFF;
        Delt = (rBL - rBLp) * (rBL - rBLm); // rBL2 + bh_spin2 - 2 * bh_mass * rBL;
        Sigm = rBL2 + bh_spin2 * costh2;
        Sigm2 = Sigm * Sigm;
        fctFF = (rBL2 + bh_spin2) * (rBL2 + bh_spin2) - Delt * bh_spin2 * sinth2; // "A" no artigo

        const CCTK_REAL psi4_2 = Sigm / rr2_2; // psi04 no codigo original
        const CCTK_REAL psi2_2 = sqrt(psi4_2);
        const CCTK_REAL psi1_2 = sqrt(psi2_2);
        const CCTK_REAL psi4_1 = exp(2. * F1_1[ind]);
        const CCTK_REAL psi2_1 = sqrt(psi4_1);
        const CCTK_REAL psi1_1 = sqrt(psi2_1);

        CCTK_REAL fctGG, fctHH;
        fctGG = rBLm / (rr2_2 * (rBL - rBLm));
        fctHH = (2.0 * bh_mass * rBL + Sigm) / (rr2_2 * Sigm2);

        // CCTK_REAL detgij;
        // detgij = psi4_2 * psi4_2 * psi4_2 * (1.0 + rr2_2 * fctGG) * (1.0 +
        // bh_spin2 * rho2_2 * fctHH);

        /*----------------------------------*/

        // const CCTK_REAL alpha02 = (4.0 * rr_2 - rBLp) * (4.0 * rr_2 - rBLp) *
        // (rBL - rBLm) / (16.0 * rr_2 * (rBL2 + bh_spin2 * (1.0 + 2.0 * bh_mass
        // * rBL * sinth2 / Sigm)));
        const CCTK_REAL alpha0 = (4.0 * rr_2 - rBLp) * sqrt(rBL - rBLm) / sqrt(16.0 * rr_2 * (rBL2 + bh_spin2 * (1.0 + 2.0 * bh_mass * rBL * sinth2 / Sigm))); // primeiro termo para schwarzschild e zero;
        const CCTK_REAL alpha02 = alpha0 * alpha0;
        // const CCTK_REAL alpha0 = sqrt(alpha02);
        // const CCTK_REAL alpha0 = sqrt(Delt * Sigm / fctFF);
        // const CCTK_REAL alpha02 = Delt * Sigm / fctFF;
        // const CCTK_REAL alpha0 = sqrt(alpha02);
        // const CCTK_REAL dalpha0_dR = 0.5 / alpha0 * (-(Delt * Sigm *
        // dfctFF_dR) + fctFF * (Sigm * dDelt_dR + Delt * dSigm_dR)) /
        // pow(fctFF, 2);
        // const CCTK_REAL dalpha0_dth = 0.5 / alpha0 * (Delt *
        // (-(Sigm * dfctFF_dth) + fctFF * dSigm_dth)) / pow(fctFF, 2);
        // const CCTK_REAL dalpha02_dR =(-(Delt * Sigm * dfctFF_dR) + fctFF * (Sigm * dDelt_dR + Delt * dSigm_dR)) /(fctFF * fctFF);
        // const CCTK_REAL dalpha02_dth = (Delt * (-(Sigm * dfctFF_dth) + fctFF * dSigm_dth)) / (fctFF * fctFF);
        const CCTK_REAL gphiphi = fctFF / Sigm * sinth2;
        const CCTK_REAL bphiup = -2.0 * bh_mass * bh_spin * rBL / fctFF;
        const CCTK_REAL bphi = bphiup * gphiphi;

        // 1st derivatives of auxiliary quantities
        CCTK_REAL drBLdR, dfctGG_dR, dfctHH_dR, dfctHH_dth, dpsi4_2_dR, dpsi4_2_dth, dDelt_dR, dSigm_dR, dSigm_dth, dfctFF_dR, dfctFF_dth;
        drBLdR = 1.0 - rBLp * rBLp / (16.0 * rr2_2);
        dDelt_dR = drBLdR * (rBL - rBLp) + drBLdR * (rBL - rBLm); //(2 * rBL - 2 * bh_mass) * drBLdR;
        dSigm_dR = 2 * rBL * drBLdR;
        dSigm_dth = -2 * bh_spin2 * costh * sinth;
        dfctFF_dR = 4 * rBL * (bh_spin2 + rBL2) * drBLdR - bh_spin2 * sinth2 * dDelt_dR;
        dfctFF_dth = -(bh_spin2 * Delt * 2 * costh * sinth);
        dfctGG_dR = -((rBLm * (-2 * rBLm + 2 * rBL + rr_2 * drBLdR)) / (rr2_2 * rr_2 * (rBLm - rBL) * (rBLm - rBL)));
        dfctHH_dR = -((2 * Sigm * (2 * bh_mass * rBL + Sigm - bh_mass * rr_2 * drBLdR) + rr_2 * (4 * bh_mass * rBL + Sigm) * dSigm_dR) / (rr2_2 * rr_2 * Sigm2 * Sigm));
        dfctHH_dth = -(((4 * bh_mass * rBL + Sigm) * dSigm_dth) / (rr2_2 * Sigm2 * Sigm));
        dpsi4_2_dR = (-2 * Sigm + rr_2 * dSigm_dR) / (rr2_2 * rr_2);
        dpsi4_2_dth = dSigm_dth / rr2_2;
        const CCTK_REAL dalpha02_dR = ((4 * rr_2 - rBLp) * (2 * bh_mass * bh_spin2 * Sigm * sinth2 * (-((4 * rr_2 + rBLp) * (rBLm - rBL) * rBL) + rr_2 * rBLm * (4 * rr_2 - rBLp) * drBLdR) + Sigm2 * (-((4 * rr_2 + rBLp) * (rBLm - rBL) * (bh_spin2 + rBL2)) + rr_2 * (4 * rr_2 - rBLp) * (bh_spin2 + 2 * rBLm * rBL - rBL2) * drBLdR) - 2 * bh_mass * bh_spin2 * rr_2 * (4 * rr_2 - rBLp) * (rBLm - rBL) * rBL * sinth2 * dSigm_dR)) / (16. * (rr_2 * (bh_spin2 + rBL2) * Sigm + 2 * bh_mass * bh_spin2 * rr_2 * rBL * sinth2) * (rr_2 * (bh_spin2 + rBL2) * Sigm + 2 * bh_mass * bh_spin2 * rr_2 * rBL * sinth2));
        const CCTK_REAL dalpha02_dth = (bh_mass * bh_spin2 * (-4 * rr_2 + rBLp) * (-4 * rr_2 + rBLp) * (rBLm - rBL) * rBL * sinth * (2 * costh * Sigm - sinth * dSigm_dth)) / (8. * rr_2 * ((bh_spin2 + rBL2) * Sigm + 2 * bh_mass * bh_spin2 * rBL * sinth2) * ((bh_spin2 + rBL2) * Sigm + 2 * bh_mass * bh_spin2 * rBL * sinth2));
        const CCTK_REAL dgphiphi_dR = (sinth2 * (Sigm * dfctFF_dR - fctFF * dSigm_dR)) / Sigm2;
        const CCTK_REAL dgphiphi_dth = (sinth * (Sigm * sinth * dfctFF_dth + fctFF * (2 * costh * Sigm - sinth * dSigm_dth))) / Sigm2;
        const CCTK_REAL dbphiup_dR = (2 * bh_mass * bh_spin * (-(fctFF * drBLdR) + rBL * dfctFF_dR)) / (fctFF * fctFF);
        const CCTK_REAL dbphiup_dth = (2 * bh_mass * bh_spin * rBL * dfctFF_dth) / (fctFF * fctFF);
        const CCTK_REAL dbphi_dR = gphiphi * dbphiup_dR + bphiup * dgphiphi_dR;
        const CCTK_REAL dbphi_dth = gphiphi * dbphiup_dth + bphiup * dgphiphi_dth;

        // 2nd derivatives
        /* 2nd derivatives (Cartesian coordinates) of auxiliary quantities.
           Skeleton for user to fill in.
           Notation suggestion:
             d2<name>_dx2  = ∂^2(<name>)/∂x^2
             d2<name>_dxy  = ∂^2(<name>)/∂x∂y
             etc.

           Chain rule pattern for a scalar F(R,th):
             Fx  = F_R * R_x + F_th * th_x
             Fxx = F_RR * R_x*R_x + 2*F_Rth * R_x*th_x + F_thth * th_x*th_x
             + F_R * R_xx + F_th * th_xx
             (analogous for mixed and other second derivatives)

           Provide (and fill) second derivatives of R and theta if needed:
             R_xx, R_xy, R_xz, R_yy, R_yz, R_zz
             th_xx, th_xy, th_xz, th_yy, th_yz, th_zz

           Below everything initialized to 0.0 for later replacement.
        */

        /* Coordinate second derivatives (fill if needed) */
        CCTK_REAL R_xx = y1_2 * z1_2 / (rr2_2 * rr_2), R_xy = -y1_2 * x1_2 * gamma / (rr2_2 * rr_2), R_xz = -x1_2 * gamma * z1_2 / (rr2_2 * rr_2);
        CCTK_REAL R_yy = x1_2 * gamma * z1_2 / (rr2_2 * rr_2), R_yz = -y1_2 * z1_2 / (rr2_2 * rr_2), R_zz = rho2_2 / (rr2_2 * rr_2);
        CCTK_REAL th_xx = z1_2 * (-2 * pow(x1_2 * gamma, 4) - pow(x1_2 * gamma * y1_2, 2) + pow(y1_2, 4) + pow(y1_2 * z1_2, 2)) / (rho3_2 * rr2_2 * rr2_2), th_xy = -x1_2 * gamma * y1_2 * z1_2 * (3 * rho2_2 + z1_2 * z1_2) / (rho3_2 * rr2_2 * rr2_2), th_xz = x1_2 * gamma * (rho2_2 - z1_2 * z1_2) / (rho3_2 * rr2_2 * rr2_2);
        CCTK_REAL th_yy = z1_2 * (pow(x1_2 * gamma, 4) - 2 * pow(y1_2, 4) + x1_2 * (-y1_2 * y1_2 + z1_2 * z1_2)) / (rho3_2 * rr2_2 * rr2_2), th_yz = y1_2 * (rho2_2 - z1_2 * z1_2) / (rho3_2 * rr2_2 * rr2_2), th_zz = 2 * z1_2 * rho_2 / rr2_2;

        /* Second derivatives in (R,theta) space (to be provided) */
        const CCTK_REAL d2rBL_dR2 = rBLp * rBLp / (8 * rr2_2 * rr_2), d2rBL_dRth = 0.0, d2rBL_dth2 = 0.0;
        const CCTK_REAL d2Delt_dR2 = 2 * (pow(drBLdR, 2) + (-bh_mass + rBL) * d2rBL_dR2), d2Delt_dRth = 0.0, d2Delt_dth2 = 0.0;
        const CCTK_REAL d2Sigm_dR2 = 2 * (pow(drBLdR, 2) + rBL * d2rBL_dR2), d2Sigm_dRth = 0.0, d2Sigm_dth2 = -2 * bh_spin2 * (costh2 - sinth2);
        const CCTK_REAL d2fctFF_dR2 = 4 * (bh_spin2 + 3 * rBL2) * pow(drBLdR, 2) - bh_spin2 * sinth2 * d2Delt_dR2 + 4 * rBL * (bh_spin2 + rBL2) * d2rBL_dR2, d2fctFF_dRth = -2 * bh_spin2 * costh * sinth * dDelt_dR, d2fctFF_dth2 = -2 * bh_spin2 * (costh2 - sinth2) * Delt;
        const CCTK_REAL d2fctGG_dR2 = (rBLm * (-6 * pow(rBLm - rBL, 2) + 4 * rr_2 * (rBLm - rBL) * drBLdR - rr2_2 * (2 * pow(drBLdR, 2) + (rBLm - rBL) * d2rBL_dR2))) / (pow(rr_2, 4) * pow(rBLm - rBL, 3)), d2fctGG_dRth = 0.0, d2fctGG_dth2 = 0.0;
        const CCTK_REAL d2fctHH_dR2 = (2 * (Sigm * (3 * Sigm2 + bh_mass * rr_2 * Sigm * (-4 * drBLdR + rr_2 * d2rBL_dR2) - 2 * bh_mass * pow(rr_2, 2) * drBLdR * dSigm_dR) + bh_mass * rBL * (6 * Sigm2 + 2 * rr2_2 * pow(dSigm_dR, 2) + rr_2 * Sigm * (4 * dSigm_dR - rr_2 * d2Sigm_dR2)))) / (pow(rr_2, 4) * pow(Sigm, 3)), d2fctHH_dRth = (2 * bh_mass * dSigm_dth * (-(rr_2 * Sigm * drBLdR) + 2 * rBL * (Sigm + rr_2 * dSigm_dR)) - 2 * bh_mass * rr_2 * rBL * Sigm * d2Sigm_dRth) / (pow(rr_2, 3) * pow(Sigm, 3)), d2fctHH_dth2 = (-2 * bh_mass * rBL * (-2 * pow(dSigm_dth, 2) + Sigm * d2Sigm_dR2)) / (rr2_2 * pow(Sigm, 3));
        const CCTK_REAL d2psi4_2_dR2 = (6 * Sigm - 4 * rr_2 * dSigm_dR + rr2_2 * d2Sigm_dR2) / pow(rr_2, 4), d2psi4_2_dRth = (-2 * dSigm_dth + rr_2 * d2Sigm_dRth) / pow(rr_2, 3), d2psi4_2_dth2 = d2Sigm_dth2 / rr2_2;
        const CCTK_REAL d2alpha02_dR2 = (2 * Delt * Sigm * pow(dfctFF_dR, 2) + fctFF * (-2 * dfctFF_dR * (Sigm * dDelt_dR + Delt * dSigm_dR) - Delt * Sigm * d2fctFF_dR2) + pow(fctFF, 2) * (Sigm * d2Delt_dR2 + 2 * dDelt_dR * dSigm_dR + Delt * d2Sigm_dR2)) / pow(fctFF, 3), d2alpha02_dRth = (2 * Delt * Sigm * dfctFF_dth * dfctFF_dR - fctFF * (Delt * dSigm_dth * dfctFF_dR + dfctFF_dth * (Sigm * dDelt_dR + Delt * dSigm_dR) + Delt * Sigm * d2fctFF_dRth) + pow(fctFF, 2) * (dDelt_dR * dSigm_dth + Delt * d2Sigm_dRth)) / pow(fctFF, 3), d2alpha02_dth2 = (Delt * (Sigm * (2 * pow(dfctFF_dth, 2) - fctFF * d2fctFF_dth2) + fctFF * (-2 * dfctFF_dth * dSigm_dth + fctFF * d2Sigm_dth2))) / pow(fctFF, 3);
        const CCTK_REAL d2gphiphi_dR2 = (sinth2 * (2 * fctFF * pow(dSigm_dR, 2) + Sigm2 * d2fctFF_dR2 - Sigm * (2 * dfctFF_dR * dSigm_dR + fctFF * d2Sigm_dR2))) / pow(Sigm, 3), d2gphiphi_dRth = (sinth * (2 * fctFF * sinth * dSigm_dth * dSigm_dR + pow(Sigm, 2) * (2 * costh * dfctFF_dR + sinth * d2fctFF_dRth) - Sigm * (sinth * dSigm_dth * dfctFF_dR + (2 * costh * fctFF + sinth * dfctFF_dth) * dSigm_dR + fctFF * sinth * d2Sigm_dRth))) / pow(Sigm, 3), d2gphiphi_dth2 = (Sigm * sinth * (2 * dfctFF_dth * (2 * costh * Sigm - sinth * dSigm_dth) + Sigm * sinth * d2fctFF_dth2) + fctFF * (2 * (costh2 - sinth2) * Sigm2 + 2 * sinth2 * pow(dSigm_dth, 2) - Sigm * sinth * (4 * costh * dSigm_dth + sinth * d2Sigm_dth2))) / pow(Sigm, 3);
        const CCTK_REAL d2bphiup_dR2 = (-2 * bh_mass * bh_spin * (pow(fctFF, 2) * d2rBL_dR2 + 2 * rBL * pow(dfctFF_dR, 2) - fctFF * (2 * drBLdR * dfctFF_dR + rBL * d2fctFF_dR2))) / pow(fctFF, 3), d2bphiup_dRth = (2 * bh_mass * bh_spin * (dfctFF_dth * (fctFF * drBLdR - 2 * rBL * dfctFF_dR) + fctFF * rBL * d2fctFF_dRth)) / pow(fctFF, 3), d2bphiup_dth2 = (2 * bh_mass * bh_spin * rBL * (-2 * pow(dfctFF_dth, 2) + fctFF * d2fctFF_dth2)) / pow(fctFF, 3);
        const CCTK_REAL d2bphi_dR2 = 2 * dbphiup_dR * dgphiphi_dR + gphiphi * d2bphiup_dR2 + bphiup * d2gphiphi_dR2, d2bphi_dRth = dgphiphi_dth * dbphiup_dR + dbphiup_dth * dgphiphi_dR + gphiphi * d2bphiup_dRth + bphiup * d2gphiphi_dRth, d2bphi_dth2 = 2 * dbphiup_dth * dgphiphi_dth + gphiphi * d2bphiup_dth2 + bphiup * d2gphiphi_dth2;

        /* Cartesian second derivatives (fill using chain rule) */
        /* rBL */
        /* rBL (depends only on R) */
        CCTK_REAL d2rBL_dx2 = d2rBL_dR2 * R_x * R_x + drBLdR * R_xx;
        CCTK_REAL d2rBL_dy2 = d2rBL_dR2 * R_y * R_y + drBLdR * R_yy;
        CCTK_REAL d2rBL_dz2 = d2rBL_dR2 * R_z * R_z + drBLdR * R_zz;
        CCTK_REAL d2rBL_dxy = d2rBL_dR2 * R_x * R_y + drBLdR * R_xy;
        CCTK_REAL d2rBL_dxz = d2rBL_dR2 * R_x * R_z + drBLdR * R_xz;
        CCTK_REAL d2rBL_dyz = d2rBL_dR2 * R_y * R_z + drBLdR * R_yz;

        /* Delt */
        CCTK_REAL dDelt_dth = 0.0; // Delt depende apenas de R
        CCTK_REAL d2Delt_dx2 = d2Delt_dR2 * R_x * R_x + 2.0 * d2Delt_dRth * R_x * th_x + d2Delt_dth2 * th_x * th_x + dDelt_dR * R_xx + dDelt_dth * th_xx;
        CCTK_REAL d2Delt_dy2 = d2Delt_dR2 * R_y * R_y + 2.0 * d2Delt_dRth * R_y * th_y + d2Delt_dth2 * th_y * th_y + dDelt_dR * R_yy + dDelt_dth * th_yy;
        CCTK_REAL d2Delt_dz2 = d2Delt_dR2 * R_z * R_z + 2.0 * d2Delt_dRth * R_z * th_z + d2Delt_dth2 * th_z * th_z + dDelt_dR * R_zz + dDelt_dth * th_zz;
        CCTK_REAL d2Delt_dxy = d2Delt_dR2 * R_x * R_y + d2Delt_dRth * (R_x * th_y + R_y * th_x) + d2Delt_dth2 * th_x * th_y + dDelt_dR * R_xy + dDelt_dth * th_xy;
        CCTK_REAL d2Delt_dxz = d2Delt_dR2 * R_x * R_z + d2Delt_dRth * (R_x * th_z + R_z * th_x) + d2Delt_dth2 * th_x * th_z + dDelt_dR * R_xz + dDelt_dth * th_xz;
        CCTK_REAL d2Delt_dyz = d2Delt_dR2 * R_y * R_z + d2Delt_dRth * (R_y * th_z + R_z * th_y) + d2Delt_dth2 * th_y * th_z + dDelt_dR * R_yz + dDelt_dth * th_yz;

        /* Sigm */
        CCTK_REAL d2Sigm_dx2 = d2Sigm_dR2 * R_x * R_x + 2.0 * d2Sigm_dRth * R_x * th_x + d2Sigm_dth2 * th_x * th_x + dSigm_dR * R_xx + dSigm_dth * th_xx;
        CCTK_REAL d2Sigm_dy2 = d2Sigm_dR2 * R_y * R_y + 2.0 * d2Sigm_dRth * R_y * th_y + d2Sigm_dth2 * th_y * th_y + dSigm_dR * R_yy + dSigm_dth * th_yy;
        CCTK_REAL d2Sigm_dz2 = d2Sigm_dR2 * R_z * R_z + 2.0 * d2Sigm_dRth * R_z * th_z + d2Sigm_dth2 * th_z * th_z + dSigm_dR * R_zz + dSigm_dth * th_zz;
        CCTK_REAL d2Sigm_dxy = d2Sigm_dR2 * R_x * R_y + d2Sigm_dRth * (R_x * th_y + R_y * th_x) + d2Sigm_dth2 * th_x * th_y + dSigm_dR * R_xy + dSigm_dth * th_xy;
        CCTK_REAL d2Sigm_dxz = d2Sigm_dR2 * R_x * R_z + d2Sigm_dRth * (R_x * th_z + R_z * th_x) + d2Sigm_dth2 * th_x * th_z + dSigm_dR * R_xz + dSigm_dth * th_xz;
        CCTK_REAL d2Sigm_dyz = d2Sigm_dR2 * R_y * R_z + d2Sigm_dRth * (R_y * th_z + R_z * th_y) + d2Sigm_dth2 * th_y * th_z + dSigm_dR * R_yz + dSigm_dth * th_yz;

        /* fctFF */
        CCTK_REAL d2fctFF_dx2 = d2fctFF_dR2 * R_x * R_x + 2.0 * d2fctFF_dRth * R_x * th_x + d2fctFF_dth2 * th_x * th_x + dfctFF_dR * R_xx + dfctFF_dth * th_xx;
        CCTK_REAL d2fctFF_dy2 = d2fctFF_dR2 * R_y * R_y + 2.0 * d2fctFF_dRth * R_y * th_y + d2fctFF_dth2 * th_y * th_y + dfctFF_dR * R_yy + dfctFF_dth * th_yy;
        CCTK_REAL d2fctFF_dz2 = d2fctFF_dR2 * R_z * R_z + 2.0 * d2fctFF_dRth * R_z * th_z + d2fctFF_dth2 * th_z * th_z + dfctFF_dR * R_zz + dfctFF_dth * th_zz;
        CCTK_REAL d2fctFF_dxy = d2fctFF_dR2 * R_x * R_y + d2fctFF_dRth * (R_x * th_y + R_y * th_x) + d2fctFF_dth2 * th_x * th_y + dfctFF_dR * R_xy + dfctFF_dth * th_xy;
        CCTK_REAL d2fctFF_dxz = d2fctFF_dR2 * R_x * R_z + d2fctFF_dRth * (R_x * th_z + R_z * th_x) + d2fctFF_dth2 * th_x * th_z + dfctFF_dR * R_xz + dfctFF_dth * th_xz;
        CCTK_REAL d2fctFF_dyz = d2fctFF_dR2 * R_y * R_z + d2fctFF_dRth * (R_y * th_z + R_z * th_y) + d2fctFF_dth2 * th_y * th_z + dfctFF_dR * R_yz + dfctFF_dth * th_yz;

        /* fctGG (only R) */
        CCTK_REAL dfctGG_dth = 0.0;
        CCTK_REAL d2fctGG_dx2 = d2fctGG_dR2 * R_x * R_x + dfctGG_dR * R_xx;
        CCTK_REAL d2fctGG_dy2 = d2fctGG_dR2 * R_y * R_y + dfctGG_dR * R_yy;
        CCTK_REAL d2fctGG_dz2 = d2fctGG_dR2 * R_z * R_z + dfctGG_dR * R_zz;
        CCTK_REAL d2fctGG_dxy = d2fctGG_dR2 * R_x * R_y + dfctGG_dR * R_xy;
        CCTK_REAL d2fctGG_dxz = d2fctGG_dR2 * R_x * R_z + dfctGG_dR * R_xz;
        CCTK_REAL d2fctGG_dyz = d2fctGG_dR2 * R_y * R_z + dfctGG_dR * R_yz;

        /* fctHH */
        CCTK_REAL d2fctHH_dx2 = d2fctHH_dR2 * R_x * R_x + 2.0 * d2fctHH_dRth * R_x * th_x + d2fctHH_dth2 * th_x * th_x + dfctHH_dR * R_xx + dfctHH_dth * th_xx;
        CCTK_REAL d2fctHH_dy2 = d2fctHH_dR2 * R_y * R_y + 2.0 * d2fctHH_dRth * R_y * th_y + d2fctHH_dth2 * th_y * th_y + dfctHH_dR * R_yy + dfctHH_dth * th_yy;
        CCTK_REAL d2fctHH_dz2 = d2fctHH_dR2 * R_z * R_z + 2.0 * d2fctHH_dRth * R_z * th_z + d2fctHH_dth2 * th_z * th_z + dfctHH_dR * R_zz + dfctHH_dth * th_zz;
        CCTK_REAL d2fctHH_dxy = d2fctHH_dR2 * R_x * R_y + d2fctHH_dRth * (R_x * th_y + R_y * th_x) + d2fctHH_dth2 * th_x * th_y + dfctHH_dR * R_xy + dfctHH_dth * th_xy;
        CCTK_REAL d2fctHH_dxz = d2fctHH_dR2 * R_x * R_z + d2fctHH_dRth * (R_x * th_z + R_z * th_x) + d2fctHH_dth2 * th_x * th_z + dfctHH_dR * R_xz + dfctHH_dth * th_xz;
        CCTK_REAL d2fctHH_dyz = d2fctHH_dR2 * R_y * R_z + d2fctHH_dRth * (R_y * th_z + R_z * th_y) + d2fctHH_dth2 * th_y * th_z + dfctHH_dR * R_yz + dfctHH_dth * th_yz;

        /* psi4_2 */
        CCTK_REAL d2psi4_2_dx2 = d2psi4_2_dR2 * R_x * R_x + 2.0 * d2psi4_2_dRth * R_x * th_x + d2psi4_2_dth2 * th_x * th_x + dpsi4_2_dR * R_xx + dpsi4_2_dth * th_xx;
        CCTK_REAL d2psi4_2_dy2 = d2psi4_2_dR2 * R_y * R_y + 2.0 * d2psi4_2_dRth * R_y * th_y + d2psi4_2_dth2 * th_y * th_y + dpsi4_2_dR * R_yy + dpsi4_2_dth * th_yy;
        CCTK_REAL d2psi4_2_dz2 = d2psi4_2_dR2 * R_z * R_z + 2.0 * d2psi4_2_dRth * R_z * th_z + d2psi4_2_dth2 * th_z * th_z + dpsi4_2_dR * R_zz + dpsi4_2_dth * th_zz;
        CCTK_REAL d2psi4_2_dxy = d2psi4_2_dR2 * R_x * R_y + d2psi4_2_dRth * (R_x * th_y + R_y * th_x) + d2psi4_2_dth2 * th_x * th_y + dpsi4_2_dR * R_xy + dpsi4_2_dth * th_xy;
        CCTK_REAL d2psi4_2_dxz = d2psi4_2_dR2 * R_x * R_z + d2psi4_2_dRth * (R_x * th_z + R_z * th_x) + d2psi4_2_dth2 * th_x * th_z + dpsi4_2_dR * R_xz + dpsi4_2_dth * th_xz;
        CCTK_REAL d2psi4_2_dyz = d2psi4_2_dR2 * R_y * R_z + d2psi4_2_dRth * (R_y * th_z + R_z * th_y) + d2psi4_2_dth2 * th_y * th_z + dpsi4_2_dR * R_yz + dpsi4_2_dth * th_yz;

        /* alpha02 */
        CCTK_REAL d2alpha02_dx2 = d2alpha02_dR2 * R_x * R_x + 2.0 * d2alpha02_dRth * R_x * th_x + d2alpha02_dth2 * th_x * th_x + dalpha02_dR * R_xx + dalpha02_dth * th_xx;
        CCTK_REAL d2alpha02_dy2 = d2alpha02_dR2 * R_y * R_y + 2.0 * d2alpha02_dRth * R_y * th_y + d2alpha02_dth2 * th_y * th_y + dalpha02_dR * R_yy + dalpha02_dth * th_yy;
        CCTK_REAL d2alpha02_dz2 = d2alpha02_dR2 * R_z * R_z + 2.0 * d2alpha02_dRth * R_z * th_z + d2alpha02_dth2 * th_z * th_z + dalpha02_dR * R_zz + dalpha02_dth * th_zz;
        CCTK_REAL d2alpha02_dxy = d2alpha02_dR2 * R_x * R_y + d2alpha02_dRth * (R_x * th_y + R_y * th_x) + d2alpha02_dth2 * th_x * th_y + dalpha02_dR * R_xy + dalpha02_dth * th_xy;
        CCTK_REAL d2alpha02_dxz = d2alpha02_dR2 * R_x * R_z + d2alpha02_dRth * (R_x * th_z + R_z * th_x) + d2alpha02_dth2 * th_x * th_z + dalpha02_dR * R_xz + dalpha02_dth * th_xz;
        CCTK_REAL d2alpha02_dyz = d2alpha02_dR2 * R_y * R_z + d2alpha02_dRth * (R_y * th_z + R_z * th_y) + d2alpha02_dth2 * th_y * th_z + dalpha02_dR * R_yz + dalpha02_dth * th_yz;

        /* gphiphi */
        CCTK_REAL d2gphiphi_dx2 = d2gphiphi_dR2 * R_x * R_x + 2.0 * d2gphiphi_dRth * R_x * th_x + d2gphiphi_dth2 * th_x * th_x + dgphiphi_dR * R_xx + dgphiphi_dth * th_xx;
        CCTK_REAL d2gphiphi_dy2 = d2gphiphi_dR2 * R_y * R_y + 2.0 * d2gphiphi_dRth * R_y * th_y + d2gphiphi_dth2 * th_y * th_y + dgphiphi_dR * R_yy + dgphiphi_dth * th_yy;
        CCTK_REAL d2gphiphi_dz2 = d2gphiphi_dR2 * R_z * R_z + 2.0 * d2gphiphi_dRth * R_z * th_z + d2gphiphi_dth2 * th_z * th_z + dgphiphi_dR * R_zz + dgphiphi_dth * th_zz;
        CCTK_REAL d2gphiphi_dxy = d2gphiphi_dR2 * R_x * R_y + d2gphiphi_dRth * (R_x * th_y + R_y * th_x) + d2gphiphi_dth2 * th_x * th_y + dgphiphi_dR * R_xy + dgphiphi_dth * th_xy;
        CCTK_REAL d2gphiphi_dxz = d2gphiphi_dR2 * R_x * R_z + d2gphiphi_dRth * (R_x * th_z + R_z * th_x) + d2gphiphi_dth2 * th_x * th_z + dgphiphi_dR * R_xz + dgphiphi_dth * th_xz;
        CCTK_REAL d2gphiphi_dyz = d2gphiphi_dR2 * R_y * R_z + d2gphiphi_dRth * (R_y * th_z + R_z * th_y) + d2gphiphi_dth2 * th_y * th_z + dgphiphi_dR * R_yz + dgphiphi_dth * th_yz;

        /* bphiup */
        CCTK_REAL d2bphiup_dx2 = d2bphiup_dR2 * R_x * R_x + 2.0 * d2bphiup_dRth * R_x * th_x + d2bphiup_dth2 * th_x * th_x + dbphiup_dR * R_xx + dbphiup_dth * th_xx;
        CCTK_REAL d2bphiup_dy2 = d2bphiup_dR2 * R_y * R_y + 2.0 * d2bphiup_dRth * R_y * th_y + d2bphiup_dth2 * th_y * th_y + dbphiup_dR * R_yy + dbphiup_dth * th_yy;
        CCTK_REAL d2bphiup_dz2 = d2bphiup_dR2 * R_z * R_z + 2.0 * d2bphiup_dRth * R_z * th_z + d2bphiup_dth2 * th_z * th_z + dbphiup_dR * R_zz + dbphiup_dth * th_zz;
        CCTK_REAL d2bphiup_dxy = d2bphiup_dR2 * R_x * R_y + d2bphiup_dRth * (R_x * th_y + R_y * th_x) + d2bphiup_dth2 * th_x * th_y + dbphiup_dR * R_xy + dbphiup_dth * th_xy;
        CCTK_REAL d2bphiup_dxz = d2bphiup_dR2 * R_x * R_z + d2bphiup_dRth * (R_x * th_z + R_z * th_x) + d2bphiup_dth2 * th_x * th_z + dbphiup_dR * R_xz + dbphiup_dth * th_xz;
        CCTK_REAL d2bphiup_dyz = d2bphiup_dR2 * R_y * R_z + d2bphiup_dRth * (R_y * th_z + R_z * th_y) + d2bphiup_dth2 * th_y * th_z + dbphiup_dR * R_yz + dbphiup_dth * th_yz;

        /* bphi */
        CCTK_REAL d2bphi_dx2 = d2bphi_dR2 * R_x * R_x + 2.0 * d2bphi_dRth * R_x * th_x + d2bphi_dth2 * th_x * th_x + dbphi_dR * R_xx + dbphi_dth * th_xx;
        CCTK_REAL d2bphi_dy2 = d2bphi_dR2 * R_y * R_y + 2.0 * d2bphi_dRth * R_y * th_y + d2bphi_dth2 * th_y * th_y + dbphi_dR * R_yy + dbphi_dth * th_yy;
        CCTK_REAL d2bphi_dz2 = d2bphi_dR2 * R_z * R_z + 2.0 * d2bphi_dRth * R_z * th_z + d2bphi_dth2 * th_z * th_z + dbphi_dR * R_zz + dbphi_dth * th_zz;
        CCTK_REAL d2bphi_dxy = d2bphi_dR2 * R_x * R_y + d2bphi_dRth * (R_x * th_y + R_y * th_x) + d2bphi_dth2 * th_x * th_y + dbphi_dR * R_xy + dbphi_dth * th_xy;
        CCTK_REAL d2bphi_dxz = d2bphi_dR2 * R_x * R_z + d2bphi_dRth * (R_x * th_z + R_z * th_x) + d2bphi_dth2 * th_x * th_z + dbphi_dR * R_xz + dbphi_dth * th_xz;
        CCTK_REAL d2bphi_dyz = d2bphi_dR2 * R_y * R_z + d2bphi_dRth * (R_y * th_z + R_z * th_y) + d2bphi_dth2 * th_y * th_z + dbphi_dR * R_yz + dbphi_dth * th_yz;

        /* Example chain rule template (replace <F> with variable, and supply needed F_RR etc):
           d2<F>_dx2 =d2<F>_dR2 * R_x* R_x +2.0 * d2<F>_dRth * R_x * th_x +d2<F>_dth2 * th_x * th_x +d<F>_dR * R_xx + d<F>_dth * th_xx;
           Similar for mixed and y,z derivatives.
        */

        CCTK_REAL G[4][4];
        // Initialize G to zero
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            G[i][j] = 0.0;

        G[0][0] = -alpha02 + bphi * bphiup;
        G[0][1] = -y1_2 / rho2_2 * bphi;
        G[0][2] = x1_2 * gamma / rho2_2 * bphi;
        G[0][3] = 0;
        G[1][0] = G[0][1];
        G[2][0] = G[0][2];
        G[3][0] = G[0][3];
        G[1][1] = psi4_2 * (1.0 + x1_2 * x1_2 * gamma2 * fctGG + bh_spin2 * y1_2 * y1_2 * fctHH);
        G[1][2] = psi4_2 * (x1_2 * gamma * y1_2 * fctGG - bh_spin2 * x1_2 * gamma * y1_2 * fctHH);
        G[1][3] = psi4_2 * (x1_2 * gamma * z1_2 * fctGG);
        G[2][1] = G[1][2];
        G[2][2] = psi4_2 * (1.0 + y1_2 * y1_2 * fctGG + bh_spin2 * x1_2 * x1_2 * gamma2 * fctHH);
        G[2][3] = psi4_2 * (y1_2 * z1_2 * fctGG);
        G[3][1] = G[1][3];
        G[3][2] = G[2][3];
        G[3][3] = psi4_2 * (1.0 + z1_2 * z1_2 * fctGG);

        /* NaN/Inf checks for metric tensor G */
        for (int aa = 0; aa < 4; ++aa) {
          for (int bb = 0; bb < 4; ++bb) {
            char name[32];
            snprintf(name, sizeof(name), "G[%d][%d]", aa, bb);
            check_nan_or_inf(name, G[aa][bb]);
          }
        }

        // Derivatives of the metric functions
        CCTK_REAL dG[4][4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              if (c == 0) {
                dG[a][b][c] = 0.0; // time derivatives of the metric are zero
              } else {
                dG[a][b][c] = NAN;
              }
            }
          }
        }

        UAv_EvalPoint P = {
            .bh_mass = bh_mass,
            .bh_spin = bh_spin,
            .bh_spin2 = bh_spin2,
            .bh_v = bh_v,
            .rBLp = rBLp,
            .rBLm = rBLm,
            .x1_2 = x1_2 * gamma,
            .y1_2 = y1_2,
            .z1_2 = z1_2,
            .rr_2 = rr_2,
            .rr2_2 = rr2_2,
            .rho_2 = rho_2,
            .rho2_2 = rho2_2,
            .R_x = R_x,
            .R_y = R_y,
            .R_z = R_z,
            .th_x = th_x,
            .th_y = th_y,
            .th_z = th_z,
            .costh = costh,
            .costh2 = costh2,
            .sinth = sinth,
            .sinth2 = sinth2,

            .R_xx = R_xx,
            .R_xy = R_xy,
            .R_xz = R_xz,
            .R_yy = R_yy,
            .R_yz = R_yz,
            .R_zz = R_zz,
            .th_xx = th_xx,
            .th_xy = th_xy,
            .th_xz = th_xz,
            .th_yy = th_yy,
            .th_yz = th_yz,
            .th_zz = th_zz,

            .rBL = rBL,
            .rBLp = rBLp,
            .rBLm = rBLm,
            .drBLdR = drBLdR,
            .d2rBL_dR2 = d2rBL_dR2,
            .d2rBL_dRth = d2rBL_dRth,
            .d2rBL_dth2 = d2rBL_dth2,

            .Delt = Delt,
            .dDelt_dR = dDelt_dR,
            .dDelt_dth = dDelt_dth,
            .d2Delt_dR2 = d2Delt_dR2,
            .d2Delt_dRth = d2Delt_dRth,
            .d2Delt_dth2 = d2Delt_dth2,

            .Sigm = Sigm,
            .Sigm2 = Sigm2,
            .dSigm_dR = dSigm_dR,
            .dSigm_dth = dSigm_dth,
            .d2Sigm_dR2 = d2Sigm_dR2,
            .d2Sigm_dRth = d2Sigm_dRth,
            .d2Sigm_dth2 = d2Sigm_dth2,

            .fctFF = fctFF,
            .dfctFF_dR = dfctFF_dR,
            .dfctFF_dth = dfctFF_dth,
            .d2fctFF_dR2 = d2fctFF_dR2,
            .d2fctFF_dRth = d2fctFF_dRth,
            .d2fctFF_dth2 = d2fctFF_dth2,

            .fctGG = fctGG,
            .dfctGG_dR = dfctGG_dR,
            .d2fctGG_dR2 = d2fctGG_dR2,
            .d2fctGG_dRth = d2fctGG_dRth,
            .d2fctGG_dth2 = d2fctGG_dth2,

            .fctHH = fctHH,
            .dfctHH_dR = dfctHH_dR,
            .dfctHH_dth = dfctHH_dth,
            .d2fctHH_dR2 = d2fctHH_dR2,
            .d2fctHH_dRth = d2fctHH_dRth,
            .d2fctHH_dth2 = d2fctHH_dth2,

            .psi4_2 = psi4_2,
            .dpsi4_2_dR = dpsi4_2_dR,
            .dpsi4_2_dth = dpsi4_2_dth,
            .d2psi4_2_dR2 = d2psi4_2_dR2,
            .d2psi4_2_dRth = d2psi4_2_dRth,
            .d2psi4_2_dth2 = d2psi4_2_dth2,

            .alpha0 = alpha0,
            .alpha02 = alpha02,
            .dalpha02_dR = dalpha02_dR,
            .dalpha02_dth = dalpha02_dth,
            .d2alpha02_dR2 = d2alpha02_dR2,
            .d2alpha02_dRth = d2alpha02_dRth,
            .d2alpha02_dth2 = d2alpha02_dth2,

            .gphiphi = gphiphi,
            .dgphiphi_dR = dgphiphi_dR,
            .dgphiphi_dth = dgphiphi_dth,
            .d2gphiphi_dR2 = d2gphiphi_dR2,
            .d2gphiphi_dRth = d2gphiphi_dRth,
            .d2gphiphi_dth2 = d2gphiphi_dth2,

            .bphiup = bphiup,
            .dbphiup_dR = dbphiup_dR,
            .dbphiup_dth = dbphiup_dth,
            .d2bphiup_dR2 = d2bphiup_dR2,
            .d2bphiup_dRth = d2bphiup_dRth,
            .d2bphiup_dth2 = d2bphiup_dth2,

            .bphi = bphi,
            .dbphi_dR = dbphi_dR,
            .dbphi_dth = dbphi_dth,
            .d2bphi_dR2 = d2bphi_dR2,
            .d2bphi_dRth = d2bphi_dRth,
            .d2bphi_dth2 = d2bphi_dth2,

            // Cartesian second derivatives
            .d2rBL_dx2 = d2rBL_dx2,
            .d2rBL_dy2 = d2rBL_dy2,
            .d2rBL_dz2 = d2rBL_dz2,
            .d2rBL_dxy = d2rBL_dxy,
            .d2rBL_dxz = d2rBL_dxz,
            .d2rBL_dyz = d2rBL_dyz,

            .d2Delt_dx2 = d2Delt_dx2,
            .d2Delt_dy2 = d2Delt_dy2,
            .d2Delt_dz2 = d2Delt_dz2,
            .d2Delt_dxy = d2Delt_dxy,
            .d2Delt_dxz = d2Delt_dxz,
            .d2Delt_dyz = d2Delt_dyz,

            .d2Sigm_dx2 = d2Sigm_dx2,
            .d2Sigm_dy2 = d2Sigm_dy2,
            .d2Sigm_dz2 = d2Sigm_dz2,
            .d2Sigm_dxy = d2Sigm_dxy,
            .d2Sigm_dxz = d2Sigm_dxz,
            .d2Sigm_dyz = d2Sigm_dyz,

            .d2fctFF_dx2 = d2fctFF_dx2,
            .d2fctFF_dy2 = d2fctFF_dy2,
            .d2fctFF_dz2 = d2fctFF_dz2,
            .d2fctFF_dxy = d2fctFF_dxy,
            .d2fctFF_dxz = d2fctFF_dxz,
            .d2fctFF_dyz = d2fctFF_dyz,

            .d2fctGG_dx2 = d2fctGG_dx2,
            .d2fctGG_dy2 = d2fctGG_dy2,
            .d2fctGG_dz2 = d2fctGG_dz2,
            .d2fctGG_dxy = d2fctGG_dxy,
            .d2fctGG_dxz = d2fctGG_dxz,
            .d2fctGG_dyz = d2fctGG_dyz,

            .d2fctHH_dx2 = d2fctHH_dx2,
            .d2fctHH_dy2 = d2fctHH_dy2,
            .d2fctHH_dz2 = d2fctHH_dz2,
            .d2fctHH_dxy = d2fctHH_dxy,
            .d2fctHH_dxz = d2fctHH_dxz,
            .d2fctHH_dyz = d2fctHH_dyz,

            .d2psi4_2_dx2 = d2psi4_2_dx2,
            .d2psi4_2_dy2 = d2psi4_2_dy2,
            .d2psi4_2_dz2 = d2psi4_2_dz2,
            .d2psi4_2_dxy = d2psi4_2_dxy,
            .d2psi4_2_dxz = d2psi4_2_dxz,
            .d2psi4_2_dyz = d2psi4_2_dyz,

            .d2alpha02_dx2 = d2alpha02_dx2,
            .d2alpha02_dy2 = d2alpha02_dy2,
            .d2alpha02_dz2 = d2alpha02_dz2,
            .d2alpha02_dxy = d2alpha02_dxy,
            .d2alpha02_dxz = d2alpha02_dxz,
            .d2alpha02_dyz = d2alpha02_dyz,

            .d2gphiphi_dx2 = d2gphiphi_dx2,
            .d2gphiphi_dy2 = d2gphiphi_dy2,
            .d2gphiphi_dz2 = d2gphiphi_dz2,
            .d2gphiphi_dxy = d2gphiphi_dxy,
            .d2gphiphi_dxz = d2gphiphi_dxz,
            .d2gphiphi_dyz = d2gphiphi_dyz,

            .d2bphiup_dx2 = d2bphiup_dx2,
            .d2bphiup_dy2 = d2bphiup_dy2,
            .d2bphiup_dz2 = d2bphiup_dz2,
            .d2bphiup_dxy = d2bphiup_dxy,
            .d2bphiup_dxz = d2bphiup_dxz,
            .d2bphiup_dyz = d2bphiup_dyz,

            .d2bphi_dx2 = d2bphi_dx2,
            .d2bphi_dy2 = d2bphi_dy2,
            .d2bphi_dz2 = d2bphi_dz2,
            .d2bphi_dxy = d2bphi_dxy,
            .d2bphi_dxz = d2bphi_dxz,
            .d2bphi_dyz = d2bphi_dyz,
        };

        UAv_MetricDerivs1 D1;
        UAv_MetricDerivs2 D2;
        UAv_ComputeMetricDerivsAtPoint(&P, &D1, &D2);

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              dG[a][b][c] = D1.dG[a][b][c];
            }
          }
        }

        CCTK_REAL ddG[4][4][4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                ddG[a][b][c][d] = D2.ddG[a][b][c][d];
              }
            }
          }
        }

        /* NaN/Inf checks for dG tensor */
        for (int aa = 0; aa < 4; ++aa) {
          for (int bb = 0; bb < 4; ++bb) {
            for (int cc = 0; cc < 4; ++cc) {
              char name[32];
              snprintf(name, sizeof(name), "dG[%d][%d][%d]", aa, bb, cc);
              check_nan_or_inf(name, dG[aa][bb][cc]);
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

        CCTK_REAL dg[4][4][4], ddg[4][4][4][4]; // dg[i][j][k] = \partial_k g_{ij} (boosted metric)

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

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                CCTK_REAL sum = 0.0;
                for (int mu = 0; mu < 4; ++mu)
                  for (int nu = 0; nu < 4; ++nu)
                    for (int lam = 0; lam < 4; ++lam)
                      for (int sig = 0; sig < 4; ++sig)
                        sum += invLambda[mu][a] * invLambda[nu][b] * invLambda[lam][c] * invLambda[sig][d] * ddG[mu][nu][lam][sig];
                ddg[a][b][c][d] = sum;
              }
            }
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

        // CCTK_REAL G00up = -((psi4_2 * (rho2_2 * rho2_2) * (1 + bh_spin2 *
        // fctHH * (rho2_2))) / (alpha02 * psi4_2 * (rho2_2 * rho2_2) * (1 +
        // bh_spin2 * fctHH * (rho2_2)) + bphi * (bphi * (rho2_2)-bphiup *
        // psi4_2 * (rho2_2 * rho2_2) * (1 + bh_spin2 * fctHH * (rho2_2))))); //
        // unboosted G^{00}

        // CCTK_REAL G0xup = -((bphi * rho2_2 * y1_2) / (alpha02 * psi4_2 *
        // (rho2_2 * rho2_2) * (1 + bh_spin2 * fctHH * (rho2_2)) + bphi * (bphi
        // * (rho2_2)-bphiup * psi4_2 * (rho2_2 * rho2_2) * (1 + bh_spin2 *
        // fctHH * (rho2_2))))); // unboosted G^{0x}

        // CCTK_REAL Gxxup = (alpha02 * psi4_2 * (rho2_2 * rho2_2) * (1 + fctGG
        // * (y1_2 * y1_2 + z1_2 * z1_2) + bh_spin2 * fctHH * x1_2 * x1_2 *
        // gamma2 * (1 + fctGG * z1_2 * z1_2)) + bphi * (bphi * x1_2 * x1_2 *
        // gamma2 * (1 + fctGG * z1_2 * z1_2) - bphiup * psi4_2 * (rho2_2 *
        // rho2_2) * (1 + fctGG * (y1_2 * y1_2 + z1_2 * z1_2) + bh_spin2 * fctHH
        // * x1_2 * x1_2 * gamma2 * (1 + fctGG * z1_2 * z1_2)))) / (psi4_2 *
        // (alpha02 * psi4_2 * (rho2_2 * rho2_2) * (1 + bh_spin2 * fctHH *
        // (rho2_2)) + bphi * (bphi * (rho2_2)-bphiup * psi4_2 * (rho2_2 *
        // rho2_2) * (1 + bh_spin2 * fctHH * (rho2_2)))) * (1 + fctGG *
        // (rr2_2))); // unboosted G^{xx}

        // CCTK_REAL Gb00up = gamma2 * bh_v2 * Gxxup + gamma2 * G00up - 2.0 *
        // gamma2 * bh_v * G0xup; // boosted Gb^{00} CCTK_REAL new_alpha = 1.0 /
        // sqrt(-Gb00up); // older version worked better

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

        // Optional checks
        check_nan_or_inf("Gb3_inv[1][1]", Gb3_inv[1][1]);
        check_nan_or_inf("Gb3_inv[1][2]", Gb3_inv[1][2]);
        check_nan_or_inf("Gb3_inv[1][3]", Gb3_inv[1][3]);
        check_nan_or_inf("Gb3_inv[2][2]", Gb3_inv[2][2]);
        check_nan_or_inf("Gb3_inv[2][3]", Gb3_inv[2][3]);
        check_nan_or_inf("Gb3_inv[3][3]", Gb3_inv[3][3]);

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

        /* Added NaN/Inf checks for beta and betaup */
        // check_nan_or_inf("beta[0]", beta[0]);
        check_nan_or_inf("beta[1]", beta[1]);
        check_nan_or_inf("beta[2]", beta[2]);
        check_nan_or_inf("beta[3]", beta[3]);
        check_nan_or_inf("betaup[0]", betaup[0]);
        check_nan_or_inf("betaup[1]", betaup[1]);
        check_nan_or_inf("betaup[2]", betaup[2]);
        check_nan_or_inf("betaup[3]", betaup[3]);

        //////////////////////////////////////////////////////////////////////////////////////////

        CCTK_REAL dW_drho_1, dW_dz_1;
        const CCTK_REAL exp_auxi_1 = exp(2. * F2_1[ind] - F0_1[ind]);

        if (rho_1 < 1e-8) {
          dW_drho_1 = 0.;
          dW_dz_1 = 0.;
        } else {
          dW_drho_1 =
              rho_1 / rr_1 * dW_dr_1[ind] + z1_1 / rr2_1 * dW_dth_1[ind];
          dW_dz_1 = z1_1 / rr_1 * dW_dr_1[ind] - rho_1 / rr2_1 * dW_dth_1[ind];
        }

        CCTK_REAL new_alpha2 = -Gb[0][0] + betaup[1] * beta[1] +
                               betaup[2] * beta[2] + betaup[3] * beta[3];
        CCTK_REAL new_alpha = sqrt(new_alpha2);

        {
          if (new_alpha2 < 0) {
            CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "alpha^2=%g at (%lf,%lf,%lf).", (double)new_alpha2, x1_2,
                       y1_2, z1_2);
            break;
          }
        }

        // check_nan_or_inf("new_alpha", new_alpha);

        /////////////////////////////////////////////////////////////////////
        // no horizonte:

        CCTK_REAL rhor = rBLp / 4.0;

        /*aux quantities at horizon*/
        CCTK_REAL hor_x1_2 = gamma * rhor * cosph_2 * sinth;
        CCTK_REAL hor_y1_2 = rhor * sinph_2 * sinth;
        CCTK_REAL hor_z1_2 = rhor * costh;
        CCTK_REAL hor_rr2_2 = hor_x1_2 * hor_x1_2 * gamma2 + hor_y1_2 * hor_y1_2 + hor_z1_2 * hor_z1_2;

        CCTK_REAL hor_rr_2 = sqrt(hor_rr2_2);

        CCTK_REAL hor_rho2_2 = hor_x1_2 * hor_x1_2 * gamma2 + hor_y1_2 * hor_y1_2;

        CCTK_REAL hor_rho_2 = sqrt(hor_rho2_2);

        CCTK_REAL hor_rho3_2 = hor_rho2_2 * hor_rho_2;

        // R0pert2 = (hor_rr_2 - R0pert)*(hor_rr_2 - R0pert) ;

        CCTK_REAL hor_R_x = hor_x1_2 * gamma / hor_rr_2;
        CCTK_REAL hor_R_y = hor_y1_2 / hor_rr_2;
        CCTK_REAL hor_R_z = hor_z1_2 / hor_rr_2;
        CCTK_REAL hor_x_R = hor_x1_2 * gamma / hor_rr_2;
        CCTK_REAL hor_y_R = hor_y1_2 / hor_rr_2;
        CCTK_REAL hor_z_R = hor_z1_2 / hor_rr_2;

        CCTK_REAL hor_costh = hor_z1_2 / hor_rr_2;
        CCTK_REAL hor_costh2 = hor_costh * hor_costh;
        CCTK_REAL hor_sinth2 = 1. - hor_costh2;
        CCTK_REAL hor_sinth = sqrt(hor_sinth2);

        CCTK_REAL hor_th_x = hor_costh * hor_R_x / hor_rho_2;
        CCTK_REAL hor_th_y = hor_costh * hor_R_y / hor_rho_2;
        CCTK_REAL hor_th_z = -hor_rho_2 / hor_rr2_2;

        // auxilliary quantities
        CCTK_REAL hor_rBL, hor_rBL2;
        hor_rBL = hor_rr_2 * (1.0 + 0.25 * rBLp / hor_rr_2) * (1.0 + 0.25 * rBLp / hor_rr_2);
        hor_rBL2 = hor_rBL * hor_rBL;
        CCTK_REAL hor_Delt, hor_Sigm, hor_Sigm2, hor_fctFF;
        hor_Delt = (hor_rBL - rBLp) * (hor_rBL - rBLm); // hor_rBL2 + bh_spin2 - 2 * bh_mass * hor_rBL;
        hor_Sigm = hor_rBL2 + bh_spin2 * hor_costh2;
        hor_Sigm2 = hor_Sigm * hor_Sigm;
        hor_fctFF = (hor_rBL2 + bh_spin2) * (hor_rBL2 + bh_spin2) - hor_Delt * bh_spin2 * hor_sinth2; // "A" no artigo

        CCTK_REAL hor_psi4_2 = hor_Sigm / hor_rr2_2; // psi04 no codigo original
        CCTK_REAL hor_psi2_2 = sqrt(hor_psi4_2);
        CCTK_REAL hor_psi1_2 = sqrt(hor_psi2_2);

        CCTK_REAL hor_fctGG, hor_fctHH;
        hor_fctGG = rBLm / (hor_rr2_2 * (hor_rBL - rBLm));
        hor_fctHH = (2.0 * bh_mass * hor_rBL + hor_Sigm) / (hor_rr2_2 * hor_Sigm2);

        // CCTK_REAL detgij;
        // detgij = hor_psi4_2 * hor_psi4_2 * hor_psi4_2 * (1.0 + hor_rr2_2 * hor_fctGG) * (1.0 +
        // bh_spin2 * hor_rho2_2 * hor_fctHH);

        /*----------------------------------*/

        //  CCTK_REAL alpha02 = (4.0 * hor_rr_2 - rBLp) * (4.0 * hor_rr_2 - rBLp) *
        // (hor_rBL - rBLm) / (16.0 * hor_rr_2 * (hor_rBL2 + bh_spin2 * (1.0 + 2.0 * bh_mass
        // * hor_rBL * hor_sinth2 / hor_Sigm)));
        CCTK_REAL hor_alpha0 = (4.0 * hor_rr_2 - rBLp) * sqrt(hor_rBL - rBLm) / sqrt(16.0 * hor_rr_2 * (hor_rBL2 + bh_spin2 * (1.0 + 2.0 * bh_mass * hor_rBL * hor_sinth2 / hor_Sigm))); // primeiro termo para schwarzschild e zero;
        CCTK_REAL hor_alpha02 = hor_alpha0 * hor_alpha0;
        //  CCTK_REAL hor_alpha0 = sqrt(alpha02);
        //  CCTK_REAL hor_alpha0 = sqrt(hor_Delt * hor_Sigm / hor_fctFF);
        //  CCTK_REAL alpha02 = hor_Delt * hor_Sigm / hor_fctFF;
        //  CCTK_REAL hor_alpha0 = sqrt(alpha02);
        //  CCTK_REAL dalpha0_dR = 0.5 / hor_alpha0 * (-(hor_Delt * hor_Sigm *
        // hor_dfctFF_dR) + hor_fctFF * (hor_Sigm * hor_dpsi4_2_dth + hor_Delt * hor_dSigm_dR)) /
        // pow(hor_fctFF, 2);
        //  CCTK_REAL dalpha0_dth = 0.5 / hor_alpha0 * (hor_Delt *
        // (-(hor_Sigm * hor_dfctFF_dth) + hor_fctFF * hor_dSigm_dth)) / pow(hor_fctFF, 2);
        //  CCTK_REAL hor_dalpha02_dR =(-(hor_Delt * hor_Sigm * hor_dfctFF_dR) + hor_fctFF * (hor_Sigm * hor_dpsi4_2_dth + hor_Delt * hor_dSigm_dR)) /(hor_fctFF * hor_fctFF);
        //  CCTK_REAL hor_dalpha02_dth = (hor_Delt * (-(hor_Sigm * hor_dfctFF_dth) + hor_fctFF * hor_dSigm_dth)) / (hor_fctFF * hor_fctFF);
        CCTK_REAL hor_gphiphi = hor_fctFF / hor_Sigm * hor_sinth2;
        CCTK_REAL hor_bphiup = -2.0 * bh_mass * bh_spin * hor_rBL / hor_fctFF;
        CCTK_REAL hor_bphi = hor_bphiup * hor_gphiphi;

        // 1st derivatives of auxiliary quantities
        CCTK_REAL hor_drBLdR, hor_dfctGG_dR, hor_dfctHH_dR, hor_dfctHH_dth, hor_dpsi4_2_dR, hor_dDelt_dR, hor_dpsi4_2_dth, hor_dSigm_dR, hor_dSigm_dth, hor_dfctFF_dR, hor_dfctFF_dth;
        hor_drBLdR = 1.0 - rBLp * rBLp / (16.0 * hor_rr2_2);
        hor_dDelt_dR = hor_drBLdR * (hor_rBL - rBLp) + hor_drBLdR * (hor_rBL - rBLm); //(2 * hor_rBL - 2 * bh_mass) * hor_drBLdR;
        hor_dSigm_dR = 2 * hor_rBL * hor_drBLdR;
        hor_dSigm_dth = -2 * bh_spin2 * hor_costh * hor_sinth;
        hor_dfctFF_dR = 4 * hor_rBL * (bh_spin2 + hor_rBL2) * hor_drBLdR - bh_spin2 * hor_sinth2 * hor_dDelt_dR;
        hor_dfctFF_dth = -(bh_spin2 * hor_Delt * 2 * hor_costh * hor_sinth);
        hor_dfctGG_dR = -((rBLm * (-2 * rBLm + 2 * hor_rBL + hor_rr_2 * hor_drBLdR)) / (hor_rr2_2 * hor_rr_2 * (rBLm - hor_rBL) * (rBLm - hor_rBL)));
        hor_dfctHH_dR = -((2 * hor_Sigm * (2 * bh_mass * hor_rBL + hor_Sigm - bh_mass * hor_rr_2 * hor_drBLdR) + hor_rr_2 * (4 * bh_mass * hor_rBL + hor_Sigm) * hor_dSigm_dR) / (hor_rr2_2 * hor_rr_2 * hor_Sigm2 * hor_Sigm));
        hor_dfctHH_dth = -(((4 * bh_mass * hor_rBL + hor_Sigm) * hor_dSigm_dth) / (hor_rr2_2 * hor_Sigm2 * hor_Sigm));
        hor_dpsi4_2_dR = (-2 * hor_Sigm + hor_rr_2 * hor_dSigm_dR) / (hor_rr2_2 * hor_rr_2);
        hor_dpsi4_2_dth = hor_dSigm_dth / hor_rr2_2;
        CCTK_REAL hor_dalpha02_dR = ((4 * hor_rr_2 - rBLp) * (2 * bh_mass * bh_spin2 * hor_Sigm * hor_sinth2 * (-((4 * hor_rr_2 + rBLp) * (rBLm - hor_rBL) * hor_rBL) + hor_rr_2 * rBLm * (4 * hor_rr_2 - rBLp) * hor_drBLdR) + hor_Sigm2 * (-((4 * hor_rr_2 + rBLp) * (rBLm - hor_rBL) * (bh_spin2 + hor_rBL2)) + hor_rr_2 * (4 * hor_rr_2 - rBLp) * (bh_spin2 + 2 * rBLm * hor_rBL - hor_rBL2) * hor_drBLdR) - 2 * bh_mass * bh_spin2 * hor_rr_2 * (4 * hor_rr_2 - rBLp) * (rBLm - hor_rBL) * hor_rBL * hor_sinth2 * hor_dSigm_dR)) / (16. * (hor_rr_2 * (bh_spin2 + hor_rBL2) * hor_Sigm + 2 * bh_mass * bh_spin2 * hor_rr_2 * hor_rBL * hor_sinth2) * (hor_rr_2 * (bh_spin2 + hor_rBL2) * hor_Sigm + 2 * bh_mass * bh_spin2 * hor_rr_2 * hor_rBL * hor_sinth2));
        CCTK_REAL hor_dalpha02_dth = (bh_mass * bh_spin2 * (-4 * hor_rr_2 + rBLp) * (-4 * hor_rr_2 + rBLp) * (rBLm - hor_rBL) * hor_rBL * hor_sinth * (2 * hor_costh * hor_Sigm - hor_sinth * hor_dSigm_dth)) / (8. * hor_rr_2 * ((bh_spin2 + hor_rBL2) * hor_Sigm + 2 * bh_mass * bh_spin2 * hor_rBL * hor_sinth2) * ((bh_spin2 + hor_rBL2) * hor_Sigm + 2 * bh_mass * bh_spin2 * hor_rBL * hor_sinth2));
        CCTK_REAL hor_dgphiphi_dR = (hor_sinth2 * (hor_Sigm * hor_dfctFF_dR - hor_fctFF * hor_dSigm_dR)) / hor_Sigm2;
        CCTK_REAL hor_dgphiphi_dth = (hor_sinth * (hor_Sigm * hor_sinth * hor_dfctFF_dth + hor_fctFF * (2 * hor_costh * hor_Sigm - hor_sinth * hor_dSigm_dth))) / hor_Sigm2;
        CCTK_REAL hor_dbphiup_dR = (2 * bh_mass * bh_spin * (-(hor_fctFF * hor_drBLdR) + hor_rBL * hor_dfctFF_dR)) / (hor_fctFF * hor_fctFF);
        CCTK_REAL hor_dbphiup_dth = (2 * bh_mass * bh_spin * hor_rBL * hor_dfctFF_dth) / (hor_fctFF * hor_fctFF);
        CCTK_REAL hor_dbphi_dR = hor_gphiphi * hor_dbphiup_dR + hor_bphiup * hor_dgphiphi_dR;
        CCTK_REAL hor_dbphi_dth = hor_gphiphi * hor_dbphiup_dth + hor_bphiup * hor_dgphiphi_dth;

        // 2nd derivatives
        /* 2nd derivatives (Cartesian coordinates) of auxiliary quantities.
           Skeleton for user to fill in.
           Notation suggestion:
             d2<name>_dx2  = ∂^2(<name>)/∂x^2
             d2<name>_dxy  = ∂^2(<name>)/∂x∂y
             etc.

           Chain rule pattern for a scalar F(R,th):
             Fx  = F_R * hor_R_x + F_th * hor_th_x
             Fxx = F_RR * hor_R_x*hor_R_x + 2*F_Rth * hor_R_x*hor_th_x + F_thth * hor_th_x*hor_th_x
             + F_R * hor_R_xx + F_th * hor_th_xx
             (analogous for mixed and other second derivatives)

           Provide (and fill) second derivatives of R and theta if needed:
             hor_R_xx, hor_R_xy, hor_R_xz, hor_R_yy, hor_R_yz, hor_R_zz
             hor_th_xx, hor_th_xy, hor_th_xz, hor_th_yy, hor_th_yz, hor_th_zz

           Below everything initialized to 0.0 for later replacement.
        */

        /* Coordinate second derivatives (fill if needed) */
        CCTK_REAL hor_R_xx = hor_y1_2 * hor_z1_2 / (hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_R_xy = -hor_y1_2 * hor_x1_2 * gamma / (hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_R_xz = -hor_x1_2 * gamma * hor_z1_2 / (hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_R_yy = hor_x1_2 * gamma * hor_z1_2 / (hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_R_yz = -hor_y1_2 * hor_z1_2 / (hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_R_zz = hor_rho2_2 / (hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_th_xx = hor_z1_2 * (-2 * pow(hor_x1_2 * gamma, 4) - pow(hor_x1_2 * gamma * hor_y1_2, 2) + pow(hor_y1_2, 4) + pow(hor_y1_2 * hor_z1_2, 2)) / (hor_rho3_2 * hor_rr2_2 * hor_rr2_2);
        CCTK_REAL hor_th_xy = -hor_x1_2 * gamma * hor_y1_2 * hor_z1_2 * (3 * hor_rho2_2 + hor_z1_2 * hor_z1_2) / (hor_rho3_2 * hor_rr2_2 * hor_rr2_2);
        CCTK_REAL hor_th_xz = hor_x1_2 * gamma * (hor_rho2_2 - hor_z1_2 * hor_z1_2) / (hor_rho3_2 * hor_rr2_2 * hor_rr2_2);
        CCTK_REAL hor_th_yy = hor_z1_2 * (pow(hor_x1_2 * gamma, 4) - 2 * pow(hor_y1_2, 4) + hor_x1_2 * (-hor_y1_2 * hor_y1_2 + hor_z1_2 * hor_z1_2)) / (hor_rho3_2 * hor_rr2_2 * hor_rr2_2);
        CCTK_REAL hor_th_yz = hor_y1_2 * (hor_rho2_2 - hor_z1_2 * hor_z1_2) / (hor_rho3_2 * hor_rr2_2 * hor_rr2_2);
        CCTK_REAL hor_th_zz = 2 * hor_z1_2 * hor_rho_2 / hor_rr2_2;

        /* Second derivatives in (R,theta) space (to be provided) */
        CCTK_REAL hor_d2rBL_dR2 = rBLp * rBLp / (8 * hor_rr2_2 * hor_rr_2);
        CCTK_REAL hor_d2rBL_dRth = 0.0;
        CCTK_REAL hor_d2rBL_dth2 = 0.0;
        CCTK_REAL hor_d2Delt_dR2 = 2 * (pow(hor_drBLdR, 2) + (-bh_mass + hor_rBL) * hor_d2rBL_dR2);
        CCTK_REAL hor_d2Delt_dRth = 0.0;
        CCTK_REAL hor_d2Delt_dth2 = 0.0;
        CCTK_REAL hor_d2Sigm_dR2 = 2 * (pow(hor_drBLdR, 2) + hor_rBL * hor_d2rBL_dR2);
        CCTK_REAL hor_d2Sigm_dRth = 0.0;
        CCTK_REAL hor_d2Sigm_dth2 = -2 * bh_spin2 * (hor_costh2 - hor_sinth2);
        CCTK_REAL hor_d2fctFF_dR2 = 4 * (bh_spin2 + 3 * hor_rBL2) * pow(hor_drBLdR, 2) - bh_spin2 * hor_sinth2 * hor_d2Delt_dR2 + 4 * hor_rBL * (bh_spin2 + hor_rBL2) * hor_d2rBL_dR2;
        CCTK_REAL hor_d2fctFF_dRth = -2 * bh_spin2 * hor_costh * hor_sinth * hor_dDelt_dR;
        CCTK_REAL hor_d2fctFF_dth2 = -2 * bh_spin2 * (hor_costh2 - hor_sinth2) * hor_Delt;
        CCTK_REAL hor_d2fctGG_dR2 = (rBLm * (-6 * pow(rBLm - hor_rBL, 2) + 4 * hor_rr_2 * (rBLm - hor_rBL) * hor_drBLdR - hor_rr2_2 * (2 * pow(hor_drBLdR, 2) + (rBLm - hor_rBL) * hor_d2rBL_dR2))) / (pow(hor_rr_2, 4) * pow(rBLm - hor_rBL, 3));
        CCTK_REAL hor_d2fctGG_dRth = 0.0;
        CCTK_REAL hor_d2fctGG_dth2 = 0.0;
        CCTK_REAL hor_d2fctHH_dR2 = (2 * (hor_Sigm * (3 * hor_Sigm2 + bh_mass * hor_rr_2 * hor_Sigm * (-4 * hor_drBLdR + hor_rr_2 * hor_d2rBL_dR2) - 2 * bh_mass * pow(hor_rr_2, 2) * hor_drBLdR * hor_dSigm_dR) + bh_mass * hor_rBL * (6 * hor_Sigm2 + 2 * hor_rr2_2 * pow(hor_dSigm_dR, 2) + hor_rr_2 * hor_Sigm * (4 * hor_dSigm_dR - hor_rr_2 * hor_d2Sigm_dR2)))) / (pow(hor_rr_2, 4) * pow(hor_Sigm, 3));
        CCTK_REAL hor_d2fctHH_dRth = (2 * bh_mass * hor_dSigm_dth * (-(hor_rr_2 * hor_Sigm * hor_drBLdR) + 2 * hor_rBL * (hor_Sigm + hor_rr_2 * hor_dSigm_dR)) - 2 * bh_mass * hor_rr_2 * hor_rBL * hor_Sigm * hor_d2Sigm_dRth) / (pow(hor_rr_2, 3) * pow(hor_Sigm, 3));
        CCTK_REAL hor_d2fctHH_dth2 = (-2 * bh_mass * hor_rBL * (-2 * pow(hor_dSigm_dth, 2) + hor_Sigm * hor_d2Sigm_dR2)) / (hor_rr2_2 * pow(hor_Sigm, 3));
        CCTK_REAL hor_d2psi4_2_dR2 = (6 * hor_Sigm - 4 * hor_rr_2 * hor_dSigm_dR + hor_rr2_2 * hor_d2Sigm_dR2) / pow(hor_rr_2, 4);
        CCTK_REAL hor_d2psi4_2_dRth = (-2 * hor_dSigm_dth + hor_rr_2 * hor_d2Sigm_dRth) / pow(hor_rr_2, 3);
        CCTK_REAL hor_d2psi4_2_dth2 = hor_d2Sigm_dth2 / hor_rr2_2;
        CCTK_REAL hor_d2alpha02_dR2 = (2 * hor_Delt * hor_Sigm * pow(hor_dfctFF_dR, 2) + hor_fctFF * (-2 * hor_dfctFF_dR * (hor_Sigm * hor_dDelt_dR + hor_Delt * hor_dSigm_dR) - hor_Delt * hor_Sigm * hor_d2fctFF_dR2) + pow(hor_fctFF, 2) * (hor_Sigm * hor_d2Delt_dR2 + 2 * hor_dDelt_dR * hor_dSigm_dR + hor_Delt * hor_d2Sigm_dR2)) / pow(hor_fctFF, 3);
        CCTK_REAL hor_d2alpha02_dRth = (2 * hor_Delt * hor_Sigm * hor_dfctFF_dth * hor_dfctFF_dR - hor_fctFF * (hor_Delt * hor_dSigm_dth * hor_dfctFF_dR + hor_dfctFF_dth * (hor_Sigm * hor_dDelt_dR + hor_Delt * hor_dSigm_dR) + hor_Delt * hor_Sigm * hor_d2fctFF_dRth) + pow(hor_fctFF, 2) * (hor_dDelt_dR * hor_dSigm_dth + hor_Delt * hor_d2Sigm_dRth)) / pow(hor_fctFF, 3);
        CCTK_REAL hor_d2alpha02_dth2 = (hor_Delt * (hor_Sigm * (2 * pow(hor_dfctFF_dth, 2) - hor_fctFF * hor_d2fctFF_dth2) + hor_fctFF * (-2 * hor_dfctFF_dth * hor_dSigm_dth + hor_fctFF * hor_d2Sigm_dth2))) / pow(hor_fctFF, 3);
        CCTK_REAL hor_d2gphiphi_dR2 = (hor_sinth2 * (2 * hor_fctFF * pow(hor_dSigm_dR, 2) + hor_Sigm2 * hor_d2fctFF_dR2 - hor_Sigm * (2 * hor_dfctFF_dR * hor_dSigm_dR + hor_fctFF * hor_d2Sigm_dR2))) / pow(hor_Sigm, 3);
        CCTK_REAL hor_d2gphiphi_dRth = (hor_sinth * (2 * hor_fctFF * hor_sinth * hor_dSigm_dth * hor_dSigm_dR + pow(hor_Sigm, 2) * (2 * hor_costh * hor_dfctFF_dR + hor_sinth * hor_d2fctFF_dRth) - hor_Sigm * (hor_sinth * hor_dSigm_dth * hor_dfctFF_dR + (2 * hor_costh * hor_fctFF + hor_sinth * hor_dfctFF_dth) * hor_dSigm_dR + hor_fctFF * hor_sinth * hor_d2Sigm_dRth))) / pow(hor_Sigm, 3);
        CCTK_REAL hor_d2gphiphi_dth2 = (hor_Sigm * hor_sinth * (2 * hor_dfctFF_dth * (2 * hor_costh * hor_Sigm - hor_sinth * hor_dSigm_dth) + hor_Sigm * hor_sinth * hor_d2fctFF_dth2) + hor_fctFF * (2 * (hor_costh2 - hor_sinth2) * hor_Sigm2 + 2 * hor_sinth2 * pow(hor_dSigm_dth, 2) - hor_Sigm * hor_sinth * (4 * hor_costh * hor_dSigm_dth + hor_sinth * hor_d2Sigm_dth2))) / pow(hor_Sigm, 3);
        CCTK_REAL hor_d2bphiup_dR2 = (-2 * bh_mass * bh_spin * (pow(hor_fctFF, 2) * hor_d2rBL_dR2 + 2 * hor_rBL * pow(hor_dfctFF_dR, 2) - hor_fctFF * (2 * hor_drBLdR * hor_dfctFF_dR + hor_rBL * hor_d2fctFF_dR2))) / pow(hor_fctFF, 3);
        CCTK_REAL hor_d2bphiup_dRth = (2 * bh_mass * bh_spin * (hor_dfctFF_dth * (hor_fctFF * hor_drBLdR - 2 * hor_rBL * hor_dfctFF_dR) + hor_fctFF * hor_rBL * hor_d2fctFF_dRth)) / pow(hor_fctFF, 3);
        CCTK_REAL hor_d2bphiup_dth2 = (2 * bh_mass * bh_spin * hor_rBL * (-2 * pow(hor_dfctFF_dth, 2) + hor_fctFF * hor_d2fctFF_dth2)) / pow(hor_fctFF, 3);
        CCTK_REAL hor_d2bphi_dR2 = 2 * hor_dbphiup_dR * hor_dgphiphi_dR + hor_gphiphi * hor_d2bphiup_dR2 + hor_bphiup * hor_d2gphiphi_dR2;
        CCTK_REAL hor_d2bphi_dRth = hor_dgphiphi_dth * hor_dbphiup_dR + hor_dbphiup_dth * hor_dgphiphi_dR + hor_gphiphi * hor_d2bphiup_dRth + hor_bphiup * hor_d2gphiphi_dRth;
        CCTK_REAL hor_d2bphi_dth2 = 2 * hor_dbphiup_dth * hor_dgphiphi_dth + hor_gphiphi * hor_d2bphiup_dth2 + hor_bphiup * hor_d2gphiphi_dth2;

        /* Cartesian second derivatives (fill using chain rule) */
        /* hor_rBL */
        /* hor_rBL (depends only on R) */
        CCTK_REAL hor_d2rBL_dx2 = hor_d2rBL_dR2 * hor_R_x * hor_R_x + hor_drBLdR * hor_R_xx;
        CCTK_REAL hor_d2rBL_dy2 = hor_d2rBL_dR2 * hor_R_y * hor_R_y + hor_drBLdR * hor_R_yy;
        CCTK_REAL hor_d2rBL_dz2 = hor_d2rBL_dR2 * hor_R_z * hor_R_z + hor_drBLdR * hor_R_zz;
        CCTK_REAL hor_d2rBL_dxy = hor_d2rBL_dR2 * hor_R_x * hor_R_y + hor_drBLdR * hor_R_xy;
        CCTK_REAL hor_d2rBL_dxz = hor_d2rBL_dR2 * hor_R_x * hor_R_z + hor_drBLdR * hor_R_xz;
        CCTK_REAL hor_d2rBL_dyz = hor_d2rBL_dR2 * hor_R_y * hor_R_z + hor_drBLdR * hor_R_yz;

        /* hor_Delt */
        CCTK_REAL hor_dDelt_dth = 0.0; // hor_Delt depende apenas de R
        CCTK_REAL hor_d2Delt_dx2 = hor_d2Delt_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2Delt_dRth * hor_R_x * hor_th_x + hor_d2Delt_dth2 * hor_th_x * hor_th_x + hor_dDelt_dR * hor_R_xx + dDelt_dth * hor_th_xx;
        CCTK_REAL hor_d2Delt_dy2 = hor_d2Delt_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2Delt_dRth * hor_R_y * hor_th_y + hor_d2Delt_dth2 * hor_th_y * hor_th_y + hor_dDelt_dR * hor_R_yy + dDelt_dth * hor_th_yy;
        CCTK_REAL hor_d2Delt_dz2 = hor_d2Delt_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2Delt_dRth * hor_R_z * hor_th_z + hor_d2Delt_dth2 * hor_th_z * hor_th_z + hor_dDelt_dR * hor_R_zz + dDelt_dth * hor_th_zz;
        CCTK_REAL hor_d2Delt_dxy = hor_d2Delt_dR2 * hor_R_x * hor_R_y + hor_d2Delt_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2Delt_dth2 * hor_th_x * hor_th_y + hor_dDelt_dR * hor_R_xy + dDelt_dth * hor_th_xy;
        CCTK_REAL hor_d2Delt_dxz = hor_d2Delt_dR2 * hor_R_x * hor_R_z + hor_d2Delt_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2Delt_dth2 * hor_th_x * hor_th_z + hor_dDelt_dR * hor_R_xz + dDelt_dth * hor_th_xz;
        CCTK_REAL hor_d2Delt_dyz = hor_d2Delt_dR2 * hor_R_y * hor_R_z + hor_d2Delt_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2Delt_dth2 * hor_th_y * hor_th_z + hor_dDelt_dR * hor_R_yz + dDelt_dth * hor_th_yz;

        /* hor_Sighor_m */
        CCTK_REAL hor_d2Sigm_dx2 = hor_d2Sigm_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2Sigm_dRth * hor_R_x * hor_th_x + hor_d2Sigm_dth2 * hor_th_x * hor_th_x + hor_dSigm_dR * hor_R_xx + hor_dSigm_dth * hor_th_xx;
        CCTK_REAL hor_d2Sigm_dy2 = hor_d2Sigm_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2Sigm_dRth * hor_R_y * hor_th_y + hor_d2Sigm_dth2 * hor_th_y * hor_th_y + hor_dSigm_dR * hor_R_yy + hor_dSigm_dth * hor_th_yy;
        CCTK_REAL hor_d2Sigm_dz2 = hor_d2Sigm_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2Sigm_dRth * hor_R_z * hor_th_z + hor_d2Sigm_dth2 * hor_th_z * hor_th_z + hor_dSigm_dR * hor_R_zz + hor_dSigm_dth * hor_th_zz;
        CCTK_REAL hor_d2Sigm_dxy = hor_d2Sigm_dR2 * hor_R_x * hor_R_y + hor_d2Sigm_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2Sigm_dth2 * hor_th_x * hor_th_y + hor_dSigm_dR * hor_R_xy + hor_dSigm_dth * hor_th_xy;
        CCTK_REAL hor_d2Sigm_dxz = hor_d2Sigm_dR2 * hor_R_x * hor_R_z + hor_d2Sigm_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2Sigm_dth2 * hor_th_x * hor_th_z + hor_dSigm_dR * hor_R_xz + hor_dSigm_dth * hor_th_xz;
        CCTK_REAL hor_d2Sigm_dyz = hor_d2Sigm_dR2 * hor_R_y * hor_R_z + hor_d2Sigm_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2Sigm_dth2 * hor_th_y * hor_th_z + hor_dSigm_dR * hor_R_yz + hor_dSigm_dth * hor_th_yz;

        /* hor_fcthor_FF */
        CCTK_REAL hor_d2fctFF_dx2 = hor_d2fctFF_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2fctFF_dRth * hor_R_x * hor_th_x + hor_d2fctFF_dth2 * hor_th_x * hor_th_x + hor_dfctFF_dR * hor_R_xx + hor_dfctFF_dth * hor_th_xx;
        CCTK_REAL hor_d2fctFF_dy2 = hor_d2fctFF_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2fctFF_dRth * hor_R_y * hor_th_y + hor_d2fctFF_dth2 * hor_th_y * hor_th_y + hor_dfctFF_dR * hor_R_yy + hor_dfctFF_dth * hor_th_yy;
        CCTK_REAL hor_d2fctFF_dz2 = hor_d2fctFF_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2fctFF_dRth * hor_R_z * hor_th_z + hor_d2fctFF_dth2 * hor_th_z * hor_th_z + hor_dfctFF_dR * hor_R_zz + hor_dfctFF_dth * hor_th_zz;
        CCTK_REAL hor_d2fctFF_dxy = hor_d2fctFF_dR2 * hor_R_x * hor_R_y + hor_d2fctFF_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2fctFF_dth2 * hor_th_x * hor_th_y + hor_dfctFF_dR * hor_R_xy + hor_dfctFF_dth * hor_th_xy;
        CCTK_REAL hor_d2fctFF_dxz = hor_d2fctFF_dR2 * hor_R_x * hor_R_z + hor_d2fctFF_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2fctFF_dth2 * hor_th_x * hor_th_z + hor_dfctFF_dR * hor_R_xz + hor_dfctFF_dth * hor_th_xz;
        CCTK_REAL hor_d2fctFF_dyz = hor_d2fctFF_dR2 * hor_R_y * hor_R_z + hor_d2fctFF_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2fctFF_dth2 * hor_th_y * hor_th_z + hor_dfctFF_dR * hor_R_yz + hor_dfctFF_dth * hor_th_yz;

        /* hor_fcthor_GG (only R) */
        CCTK_REAL hor_dfctGG_dth = 0.0;
        CCTK_REAL hor_d2fctGG_dx2 = hor_d2fctGG_dR2 * hor_R_x * hor_R_x + hor_dfctGG_dR * hor_R_xx;
        CCTK_REAL hor_d2fctGG_dy2 = hor_d2fctGG_dR2 * hor_R_y * hor_R_y + hor_dfctGG_dR * hor_R_yy;
        CCTK_REAL hor_d2fctGG_dz2 = hor_d2fctGG_dR2 * hor_R_z * hor_R_z + hor_dfctGG_dR * hor_R_zz;
        CCTK_REAL hor_d2fctGG_dxy = hor_d2fctGG_dR2 * hor_R_x * hor_R_y + hor_dfctGG_dR * hor_R_xy;
        CCTK_REAL hor_d2fctGG_dxz = hor_d2fctGG_dR2 * hor_R_x * hor_R_z + hor_dfctGG_dR * hor_R_xz;
        CCTK_REAL hor_d2fctGG_dyz = hor_d2fctGG_dR2 * hor_R_y * hor_R_z + hor_dfctGG_dR * hor_R_yz;

        /* hor_fcthor_HH */
        CCTK_REAL hor_d2fctHH_dx2 = hor_d2fctHH_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2fctHH_dRth * hor_R_x * hor_th_x + hor_d2fctHH_dth2 * hor_th_x * hor_th_x + hor_dfctHH_dR * hor_R_xx + hor_dfctHH_dth * hor_th_xx;
        CCTK_REAL hor_d2fctHH_dy2 = hor_d2fctHH_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2fctHH_dRth * hor_R_y * hor_th_y + hor_d2fctHH_dth2 * hor_th_y * hor_th_y + hor_dfctHH_dR * hor_R_yy + hor_dfctHH_dth * hor_th_yy;
        CCTK_REAL hor_d2fctHH_dz2 = hor_d2fctHH_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2fctHH_dRth * hor_R_z * hor_th_z + hor_d2fctHH_dth2 * hor_th_z * hor_th_z + hor_dfctHH_dR * hor_R_zz + hor_dfctHH_dth * hor_th_zz;
        CCTK_REAL hor_d2fctHH_dxy = hor_d2fctHH_dR2 * hor_R_x * hor_R_y + hor_d2fctHH_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2fctHH_dth2 * hor_th_x * hor_th_y + hor_dfctHH_dR * hor_R_xy + hor_dfctHH_dth * hor_th_xy;
        CCTK_REAL hor_d2fctHH_dxz = hor_d2fctHH_dR2 * hor_R_x * hor_R_z + hor_d2fctHH_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2fctHH_dth2 * hor_th_x * hor_th_z + hor_dfctHH_dR * hor_R_xz + hor_dfctHH_dth * hor_th_xz;
        CCTK_REAL hor_d2fctHH_dyz = hor_d2fctHH_dR2 * hor_R_y * hor_R_z + hor_d2fctHH_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2fctHH_dth2 * hor_th_y * hor_th_z + hor_dfctHH_dR * hor_R_yz + hor_dfctHH_dth * hor_th_yz;

        /* hor_psihor_4_2 */
        CCTK_REAL hor_d2psi4_2_dx2 = hor_d2psi4_2_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2psi4_2_dRth * hor_R_x * hor_th_x + hor_d2psi4_2_dth2 * hor_th_x * hor_th_x + hor_dpsi4_2_dR * hor_R_xx + hor_dpsi4_2_dth * hor_th_xx;
        CCTK_REAL hor_d2psi4_2_dy2 = hor_d2psi4_2_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2psi4_2_dRth * hor_R_y * hor_th_y + hor_d2psi4_2_dth2 * hor_th_y * hor_th_y + hor_dpsi4_2_dR * hor_R_yy + hor_dpsi4_2_dth * hor_th_yy;
        CCTK_REAL hor_d2psi4_2_dz2 = hor_d2psi4_2_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2psi4_2_dRth * hor_R_z * hor_th_z + hor_d2psi4_2_dth2 * hor_th_z * hor_th_z + hor_dpsi4_2_dR * hor_R_zz + hor_dpsi4_2_dth * hor_th_zz;
        CCTK_REAL hor_d2psi4_2_dxy = hor_d2psi4_2_dR2 * hor_R_x * hor_R_y + hor_d2psi4_2_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2psi4_2_dth2 * hor_th_x * hor_th_y + hor_dpsi4_2_dR * hor_R_xy + hor_dpsi4_2_dth * hor_th_xy;
        CCTK_REAL hor_d2psi4_2_dxz = hor_d2psi4_2_dR2 * hor_R_x * hor_R_z + hor_d2psi4_2_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2psi4_2_dth2 * hor_th_x * hor_th_z + hor_dpsi4_2_dR * hor_R_xz + hor_dpsi4_2_dth * hor_th_xz;
        CCTK_REAL hor_d2psi4_2_dyz = hor_d2psi4_2_dR2 * hor_R_y * hor_R_z + hor_d2psi4_2_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2psi4_2_dth2 * hor_th_y * hor_th_z + hor_dpsi4_2_dR * hor_R_yz + hor_dpsi4_2_dth * hor_th_yz;

        /* alpha02hor_ */
        CCTK_REAL hor_d2alpha02_dx2 = hor_d2alpha02_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2alpha02_dRth * hor_R_x * hor_th_x + hor_d2alpha02_dth2 * hor_th_x * hor_th_x + hor_dalpha02_dR * hor_R_xx + hor_dalpha02_dth * hor_th_xx;
        CCTK_REAL hor_d2alpha02_dy2 = hor_d2alpha02_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2alpha02_dRth * hor_R_y * hor_th_y + hor_d2alpha02_dth2 * hor_th_y * hor_th_y + hor_dalpha02_dR * hor_R_yy + hor_dalpha02_dth * hor_th_yy;
        CCTK_REAL hor_d2alpha02_dz2 = hor_d2alpha02_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2alpha02_dRth * hor_R_z * hor_th_z + hor_d2alpha02_dth2 * hor_th_z * hor_th_z + hor_dalpha02_dR * hor_R_zz + hor_dalpha02_dth * hor_th_zz;
        CCTK_REAL hor_d2alpha02_dxy = hor_d2alpha02_dR2 * hor_R_x * hor_R_y + hor_d2alpha02_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2alpha02_dth2 * hor_th_x * hor_th_y + hor_dalpha02_dR * hor_R_xy + hor_dalpha02_dth * hor_th_xy;
        CCTK_REAL hor_d2alpha02_dxz = hor_d2alpha02_dR2 * hor_R_x * hor_R_z + hor_d2alpha02_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2alpha02_dth2 * hor_th_x * hor_th_z + hor_dalpha02_dR * hor_R_xz + hor_dalpha02_dth * hor_th_xz;
        CCTK_REAL hor_d2alpha02_dyz = hor_d2alpha02_dR2 * hor_R_y * hor_R_z + hor_d2alpha02_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2alpha02_dth2 * hor_th_y * hor_th_z + hor_dalpha02_dR * hor_R_yz + hor_dalpha02_dth * hor_th_yz;

        /* hor_gphhor_iphi */
        CCTK_REAL hor_d2gphiphi_dx2 = hor_d2gphiphi_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2gphiphi_dRth * hor_R_x * hor_th_x + hor_d2gphiphi_dth2 * hor_th_x * hor_th_x + hor_dgphiphi_dR * hor_R_xx + hor_dgphiphi_dth * hor_th_xx;
        CCTK_REAL hor_d2gphiphi_dy2 = hor_d2gphiphi_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2gphiphi_dRth * hor_R_y * hor_th_y + hor_d2gphiphi_dth2 * hor_th_y * hor_th_y + hor_dgphiphi_dR * hor_R_yy + hor_dgphiphi_dth * hor_th_yy;
        CCTK_REAL hor_d2gphiphi_dz2 = hor_d2gphiphi_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2gphiphi_dRth * hor_R_z * hor_th_z + hor_d2gphiphi_dth2 * hor_th_z * hor_th_z + hor_dgphiphi_dR * hor_R_zz + hor_dgphiphi_dth * hor_th_zz;
        CCTK_REAL hor_d2gphiphi_dxy = hor_d2gphiphi_dR2 * hor_R_x * hor_R_y + hor_d2gphiphi_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2gphiphi_dth2 * hor_th_x * hor_th_y + hor_dgphiphi_dR * hor_R_xy + hor_dgphiphi_dth * hor_th_xy;
        CCTK_REAL hor_d2gphiphi_dxz = hor_d2gphiphi_dR2 * hor_R_x * hor_R_z + hor_d2gphiphi_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2gphiphi_dth2 * hor_th_x * hor_th_z + hor_dgphiphi_dR * hor_R_xz + hor_dgphiphi_dth * hor_th_xz;
        CCTK_REAL hor_d2gphiphi_dyz = hor_d2gphiphi_dR2 * hor_R_y * hor_R_z + hor_d2gphiphi_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2gphiphi_dth2 * hor_th_y * hor_th_z + hor_dgphiphi_dR * hor_R_yz + hor_dgphiphi_dth * hor_th_yz;

        /* hor_bphhor_iup */
        CCTK_REAL hor_d2bphiup_dx2 = hor_d2bphiup_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2bphiup_dRth * hor_R_x * hor_th_x + hor_d2bphiup_dth2 * hor_th_x * hor_th_x + hor_dbphiup_dR * hor_R_xx + hor_dbphiup_dth * hor_th_xx;
        CCTK_REAL hor_d2bphiup_dy2 = hor_d2bphiup_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2bphiup_dRth * hor_R_y * hor_th_y + hor_d2bphiup_dth2 * hor_th_y * hor_th_y + hor_dbphiup_dR * hor_R_yy + hor_dbphiup_dth * hor_th_yy;
        CCTK_REAL hor_d2bphiup_dz2 = hor_d2bphiup_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2bphiup_dRth * hor_R_z * hor_th_z + hor_d2bphiup_dth2 * hor_th_z * hor_th_z + hor_dbphiup_dR * hor_R_zz + hor_dbphiup_dth * hor_th_zz;
        CCTK_REAL hor_d2bphiup_dxy = hor_d2bphiup_dR2 * hor_R_x * hor_R_y + hor_d2bphiup_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2bphiup_dth2 * hor_th_x * hor_th_y + hor_dbphiup_dR * hor_R_xy + hor_dbphiup_dth * hor_th_xy;
        CCTK_REAL hor_d2bphiup_dxz = hor_d2bphiup_dR2 * hor_R_x * hor_R_z + hor_d2bphiup_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2bphiup_dth2 * hor_th_x * hor_th_z + hor_dbphiup_dR * hor_R_xz + hor_dbphiup_dth * hor_th_xz;
        CCTK_REAL hor_d2bphiup_dyz = hor_d2bphiup_dR2 * hor_R_y * hor_R_z + hor_d2bphiup_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2bphiup_dth2 * hor_th_y * hor_th_z + hor_dbphiup_dR * hor_R_yz + hor_dbphiup_dth * hor_th_yz;

        /* hor_bphhor_i */
        CCTK_REAL hor_d2bphi_dx2 = hor_d2bphi_dR2 * hor_R_x * hor_R_x + 2.0 * hor_d2bphi_dRth * hor_R_x * hor_th_x + hor_d2bphi_dth2 * hor_th_x * hor_th_x + hor_dbphi_dR * hor_R_xx + hor_dbphi_dth * hor_th_xx;
        CCTK_REAL hor_d2bphi_dy2 = hor_d2bphi_dR2 * hor_R_y * hor_R_y + 2.0 * hor_d2bphi_dRth * hor_R_y * hor_th_y + hor_d2bphi_dth2 * hor_th_y * hor_th_y + hor_dbphi_dR * hor_R_yy + hor_dbphi_dth * hor_th_yy;
        CCTK_REAL hor_d2bphi_dz2 = hor_d2bphi_dR2 * hor_R_z * hor_R_z + 2.0 * hor_d2bphi_dRth * hor_R_z * hor_th_z + hor_d2bphi_dth2 * hor_th_z * hor_th_z + hor_dbphi_dR * hor_R_zz + hor_dbphi_dth * hor_th_zz;
        CCTK_REAL hor_d2bphi_dxy = hor_d2bphi_dR2 * hor_R_x * hor_R_y + hor_d2bphi_dRth * (hor_R_x * hor_th_y + hor_R_y * hor_th_x) + hor_d2bphi_dth2 * hor_th_x * hor_th_y + hor_dbphi_dR * hor_R_xy + hor_dbphi_dth * hor_th_xy;
        CCTK_REAL hor_d2bphi_dxz = hor_d2bphi_dR2 * hor_R_x * hor_R_z + hor_d2bphi_dRth * (hor_R_x * hor_th_z + hor_R_z * hor_th_x) + hor_d2bphi_dth2 * hor_th_x * hor_th_z + hor_dbphi_dR * hor_R_xz + hor_dbphi_dth * hor_th_xz;
        CCTK_REAL hor_d2bphi_dyz = hor_d2bphi_dR2 * hor_R_y * hor_R_z + hor_d2bphi_dRth * (hor_R_y * hor_th_z + hor_R_z * hor_th_y) + hor_d2bphi_dth2 * hor_th_y * hor_th_z + hor_dbphi_dR * hor_R_yz + hor_dbphi_dth * hor_th_yz;

        /* horizon shift derivative*/
        CCTK_REAL hor_beta[4];
        hor_beta[0] = 0;
        hor_beta[1] = -hor_y1_2*hor_bphi/hor_rho2_2;
        hor_beta[2] = hor_bphi * hor_x1_2 * gamma/hor_rho2_2;
        hor_beta[3] = 0.0; // because b^z=0
                                                          //estas partes estão mal!!!!!!!!! tem de ser o shift da metrica boosted.

        CCTK_REAL hor_betaup[4];
        hor_betaup[0] = 0;
        hor_betaup[1] = -hor_y1_2*hor_bphiup;
        hor_betaup[2] = hor_x1_2*gamma*hor_bphiup;
        hor_betaup[3] = 0.0; // because b^z=0


        CCTK_REAL hor_dbetaup[4][4];
        hor_dbetaup[1][0] = gamma * bh_v * (-hor_y1_2) * (hor_dbphiup_dR * hor_R_x + hor_dbphiup_dth * hor_th_x);
        hor_dbetaup[2][0] = gamma * bh_v * (hor_bphiup + hor_x1_2 * gamma * (hor_dbphiup_dR * hor_R_x + hor_dbphiup_dth * hor_th_x));
        hor_dbetaup[3][0] = 0.0; // because b^z=0
        hor_dbetaup[1][1] = (-hor_y1_2) * (hor_dbphiup_dR * hor_R_x + hor_dbphiup_dth * hor_th_x);
        hor_dbetaup[2][1] = (hor_bphiup * hor_x1_2 * gamma * (hor_dbphiup_dR * hor_R_x + hor_dbphiup_dth * hor_th_x));
        hor_dbetaup[3][1] = 0.0; // because b^z=0
        hor_dbetaup[1][2] = (-hor_bphiup - hor_y1_2 * (hor_dbphiup_dR * hor_R_y + hor_dbphiup_dth * hor_th_y));
        hor_dbetaup[2][2] = hor_x1_2 * gamma * (hor_dbphiup_dR * hor_R_y + hor_dbphiup_dth * hor_th_y);
        hor_dbetaup[3][2] = 0.0; // because b^z=0
        hor_dbetaup[1][3] = (-hor_y1_2) * (hor_dbphiup_dR * hor_R_z + hor_dbphiup_dth * hor_th_z);
        hor_dbetaup[2][3] = hor_x1_2 * gamma * (hor_dbphiup_dR * hor_R_z + hor_dbphiup_dth * hor_th_z);
        hor_dbetaup[3][3] = 0.0; // because b^z=0

        CCTK_REAL hor_dbetaup_dR[4];
        hor_dbetaup_dR[0] = 0.0;
        hor_dbetaup_dR[1] = hor_dbetaup[1][1] * hor_x_R + hor_dbetaup[1][2] * hor_y_R + hor_dbetaup[1][3] * hor_z_R;
        hor_dbetaup_dR[2] = hor_dbetaup[2][1] * hor_x_R + hor_dbetaup[2][2] * hor_y_R + hor_dbetaup[2][3] * hor_z_R;
        hor_dbetaup_dR[3] = hor_dbetaup[3][1] * hor_x_R + hor_dbetaup[3][2] * hor_y_R + hor_dbetaup[3][3] * hor_z_R;

        UAv_EvalPoint hor_P = {
            .bh_mass = bh_mass,
            .bh_spin = bh_spin,
            .bh_spin2 = bh_spin2,
            .bh_v = bh_v,
            .rBLp = rBLp,
            .rBLm = rBLm,
            .x1_2 = gamma * rhor * cosph_2 * sinth,
            .y1_2 = rhor * sinph_2 * sinth,
            .z1_2 = rhor * costh,
            .rr_2 = rhor,
            .rr2_2 = rhor * rhor,
            .rho_2 = rhor * sinth,
            .rho2_2 = rhor * rhor * sinth * sinth,

            .R_x = hor_R_x,
            .R_y = hor_R_y,
            .R_z = hor_R_z,

            .th_x = hor_th_x,
            .th_y = hor_th_y,
            .th_z = hor_th_z,

            .costh = hor_costh,
            .costh2 = hor_costh2,
            .sinth2 = hor_sinth2,
            .sinth = hor_sinth,
            // auxilliary quantities

            .rBL = hor_rBL,
            .Delt = hor_Delt,
            .Sigm = hor_Sigm,
            .Sigm2 = hor_Sigm2,
            .fctFF = hor_fctFF,

            .psi4_2 = hor_psi4_2, // psi04 no codigo original

            .fctGG = hor_fctGG,
            .fctHH = hor_fctHH,

            .alpha0 = hor_alpha0,
            .alpha02 = alpha02,

            .gphiphi = hor_gphiphi,
            .bphiup = hor_bphiup,
            .bphi = hor_bphi,

            // 1st derivatives of auxiliary quantities

            .drBLdR = hor_drBLdR,
            .dDelt_dR = hor_dDelt_dR,
            .dSigm_dR = hor_dSigm_dR,
            .dSigm_dth = hor_dSigm_dth,
            .dfctFF_dR = hor_dfctFF_dR,
            .dfctFF_dth = hor_dfctFF_dth,
            .dfctGG_dR = hor_dfctGG_dR,
            .dfctHH_dR = hor_dfctHH_dR,
            .dfctHH_dth = hor_dfctHH_dth,
            .dpsi4_2_dR = hor_dpsi4_2_dR,
            .dpsi4_2_dth = hor_dpsi4_2_dth,
            .dalpha02_dR = hor_dalpha02_dR,
            .dalpha02_dth = hor_dalpha02_dth,
            .dgphiphi_dR = hor_dgphiphi_dR,
            .dgphiphi_dth = hor_dgphiphi_dth,
            .dbphiup_dR = hor_dbphiup_dR,
            .dbphiup_dth = hor_dbphiup_dth,
            .dbphi_dR = hor_dbphi_dR,
            .dbphi_dth = hor_dbphi_dth,

            // 2nd derivatives

            /* Coordinate second derivatives (fill if needed) */
            .R_xx = hor_R_xx,
            .R_xy = hor_R_xy,
            .R_xz = hor_R_xz,
            .R_yy = hor_R_yy,
            .R_yz = hor_R_yz,
            .R_zz = hor_R_zz,
            .th_xx = hor_th_xx,
            .th_xy = hor_th_xy,
            .th_xz = hor_th_xz,
            .th_yy = hor_th_yy,
            .th_yz = hor_th_yz,
            .th_zz = hor_th_zz,

            /* Second derivatives in (R,theta) space (to be provided) */
            .d2rBL_dR2 = hor_d2rBL_dR2,
            .d2rBL_dRth = hor_d2rBL_dRth,
            .d2rBL_dth2 = hor_d2rBL_dth2, /* note: duplicated symbol in original */
            .d2Delt_dR2 = hor_d2Delt_dR2,
            .d2Delt_dRth = hor_d2Delt_dRth,
            .d2Delt_dth2 = hor_d2Delt_dth2,
            .d2Sigm_dR2 = hor_d2Sigm_dR2,
            .d2Sigm_dRth = hor_d2Sigm_dRth,
            .d2Sigm_dth2 = hor_d2Sigm_dth2,
            .d2fctFF_dR2 = hor_d2fctFF_dR2,
            .d2fctFF_dRth = hor_d2fctFF_dRth,
            .d2fctFF_dth2 = hor_d2fctFF_dth2,
            .d2fctGG_dR2 = hor_d2fctGG_dR2,
            .d2fctGG_dRth = hor_d2fctGG_dRth,
            .d2fctGG_dth2 = hor_d2fctGG_dth2,
            .d2fctHH_dR2 = hor_d2fctHH_dR2,
            .d2fctHH_dRth = hor_d2fctHH_dRth,
            .d2fctHH_dth2 = hor_d2fctHH_dth2,
            .d2psi4_2_dR2 = hor_d2psi4_2_dR2,
            .d2psi4_2_dRth = hor_d2psi4_2_dRth,
            .d2psi4_2_dth2 = hor_d2psi4_2_dth2,
            .d2alpha02_dR2 = hor_d2alpha02_dR2,
            .d2alpha02_dRth = hor_d2alpha02_dRth,
            .d2alpha02_dth2 = hor_d2alpha02_dth2,
            .d2gphiphi_dR2 = hor_d2gphiphi_dR2,
            .d2gphiphi_dRth = hor_d2gphiphi_dRth,
            .d2gphiphi_dth2 = hor_d2gphiphi_dth2,
            .d2bphiup_dR2 = hor_d2bphiup_dR2,
            .d2bphiup_dRth = hor_d2bphiup_dRth,
            .d2bphiup_dth2 = hor_d2bphiup_dth2,
            .d2bphi_dR2 = hor_d2bphi_dR2,
            .d2bphi_dRth = hor_d2bphi_dRth,
            .d2bphi_dth2 = hor_d2bphi_dth2,

            /* Cartesian second derivatives (fill using chain rule) */
            /* rBL (depends only on R) */
            .d2rBL_dx2 = hor_d2rBL_dx2,
            .d2rBL_dy2 = hor_d2rBL_dy2,
            .d2rBL_dz2 = hor_d2rBL_dz2,
            .d2rBL_dxy = hor_d2rBL_dxy,
            .d2rBL_dxz = hor_d2rBL_dxz,
            .d2rBL_dyz = hor_d2rBL_dyz,

            /* Delt */
            .dDelt_dth = hor_dDelt_dth, // Delt depende apenas de R
            .d2Delt_dx2 = hor_d2Delt_dx2,
            .d2Delt_dy2 = hor_d2Delt_dy2,
            .d2Delt_dz2 = hor_d2Delt_dz2,
            .d2Delt_dxy = hor_d2Delt_dxy,
            .d2Delt_dxz = hor_d2Delt_dxz,
            .d2Delt_dyz = hor_d2Delt_dyz,

            /* Sigm */
            .d2Sigm_dx2 = hor_d2Sigm_dx2,
            .d2Sigm_dy2 = hor_d2Sigm_dy2,
            .d2Sigm_dz2 = hor_d2Sigm_dz2,
            .d2Sigm_dxy = hor_d2Sigm_dxy,
            .d2Sigm_dxz = hor_d2Sigm_dxz,
            .d2Sigm_dyz = hor_d2Sigm_dyz,

            /* fctFF */
            .d2fctFF_dx2 = hor_d2fctFF_dx2,
            .d2fctFF_dy2 = hor_d2fctFF_dy2,
            .d2fctFF_dz2 = hor_d2fctFF_dz2,
            .d2fctFF_dxy = hor_d2fctFF_dxy,
            .d2fctFF_dxz = hor_d2fctFF_dxz,
            .d2fctFF_dyz = hor_d2fctFF_dyz,

            /* fctGG (only R) */
            .d2fctGG_dx2 = hor_d2fctGG_dx2,
            .d2fctGG_dy2 = hor_d2fctGG_dy2,
            .d2fctGG_dz2 = hor_d2fctGG_dz2,
            .d2fctGG_dxy = hor_d2fctGG_dxy,
            .d2fctGG_dxz = hor_d2fctGG_dxz,
            .d2fctGG_dyz = hor_d2fctGG_dyz,

            /* fctHH */
            .d2fctHH_dx2 = hor_d2fctHH_dx2,
            .d2fctHH_dy2 = hor_d2fctHH_dy2,
            .d2fctHH_dz2 = hor_d2fctHH_dz2,
            .d2fctHH_dxy = hor_d2fctHH_dxy,
            .d2fctHH_dxz = hor_d2fctHH_dxz,
            .d2fctHH_dyz = hor_d2fctHH_dyz,

            /* psi4_2 */
            .d2psi4_2_dx2 = hor_d2psi4_2_dx2,
            .d2psi4_2_dy2 = hor_d2psi4_2_dy2,
            .d2psi4_2_dz2 = hor_d2psi4_2_dz2,
            .d2psi4_2_dxy = hor_d2psi4_2_dxy,
            .d2psi4_2_dxz = hor_d2psi4_2_dxz,
            .d2psi4_2_dyz = hor_d2psi4_2_dyz,

            /* alpha02 */
            .d2alpha02_dx2 = hor_d2alpha02_dx2,
            .d2alpha02_dy2 = hor_d2alpha02_dy2,
            .d2alpha02_dz2 = hor_d2alpha02_dz2,
            .d2alpha02_dxy = hor_d2alpha02_dxy,
            .d2alpha02_dxz = hor_d2alpha02_dxz,
            .d2alpha02_dyz = hor_d2alpha02_dyz,

            /* gphiphi */
            .d2gphiphi_dx2 = hor_d2gphiphi_dx2,
            .d2gphiphi_dy2 = hor_d2gphiphi_dy2,
            .d2gphiphi_dz2 = hor_d2gphiphi_dz2,
            .d2gphiphi_dxy = hor_d2gphiphi_dxy,
            .d2gphiphi_dxz = hor_d2gphiphi_dxz,
            .d2gphiphi_dyz = hor_d2gphiphi_dyz,

            /* bphiup */
            .d2bphiup_dx2 = hor_d2bphiup_dx2,
            .d2bphiup_dy2 = hor_d2bphiup_dy2,
            .d2bphiup_dz2 = hor_d2bphiup_dz2,
            .d2bphiup_dxy = hor_d2bphiup_dxy,
            .d2bphiup_dxz = hor_d2bphiup_dxz,
            .d2bphiup_dyz = hor_d2bphiup_dyz,

            /* bphi */
            .d2bphi_dx2 = hor_d2bphi_dx2,
            .d2bphi_dy2 = hor_d2bphi_dy2,
            .d2bphi_dz2 = hor_d2bphi_dz2,
            .d2bphi_dxy = hor_d2bphi_dxy,
            .d2bphi_dxz = hor_d2bphi_dxz,
            .d2bphi_dyz = hor_d2bphi_dyz,
        };

        UAv_MetricDerivs1 hor_D1;
        UAv_MetricDerivs2 hor_D2;
        UAv_ComputeMetricDerivsAtPoint(&hor_P, &hor_D1, &hor_D2);

        // make if statement at the horizon for betaup and dg.

        CCTK_REAL hor_dg[4][4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              CCTK_REAL sum = 0.0;
              for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                  for (int lam = 0; lam < 4; ++lam)
                    sum += invLambda[mu][a] * invLambda[nu][b] * invLambda[lam][c] * hor_D1.dG[mu][nu][lam];
              hor_dg[a][b][c] = sum;
            }
          }
        }
        
        
        
        CCTK_REAL hor_ddg[4][4][4][4];

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                CCTK_REAL sum = 0.0;
                for (int mu = 0; mu < 4; ++mu)
                  for (int nu = 0; nu < 4; ++nu)
                    for (int lam = 0; lam < 4; ++lam)
                      for (int sig = 0; sig < 4; ++sig)
                        sum += invLambda[mu][a] * invLambda[nu][b] * invLambda[lam][c] * invLambda[sig][d] * hor_D2.ddG[mu][nu][lam][sig];
                hor_ddg[a][b][c][d] = sum;
              }
            }
          }
        }

        CCTK_REAL K_B[4][4]; // extrinsic curvature
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            K_B[a][b] = 0.0;
          }
        } // K_0\mu might not be zero but irrelevant for what i want to compute

        CCTK_REAL hor_Gb00 = gamma2* (-hor_alpha02+hor_bphi*hor_bphiup) + gamma2*bh_v2*(hor_psi4_2 * (1.0 + hor_x1_2 * hor_x1_2 * gamma2 * hor_fctGG + bh_spin2 * hor_y1_2 * hor_y1_2 * hor_fctHH)) + 2.0*gamma2*bh_v*(-hor_y1_2 / hor_rho2_2 * hor_bphi);
        CCTK_REAL hor_new_lapse = -hor_Gb00 + hor_betaup[1] * hor_beta[1] + hor_betaup[2] * hor_beta[2] + hor_betaup[3] * hor_beta[3];
        // CCTK_REAL hor_dnew_lapse_dR = hor_ddg[0][0][1] * hor_x_R + hor_ddg[0][0][2] * hor_y_R + hor_ddg[0][0][3] * hor_z_R + hor_dbetaup_dR[1] * hor_beta[1] + hor_dbetaup_dR[2] * hor_beta[2] + hor_dbetaup_dR[3] * hor_beta[3];

        if (fabs(rr_2-rBLp/4) > 1e-8) {
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
        } else if (fabs(rr_2-rBLp/4) <= 1e-8) {
          // inside the BH horizon
          for (int a = 1; a < 4; ++a) {
            for (int b = 1; b < 4; ++b) {
              CCTK_REAL sum1 = 0.0;
              CCTK_REAL sum2 = 0.0;
              CCTK_REAL sum3 = 0.0;
              for (int c = 1; c < 4; ++c) {
                sum1 += hor_dbetaup_dR[c] * hor_dg[a][b][c] + hor_betaup[c] * (hor_ddg[a][b][c][1] * hor_x_R + hor_ddg[a][b][c][2] * hor_y_R + hor_ddg[a][b][c][3] * hor_z_R);
                sum2 += hor_dbetaup_dR[c] * hor_dg[b][c][a] + hor_betaup[c] * (hor_ddg[b][c][a][1] * hor_x_R + hor_ddg[b][c][a][2] * hor_y_R + hor_ddg[b][c][a][3] * hor_z_R);
                sum3 += hor_dbetaup_dR[c] * hor_dg[a][c][b] + hor_betaup[c] * (hor_ddg[a][c][b][1] * hor_x_R + hor_ddg[a][c][b][2] * hor_y_R + hor_ddg[a][c][b][3] * hor_z_R);
              }
              K_B[a][b] = -0.5 * ((hor_ddg[a][b][0][1] * hor_x_R + hor_ddg[a][b][0][2] * hor_y_R + hor_ddg[a][b][0][3]*hor_z_R) - sum1 - ((hor_ddg[0][b][a][1] * hor_x_R + hor_ddg[0][b][a][2] * hor_y_R + hor_ddg[0][b][a][3] * hor_z_R) - sum2) - ((hor_ddg[0][a][b][1] * hor_x_R + hor_ddg[0][a][b][2] * hor_y_R + hor_ddg[0][a][b][3] * hor_z_R) - sum3)) / (hor_dalpha02_dR); //por o alpha boosted
            }
          }
        }

        kxx[ind] = K_B[1][1];
        kxy[ind] = K_B[1][2];
        kxz[ind] = K_B[1][3];
        kyy[ind] = K_B[2][2];
        kyz[ind] = K_B[2][3];
        kzz[ind] = K_B[3][3];

        // if (bh_mass == 0.0) {
        //   CCTK_REAL maxK = 0.0;
        //   for (int i3 = 1; i3 <= 3; ++i3)
        //     for (int j3 = 1; j3 <= 3; ++j3)
        //       maxK = fmax(maxK, fabs(K_B[i3][j3]));
        //   if (maxK > 1e-12) {
        //     CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
        //                "Boosted Minkowski test: |K|=%g", maxK);
        //     break;
        //   }
        // }

        // CCTK_REAL auxKij, facKij, facKijRho, facKijZ;
        // auxKij = 2.0 * rBL2 * (rBL2 + bh_spin2) + Sigm * (rBL2 - bh_spin2);
        // facKij = alpha0 * bh_spin * bh_mass * sinth2 / (rr2_2 * pow(rho_2, 3)
        // * Delt * Sigm2); facKijRho = 2.0 * z1_2 * bh_spin2 * rBL * Delt *
        // costh * sinth - rho_2 * rr_2 * drBLdR * auxKij; facKijZ = 2.0 * rho_2
        // * bh_spin2 * rBL * Delt * costh * sinth + z1_2 * rr_2 * drBLdR *
        // auxKij;

        // kxx[ind] = 2.0 * x1_2*gamma  * y1_2 * facKij * facKijRho;
        // kxy[ind] = (y1_2 * y1_2 - x1_2 * x1_2*gamma2) * facKij * facKijRho;
        // kxz[ind] = -y1_2 * rho_2 * facKij * facKijZ;
        // kyy[ind] = -2.0 * x1_2 *gamma * y1_2 * facKij * facKijRho;
        // kyz[ind] = x1_2 * gamma * rho_2 * facKij * facKijZ;
        // kzz[ind] = 0.0;

        // }

        // check_nan_or_inf("kxx", kxx[ind]);
        // check_nan_or_inf("kxy", kxy[ind]);
        // check_nan_or_inf("kxz", kxz[ind]);
        // check_nan_or_inf("kyy", kyy[ind]);
        // check_nan_or_inf("kyz", kyz[ind]);
        // check_nan_or_inf("kzz", kzz[ind]);

        const CCTK_REAL phi0_l_1 = phi0_1[ind]; // * pert_phi_1;
        const CCTK_REAL phi0_l_2 = 0.;

        // scalar fields

        // star 1
        CCTK_REAL phi1_1 = phi0_l_1 * (coswt * cosmph_1 + sinwt * sinmph_1);
        CCTK_REAL phi2_1 = phi0_l_1 * (coswt * sinmph_1 - sinwt * cosmph_1);

        // BH 2
        CCTK_REAL phi1_2 =
            0; // phi0_l_2 * (coswt * cosmph_2 + sinwt * sinmph_2);
        CCTK_REAL phi2_2 =
            0; // phi0_l_2 * (coswt * sinmph_2 - sinwt * cosmph_2);
        /////////////////////////////

        phi1[ind] = phi1_1 + phi1_2;
        phi2[ind] = phi2_1 + phi2_2;

        const CCTK_REAL alph_1 = exp(F0_1[ind]);
        const CCTK_REAL alph_2 = new_alpha;

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
        // CCTK_REAL fctFF = ( rBL*rBL + bh_spin2 ) * ( rBL*rBL + bh_spin2 ) -
        // Delt * bh_spin2 * sinth2; bphi = 2.0 * bh_spin * bh_mass * rBL /
        // fctFF;

        // shift
        if (CCTK_EQUALS(initial_shift, "Kerr_BS")) {
          betax[ind] = betaup[1];
          betay[ind] = betaup[2];
          betaz[ind] = betaup[3];
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
