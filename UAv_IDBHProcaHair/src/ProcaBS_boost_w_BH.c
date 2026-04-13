
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

void UAv_IDBHProcaHair_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL[], CCTK_REAL[],
                                 CCTK_REAL[], CCTK_REAL[], CCTK_REAL[], CCTK_REAL[],
                                 CCTK_REAL[], CCTK_REAL[], CCTK_REAL[], CCTK_REAL[]);

void check_nan_or_inf(const char *var_name, double value);

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

void UAv_IDProcaBSboostBH(CCTK_ARGUMENTS) {
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

  CCTK_INT NF;     // NF will be the actual size of the arrays
  CCTK_INT NX;     // NX will be the number of X points
  CCTK_INT Ntheta; // Ntheta will be the number of theta points

  CCTK_REAL *Xtmp, *thtmp, *F1_in, *F2_in, *F0_in, *Wbar_in, *H1_in, *H2_in, *H3_in, *V_in;
  Xtmp = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  thtmp = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  F1_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  F2_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  F0_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  Wbar_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  H1_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  H2_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  H3_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));
  V_in = (CCTK_REAL *)malloc(maxNF * sizeof(CCTK_REAL));

  // we get the data from the input file
  UAv_IDBHProcaHair_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, Wbar_in, H1_in, H2_in, H3_in, V_in);

  Ntheta = NF / NX;

  CCTK_VInfo(CCTK_THORNSTRING, "NX     = %d", NX);
  CCTK_VInfo(CCTK_THORNSTRING, "Ntheta = %d", Ntheta);
  CCTK_VInfo(CCTK_THORNSTRING, "NF     = %d", NF);

  // Minimum number of theta points required by the FD stencils, and the j index assignment below
  const CCTK_INT min_theta_pts = 5;
  if (Ntheta < min_theta_pts) {
    CCTK_VERROR("The initial data file doesn't have enough points in theta to be consistent with the implementation.\n "
                "Ntheta = %d. min_theta_points = %d.",
                Ntheta, min_theta_pts);
  }

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

  // To take care properly of z=0 symmetry (i.e. theta <-> pi-theta) in interpolation
  // we need to extend the arrays to z<0 values.
  // For convenience, we keep Ntheta as the number of points in the input half-space.

  NF = NX * (2 * Ntheta - 1);

  CCTK_REAL *F1_extd, *F2_extd, *F0_extd, *H2_extd, *H3_extd, *V_extd;
  F1_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  F2_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  F0_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  H2_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  H3_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  V_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  // We'll use A_r (H1r) rather than the input H1_in
  // Notation here: A ~ H1r dr + ... = H1_in/r dr + ...
  CCTK_REAL *H1r_extd;
  H1r_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  // now we need to take the derivatives of the Wbar function
  // Then we convert to W and store the values

  CCTK_REAL *W_extd, *dW_dr_extd, *dW_dth_extd;
  W_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dW_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dW_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  // Same for H3 and V

  CCTK_REAL *dH3_dr_extd, *dH3_dth_extd;
  dH3_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dH3_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  CCTK_REAL *dV_dr_extd, *dV_dth_extd, *dH1_dr_extd, *dH1_dth_extd, *dH2_dr_extd, *dH2_dth_extd;
  dV_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dV_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dH1_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dH1_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dH2_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dH2_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  // New: same for F0, F1, F2
  CCTK_REAL *dF0_dr_extd, *dF1_dr_extd, *dF2_dr_extd, *dF0_dth_extd, *dF1_dth_extd, *dF2_dth_extd;
  dF0_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dF1_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dF2_dr_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dF0_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dF1_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));
  dF2_dth_extd = (CCTK_REAL *)malloc(NF * sizeof(CCTK_REAL));

  // // Some auxi file for debug
  // FILE* debugfile = fopen ("testdebug.txt", "w");
  // if (debugfile == NULL) {
  //   CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
  //   "Unable to open file %s\n", "testdebug.txt");
  // } else {
  //   CCTK_VInfo(CCTK_THORNSTRING, "Write test file %s", "testdebug.txt");
  // }

  const CCTK_REAL oodX = 1. / dX;
  const CCTK_REAL oodX12 = 1. / (12. * dX);
  const CCTK_REAL oodXsq12 = oodX * oodX12;
  const CCTK_REAL oodth12 = 1. / (12. * dtheta);

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
      const CCTK_INT j = jj;
      const CCTK_INT jm1 = abs(jj - 1);
      const CCTK_INT jm2 = abs(jj - 2);
      const CCTK_INT jm3 = abs(jj - 3);
      const CCTK_INT jm4 = abs(jj - 4);
      const CCTK_INT jp1 = Ntheta - 1 - abs(Ntheta - 1 - (jj + 1));
      const CCTK_INT jp2 = Ntheta - 1 - abs(Ntheta - 1 - (jj + 2));
      const CCTK_INT jp3 = Ntheta - 1 - abs(Ntheta - 1 - (jj + 3));
      const CCTK_INT jp4 = Ntheta - 1 - abs(Ntheta - 1 - (jj + 4));

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
      const CCTK_INT indjm3 = i + jm3 * NX;
      const CCTK_INT indjm4 = i + jm4 * NX;
      const CCTK_INT indjp1 = i + jp1 * NX;
      const CCTK_INT indjp2 = i + jp2 * NX;
      const CCTK_INT indjp3 = i + jp3 * NX;
      const CCTK_INT indjp4 = i + jp4 * NX;

      // Just copy input values of ansatz functions
      F1_extd[ind] = F1_in[ind];
      F2_extd[ind] = F2_in[ind];
      F0_extd[ind] = F0_in[ind];
      H2_extd[ind] = H2_in[ind];
      H3_extd[ind] = H3_in[ind];
      V_extd[ind] = V_in[ind];

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

      const CCTK_REAL F0_th = (-F0_in[indjp2] + 8 * F0_in[indjp1] - 8 * F0_in[indjm1] + F0_in[indjm2]) *
                              oodth12;
      const CCTK_REAL F1_th = (-F1_in[indjp2] + 8 * F1_in[indjp1] - 8 * F1_in[indjm1] + F1_in[indjm2]) *
                              oodth12;
      const CCTK_REAL F2_th = (-F2_in[indjp2] + 8 * F2_in[indjp1] - 8 * F2_in[indjm1] + F2_in[indjm2]) *
                              oodth12;

      // Symmetries of V on the axis and/or the equator can vary (theta = 0, pi/2 resp.).
      // In particular, it can occur that dV/dth != 0, which can't be captured by centered finite differences and theta symmetry.
      // Since different systems have different symmetries, we resort to non-symmetric stencils in any case
      CCTK_REAL V_th, H1_th, H2_th;
      if (jj == 0) {
        // 1st derivative with 4th order accuracy (forward stencils)
        V_th = (-25 * V_in[ind] + 48 * V_in[indjp1] - 36 * V_in[indjp2] + 16 * V_in[indjp3] - 3 * V_in[indjp4]) *
               oodth12;
        H1_th = (-25 * H1_in[ind] + 48 * H1_in[indjp1] - 36 * H1_in[indjp2] + 16 * H1_in[indjp3] - 3 * H1_in[indjp4]) *
                oodth12;
        H2_th = (-25 * H2_in[ind] + 48 * H2_in[indjp1] - 36 * H2_in[indjp2] + 16 * H2_in[indjp3] - 3 * H2_in[indjp4]) *
                oodth12;
      } else if (jj == 1) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (-3 * V_in[indjm1] - 10 * V_in[ind] + 18 * V_in[indjp1] - 6 * V_in[indjp2] + V_in[indjp3]) *
               oodth12;
        H1_th = (-3 * H1_in[indjm1] - 10 * H1_in[ind] + 18 * H1_in[indjp1] - 6 * H1_in[indjp2] + H1_in[indjp3]) *
                oodth12;
        H2_th = (-3 * H2_in[indjm1] - 10 * H2_in[ind] + 18 * H2_in[indjp1] - 6 * H2_in[indjp2] + H2_in[indjp3]) *
                oodth12;
      } else if (jj == Ntheta - 2) {
        // 1st derivative with 4th order accuracy (mixed stencils)
        V_th = (3 * V_in[indjp1] + 10 * V_in[ind] - 18 * V_in[indjm1] + 6 * V_in[indjm2] - V_in[indjm3]) *
               oodth12;
        H1_th = (3 * H1_in[indjp1] + 10 * H1_in[ind] - 18 * H1_in[indjm1] + 6 * H1_in[indjm2] - H1_in[indjm3]) *
                oodth12;
        H2_th = (3 * H2_in[indjp1] + 10 * H2_in[ind] - 18 * H2_in[indjm1] + 6 * H2_in[indjm2] - H2_in[indjm3]) *
                oodth12;
      } else if (jj == Ntheta - 1) {
        // 1st derivative with 4th order accuracy (backward stencils)
        V_th = (25 * V_in[ind] - 48 * V_in[indjm1] + 36 * V_in[indjm2] - 16 * V_in[indjm3] + 3 * V_in[indjm4]) *
               oodth12;
        H1_th = (25 * H1_in[ind] - 48 * H1_in[indjm1] + 36 * H1_in[indjm2] - 16 * H1_in[indjm3] + 3 * H1_in[indjm4]) *
                oodth12;
        H2_th = (25 * H2_in[ind] - 48 * H2_in[indjm1] + 36 * H2_in[indjm2] - 16 * H2_in[indjm3] + 3 * H2_in[indjm4]) *
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

      CCTK_REAL Wbar_X, H3_X, V_X, F0_X, F1_X, F2_X, H2_X; // radial derivatives
      CCTK_REAL Wbar_XX = 0.;                              // Used for r=0 (i==0), if Wbar_r_power == 2.
      CCTK_REAL H1_X = 0.;                                 // Used for r=0 (i==0), due to H1_in/r.

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
        Wbar_X = (-25 * Wbar_in[ind] + 48 * Wbar_in[indip1] - 36 * Wbar_in[indip2] + 16 * Wbar_in[indip3] - 3 * Wbar_in[indip4]) * oodX12;

        // 1st derivative with 4th order accuracy (forward stencils)
        H3_X = (-25 * H3_in[ind] + 48 * H3_in[indip1] - 36 * H3_in[indip2] + 16 * H3_in[indip3] - 3 * H3_in[indip4]) * oodX12;

        // 1st derivative with 4th order accuracy (forward stencils)
        V_X = (-25 * V_in[ind] + 48 * V_in[indip1] - 36 * V_in[indip2] + 16 * V_in[indip3] - 3 * V_in[indip4]) * oodX12;

        H2_X = (-25 * H2_in[ind] + 48 * H2_in[indip1] - 36 * H2_in[indip2] + 16 * H2_in[indip3] - 3 * H2_in[indip4]) * oodX12;

        // 1st derivative with 4th order accuracy (forward stencils)
        F0_X = (-25 * F0_in[ind] + 48 * F0_in[indip1] - 36 * F0_in[indip2] + 16 * F0_in[indip3] - 3 * F0_in[indip4]) * oodX12;
        // 1st derivative with 4th order accuracy (forward stencils)
        F1_X = (-25 * F1_in[ind] + 48 * F1_in[indip1] - 36 * F1_in[indip2] + 16 * F1_in[indip3] - 3 * F1_in[indip4]) * oodX12;
        // 1st derivative with 4th order accuracy (forward stencils)
        F2_X = (-25 * F2_in[ind] + 48 * F2_in[indip1] - 36 * F2_in[indip2] + 16 * F2_in[indip3] - 3 * F2_in[indip4]) * oodX12;

        // Special care at r=0

        // H1_X required for H1r computed from H1_in/r
        // 1st derivative with 4th order accuracy (forward stencils)
        H1_X = (-25 * H1_in[ind] + 48 * H1_in[indip1] - 36 * H1_in[indip2] + 16 * H1_in[indip3] - 3 * H1_in[indip4]) * oodX12;

        if (Wbar_r_power == 2) {
          // If Wbar = r^2 * W, to compute W(r=0), we need to compute Wbar_XX.
          // 2nd derivative with 4th order accuracy (forward stencils)
          Wbar_XX = (45 * Wbar_in[ind] - 154 * Wbar_in[indip1] + 214 * Wbar_in[indip2] - 156 * Wbar_in[indip3] + 61 * Wbar_in[indip4] - 10 * Wbar_in[indip5]) * oodXsq12;
        }

      } else if (i == 1) {
        // 1st derivative, 4th order accuracy
        Wbar_X = (-3 * Wbar_in[indim1] - 10 * Wbar_in[ind] + 18 * Wbar_in[indip1] - 6 * Wbar_in[indip2] + Wbar_in[indip3]) * oodX12;

        // 1st derivative, 4th order accuracy
        H3_X = (-3 * H3_in[indim1] - 10 * H3_in[ind] + 18 * H3_in[indip1] - 6 * H3_in[indip2] + H3_in[indip3]) * oodX12;

        // 1st derivative, 4th order accuracy
        V_X = (-3 * V_in[indim1] - 10 * V_in[ind] + 18 * V_in[indip1] - 6 * V_in[indip2] + V_in[indip3]) * oodX12;

        H2_X = (-3 * H2_in[indim1] - 10 * H2_in[ind] + 18 * H2_in[indip1] - 6 * H2_in[indip2] + H2_in[indip3]) * oodX12;

        H1_X = (-3 * H1_in[indim1] - 10 * H1_in[ind] + 18 * H1_in[indip1] - 6 * H1_in[indip2] + H1_in[indip3]) * oodX12;

        // 1st derivative, 4th order accuracy
        F0_X = (-3 * F0_in[indim1] - 10 * F0_in[ind] + 18 * F0_in[indip1] - 6 * F0_in[indip2] + F0_in[indip3]) * oodX12;
        // 1st derivative, 4th order accuracy
        F1_X = (-3 * F1_in[indim1] - 10 * F1_in[ind] + 18 * F1_in[indip1] - 6 * F1_in[indip2] + F1_in[indip3]) * oodX12;
        // 1st derivative, 4th order accuracy
        F2_X = (-3 * F2_in[indim1] - 10 * F2_in[ind] + 18 * F2_in[indip1] - 6 * F2_in[indip2] + F2_in[indip3]) * oodX12;

      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        Wbar_X = (Wbar_in[indim2] - 4 * Wbar_in[indim1] + 3 * Wbar_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        H3_X = (H3_in[indim2] - 4 * H3_in[indim1] + 3 * H3_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        V_X = (V_in[indim2] - 4 * V_in[indim1] + 3 * V_in[ind]) * 0.5 * oodX;

        H2_X = (H2_in[indim2] - 4 * H2_in[indim1] + 3 * H2_in[ind]) * 0.5 * oodX;

        H1_X = (H1_in[indim2] - 4 * H1_in[indim1] + 3 * H1_in[ind]) * 0.5 * oodX;

        // 1st derivative with 2nd order accuracy (backward stencils)
        F0_X = (F0_in[indim2] - 4 * F0_in[indim1] + 3 * F0_in[ind]) * 0.5 * oodX;
        // 1st derivative with 2nd order accuracy (backward stencils)
        F1_X = (F1_in[indim2] - 4 * F1_in[indim1] + 3 * F1_in[ind]) * 0.5 * oodX;
        // 1st derivative with 2nd order accuracy (backward stencils)
        F2_X = (F2_in[indim2] - 4 * F2_in[indim1] + 3 * F2_in[ind]) * 0.5 * oodX;

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
        Wbar_X = (-Wbar_in[indip2] + 8 * Wbar_in[indip1] - 8 * Wbar_in[indim1] + Wbar_in[indim2]) * oodX12;

        // 4th order accurate stencils
        H3_X = (-H3_in[indip2] + 8 * H3_in[indip1] - 8 * H3_in[indim1] + H3_in[indim2]) * oodX12;

        // 4th order accurate stencils
        V_X = (-V_in[indip2] + 8 * V_in[indip1] - 8 * V_in[indim1] + V_in[indim2]) * oodX12;

        H2_X = (-H2_in[indip2] + 8 * H2_in[indip1] - 8 * H2_in[indim1] + H2_in[indim2]) * oodX12;

        H1_X = (-H1_in[indip2] + 8 * H1_in[indip1] - 8 * H1_in[indim1] + H1_in[indim2]) * oodX12;

        // 4th order accurate stencils
        F0_X = (-F0_in[indip2] + 8 * F0_in[indip1] - 8 * F0_in[indim1] + F0_in[indim2]) * oodX12;
        // 4th order accurate stencils
        F1_X = (-F1_in[indip2] + 8 * F1_in[indip1] - 8 * F1_in[indim1] + F1_in[indim2]) * oodX12;
        // 4th order accurate stencils
        F2_X = (-F2_in[indip2] + 8 * F2_in[indip1] - 8 * F2_in[indim1] + F2_in[indim2]) * oodX12;
      }

      // From the X coordinate used in the input files to the r coordinate (coincides with x1 for the Boson Star, rH=0).
      // We also do the conversion from Wbar to W here, and H1_in to H1r, to tackle r = 0 (X = 0).

      // i == 0  <=>  X == 0  <=>  r == 0
      if (i == 0) {

        // At r=0 we have dW/dr = 0 and dW/dth = 0
        dW_dr_extd[ind] = 0.;
        dW_dth_extd[ind] = 0.;

        // At X==0 (rr==0), dXdr = 1/C0
        dH3_dr_extd[ind] = H3_X / C0;
        dV_dr_extd[ind] = V_X / C0;
        dH1_dr_extd[ind] = H1_X / C0;
        dH2_dr_extd[ind] = H2_X / C0;
        dF0_dr_extd[ind] = F0_X / C0;
        dF1_dr_extd[ind] = F1_X / C0;
        dF2_dr_extd[ind] = F2_X / C0;

        // For W we need more care depending on the power
        switch (Wbar_r_power) {
        case 0: // Wbar = W
          W_extd[ind] = Wbar_in[ind];
          break;

        case 1: // Wbar = r * W
          /*
          dWbar/dr = W + r * dW/dr
                   = W + 0 * 0     at r=0

          dWbar/dr = dWbar/dX * dX/dr
          dX/dr = C/(C+r)^2 = 1/C  at r=0
          */
          W_extd[ind] = Wbar_X / C0;
          break;

        case 2: // Wbar = r^2 * W
          /*
          dWbar/dr   = 2r * W + r^2 * dW/dr
          d2Wbar/dr2 = 2  * W + 4r  * dW/dr + r^2 * d2W/dr2
                     = 2  * W + 0 + 0        at r=0

          d2Wbar/dr2 = d2Wbar/dX2 * (dX/dr)^2 + dWbar/dX * d2X/dr2
          dX/dr   =   C/(C+r)^2 =  1/C      at r=0
          d2X/dr2 = -2C/(C+r)^3 = -2/C^2    at r=0

          W (r=0) = 1/C^2 * [1/2 * d2Wbar/dX^2 (X=0)  -  dWbar/dX (X=0)]
          */
          W_extd[ind] = (0.5 * Wbar_XX - Wbar_X) / (C0 * C0);
          break;

        default: // As of writing, this should be prevented by the scope of the parameter anyway
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
        W_extd[ind] = 0.;

        // Actually, the asymptotic expansion (Appendix B in the construction paper) also gives:
        dW_dr_extd[ind] = 0.;
        dW_dth_extd[ind] = 0.;

        dH3_dr_extd[ind] = 0.;
        dV_dr_extd[ind] = 0.;
        dH1_dr_extd[ind] = 0.;
        dH2_dr_extd[ind] = 0.;
        dF0_dr_extd[ind] = 0.;
        dF1_dr_extd[ind] = 0.;
        dF2_dr_extd[ind] = 0.;

        H1r_extd[ind] = 0.; // A_r = 0 at infinity

      } else {

        const CCTK_REAL rr = C0 * lX / (1. - lX);

        // corresponding derivatives
        // const CCTK_REAL dXdr = 1./(C0 + rr) - rr/((C0 + rr)*(C0 + rr));
        const CCTK_REAL dXdr = C0 / ((C0 + rr) * (C0 + rr));

        const CCTK_REAL Wbar_r = dXdr * Wbar_X;

        dH3_dr_extd[ind] = dXdr * H3_X;
        dV_dr_extd[ind] = dXdr * V_X;
        dH1_dr_extd[ind] = dXdr * H1_X;
        dH2_dr_extd[ind] = dXdr * H2_X;
        dF0_dr_extd[ind] = dXdr * F0_X;
        dF1_dr_extd[ind] = dXdr * F1_X;
        dF2_dr_extd[ind] = dXdr * F2_X;

        // Now translate from Wbar to W
        switch (Wbar_r_power) // We could put a generic power for the computation here I guess...
        {
        case 0: // Wbar = W
          W_extd[ind] = Wbar_in[ind];
          dW_dr_extd[ind] = Wbar_r;
          dW_dth_extd[ind] = Wbar_th;
          break;

        case 1: // Wbar = r * W
          W_extd[ind] = Wbar_in[ind] / rr;
          dW_dr_extd[ind] = (Wbar_r - W_extd[ind]) / rr; // dW/dr  =  1/r * dWbar/dr - Wbar / r^2  =  (dWbar/dr - W) / r
          dW_dth_extd[ind] = Wbar_th / rr;
          break;

        case 2:; // Wbar = r^2 * W
          // empty statement after case to prevent compilation error on some gcc versions...
          const CCTK_REAL rr2 = rr * rr;
          W_extd[ind] = Wbar_in[ind] / rr2;
          dW_dr_extd[ind] = Wbar_r / rr2 - 2 * W_extd[ind] / rr; // dW/dr  =  1/r^2 * dWbar/dr - 2 * Wbar / r^3  =  1/r^2 * dWbar/dr - 2 * W / r
          dW_dth_extd[ind] = Wbar_th / rr2;
          break;

        default: // As of writing, this should be prevented by the scope of the parameter anyway
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Unknown value of Wbar_r_power: %d. Aborting.", Wbar_r_power);
          break;
        }

        // From H1_in/r to H1r
        H1r_extd[ind] = H1_in[ind] / rr;

      } // if/else i==...

      dH3_dth_extd[ind] = H3_th;
      dV_dth_extd[ind] = V_th;
      dH1_dth_extd[ind] = H1_th;
      dH2_dth_extd[ind] = H2_th;
      dF0_dth_extd[ind] = F0_th;
      dF1_dth_extd[ind] = F1_th;
      dF2_dth_extd[ind] = F2_th;

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
  const CCTK_INT V_z_sign = Amu_z_sym_is_odd ? -1 : +1;

  for (int jj = 1; jj < Ntheta; jj++) { // don't repeat theta == pi/2
    for (int i = 0; i < NX; i++) {

      // j or jsym == Ntheta - 1  is theta == pi/2
      const CCTK_INT j = Ntheta - 1 + jj;
      const CCTK_INT jsym = Ntheta - 1 - jj;

      const CCTK_INT ind = i + j * NX;
      const CCTK_INT indsym = i + jsym * NX;

      // Even
      F1_extd[ind] = F1_extd[indsym];
      F2_extd[ind] = F2_extd[indsym];
      F0_extd[ind] = F0_extd[indsym];

      W_extd[ind] = W_extd[indsym];
      dW_dr_extd[ind] = dW_dr_extd[indsym];
      dF0_dr_extd[ind] = dF0_dr_extd[indsym];
      dF1_dr_extd[ind] = dF1_dr_extd[indsym];
      dF2_dr_extd[ind] = dF2_dr_extd[indsym];

      // Odd
      dW_dth_extd[ind] = -dW_dth_extd[indsym];
      dF0_dth_extd[ind] = -dF0_dth_extd[indsym];
      dF1_dth_extd[ind] = -dF1_dth_extd[indsym];
      dF2_dth_extd[ind] = -dF2_dth_extd[indsym];

      // Vector potential

      H1r_extd[ind] = H1_z_sign * H1r_extd[indsym];
      H2_extd[ind] = H2_z_sign * H2_extd[indsym];
      H3_extd[ind] = H3_z_sign * H3_extd[indsym];
      V_extd[ind] = V_z_sign * V_extd[indsym];

      dH3_dr_extd[ind] = H3_z_sign * dH3_dr_extd[indsym];
      dV_dr_extd[ind] = V_z_sign * dV_dr_extd[indsym];
      dH1_dr_extd[ind] = H1_z_sign * dH1_dr_extd[indsym];
      dH2_dr_extd[ind] = H2_z_sign * dH2_dr_extd[indsym];

      dH3_dth_extd[ind] = -H3_z_sign * dH3_dth_extd[indsym];
      dV_dth_extd[ind] = -V_z_sign * dV_dth_extd[indsym];
      dH1_dth_extd[ind] = -H1_z_sign * dH1_dth_extd[indsym];
      dH2_dth_extd[ind] = -H2_z_sign * dH2_dth_extd[indsym];

    } // for i
  } // for jj

  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]; // total points

  CCTK_REAL *X_g, *theta_g;
  X_g = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  theta_g = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));

  // const CCTK_REAL bs_v2 = bs_v * bs_v;
  const CCTK_REAL bs_v2 = bs_vx * bs_vx + bs_vy * bs_vy + bs_vz * bs_vz;
  const CCTK_REAL bs_v = sqrt(bs_v2);
  const CCTK_REAL gamma2 = 1. / (1. - bs_v2);
  const CCTK_REAL gamma = sqrt(gamma2);

  for (int k = 0; k < cctk_lsh[2]; ++k) { // code is in lab-frame coordinates. need to write functions as functions of rest frame coords
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_REAL x1 = x[ind] - x0;
        const CCTK_REAL y1 = y[ind] - y0;
        const CCTK_REAL z1 = z[ind] - z0;

        const CCTK_REAL x1rest = x1 + gamma2 / (gamma + 1.) * bs_vx * (bs_vx * x1 + bs_vy * y1 + bs_vz * z1);
        const CCTK_REAL y1rest = y1 + gamma2 / (gamma + 1.) * bs_vy * (bs_vx * x1 + bs_vy * y1 + bs_vz * z1);
        const CCTK_REAL z1rest = z1 + gamma2 / (gamma + 1.) * bs_vz * (bs_vx * x1 + bs_vy * y1 + bs_vz * z1);

        const CCTK_REAL hx = mu * x1rest;
        const CCTK_REAL hy = mu * y1rest;
        const CCTK_REAL hz = mu * z1rest;

        const CCTK_REAL rr2 = hx * hx + hy * hy + hz * hz;

        CCTK_REAL rr = sqrt(rr2);
        /* For the Boson Star, x, r and R coordinates coincide (rH=0). */
        /* note that there are divisions by rr in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        // From r to the X radial coordinate (used in input files)
        const CCTK_REAL lX = rr / (C0 + rr);

        const CCTK_REAL ltheta = rr < 1e-16 ? 0 : acos(hz / rr); // There should be at most one point in the grid with rr~0. Not sure about the threshold.

        X_g[ind] = lX;
        theta_g[ind] = ltheta;
      }
    }
  }

  /* now for the interpolation */

  const CCTK_INT N_dims = 2; // 2-D interpolation

  const CCTK_INT N_input_arrays = 24;
  const CCTK_INT N_output_arrays = 24;

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
  const void *interp_coords[N_dims];
  interp_coords[0] = (const void *)X_g;
  interp_coords[1] = (const void *)theta_g;

  /* input arrays */
  const void *input_arrays[N_input_arrays];
  CCTK_INT input_array_type_codes[N_input_arrays];
  CCTK_INT input_array_dims[N_dims];
  input_array_dims[0] = NX;
  input_array_dims[1] = 2 * Ntheta - 1;

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
  input_array_type_codes[10] = CCTK_VARIABLE_REAL;
  input_array_type_codes[11] = CCTK_VARIABLE_REAL;
  input_array_type_codes[12] = CCTK_VARIABLE_REAL;
  input_array_type_codes[13] = CCTK_VARIABLE_REAL;
  input_array_type_codes[14] = CCTK_VARIABLE_REAL;
  input_array_type_codes[15] = CCTK_VARIABLE_REAL;
  input_array_type_codes[16] = CCTK_VARIABLE_REAL;
  input_array_type_codes[17] = CCTK_VARIABLE_REAL;
  input_array_type_codes[18] = CCTK_VARIABLE_REAL;
  input_array_type_codes[19] = CCTK_VARIABLE_REAL;
  input_array_type_codes[20] = CCTK_VARIABLE_REAL;
  input_array_type_codes[21] = CCTK_VARIABLE_REAL;
  input_array_type_codes[22] = CCTK_VARIABLE_REAL;
  input_array_type_codes[23] = CCTK_VARIABLE_REAL;

  /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
  input_arrays[0] = (const void *)F1_extd;
  input_arrays[1] = (const void *)F2_extd;
  input_arrays[2] = (const void *)F0_extd;
  input_arrays[3] = (const void *)W_extd;
  input_arrays[4] = (const void *)dW_dr_extd;
  input_arrays[5] = (const void *)dW_dth_extd;
  input_arrays[6] = (const void *)H1r_extd;
  input_arrays[7] = (const void *)H2_extd;
  input_arrays[8] = (const void *)H3_extd;
  input_arrays[9] = (const void *)V_extd;
  input_arrays[10] = (const void *)dH3_dr_extd;
  input_arrays[11] = (const void *)dH3_dth_extd;
  input_arrays[12] = (const void *)dV_dr_extd;
  input_arrays[13] = (const void *)dV_dth_extd;
  input_arrays[14] = (const void *)dF0_dr_extd;
  input_arrays[15] = (const void *)dF1_dr_extd;
  input_arrays[16] = (const void *)dF2_dr_extd;
  input_arrays[17] = (const void *)dF0_dth_extd;
  input_arrays[18] = (const void *)dF1_dth_extd;
  input_arrays[19] = (const void *)dF2_dth_extd;
  input_arrays[20] = (const void *)dH1_dr_extd;
  input_arrays[21] = (const void *)dH1_dth_extd;
  input_arrays[22] = (const void *)dH2_dr_extd;
  input_arrays[23] = (const void *)dH2_dth_extd;

  /* output arrays */
  void *output_arrays[N_output_arrays];
  CCTK_INT output_array_type_codes[N_output_arrays];
  CCTK_REAL *F1, *F2, *F0, *W, *H1r, *H2, *H3, *V;
  CCTK_REAL *dW_dr, *dW_dth;
  CCTK_REAL *dH3_dr, *dH3_dth;
  CCTK_REAL *dV_dr, *dV_dth, *dH1_dr, *dH1_dth, *dH2_dr, *dH2_dth;
  CCTK_REAL *dF0_dr, *dF1_dr, *dF2_dr, *dF0_dth, *dF1_dth, *dF2_dth;

  F1 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  F2 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  F0 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  W = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  H1r = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  H2 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  H3 = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  V = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dH3_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dV_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dH1_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dH1_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dF0_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dF1_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dF2_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dF0_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dF1_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dF2_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dH2_dr = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));
  dH2_dth = (CCTK_REAL *)malloc(N_interp_points * sizeof(CCTK_REAL));

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
  output_array_type_codes[10] = CCTK_VARIABLE_REAL;
  output_array_type_codes[11] = CCTK_VARIABLE_REAL;
  output_array_type_codes[12] = CCTK_VARIABLE_REAL;
  output_array_type_codes[13] = CCTK_VARIABLE_REAL;
  output_array_type_codes[14] = CCTK_VARIABLE_REAL;
  output_array_type_codes[15] = CCTK_VARIABLE_REAL;
  output_array_type_codes[16] = CCTK_VARIABLE_REAL;
  output_array_type_codes[17] = CCTK_VARIABLE_REAL;
  output_array_type_codes[18] = CCTK_VARIABLE_REAL;
  output_array_type_codes[19] = CCTK_VARIABLE_REAL;
  output_array_type_codes[20] = CCTK_VARIABLE_REAL;
  output_array_type_codes[21] = CCTK_VARIABLE_REAL;
  output_array_type_codes[22] = CCTK_VARIABLE_REAL;
  output_array_type_codes[23] = CCTK_VARIABLE_REAL;

  output_arrays[0] = (void *)F1;
  output_arrays[1] = (void *)F2;
  output_arrays[2] = (void *)F0;
  output_arrays[3] = (void *)W;
  output_arrays[4] = (void *)dW_dr;
  output_arrays[5] = (void *)dW_dth;
  output_arrays[6] = (void *)H1r;
  output_arrays[7] = (void *)H2;
  output_arrays[8] = (void *)H3;
  output_arrays[9] = (void *)V;
  output_arrays[10] = (void *)dH3_dr;
  output_arrays[11] = (void *)dH3_dth;
  output_arrays[12] = (void *)dV_dr;
  output_arrays[13] = (void *)dV_dth;
  output_arrays[14] = (void *)dF0_dr;
  output_arrays[15] = (void *)dF1_dr;
  output_arrays[16] = (void *)dF2_dr;
  output_arrays[17] = (void *)dF0_dth;
  output_arrays[18] = (void *)dF1_dth;
  output_arrays[19] = (void *)dF2_dth;
  output_arrays[20] = (void *)dH1_dr;
  output_arrays[21] = (void *)dH1_dth;
  output_arrays[22] = (void *)dH2_dr;
  output_arrays[23] = (void *)dH2_dth;

  /* handle and settings for the interpolation routine */
  int operator_handle, param_table_handle;
  operator_handle = CCTK_InterpHandle("Lagrange polynomial interpolation");
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

  free(X_g);
  free(theta_g);
  free(Xtmp);
  free(thtmp);
  free(F1_in);
  free(F2_in);
  free(F0_in);
  free(Wbar_in);
  free(H1_in);
  free(H2_in);
  free(H3_in);
  free(V_in);

  free(F1_extd);
  free(F2_extd);
  free(F0_extd);
  free(W_extd);
  free(H1r_extd);
  free(H2_extd);
  free(H3_extd);
  free(V_extd);
  free(dW_dr_extd);
  free(dW_dth_extd);
  free(dH3_dr_extd);
  free(dH3_dth_extd);
  free(dV_dr_extd);
  free(dV_dth_extd);
  free(dH1_dr_extd);
  free(dH1_dth_extd);
  free(dF0_dr_extd);
  free(dF1_dr_extd);
  free(dF2_dr_extd);
  free(dF0_dth_extd);
  free(dF1_dth_extd);
  free(dF2_dth_extd);

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

  // const CCTK_REAL coswt = cos(omega_BS * tt * gamma);
  // const CCTK_REAL sinwt = sin(omega_BS * tt * gamma);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        // Proca Star A

        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_REAL x1 = x[ind] - x0;
        const CCTK_REAL y1 = y[ind] - y0;
        const CCTK_REAL z1 = z[ind] - z0;

        const CCTK_REAL x1rest = x1 + gamma2 / (gamma + 1.) * bs_vx * (bs_vx * x1 + bs_vy * y1 + bs_vz * z1);
        const CCTK_REAL y1rest = y1 + gamma2 / (gamma + 1.) * bs_vy * (bs_vx * x1 + bs_vy * y1 + bs_vz * z1);
        const CCTK_REAL z1rest = z1 + gamma2 / (gamma + 1.) * bs_vz * (bs_vx * x1 + bs_vy * y1 + bs_vz * z1);
        const CCTK_REAL ttrest = gamma * (tt + (bs_vx * x1 + bs_vy * y1 + bs_vz * z1));

        CCTK_REAL hx = x1rest * mu, hy = y1rest * mu, hz = z1rest * mu, ht = mu * ttrest;

        // For the Boson Star, r = R, no coordinate change needed.
        CCTK_REAL rr2 = hx * hx + hy * hy + hz * hz;
        // if( rr2 < pow( eps_r, 2 ) ) {
        // rr2 = pow( eps_r, 2 );
        // }
        const CCTK_REAL rr = sqrt(rr2);
        /* note that there are divisions by rr in the following expressions.
          divisions by zero should be avoided by choosing a non-zero value for
          z0 (for instance)
        */

        CCTK_REAL rho2 = hx * hx + hy * hy;
        // if( rho2 < pow( eps_r, 2 ) ){
        // rho2 = pow( eps_r, 2 );
        // }
        const CCTK_REAL rho = sqrt(rho2);

        const CCTK_REAL coswt = cos(omega_BS * ht); // cos( omega_BS * ttrest);
        const CCTK_REAL sinwt = sin(omega_BS * ht);

        const CCTK_REAL costh = hz / rr;
        const CCTK_REAL costh2 = costh * costh;
        /*
          For some grid points actually on the axis, it occurred that costh = 1-1e-16, resulting in sinth ~ 1.5e-8 instead of 0.
          Thus we force it in that case.
          Even if there is a legit grid point such that theta ~ a few 1e-8, it should mean RR >> rho and the axis treatment should be fine.
        */
        CCTK_REAL sinth, sinth2;
        if (1 - costh2 < 1e-15) {
          sinth2 = 0.;
          sinth = 0.;
        } else {
          sinth2 = 1. - costh2;
          sinth = sqrt(sinth2);
        }

        const CCTK_REAL ph = atan2(hy, hx);
        // If hx=hy=0, should return 0? The other metric functions should vanish anyway to make sure that this doesn't matter,
        // but can this lead to nan depending on the C implementation?

        const CCTK_REAL cosph = cos(ph);
        const CCTK_REAL sinph = sin(ph);

        const CCTK_REAL cosmph = cos(mm * ph);
        const CCTK_REAL sinmph = sin(mm * ph);

        const CCTK_REAL d_sinph_dx = -hx * hy / (rho2 * rho);
        const CCTK_REAL d_sinph_dy = (hx * hx) / (rho2 * rho);
        const CCTK_REAL d_sinph_dz = 0;

        const CCTK_REAL d_cosph_dx = (hy * hy) / (rho2 * rho);
        const CCTK_REAL d_cosph_dy = -(hx * hy) / (rho2 * rho);
        const CCTK_REAL d_cosph_dz = 0;

        const CCTK_REAL h_rho2 = exp(2. * (F2[ind] - F1[ind])) - 1.;

        const CCTK_REAL R_x = hx / rr;
        const CCTK_REAL R_y = hy / rr;
        const CCTK_REAL R_z = hz / rr;

        const CCTK_REAL th_x = costh * R_x / rho;
        const CCTK_REAL th_y = costh * R_y / rho;
        const CCTK_REAL th_z = -rho / rr2;

        const CCTK_REAL psi4 = exp(2. * F1[ind]);
        const CCTK_REAL psi2 = sqrt(psi4);
        const CCTK_REAL psi1 = sqrt(psi2);

        const CCTK_REAL dF1_dx = dF1_dr[ind] * R_x + dF1_dth[ind] * th_x;
        const CCTK_REAL dF1_dy = dF1_dr[ind] * R_y + dF1_dth[ind] * th_y;
        const CCTK_REAL dF1_dz = dF1_dr[ind] * R_z + dF1_dth[ind] * th_z;
        const CCTK_REAL dF2_dx = dF2_dr[ind] * R_x + dF2_dth[ind] * th_x;
        const CCTK_REAL dF2_dy = dF2_dr[ind] * R_y + dF2_dth[ind] * th_y;
        const CCTK_REAL dF2_dz = dF2_dr[ind] * R_z + dF2_dth[ind] * th_z;
        const CCTK_REAL dF0_dx = dF0_dr[ind] * R_x + dF0_dth[ind] * th_x;
        const CCTK_REAL dF0_dy = dF0_dr[ind] * R_y + dF0_dth[ind] * th_y;
        const CCTK_REAL dF0_dz = dF0_dr[ind] * R_z + dF0_dth[ind] * th_z;

        CCTK_REAL G[4][4]; // temporary storage for the 4-metric
        // CCTK_REAL G3_inv[4][4]; // temporary storage for the inverse of the 3-metric
        CCTK_REAL Gb[4][4];         // temporary storage for the boosted metric
        CCTK_REAL gammaA_inv[4][4]; // temporary storage for the inverse of the boosted 3-metric

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            G[a][b] = 0.0;
            gammaA_inv[a][b] = 0.0;
            // G3_inv[a][b] = 0.0;
            Gb[a][b] = 0.0;
          }
        }

        G[0][0] = -exp(2. * F0[ind]);
        G[1][1] = psi4 * (1. + h_rho2 * sinph * sinph);
        G[1][2] = -psi4 * h_rho2 * sinph * cosph;
        G[2][1] = G[1][2];
        G[2][2] = psi4 * (1. + h_rho2 * cosph * cosph);
        G[3][3] = psi4;

        // Derivatives of the metric functions
        CCTK_REAL dG[4][4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              dG[a][b][c] = 0.0;
            }
          }
        }
        // dG[a][b][c] = dG_ab/dx^c at rest

        dG[1][1][1] = (2 * exp(2. * F1[ind]) * hx * (pow(hy, 2) + hx * (pow(hx, 2) + pow(hy, 2)) * dF1_dx) + 2 * exp(2. * F2[ind]) * pow(hy, 2) * (-hx + (pow(hx, 2) + pow(hy, 2)) * dF2_dx)) / pow(pow(hx, 2) + pow(hy, 2), 2);

        dG[1][1][2] = (2 * exp(2. * F1[ind]) * pow(hx, 2) * (-hy + (pow(hx, 2) + pow(hy, 2)) * dF1_dy) + 2 * exp(2. * F2[ind]) * hy * (pow(hx, 2) + hy * (pow(hx, 2) + pow(hy, 2)) * dF2_dy)) / pow(pow(hx, 2) + pow(hy, 2), 2);

        dG[1][1][3] = (2 * (exp(2. * F1[ind]) * pow(hx, 2) * dF1_dz + exp(2. * F2[ind]) * pow(hy, 2) * dF2_dz)) / (pow(hx, 2) + pow(hy, 2));

        dG[1][2][1] = (exp(2. * F1[ind]) * hy * (-pow(hx, 2) + pow(hy, 2) + 2 * hx * (pow(hx, 2) + pow(hy, 2)) * dF1_dx) - exp(2. * F2[ind]) * hy * (-pow(hx, 2) + pow(hy, 2) + 2 * hx * (pow(hx, 2) + pow(hy, 2)) * dF2_dx)) / pow(pow(hx, 2) + pow(hy, 2), 2);

        dG[1][2][2] = (exp(2. * F1[ind]) * hx * (pow(hx, 2) - pow(hy, 2) + 2 * hy * (pow(hx, 2) + pow(hy, 2)) * dF1_dy) - exp(2. * F2[ind]) * hx * (pow(hx, 2) - pow(hy, 2) + 2 * hy * (pow(hx, 2) + pow(hy, 2)) * dF2_dy)) / pow(pow(hx, 2) + pow(hy, 2), 2);

        dG[1][2][3] = (2 * hx * hy * (exp(2. * F1[ind]) * dF1_dz - exp(2. * F2[ind]) * dF2_dz)) / rho2;

        dG[2][2][1] = (2 * exp(2. * F1[ind]) * pow(hy, 2) * (-hx + (pow(hx, 2) + pow(hy, 2)) * dF1_dx) + 2 * exp(2. * F2[ind]) * hx * (pow(hy, 2) + hx * (pow(hx, 2) + pow(hy, 2)) * dF2_dx)) / pow(pow(hx, 2) + pow(hy, 2), 2);

        dG[2][2][2] = (2 * exp(2. * F1[ind]) * hy * (pow(hx, 2) + hy * (pow(hx, 2) + pow(hy, 2)) * dF1_dy) + 2 * exp(2. * F2[ind]) * pow(hx, 2) * (-hy + (pow(hx, 2) + pow(hy, 2)) * dF2_dy)) / pow(pow(hx, 2) + pow(hy, 2), 2);

        dG[2][2][3] = (2 * (exp(2. * F1[ind]) * pow(hy, 2) * dF1_dz + exp(2. * F2[ind]) * pow(hx, 2) * dF2_dz)) / (pow(hx, 2) + pow(hy, 2));

        // dG23_dx^i = 0

        dG[3][3][1] = 2 * exp(2. * F1[ind]) * dF1_dx;
        dG[3][3][2] = 2 * exp(2. * F1[ind]) * dF1_dy;
        dG[3][3][3] = 2 * exp(2. * F1[ind]) * dF1_dz;

        dG[0][0][1] = -2 * exp(2. * F0[ind]) * dF0_dx;
        dG[0][0][2] = -2 * exp(2. * F0[ind]) * dF0_dy;
        dG[0][0][3] = -2 * exp(2. * F0[ind]) * dF0_dz;

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

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              if (isnan(dG[a][b][c]) || isinf(dG[a][b][c])) {
                fprintf(stderr, "Error: dG[%d][%d][%d] is nan or inf at grid point (%lf,%lf,%lf)\n", a, b, c, hx, hy, hz);
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
        invLambda[0][1] = gamma * bs_vx;
        invLambda[0][2] = gamma * bs_vy;
        invLambda[0][3] = gamma * bs_vz;
        invLambda[1][0] = gamma * bs_vx;
        invLambda[1][1] = 1. + gamma2 / (gamma + 1.) * (bs_vx * bs_vx);
        invLambda[1][2] = gamma2 / (gamma + 1.) * (bs_vx * bs_vy);
        invLambda[1][3] = gamma2 / (gamma + 1.) * (bs_vx * bs_vz);
        invLambda[2][0] = gamma * bs_vy;
        invLambda[2][1] = gamma2 / (gamma + 1.) * (bs_vy * bs_vx);
        invLambda[2][2] = 1. + gamma2 / (gamma + 1.) * (bs_vy * bs_vy);
        invLambda[2][3] = gamma2 / (gamma + 1.) * (bs_vy * bs_vz);
        invLambda[3][0] = gamma * bs_vz;
        invLambda[3][1] = gamma2 / (gamma + 1.) * (bs_vz * bs_vx);
        invLambda[3][2] = gamma2 / (gamma + 1.) * (bs_vz * bs_vy);
        invLambda[3][3] = 1. + gamma2 / (gamma + 1.) * (bs_vz * bs_vz);

        CCTK_REAL Lambda[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            Lambda[a][b] = 0.0;
          }
        }
        Lambda[0][0] = gamma;
        Lambda[0][1] = -gamma * bs_vx;
        Lambda[0][2] = -gamma * bs_vy;
        Lambda[0][3] = -gamma * bs_vz;
        Lambda[1][0] = -gamma * bs_vx;
        Lambda[1][1] = 1. + gamma2 / (gamma + 1.) * (bs_vx * bs_vx);
        Lambda[1][2] = gamma2 / (gamma + 1.) * (bs_vx * bs_vy);
        Lambda[1][3] = gamma2 / (gamma + 1.) * (bs_vx * bs_vz);
        Lambda[2][0] = -gamma * bs_vy;
        Lambda[2][1] = gamma2 / (gamma + 1.) * (bs_vy * bs_vx);
        Lambda[2][2] = 1. + gamma2 / (gamma + 1.) * (bs_vy * bs_vy);
        Lambda[2][3] = gamma2 / (gamma + 1.) * (bs_vy * bs_vz);
        Lambda[3][0] = -gamma * bs_vz;
        Lambda[3][1] = gamma2 / (gamma + 1.) * (bs_vz * bs_vx);
        Lambda[3][2] = gamma2 / (gamma + 1.) * (bs_vz * bs_vy);
        Lambda[3][3] = 1. + gamma2 / (gamma + 1.) * (bs_vz * bs_vz);

        // Boosted metric

        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            CCTK_REAL sum = 0.0;
            for (int chi = 0; chi < 4; ++chi) // making \mu\ -> \chi to avoid conflict with field mass \mu
              for (int nu = 0; nu < 4; ++nu)
                sum += invLambda[chi][a] * invLambda[nu][b] * G[chi][nu];
            Gb[a][b] = sum;
          }
        }

        gammaA_inv[1][1] = (-(exp(2 * F0[ind] + 2 * F2[ind]) * pow(bs_vz, 2) * pow(hx, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F1[ind]) * (pow(bs_vz, 2) * pow(hy, 2) + pow(bs_vy, 2) * (pow(hx, 2) + pow(hy, 2))) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(4 * F1[ind]) * pow(hy + hy * gamma + (bs_vx * bs_vy * hx + (pow(bs_vy, 2) + pow(bs_vz, 2)) * hy) * pow(gamma, 2), 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(bs_vx, 2) * (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2) * pow(gamma, 4) - 2 * bs_vx * bs_vy * hx * hy * pow(gamma, 2) * (1 + gamma + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2)) + pow(hx, 2) * (1 + gamma * (2 + gamma + pow(bs_vy, 4) * pow(gamma, 3) + 2 * pow(bs_vy, 2) * gamma * (1 + gamma + pow(bs_vz, 2) * pow(gamma, 2)) + pow(bs_vz, 2) * gamma * (2 + gamma * (2 + (pow(bs_vx, 2) + pow(bs_vz, 2)) * gamma)))))) / (exp(2. * F1[ind]) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * ((pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * hx * hy + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2)) * pow(1 + gamma + (pow(bs_vx, 2) + pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2), 2)));

        gammaA_inv[1][2] = (-(exp(2 * F0[ind] + 2 * F2[ind]) * pow(bs_vz, 2) * hx * hy * pow(gamma, 2) * pow(1 + gamma, 2)) + exp(2 * F0[ind] + 2 * F1[ind]) * (pow(bs_vz, 2) * hx * hy + bs_vx * bs_vy * (pow(hx, 2) + pow(hy, 2))) * pow(gamma, 2) * pow(1 + gamma, 2) - exp(4 * F1[ind]) * (hx + hx * gamma + ((pow(bs_vx, 2) + pow(bs_vz, 2)) * hx + bs_vx * bs_vy * hy) * pow(gamma, 2)) * (hy + hy * gamma + (bs_vx * bs_vy * hx + (pow(bs_vy, 2) + pow(bs_vz, 2)) * hy) * pow(gamma, 2)) + exp(2 * F1[ind] + 2 * F2[ind]) * (-(bs_vx * bs_vy * pow(hy, 2) * pow(gamma, 2) * (1 + gamma + pow(bs_vx, 2) * pow(gamma, 2))) - bs_vx * bs_vy * pow(hx, 2) * pow(gamma, 2) * (1 + gamma + pow(bs_vy, 2) * pow(gamma, 2)) + hx * hy * (1 + gamma * (2 + gamma * (1 + pow(bs_vz, 4) * pow(gamma, 2) + 2 * pow(bs_vz, 2) * (1 + gamma) + pow(bs_vy, 2) * (1 + gamma + pow(bs_vz, 2) * pow(gamma, 2)) + pow(bs_vx, 2) * (1 + gamma + (2 * pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2))))))) / (exp(2. * F1[ind]) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * ((pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * hx * hy + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2)) * pow(1 + gamma + (pow(bs_vx, 2) + pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2), 2)));

        gammaA_inv[1][3] = (bs_vz * pow(gamma, 2) * (exp(2 * F0[ind] + 2 * F1[ind]) * hy * (-(bs_vy * hx) + bs_vx * hy) * pow(1 + gamma, 2) + exp(2 * F0[ind] + 2 * F2[ind]) * hx * (bs_vx * hx + bs_vy * hy) * pow(1 + gamma, 2) + exp(4 * F1[ind]) * (bs_vy * hx - bs_vx * hy) * (hy + hy * gamma + (bs_vx * bs_vy * hx + (pow(bs_vy, 2) + pow(bs_vz, 2)) * hy) * pow(gamma, 2)) - exp(2 * F1[ind] + 2 * F2[ind]) * (-(pow(bs_vx, 2) * bs_vy * hx * hy * pow(gamma, 2)) + pow(bs_vx, 3) * (pow(hx, 2) + pow(hy, 2)) * pow(gamma, 2) + bs_vy * hx * hy * (1 + gamma + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2)) + bs_vx * (pow(hy, 2) * (1 + gamma) + pow(hx, 2) * (2 + gamma * (2 + 2 * pow(bs_vy, 2) * gamma + pow(bs_vz, 2) * gamma)))))) / (exp(2. * F1[ind]) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * ((pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * hx * hy + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2)) * pow(1 + gamma + (pow(bs_vx, 2) + pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2), 2)));

        gammaA_inv[2][1] = gammaA_inv[1][2];

        gammaA_inv[2][2] = (-(exp(2 * F0[ind] + 2 * F2[ind]) * pow(bs_vz, 2) * pow(hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) -
                            exp(2 * F0[ind] + 2 * F1[ind]) * (pow(bs_vz, 2) * pow(hx, 2) + pow(bs_vx, 2) * (pow(hx, 2) + pow(hy, 2))) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(4 * F1[ind]) * pow(hx + hx * gamma + ((pow(bs_vx, 2) + pow(bs_vz, 2)) * hx + bs_vx * bs_vy * hy) * pow(gamma, 2), 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(bs_vy, 2) * (pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) * pow(gamma, 4) - 2 * bs_vx * bs_vy * hx * hy * pow(gamma, 2) * (1 + gamma + (pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(gamma, 2)) + pow(hy, 2) * (1 + gamma * (2 + gamma + pow(bs_vx, 4) * pow(gamma, 3) + 2 * pow(bs_vx, 2) * gamma * (1 + gamma + pow(bs_vz, 2) * pow(gamma, 2)) + pow(bs_vz, 2) * gamma * (2 + gamma * (2 + (pow(bs_vy, 2) + pow(bs_vz, 2)) * gamma)))))) /
                           (exp(2. * F1[ind]) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * ((pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * hx * hy + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2)) * pow(1 + gamma + (pow(bs_vx, 2) + pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2), 2)));

        gammaA_inv[2][3] = -((bs_vz * pow(gamma, 2) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * hx * (bs_vy * hx - bs_vx * hy) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * hy * (bs_vx * hx + bs_vy * hy) * pow(1 + gamma, 2) + exp(4 * F1[ind]) * (bs_vy * hx - bs_vx * hy) * (hx + hx * gamma + ((pow(bs_vx, 2) + pow(bs_vz, 2)) * hx + bs_vx * bs_vy * hy) * pow(gamma, 2)) + exp(2 * F1[ind] + 2 * F2[ind]) * (-(bs_vx * pow(bs_vy, 2) * hx * hy * pow(gamma, 2)) + pow(bs_vy, 3) * (pow(hx, 2) + pow(hy, 2)) * pow(gamma, 2) + bs_vx * hx * hy * (1 + gamma + (pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(gamma, 2)) + bs_vy * (pow(hx, 2) * (1 + gamma) + pow(hy, 2) * (2 + gamma * (2 + 2 * pow(bs_vx, 2) * gamma + pow(bs_vz, 2) * gamma)))))) / (exp(2. * F1[ind]) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * ((pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * hx * hy + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2)) * pow(1 + gamma + (pow(bs_vx, 2) + pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2), 2))));

        gammaA_inv[3][1] = gammaA_inv[1][3];

        gammaA_inv[3][2] = gammaA_inv[2][3];

        gammaA_inv[3][3] = (exp(4 * F1[ind]) * pow(bs_vz, 2) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 4) - exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2) - exp(2 * F0[ind] + 2 * F2[ind]) * pow(bs_vx * hx + bs_vy * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2) + 2 * (pow(hx, 2) + pow(hy, 2)) * gamma + (1 + 2 * pow(bs_vx, 2) + 2 * pow(bs_vy, 2)) * (pow(hx, 2) + pow(hy, 2)) * pow(gamma, 2) + 2 * (pow(bs_vx, 2) + pow(bs_vy, 2)) * (pow(hx, 2) + pow(hy, 2)) * pow(gamma, 3) + ((pow(pow(bs_vx, 2) + pow(bs_vy, 2), 2) + pow(bs_vx, 2) * pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * pow(bs_vz, 2) * hx * hy + (pow(pow(bs_vx, 2) + pow(bs_vy, 2), 2) + pow(bs_vy, 2) * pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 4))) / (exp(2. * F1[ind]) * (-(exp(2 * F0[ind] + 2 * F1[ind]) * pow(bs_vy * hx - bs_vx * hy, 2) * pow(gamma, 2) * pow(1 + gamma, 2)) - exp(2 * F0[ind] + 2 * F2[ind]) * ((pow(bs_vx, 2) + pow(bs_vz, 2)) * pow(hx, 2) + 2 * bs_vx * bs_vy * hx * hy + (pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(hy, 2)) * pow(gamma, 2) * pow(1 + gamma, 2) + exp(2 * F1[ind] + 2 * F2[ind]) * (pow(hx, 2) + pow(hy, 2)) * pow(1 + gamma + (pow(bs_vx, 2) + pow(bs_vy, 2) + pow(bs_vz, 2)) * pow(gamma, 2), 2)));

        // // Build spatial metric from boosted 4-metric Gb
        // CCTK_REAL gammaA[3][3] = {
        //   { Gb[1][1], Gb[1][2], Gb[1][3] },
        //   { Gb[2][1], Gb[2][2], Gb[2][3] },
        //   { Gb[3][1], Gb[3][2], Gb[3][3] }
        // };

        // // Invert 3x3 gammaA numerically (cofactor formula)
        // CCTK_REAL det =
        //     gammaA[0][0]*(gammaA[1][1]*gammaA[2][2]-gammaA[1][2]*gammaA[2][1])
        //   - gammaA[0][1]*(gammaA[1][0]*gammaA[2][2]-gammaA[1][2]*gammaA[2][0])
        //   + gammaA[0][2]*(gammaA[1][0]*gammaA[2][1]-gammaA[1][1]*gammaA[2][0]);

        // if (fabs(det) < 1e-30) {
        //   fprintf(stderr,"Error: det(gammaA) ~ 0 at (%lf,%lf,%lf)\n", hx,hy,hz);
        // }

        // CCTK_REAL invdet = 1.0/det;
        // gammaA_inv[1][1] =  (gammaA[1][1]*gammaA[2][2]-gammaA[1][2]*gammaA[2][1])*invdet;
        // gammaA_inv[1][2] = -(gammaA[0][1]*gammaA[2][2]-gammaA[0][2]*gammaA[2][1])*invdet;
        // gammaA_inv[1][3] =  (gammaA[0][1]*gammaA[1][2]-gammaA[0][2]*gammaA[1][1])*invdet;
        // gammaA_inv[2][1] = -(gammaA[1][0]*gammaA[2][2]-gammaA[1][2]*gammaA[2][0])*invdet;
        // gammaA_inv[2][2] =  (gammaA[0][0]*gammaA[2][2]-gammaA[0][2]*gammaA[2][0])*invdet;
        // gammaA_inv[2][3] = -(gammaA[0][0]*gammaA[1][2]-gammaA[0][2]*gammaA[1][0])*invdet;
        // gammaA_inv[3][1] =  (gammaA[1][0]*gammaA[2][1]-gammaA[1][1]*gammaA[2][0])*invdet;
        // gammaA_inv[3][2] = -(gammaA[0][0]*gammaA[2][1]-gammaA[0][1]*gammaA[2][0])*invdet;
        // gammaA_inv[3][3] =  (gammaA[0][0]*gammaA[1][1]-gammaA[0][1]*gammaA[1][0])*invdet;

        // Check for NaN or Inf in gammaA_inv
        for (int a = 1; a < 4; ++a) {
          for (int b = 1; b < 4; ++b) {
            if (isnan(gammaA_inv[a][b]) || isinf(gammaA_inv[a][b])) {
              fprintf(stderr, "Error: gammaA_inv[%d][%d] is nan or inf at grid point (%lf,%lf,%lf)\n", a, b, hx, hy, hz);
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
              for (int chi = 0; chi < 4; ++chi)
                for (int nu = 0; nu < 4; ++nu)
                  for (int lam = 0; lam < 4; ++lam)
                    sum += invLambda[chi][a] * invLambda[nu][b] * invLambda[lam][c] * dG[chi][nu][lam];
              dGb[a][b][c] = sum;
            }
          }
        }
        // Check for NaN or Inf in dGb
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
              if (isnan(dGb[a][b][c]) || isinf(dGb[a][b][c])) {
                fprintf(stderr, "Error: dGb[%d][%d][%d] is nan or inf at grid point (%lf,%lf,%lf)\n", a, b, c, hx, hy, hz);
              }
            }
          }
        }

        // Now we compute the 3+1 quantities

        // Shift
        CCTK_REAL betad[4], betaup[4];
        betad[0] = 0;
        betad[1] = Gb[0][1];
        betad[2] = Gb[0][2];
        betad[3] = Gb[0][3];
        betaup[0] = 0;
        betaup[1] = gammaA_inv[1][1] * betad[1] + gammaA_inv[1][2] * betad[2] + gammaA_inv[1][3] * betad[3];
        betaup[2] = gammaA_inv[2][1] * betad[1] + gammaA_inv[2][2] * betad[2] + gammaA_inv[2][3] * betad[3];
        betaup[3] = gammaA_inv[3][1] * betad[1] + gammaA_inv[3][2] * betad[2] + gammaA_inv[3][3] * betad[3];

        // Lapse
        const CCTK_REAL alpha2 = -Gb[0][0] + betaup[1] * betad[1] + betaup[2] * betad[2] + betaup[3] * betad[3];
        if (alpha2 < 0) {
          fprintf(stderr, "Error: negative argument in sqrt for alpha, alpha2=%lf at grid point (%lf,%lf,%lf)\n", alpha2, hx, hy, hz);
        }
        const CCTK_REAL alpha = sqrt(alpha2);

        // Check for NaN in betad and betaup
        for (int idx = 0; idx < 4; ++idx) {
          if (isnan(betad[idx]) || isinf(betad[idx])) {
            fprintf(stderr, "Error: betad[%d] is NaN at grid point (%lf,%lf,%lf)\n", idx, hx, hy, hz);
          }
          if (isnan(betaup[idx]) || isinf(betaup[idx])) {
            fprintf(stderr, "Error: betaup[%d] is NaN at grid point (%lf,%lf,%lf)\n", idx, hx, hy, hz);
          }
        }

        check_nan_or_inf("1/alpha", 1 / alpha);

        CCTK_REAL K_A[4][4]; // extrinsic curvature
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            K_A[a][b] = 0.0;
          }
        } // K_0\chi might not be zero but irrelevant for what i want to compute

        for (int a = 1; a < 4; ++a) {
          for (int b = 1; b < 4; ++b) {
            CCTK_REAL sum1 = 0.0;
            CCTK_REAL sum2 = 0.0;
            CCTK_REAL sum3 = 0.0;
            for (int c = 1; c < 4; ++c) {
              sum1 += betaup[c] * dGb[a][b][c];
              sum2 += betaup[c] * dGb[b][c][a];
              sum3 += betaup[c] * dGb[a][c][b];
            }
            K_A[a][b] = (-1 / (2. * alpha) * (dGb[a][b][0] - sum1 - (dGb[0][b][a] - sum2) - (dGb[0][a][b] - sum3))) * mu;
          }
        }

        for (int a = 1; a < 4; ++a) {
          for (int b = 1; b < 4; ++b) {
            if (isnan(K_A[a][b]) || isinf(K_A[a][b])) {
              fprintf(stderr, "Error: K_{%d,%d} is nan at grid point (%lf,%lf,%lf)\n", a, b, hx, hy, hz);
            }
          }
        }

        CCTK_REAL gammaB[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            gammaB[a][b] = 0.0;
          }
        }

        gammaB[1][1] = gxx[ind];
        gammaB[1][2] = gxy[ind];
        gammaB[1][3] = gxz[ind];
        gammaB[2][1] = gammaB[1][2];
        gammaB[2][2] = gyy[ind];
        gammaB[2][3] = gyz[ind];
        gammaB[3][1] = gammaB[1][3];
        gammaB[3][2] = gammaB[2][3];
        gammaB[3][3] = gzz[ind];

        CCTK_REAL gammaB_inv[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            gammaB_inv[a][b] = 0.0;
          }
        }

        // FAZER A INVERSÃO NUMÉRICA!!!!

        // Build symmetric 3x3 block M from Gb
        CCTK_REAL M[4][4] = {{0}};
        M[1][1] = gammaB[1][1];
        M[1][2] = 0.5 * (gammaB[1][2] + gammaB[2][1]);
        M[1][3] = 0.5 * (gammaB[1][3] + gammaB[3][1]);
        M[2][1] = M[1][2];
        M[2][2] = gammaB[2][2];
        M[2][3] = 0.5 * (gammaB[2][3] + gammaB[3][2]);
        M[3][1] = M[1][3];
        M[3][2] = M[2][3];
        M[3][3] = gammaB[3][3];

        if (!invert_spd3x3(M, gammaB_inv)) {
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "gammaB not SPD or ill-conditioned at (%d,%d,%d).",
                     i, j, k);
        }

        // Optional checks
        check_nan_or_inf("gammaB_inv[1][1]", gammaB_inv[1][1]);
        check_nan_or_inf("gammaB_inv[1][2]", gammaB_inv[1][2]);
        check_nan_or_inf("gammaB_inv[1][3]", gammaB_inv[1][3]);
        check_nan_or_inf("gammaB_inv[2][2]", gammaB_inv[2][2]);
        check_nan_or_inf("gammaB_inv[2][3]", gammaB_inv[2][3]);
        check_nan_or_inf("gammaB_inv[3][3]", gammaB_inv[3][3]);

        // Verify gammaB_inv * gammaB ≈ I (indices 1..3)
        {
          const CCTK_REAL tol = SMALL;
          for (int i3 = 1; i3 <= 3; ++i3) {
            for (int j3 = 1; j3 <= 3; ++j3) {
              CCTK_REAL s = 0.0;
              for (int k3 = 1; k3 <= 3; ++k3)
                s += gammaB_inv[i3][k3] * gammaB[k3][j3];
              const CCTK_REAL delta = (i3 == j3) ? 1.0 : 0.0;
              if (fabs(s - delta) > tol) {
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "gammaB_inv check failed at (%d,%d): %g (tol=%g)", i3,
                           j3, (double)s, (double)tol);
                break;
              }
            }
          }
        }

        CCTK_REAL separation = (center_offset[0] + 1) - x0; // only for separations along the x-axis, need to be modified for general case

        // 3-metric (added Bowen-York 3-metric)
        gxx[ind] = gammaB[1][1] + Gb[1][1] - pow(1 + par_m_plus / (2 * separation), 4);
        gxy[ind] = gammaB[1][2] + Gb[1][2];
        gxz[ind] = gammaB[1][3] + Gb[1][3];
        gyy[ind] = gammaB[2][2] + Gb[2][2] - pow(1 + par_m_plus / (2 * separation), 4);
        gyz[ind] = gammaB[2][3] + Gb[2][3];
        gzz[ind] = gammaB[3][3] + Gb[3][3] - pow(1 + par_m_plus / (2 * separation), 4);

        check_nan_or_inf("gxx", gxx[ind]);
        check_nan_or_inf("gxy", gxy[ind]);
        check_nan_or_inf("gxz", gxz[ind]);
        check_nan_or_inf("gyy", gyy[ind]);
        check_nan_or_inf("gyz", gyz[ind]);
        check_nan_or_inf("gzz", gzz[ind]);

        ///////////////////////////////////////////
        //aqui ate deveria ser do 2º buraco negro dado pelo two punctures que tem de ser adicionado ao do buraco negro original.

        static inline double delta_ij(int i, int j) { return (i == j) ? 1.0 : 0.0; }
        
        // n_i = x_i / r
        CCTK_REAL n[3] = {hx / rr, hy / rr, hz / rr};

        // P·n
        const CCTK_REAL Pdotn = par_P_plus[0] * n[0] + par_P_plus[1] * n[1] + par_P_plus[2] * n[2];

        // w = n × S  (components w_i = eps_{ijk} n_j S_k)
        CCTK_REAL w[3] = {
            n[1] * par_S_plus[2] - n[2] * par_S_plus[1],
            n[2] * par_S_plus[0] - n[0] * par_S_plus[2],
            n[0] * par_S_plus[1] - n[1] * par_S_plus[0]};

        const CCTK_REAL cP = 3.0 / (2.0 * rr2);
        const CCTK_REAL cS = -3.0 / (rr * rr2);

        CCTK_REAL A[4][4];

        //aqui os indices i,j vão de 0 a 2
        for (int i = 0; i < 3; ++i) {
          for (int j = i; j < 3; ++j) {
            const CCTK_REAL mom =
                cP * (n[i] * par_P_plus[j] + n[j] * par_P_plus[i] + Pdotn * (n[i] * n[j] - delta_ij(i, j)));

            const CCTK_REAL spin =
                cS * (n[i] * w[j] + n[j] * w[i]);

            A[i][j] = mom + spin;
            A[j][i] = A[i][j];
          }
        }

        CCTK_REAL Ktp[4][4];
        for (int i = 1; i < 4; ++i)
        {
          for (int j = 1; j < 4; j++)
          {
            Ktp[i][j] = 1/psi2 * A[i-1][j-1];
          }
        }
        
        ////////////////////////////////////////////

        CCTK_REAL dW_drho, dW_dz;
        const CCTK_REAL exp_auxi = exp(2. * F2[ind] - F0[ind]);

        if (rho < 1e-8) {
          dW_drho = 0.;
          dW_dz = 0.;
        } else {
          dW_drho = rho / rr * dW_dr[ind] + hz / rr2 * dW_dth[ind];
          dW_dz = hz / rr * dW_dr[ind] - rho / rr2 * dW_dth[ind];
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
        } // K_0\chi might not be zero but irrelevant for what i want to compute

        K_B[1][1] = kxx[ind];
        K_B[1][2] = kxy[ind];
        K_B[1][3] = kxz[ind];
        K_B[2][1] = K_B[1][2];
        K_B[2][2] = kyy[ind];
        K_B[2][3] = kyz[ind];
        K_B[3][1] = K_B[1][3];
        K_B[3][2] = K_B[2][3];
        K_B[3][3] = kzz[ind];

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

        check_nan_or_inf("kxx", kxx[ind]);
        check_nan_or_inf("kxy", kxy[ind]);
        check_nan_or_inf("kxz", kxz[ind]);
        check_nan_or_inf("kyy", kyy[ind]);
        check_nan_or_inf("kyz", kyz[ind]);
        check_nan_or_inf("kzz", kzz[ind]);

        // lapse value (field initialization below)
        // No lapse regularization needed for the BS, the lapse is non-zero
        // pre boost since i want to compute unboosted quantities
        const CCTK_REAL alph = exp(F0[ind]); // + alpha0 - 1;

        // let's add a perturbation to the Proca field as well
        // NOTE: the perturbation is added directed to every instance of e^{i m \varphi}, hence its derivatives are not taken into account
        // TODO (?): Design perturbation more generically as ~ cos((m+1)\varphi)
        const CCTK_REAL argpert_Proca = (rr - R0pert_Proca) / Sigmapert_Proca;
        const CCTK_REAL pert_Proca = 1. + Apert_Proca * (hx * hx - hy * hy) * mu * mu * exp(-0.5 * argpert_Proca * argpert_Proca); // ignorar por agora

        // ----- Proca fields -----

        // TODO: check what happens with divisions by rr and sinth, can we work around them?

        // Real and imaginay part of the harmonic dependence: exp[i(m\varphi - \omega t)] rotation not implemented as of yet
        const CCTK_REAL harm_re = (coswt * cosmph + sinwt * sinmph) * pert_Proca;
        const CCTK_REAL harm_im = (coswt * sinmph - sinwt * cosmph) * pert_Proca;

        // No need to change the radial component, R and r coincide

        CCTK_REAL A1_unboosted[4]; // A_\chi real part
        CCTK_REAL A2_unboosted[4]; // A_\chi imag part

        // A_t
        A1_unboosted[0] = V[ind] * sinwt;
        A2_unboosted[0] = V[ind] * coswt;

        // A_x
        A1_unboosted[1] = hx / rr * H1r[ind] * harm_re + costh * cosph / rr * H2[ind] * harm_re + sinph / rr * H3[ind] * harm_im;
        A2_unboosted[1] = hx / rr * H1r[ind] * harm_im + costh * cosph / rr * H2[ind] * harm_im - sinph / rr * H3[ind] * harm_re;

        // A_y
        A1_unboosted[2] = hy / rr * H1r[ind] * harm_re + costh * sinph / rr * H2[ind] * harm_re - cosph / rr * H3[ind] * harm_im;
        A2_unboosted[2] = hy / rr * H1r[ind] * harm_im + costh * sinph / rr * H2[ind] * harm_im + cosph / rr * H3[ind] * harm_re;

        // A_z
        A1_unboosted[3] = (hz / rr * H1r[ind] - sinth / rr * H2[ind]) * harm_re;
        A2_unboosted[3] = (hz / rr * H1r[ind] - sinth / rr * H2[ind]) * harm_im;

        const CCTK_REAL dH1r_dr = dH1_dr[ind] / rr - H1r[ind] / rr;
        const CCTK_REAL dH1r_dth = dH1_dth[ind] / rr;

        // Build unboosted field-strength tensor F_{chi nu} (only 0i components from time/spatial derivatives)
        CCTK_REAL F1_unb[4][4], F2_unb[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            F1_unb[a][b] = 0.0;
            F2_unb[a][b] = 0.0;
          }
        }

        const CCTK_REAL dV_dx = dV_dr[ind] * R_x + dV_dth[ind] * th_x;
        const CCTK_REAL dV_dy = dV_dr[ind] * R_y + dV_dth[ind] * th_y;
        const CCTK_REAL dV_dz = dV_dr[ind] * R_z + dV_dth[ind] * th_z;

        const CCTK_REAL dH1r_dx = dH1r_dr * R_x + dH1r_dth * th_x;
        const CCTK_REAL dH1r_dy = dH1r_dr * R_y + dH1r_dth * th_y;
        const CCTK_REAL dH1r_dz = dH1r_dr * R_z + dH1r_dth * th_z;

        const CCTK_REAL dH2_dx = dH2_dr[ind] * R_x + dH2_dth[ind] * th_x;
        const CCTK_REAL dH2_dy = dH2_dr[ind] * R_y + dH2_dth[ind] * th_y;
        const CCTK_REAL dH2_dz = dH2_dr[ind] * R_z + dH2_dth[ind] * th_z;

        const CCTK_REAL dH3_dx = dH3_dr[ind] * R_x + dH3_dth[ind] * th_x;
        const CCTK_REAL dH3_dy = dH3_dr[ind] * R_y + dH3_dth[ind] * th_y;
        const CCTK_REAL dH3_dz = dH3_dr[ind] * R_z + dH3_dth[ind] * th_z;

        const CCTK_REAL dA1x_dt = omega_BS * (hx / rr * H1r[ind] * harm_im + costh * cosph / rr * H2[ind] * harm_im - sinph / rr * H3[ind] * harm_re);
        const CCTK_REAL dA1y_dt = omega_BS * (hy / rr * H1r[ind] * harm_im + costh * sinph / rr * H2[ind] * harm_im + cosph / rr * H3[ind] * harm_re);
        const CCTK_REAL dA1z_dt = omega_BS * ((hz / rr * H1r[ind] - sinth / rr * H2[ind]) * harm_im);

        const CCTK_REAL dA2x_dt = -omega_BS * (hx / rr * H1r[ind] * harm_re + costh * cosph / rr * H2[ind] * harm_re + sinph / rr * H3[ind] * harm_im);
        const CCTK_REAL dA2y_dt = -omega_BS * (hy / rr * H1r[ind] * harm_re + costh * sinph / rr * H2[ind] * harm_re - cosph / rr * H3[ind] * harm_im);
        const CCTK_REAL dA2z_dt = -omega_BS * ((hz / rr * H1r[ind] - sinth / rr * H2[ind]) * harm_re);

        const CCTK_REAL dA1t_dx = dV_dx * sinwt;
        const CCTK_REAL dA1t_dy = dV_dy * sinwt;
        const CCTK_REAL dA1t_dz = dV_dz * sinwt;

        const CCTK_REAL dA2t_dx = dV_dx * coswt;
        const CCTK_REAL dA2t_dy = dV_dy * coswt;
        const CCTK_REAL dA2t_dz = dV_dz * coswt;

        const CCTK_REAL dA1x_dx = (-(R_x * ((cosph * costh * H2[ind] + H1r[ind] * hx) * coswt -
                                            H3[ind] * sinph * sinwt)) +
                                   rr * (-((d_sinph_dx * H3[ind] +
                                            dH3_dx * sinph) *
                                           sinwt) +
                                         coswt * (H1r[ind] + costh * (cosph * dH2_dx + d_cosph_dx * H2[ind]) + dH1r_dx * hx + cosph * H2[ind] * (-sinth * th_x)))) /
                                  rr2;
        const CCTK_REAL dA1x_dy = (-(R_y * ((cosph * costh * H2[ind] + H1r[ind] * hx) * coswt -
                                            H3[ind] * sinph * sinwt)) +
                                   rr * (-((d_sinph_dy * H3[ind] +
                                            dH3_dy * sinph) *
                                           sinwt) +
                                         coswt * (costh * (cosph * dH2_dy +
                                                           d_cosph_dy * H2[ind]) +
                                                  dH1r_dy * hx + cosph * H2[ind] * (-sinth * th_y)))) /
                                  rr2;
        const CCTK_REAL dA1x_dz = (-(R_z * ((cosph * costh * H2[ind] + H1r[ind] * hx) * coswt -
                                            H3[ind] * sinph * sinwt)) +
                                   rr * (-(dH3_dz * sinph * sinwt) +
                                         coswt * (dH1r_dz * hx + cosph * (costh * dH2_dz + H2[ind] * (-sinth *
                                                                                                      th_z))))) /
                                  rr2;

        const CCTK_REAL dA1y_dx = (-(R_x * ((costh * H2[ind] * sinph + H1r[ind] * hy) * coswt +
                                            cosph * H3[ind] * sinwt)) +
                                   rr * ((cosph * dH3_dx +
                                          d_cosph_dx * H3[ind]) *
                                             sinwt +
                                         coswt * (costh * (d_sinph_dx * H2[ind] +
                                                           dH2_dx * sinph) +
                                                  dH1r_dx * hy + H2[ind] * sinph * (-sinth * th_x)))) /
                                  rr2;
        const CCTK_REAL dA1y_dy = (-(R_y * ((costh * H2[ind] * sinph + H1r[ind] * hy) * coswt +
                                            cosph * H3[ind] * sinwt)) +
                                   rr * ((cosph * dH3_dy +
                                          d_cosph_dy * H3[ind]) *
                                             sinwt +
                                         coswt * (H1r[ind] +
                                                  costh * (d_sinph_dy * H2[ind] + dH2_dy * sinph) + dH1r_dy * hy +
                                                  H2[ind] * sinph * (-sinth * th_y)))) /
                                  rr2;
        const CCTK_REAL dA1y_dz = (-(R_z * ((costh * H2[ind] * sinph + H1r[ind] * hy) * coswt +
                                            cosph * H3[ind] * sinwt)) +
                                   rr * (cosph * dH3_dz * sinwt +
                                         coswt * (dH1r_dz * hy + sinph * (costh * dH2_dz + H2[ind] * (-sinth *
                                                                                                      th_z))))) /
                                  rr2;

        const CCTK_REAL dA1z_dx = (coswt * (R_x * (H2[ind] * sinth - H1r[ind] * hz) +
                                            rr * (-(dH2_dx * sinth) + dH1r_dx * hz - H2[ind] * costh * th_x))) /
                                  rr2;
        const CCTK_REAL dA1z_dy = (coswt * (R_y * (H2[ind] * sinth - H1r[ind] * hz) +
                                            rr * (-(dH2_dy * sinth) + dH1r_dy * hz - H2[ind] * costh * th_y))) /
                                  rr2;
        const CCTK_REAL dA1z_dz = (coswt * (H2[ind] * R_z * sinth + H1r[ind] * (rr - R_z * hz) +
                                            rr * (-(dH2_dz * sinth) + dH1r_dz * hz - H2[ind] * costh * th_z))) /
                                  rr2;

        const CCTK_REAL dA2x_dx = (R_x * (H3[ind] * sinph * coswt + (cosph * costh * H2[ind] +
                                                                     H1r[ind] * hx) *
                                                                        sinwt) +
                                   rr * (-((d_sinph_dx * H3[ind] +
                                            dH3_dx * sinph) *
                                           coswt) -
                                         sinwt * (H1r[ind] + costh * (cosph * dH2_dx + d_cosph_dx * H2[ind]) + dH1r_dx * hx + cosph * H2[ind] * (-sinth * th_x)))) /
                                  rr2;
        const CCTK_REAL dA2x_dy = (R_y * (H3[ind] * sinph * coswt + (cosph * costh * H2[ind] +
                                                                     H1r[ind] * hx) *
                                                                        sinwt) +
                                   rr * (-((d_sinph_dy * H3[ind] +
                                            dH3_dy * sinph) *
                                           coswt) -
                                         sinwt * (costh * (cosph * dH2_dy +
                                                           d_cosph_dy * H2[ind]) +
                                                  dH1r_dy * hx + cosph * H2[ind] * (-sinth * th_y)))) /
                                  rr2;
        const CCTK_REAL dA2x_dz = (R_z * (H3[ind] * sinph * coswt + (cosph * costh * H2[ind] +
                                                                     H1r[ind] * hx) *
                                                                        sinwt) +
                                   rr * (-(dH3_dz * sinph * coswt) -
                                         sinwt * (dH1r_dz * hx + cosph * (costh * dH2_dz + H2[ind] * (-sinth *
                                                                                                      th_z))))) /
                                  rr2;

        const CCTK_REAL dA2y_dx = (R_x * (-(cosph * H3[ind] * coswt) + (costh * H2[ind] * sinph +
                                                                        H1r[ind] * hy) *
                                                                           sinwt) +
                                   rr * ((cosph * dH3_dx +
                                          d_cosph_dx * H3[ind]) *
                                             coswt -
                                         sinwt * (costh * (d_sinph_dx * H2[ind] +
                                                           dH2_dx * sinph) +
                                                  dH1r_dx * hy + H2[ind] * sinph * (-sinth * th_x)))) /
                                  rr2;
        const CCTK_REAL dA2y_dy = (R_y * (-(cosph * H3[ind] * coswt) + (costh * H2[ind] * sinph +
                                                                        H1r[ind] * hy) *
                                                                           sinwt) +
                                   rr * ((cosph * dH3_dy +
                                          d_cosph_dy * H3[ind]) *
                                             coswt -
                                         sinwt * (H1r[ind] +
                                                  costh * (d_sinph_dy * H2[ind] + dH2_dy * sinph) + dH1r_dy * hy +
                                                  H2[ind] * sinph * (-sinth * th_y)))) /
                                  rr2;
        const CCTK_REAL dA2y_dz = (R_z * (-(cosph * H3[ind] * coswt) + (costh * H2[ind] * sinph +
                                                                        H1r[ind] * hy) *
                                                                           sinwt) +
                                   rr * (cosph * dH3_dz * coswt -
                                         sinwt * (dH1r_dz * hy + sinph * (costh * dH2_dz + H2[ind] * (-sinth *
                                                                                                      th_z))))) /
                                  rr2;

        const CCTK_REAL dA2z_dx = (sinwt * (R_x * (-(H2[ind] * sinth) + H1r[ind] * hz) +
                                            rr * (dH2_dx * sinth - dH1r_dx * hz + H2[ind] * costh * th_x))) /
                                  rr2;
        const CCTK_REAL dA2z_dy = (sinwt * (R_y * (-(H2[ind] * sinth) + H1r[ind] * hz) +
                                            rr * (dH2_dy * sinth - dH1r_dy * hz + H2[ind] * costh * th_y))) /
                                  rr2;
        const CCTK_REAL dA2z_dz = (sinwt * (-(H2[ind] * R_z * sinth) + H1r[ind] * (-rr + R_z * hz) +
                                            rr * (dH2_dz * sinth - dH1r_dz * hz + H2[ind] * costh * th_z))) /
                                  rr2;

        // // Spatial derivatives of A_0 = V * {sinwt, coswt}
        // const CCTK_REAL dV_dx = dV_dr[ind]*R_x + dV_dth[ind]*th_x;
        // const CCTK_REAL dV_dy = dV_dr[ind]*R_y + dV_dth[ind]*th_y;
        // const CCTK_REAL dV_dz = dV_dr[ind]*R_z + dV_dth[ind]*th_z;

        // F_{0i} = d_t A_i - d_i A_0  (covariant indices)
        F1_unb[0][1] = dA1x_dt - dA1t_dx;
        F1_unb[1][0] = -F1_unb[0][1];
        F1_unb[0][2] = dA1y_dt - dA1t_dy;
        F1_unb[2][0] = -F1_unb[0][2];
        F1_unb[0][3] = dA1z_dt - dA1t_dz;
        F1_unb[3][0] = -F1_unb[0][3];
        F1_unb[1][2] = dA1y_dx - dA1x_dy;
        F1_unb[2][1] = -F1_unb[1][2];
        F1_unb[1][3] = dA1z_dx - dA1x_dz;
        F1_unb[3][1] = -F1_unb[1][3];
        F1_unb[2][3] = dA1z_dy - dA1y_dz;
        F1_unb[3][2] = -F1_unb[2][3];

        F2_unb[0][1] = dA2x_dt - dA2t_dx;
        F2_unb[1][0] = -F2_unb[0][1];
        F2_unb[0][2] = dA2y_dt - dA2t_dy;
        F2_unb[2][0] = -F2_unb[0][2];
        F2_unb[0][3] = dA2z_dt - dA2t_dz;
        F2_unb[3][0] = -F2_unb[0][3];
        F2_unb[1][2] = dA2y_dx - dA2x_dy;
        F2_unb[2][1] = -F2_unb[1][2];
        F2_unb[1][3] = dA2z_dx - dA2x_dz;
        F2_unb[3][1] = -F2_unb[1][3];
        F2_unb[2][3] = dA2z_dy - dA2y_dz;
        F2_unb[3][2] = -F2_unb[2][3];

        CCTK_REAL A1_boosted[4]; // A_\chi real part
        CCTK_REAL A2_boosted[4]; // A_\chi imag part
        // Boosted components
        for (int a = 0; a < 4; ++a) {
          A1_boosted[a] = 0.0;
          A2_boosted[a] = 0.0;
          for (int chi = 0; chi < 4; ++chi) {
            A1_boosted[a] += invLambda[chi][a] * A1_unboosted[chi];
            A2_boosted[a] += invLambda[chi][a] * A2_unboosted[chi];
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
          A_\phi = -n^\chi A_\chi = - (A_t + W*A_ph)/alpha
                = -i * e^{i (m ph - w t)} * (V + W H3 sinth) / alpha
        */
        // Aphi1[ind] = (V[ind] + W[ind] * sinth * H3[ind]) / alph * harm_im;
        // Aphi2[ind] =-(V[ind] + W[ind] * sinth * H3[ind]) / alph * harm_re;

        Aphi1[ind] = 1 / alpha * (-A1_boosted[0] + betaup[1] * A1_boosted[1] + betaup[2] * A1_boosted[2] + betaup[3] * A1_boosted[3]);
        Aphi2[ind] = 1 / alpha * (-A2_boosted[0] + betaup[1] * A2_boosted[1] + betaup[2] * A2_boosted[2] + betaup[3] * A2_boosted[3]);

        //------ Boosted field-strength tensor F_{chi nu} ------
        CCTK_REAL F1_boosted[4][4], F2_boosted[4][4];
        for (int a = 0; a < 4; ++a) {
          for (int b = 0; b < 4; ++b) {
            F1_boosted[a][b] = 0.0;
            F2_boosted[a][b] = 0.0;
            for (int chi = 0; chi < 4; ++chi) {
              for (int nu = 0; nu < 4; ++nu) {
                F1_boosted[a][b] += invLambda[chi][a] * invLambda[nu][b] * F1_unb[chi][nu];
                F2_boosted[a][b] += invLambda[chi][a] * invLambda[nu][b] * F2_unb[chi][nu];
              }
            }
          }
        }

        // Boosted components

        // cannot be boosted since its foliation dependent
        CCTK_REAL E1_boosted[4]; // E_\chi real part
        CCTK_REAL E2_boosted[4]; // E_\chi imag part
        for (int a = 0; a < 4; ++a) {
          E1_boosted[a] = 0.0;
          E2_boosted[a] = 0.0;
          // for (int chi = 0; chi < 4; ++chi) {
          //   E1_boosted[a] += Lambda[a][chi] * E1_unboosted[chi];
          //   E2_boosted[a] += Lambda[a][chi] * E2_unboosted[chi];
          // }
        }

        // E_\chi
        E1_boosted[1] = 1 / alpha * (F1_boosted[1][0] - betaup[2] * F1_boosted[1][2] - betaup[3] * F1_boosted[1][3]);
        E1_boosted[2] = 1 / alpha * (F1_boosted[2][0] - betaup[1] * F1_boosted[2][1] - betaup[3] * F1_boosted[2][3]);
        E1_boosted[3] = 1 / alpha * (F1_boosted[3][0] - betaup[1] * F1_boosted[3][1] - betaup[2] * F1_boosted[3][2]);

        E2_boosted[1] = 1 / alpha * (F2_boosted[1][0] - betaup[2] * F2_boosted[1][2] - betaup[3] * F2_boosted[1][3]);
        E2_boosted[2] = 1 / alpha * (F2_boosted[2][0] - betaup[1] * F2_boosted[2][1] - betaup[3] * F2_boosted[2][3]);
        E2_boosted[3] = 1 / alpha * (F2_boosted[3][0] - betaup[1] * F2_boosted[3][1] - betaup[2] * F2_boosted[3][2]);

        CCTK_REAL E1up_boosted[4]; // E^\chi real part
        CCTK_REAL E2up_boosted[4]; // E^\chi imag part
        for (int a = 0; a < 4; ++a) {
          E1up_boosted[a] = 0.0;
          E2up_boosted[a] = 0.0;
        }
        // E^\chi /can be raised using the 3 metric since it is a spatial vector/
        E1up_boosted[1] = gammaA_inv[1][1] * E1_boosted[1] + gammaA_inv[1][2] * E1_boosted[2] + gammaA_inv[1][3] * E1_boosted[3];
        E1up_boosted[2] = gammaA_inv[2][1] * E1_boosted[1] + gammaA_inv[2][2] * E1_boosted[2] + gammaA_inv[2][3] * E1_boosted[3];
        E1up_boosted[3] = gammaA_inv[3][1] * E1_boosted[1] + gammaA_inv[3][2] * E1_boosted[2] + gammaA_inv[3][3] * E1_boosted[3];

        E2up_boosted[1] = gammaA_inv[1][1] * E2_boosted[1] + gammaA_inv[1][2] * E2_boosted[2] + gammaA_inv[1][3] * E2_boosted[3];
        E2up_boosted[2] = gammaA_inv[2][1] * E2_boosted[1] + gammaA_inv[2][2] * E2_boosted[2] + gammaA_inv[2][3] * E2_boosted[3];
        E2up_boosted[3] = gammaA_inv[3][1] * E2_boosted[1] + gammaA_inv[3][2] * E2_boosted[2] + gammaA_inv[3][3] * E2_boosted[3];

        /* store spatial components E^\chi */
        E1x[ind] = E1up_boosted[1] * mu;
        E1y[ind] = E1up_boosted[2] * mu;
        E1z[ind] = E1up_boosted[3] * mu;

        E2x[ind] = E2up_boosted[1] * mu;
        E2y[ind] = E2up_boosted[2] * mu;
        E2z[ind] = E2up_boosted[3] * mu;

        check_nan_or_inf("E1x", E1x[ind]);
        check_nan_or_inf("E1y", E1y[ind]);
        check_nan_or_inf("E1z", E1z[ind]);
        check_nan_or_inf("E2x", E2x[ind]);
        check_nan_or_inf("E2y", E2y[ind]);
        check_nan_or_inf("E2z", E2z[ind]);

        // zero-initialize constraint damping variable Z
        Zeta1[ind] = 0;
        Zeta2[ind] = 0;

        // lapse
        if (CCTK_EQUALS(post_initial_lapse, "PS_single"))
          alp[ind] = alp[ind] + pow(psi1, initial_lapse_psi_exponent) - 1.0;
        else if (CCTK_EQUALS(post_initial_lapse, "PS_lapse_single")) {
          alp[ind] = alp[ind] + alpha - 1.0;
          if (alp[ind] < SMALL)
            alp[ind] = SMALL;
        }

        // shift
        if (CCTK_EQUALS(post_initial_shift, "PS_single")) {
          betax[ind] = betax[ind] + betaup[1];
          betay[ind] = betay[ind] + betaup[2];
          betaz[ind] = betaz[ind] + betaup[3];
        } else if (CCTK_EQUALS(post_initial_shift, "PS_zero")) {
          betax[ind] = betax[ind] + 0.0;
          betay[ind] = betay[ind] + 0.0;
          betaz[ind] = betaz[ind] + 0.0;
        }

      } /* for i */
    } /* for j */
  } /* for k */

  free(F1);
  free(F2);
  free(F0);
  free(W);
  free(H1r);
  free(H2);
  free(H3);
  free(V);
  free(dW_dr);
  free(dW_dth);
  free(dH3_dr);
  free(dH3_dth);
  free(dV_dr);
  free(dV_dth);

  return;
}