#pragma once
#include "cctk.h"

// Point where you evaluate (fill these from Kerr_BS.c)
typedef struct {
  CCTK_REAL bh_mass;
  CCTK_REAL bh_spin;
  CCTK_REAL bh_spin2;
  CCTK_REAL bh_v;
  CCTK_REAL x1_2, y1_2, z1_2;
  CCTK_REAL rr_2, rr2_2, rho_2, rho2_2;

  // First derivatives of coordinates
  CCTK_REAL R_x, R_y, R_z;
  CCTK_REAL th_x, th_y, th_z;
  CCTK_REAL costh, costh2, sinth, sinth2;

  // Coordinate second derivatives
  CCTK_REAL R_xx, R_xy, R_xz, R_yy, R_yz, R_zz;
  CCTK_REAL th_xx, th_xy, th_xz, th_yy, th_yz, th_zz;

  // rBL and (R,th) derivatives
  CCTK_REAL rBL, rBLp, rBLm, drBLdR;
  CCTK_REAL d2rBL_dR2, d2rBL_dRth, d2rBL_dth2;

  // Delt
  CCTK_REAL Delt, dDelt_dR, dDelt_dth;
  CCTK_REAL d2Delt_dR2, d2Delt_dRth, d2Delt_dth2;

  // Sigm
  CCTK_REAL Sigm, Sigm2, dSigm_dR, dSigm_dth;
  CCTK_REAL d2Sigm_dR2, d2Sigm_dRth, d2Sigm_dth2;

  // fctFF
  CCTK_REAL fctFF, dfctFF_dR, dfctFF_dth;
  CCTK_REAL d2fctFF_dR2, d2fctFF_dRth, d2fctFF_dth2;

  // fctGG
  CCTK_REAL fctGG, dfctGG_dR;
  CCTK_REAL d2fctGG_dR2, d2fctGG_dRth, d2fctGG_dth2;

  // fctHH
  CCTK_REAL fctHH, dfctHH_dR, dfctHH_dth;
  CCTK_REAL d2fctHH_dR2, d2fctHH_dRth, d2fctHH_dth2;

  // psi4_2
  CCTK_REAL psi4_2, dpsi4_2_dR, dpsi4_2_dth;
  CCTK_REAL d2psi4_2_dR2, d2psi4_2_dRth, d2psi4_2_dth2;

  // alpha0/alpha02
  CCTK_REAL alpha0, alpha02, dalpha02_dR, dalpha02_dth;
  CCTK_REAL d2alpha02_dR2, d2alpha02_dRth, d2alpha02_dth2;

  // gphiphi
  CCTK_REAL gphiphi, dgphiphi_dR, dgphiphi_dth;
  CCTK_REAL d2gphiphi_dR2, d2gphiphi_dRth, d2gphiphi_dth2;

  // bphiup
  CCTK_REAL bphiup, dbphiup_dR, dbphiup_dth;
  CCTK_REAL d2bphiup_dR2, d2bphiup_dRth, d2bphiup_dth2;

  // bphi
  CCTK_REAL bphi, dbphi_dR, dbphi_dth;
  CCTK_REAL d2bphi_dR2, d2bphi_dRth, d2bphi_dth2;

  // Cartesian second derivatives (six each)
  CCTK_REAL d2rBL_dx2, d2rBL_dy2, d2rBL_dz2, d2rBL_dxy, d2rBL_dxz, d2rBL_dyz;
  CCTK_REAL d2Delt_dx2, d2Delt_dy2, d2Delt_dz2, d2Delt_dxy, d2Delt_dxz, d2Delt_dyz;
  CCTK_REAL d2Sigm_dx2, d2Sigm_dy2, d2Sigm_dz2, d2Sigm_dxy, d2Sigm_dxz, d2Sigm_dyz;
  CCTK_REAL d2fctFF_dx2, d2fctFF_dy2, d2fctFF_dz2, d2fctFF_dxy, d2fctFF_dxz, d2fctFF_dyz;
  CCTK_REAL d2fctGG_dx2, d2fctGG_dy2, d2fctGG_dz2, d2fctGG_dxy, d2fctGG_dxz, d2fctGG_dyz;
  CCTK_REAL d2fctHH_dx2, d2fctHH_dy2, d2fctHH_dz2, d2fctHH_dxy, d2fctHH_dxz, d2fctHH_dyz;
  CCTK_REAL d2psi4_2_dx2, d2psi4_2_dy2, d2psi4_2_dz2, d2psi4_2_dxy, d2psi4_2_dxz, d2psi4_2_dyz;
  CCTK_REAL d2alpha02_dx2, d2alpha02_dy2, d2alpha02_dz2, d2alpha02_dxy, d2alpha02_dxz, d2alpha02_dyz;
  CCTK_REAL d2gphiphi_dx2, d2gphiphi_dy2, d2gphiphi_dz2, d2gphiphi_dxy, d2gphiphi_dxz, d2gphiphi_dyz;
  CCTK_REAL d2bphiup_dx2, d2bphiup_dy2, d2bphiup_dz2, d2bphiup_dxy, d2bphiup_dxz, d2bphiup_dyz;
  CCTK_REAL d2bphi_dx2, d2bphi_dy2, d2bphi_dz2, d2bphi_dxy, d2bphi_dxz, d2bphi_dyz;

} UAv_EvalPoint;

// Outputs
typedef struct {
  CCTK_REAL dG[4][4][4];
} UAv_MetricDerivs1;
typedef struct {
  CCTK_REAL ddG[4][4][4][4];
} UAv_MetricDerivs2;

// Zero-initializers
void UAv_InitMetricDerivs1(UAv_MetricDerivs1 *D1);
void UAv_InitMetricDerivs2(UAv_MetricDerivs2 *D2);

// Stubs you will fill later (currently zero)
void UAv_ComputeMetricDerivsAtPoint(const UAv_EvalPoint *P,
                                    UAv_MetricDerivs1 *D1,
                                    UAv_MetricDerivs2 *D2);