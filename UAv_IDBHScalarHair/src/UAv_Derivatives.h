#pragma once
#include "cctk.h"

// Point where you evaluate (fill these from Kerr_BS.c)
typedef struct {
  CCTK_REAL bh_mass;
  CCTK_REAL bh_spin;
  CCTK_REAL bh_spin2;
  CCTK_REAL bh_v;
  CCTK_REAL x1_2, y1_2, z1_2; //começo a achar que nao preciso de as escrever aqui. basta no codigo principal. as derivadas e outra conversa. so que tenho de por as varaiveis auxiliares como .psi4_2 etc.
  CCTK_REAL rr_2, rr2_2, rho_2, rho2_2;
  CCTK_REAL R_x, R_y, R_z;
  CCTK_REAL th_x, th_y, th_z;

  CCTK_REAL rBL, rBLp, rBLm, drBLdR;
  CCTK_REAL Delt, dDelt_dR;
  CCTK_REAL Sigm, Sigm2, dSigm_dR, dSigm_dth;
  CCTK_REAL fctFF, dfctFF_dR, dfctFF_dth;
  CCTK_REAL fctGG, dfctGG_dR;
  CCTK_REAL fctHH, dfctHH_dR, dfctHH_dth;
  CCTK_REAL psi4_2, dpsi4_2_dR, dpsi4_2_dth;
  CCTK_REAL alpha0, alpha02, dalpha02_dR, dalpha02_dth;
  CCTK_REAL gphiphi, dgphiphi_dR, dgphiphi_dth;
  CCTK_REAL bphiup, dbphiup_dR, dbphiup_dth;
  CCTK_REAL bphi, dbphi_dR, dbphi_dth;

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