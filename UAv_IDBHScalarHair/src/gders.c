#include "UAv_Derivatives.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void UAv_InitMetricDerivs1(UAv_MetricDerivs1 *D1) {
  if (D1)
    memset(D1, 0, sizeof(*D1));
}
void UAv_InitMetricDerivs2(UAv_MetricDerivs2 *D2) {
  if (D2)
    memset(D2, 0, sizeof(*D2));
}

void UAv_ComputeMetricDerivsAtPoint(const UAv_EvalPoint *P,
                                    UAv_MetricDerivs1 *D1,
                                    UAv_MetricDerivs2 *D2) {
  if (D1)
    UAv_InitMetricDerivs1(D1);
  if (D2)
    UAv_InitMetricDerivs2(D2);
  if (!P || !D1)
    return;

  D1->dG[0][0][1] =
      -(P->dalpha02_dR * P->R_x + P->dalpha02_dth * P->th_x) +
      P->bphiup * (P->dbphi_dR * P->R_x + P->dbphi_dth * P->th_x) +
      P->bphi * (P->dbphiup_dR * P->R_x + P->dbphiup_dth * P->th_x);

  D1->dG[0][0][2] =
      -(P->dalpha02_dR * P->R_y + P->dalpha02_dth * P->th_y) +
      P->bphiup * (P->dbphi_dR * P->R_y + P->dbphi_dth * P->th_y) +
      P->bphi * (P->dbphiup_dR * P->R_y + P->dbphiup_dth * P->th_y);

  D1->dG[0][0][3] =
      -(P->dalpha02_dR * P->R_z + P->dalpha02_dth * P->th_z) +
      P->bphiup * (P->dbphi_dR * P->R_z + P->dbphi_dth * P->th_z) +
      P->bphi * (P->dbphiup_dR * P->R_z + P->dbphiup_dth * P->th_z);

  D1->dG[0][1][1] =
      (P->y1_2 *
       (P->bphi * 2.0 * P->x1_2 -
        P->rho2_2 * (P->dbphi_dR * P->R_x + P->dbphi_dth * P->th_x))) /
      pow(P->rho2_2, 2);
  D1->dG[0][1][2] =
      (P->bphi * (-P->rho2_2 + P->y1_2 * 2.0 * P->y1_2) -
       P->y1_2 * P->rho2_2 * (P->dbphi_dR * P->R_y + P->dbphi_dth * P->th_y)) /
      pow(P->rho2_2, 2);
  D1->dG[0][1][3] = -(
      (P->y1_2 * (P->dbphi_dR * P->R_z + P->dbphi_dth * P->th_z)) / P->rho2_2);

  D1->dG[0][2][1] =
      (P->bphi * (P->rho2_2 - P->x1_2 * 2.0 * P->x1_2) +
       P->x1_2 * P->rho2_2 * (P->dbphi_dR * P->R_x + P->dbphi_dth * P->th_x)) /
      pow(P->rho2_2, 2);
  D1->dG[0][2][2] =
      (P->x1_2 *
       (-(P->bphi * 2.0 * P->y1_2) +
        P->rho2_2 * (P->dbphi_dR * P->R_y + P->dbphi_dth * P->th_y))) /
      pow(P->rho2_2, 2);
  D1->dG[0][2][3] =
      (P->x1_2 * (P->dbphi_dR * P->R_z + P->dbphi_dth * P->th_z)) / P->rho2_2;

  D1->dG[0][3][1] = 0;
  D1->dG[0][3][2] = 0;
  D1->dG[0][3][3] = 0;

  D1->dG[1][1][1] =
      P->psi4_2 *
          (2 * P->x1_2 * P->fctGG + pow(P->x1_2, 2) * (P->dfctGG_dR * P->R_x) +
           P->bh_spin2 * pow(P->y1_2, 2) *
               (P->dfctHH_dR * P->R_x + P->dfctHH_dth * P->th_x)) +
      (1 + pow(P->x1_2, 2) * P->fctGG +
       P->bh_spin2 * pow(P->y1_2, 2) * P->fctHH) *
          (P->dpsi4_2_dR * P->R_x + P->dpsi4_2_dth * P->th_x);
  D1->dG[1][1][2] =
      P->psi4_2 * (pow(P->x1_2, 2) * (P->dfctGG_dR * P->R_y) +
                   P->bh_spin2 * P->y1_2 *
                       (2 * P->fctHH + P->y1_2 * (P->dfctHH_dR * P->R_y +
                                                  P->dfctHH_dth * P->th_y))) +
      (1 + pow(P->x1_2, 2) * P->fctGG +
       P->bh_spin2 * pow(P->y1_2, 2) * P->fctHH) *
          (P->dpsi4_2_dR * P->R_y + P->dpsi4_2_dth * P->th_y);
  D1->dG[1][1][3] =
      P->psi4_2 * (pow(P->x1_2, 2) * (P->dfctGG_dR * P->R_z) +
                   P->bh_spin2 * pow(P->y1_2, 2) *
                       (P->dfctHH_dR * P->R_z + P->dfctHH_dth * P->th_z)) +
      (1 + pow(P->x1_2, 2) * P->fctGG +
       P->bh_spin2 * pow(P->y1_2, 2) * P->fctHH) *
          (P->dpsi4_2_dR * P->R_z + P->dpsi4_2_dth * P->th_z);

  D1->dG[1][2][1] =
      P->y1_2 * P->psi4_2 *
          (P->fctGG + P->x1_2 * (P->dfctGG_dR * P->R_x) -
           P->bh_spin2 * (P->fctHH + P->x1_2 * (P->dfctHH_dR * P->R_x +
                                                P->dfctHH_dth * P->th_x))) +
      P->x1_2 * P->y1_2 * (P->fctGG - P->bh_spin2 * P->fctHH) *
          (P->dpsi4_2_dR * P->R_x + P->dpsi4_2_dth * P->th_x);
  D1->dG[1][2][2] =
      P->x1_2 * P->psi4_2 *
          (P->fctGG + P->y1_2 * (P->dfctGG_dR * P->R_y) -
           P->bh_spin2 * (P->fctHH + P->y1_2 * (P->dfctHH_dR * P->R_y +
                                                P->dfctHH_dth * P->th_y))) +
      P->x1_2 * P->y1_2 * (P->fctGG - P->bh_spin2 * P->fctHH) *
          (P->dpsi4_2_dR * P->R_y + P->dpsi4_2_dth * P->th_y);
  D1->dG[1][2][3] =
      P->x1_2 * P->y1_2 *
      (P->psi4_2 *
           ((P->dfctGG_dR * P->R_z) -
            P->bh_spin2 * (P->dfctHH_dR * P->R_z + P->dfctHH_dth * P->th_z)) +
       (P->fctGG - P->bh_spin2 * P->fctHH) *
           (P->dpsi4_2_dR * P->R_z + P->dpsi4_2_dth * P->th_z));

  D1->dG[1][3][1] =
      P->z1_2 * (P->x1_2 * P->psi4_2 * (P->dfctGG_dR * P->R_x) +
                 P->fctGG * (P->psi4_2 + P->x1_2 * (P->dpsi4_2_dR * P->R_x +
                                                    P->dpsi4_2_dth * P->th_x)));
  D1->dG[1][3][2] =
      P->x1_2 * P->z1_2 *
      (P->psi4_2 * (P->dfctGG_dR * P->R_y) +
       P->fctGG * (P->dpsi4_2_dR * P->R_y + P->dpsi4_2_dth * P->th_y));
  D1->dG[1][3][3] =
      P->x1_2 * (P->z1_2 * P->psi4_2 * (P->dfctGG_dR * P->R_z) +
                 P->fctGG * (P->psi4_2 + P->z1_2 * (P->dpsi4_2_dR * P->R_z +
                                                    P->dpsi4_2_dth * P->th_z)));

  D1->dG[2][2][1] =
      P->psi4_2 * (pow(P->y1_2, 2) * (P->dfctGG_dR * P->R_x) +
                   P->bh_spin2 * P->x1_2 *
                       (2 * P->fctHH + P->x1_2 * (P->dfctHH_dR * P->R_x +
                                                  P->dfctHH_dth * P->th_x))) +
      (1 + pow(P->y1_2, 2) * P->fctGG +
       P->bh_spin2 * pow(P->x1_2, 2) * P->fctHH) *
          (P->dpsi4_2_dR * P->R_x + P->dpsi4_2_dth * P->th_x);
  D1->dG[2][2][2] =
      P->psi4_2 *
          (2 * P->y1_2 * P->fctGG + pow(P->y1_2, 2) * (P->dfctGG_dR * P->R_y) +
           P->bh_spin2 * pow(P->x1_2, 2) *
               (P->dfctHH_dR * P->R_y + P->dfctHH_dth * P->th_y)) +
      (1 + pow(P->y1_2, 2) * P->fctGG +
       P->bh_spin2 * pow(P->x1_2, 2) * P->fctHH) *
          (P->dpsi4_2_dR * P->R_y + P->dpsi4_2_dth * P->th_y);
  D1->dG[2][2][3] =
      P->psi4_2 * (pow(P->y1_2, 2) * (P->dfctGG_dR * P->R_z) +
                   P->bh_spin2 * pow(P->x1_2, 2) *
                       (P->dfctHH_dR * P->R_z + P->dfctHH_dth * P->th_z)) +
      (1 + pow(P->y1_2, 2) * P->fctGG +
       P->bh_spin2 * pow(P->x1_2, 2) * P->fctHH) *
          (P->dpsi4_2_dR * P->R_z + P->dpsi4_2_dth * P->th_z);

  D1->dG[2][3][1] =
      P->y1_2 * P->z1_2 *
      (P->psi4_2 * (P->dfctGG_dR * P->R_x) +
       P->fctGG * (P->dpsi4_2_dR * P->R_x + P->dpsi4_2_dth * P->th_x));
  D1->dG[2][3][2] =
      P->z1_2 * (P->y1_2 * P->psi4_2 * (P->dfctGG_dR * P->R_y) +
                 P->fctGG * (P->psi4_2 + P->y1_2 * (P->dpsi4_2_dR * P->R_y +
                                                    P->dpsi4_2_dth * P->th_y)));
  D1->dG[2][3][3] =
      P->y1_2 * (P->z1_2 * P->psi4_2 * (P->dfctGG_dR * P->R_z) +
                 P->fctGG * (P->psi4_2 + P->z1_2 * (P->dpsi4_2_dR * P->R_z +
                                                    P->dpsi4_2_dth * P->th_z)));

  D1->dG[3][3][1] = pow(P->z1_2, 2) * P->psi4_2 * (P->dfctGG_dR * P->R_x) +
                    (1 + pow(P->z1_2, 2) * P->fctGG) *
                        (P->dpsi4_2_dR * P->R_x + P->dpsi4_2_dth * P->th_x);
  D1->dG[3][3][2] = pow(P->z1_2, 2) * P->psi4_2 * (P->dfctGG_dR * P->R_y) +
                    (1 + pow(P->z1_2, 2) * P->fctGG) *
                        (P->dpsi4_2_dR * P->R_y + P->dpsi4_2_dth * P->th_y);
  D1->dG[3][3][3] =
      P->z1_2 * P->psi4_2 * (2 * P->fctGG + P->z1_2 * (P->dfctGG_dR * P->R_z)) +
      (1 + pow(P->z1_2, 2) * P->fctGG) *
          (P->dpsi4_2_dR * P->R_z + P->dpsi4_2_dth * P->th_z);

  // symmetries
  D1->dG[2][1][1] = D1->dG[1][2][1];
  D1->dG[2][1][2] = D1->dG[1][2][2];
  D1->dG[2][1][3] = D1->dG[1][2][3];
  D1->dG[3][1][1] = D1->dG[1][3][1];
  D1->dG[3][1][2] = D1->dG[1][3][2];
  D1->dG[3][1][3] = D1->dG[1][3][3];
  D1->dG[3][2][1] = D1->dG[2][3][1];
  D1->dG[3][2][2] = D1->dG[2][3][2];
  D1->dG[3][2][3] = D1->dG[2][3][3];
  D1->dG[1][0][1] = D1->dG[0][1][1];
  D1->dG[1][0][2] = D1->dG[0][1][2];
  D1->dG[1][0][3] = D1->dG[0][1][3];
  D1->dG[2][0][1] = D1->dG[0][2][1];
  D1->dG[2][0][2] = D1->dG[0][2][2];
  D1->dG[2][0][3] = D1->dG[0][2][3];
  D1->dG[3][0][1] = D1->dG[0][3][1];
  D1->dG[3][0][2] = D1->dG[0][3][2];
  D1->dG[3][0][3] = D1->dG[0][3][3];

  // second derivatives (retirar os gammas)
    D2->ddG[0][0][1][1] = 2*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x)*(P->dbphiup_dR*P->R_x+P->dbphiup_dth*P->th_x) - P->d2alpha02_dx2 + P->bphiup*P->d2bphi_dx2 + P->bphi*P->d2bphiup_dx2;
    D2->ddG[0][0][1][2] = (P->dbphiup_dR*P->R_y+P->dbphiup_dth*P->th_y)*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x) + (P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y)*(P->dbphiup_dR*P->R_x+P->dbphiup_dth*P->th_x) - P->d2alpha02_dxy + P->bphiup*P->d2bphi_dxy + P->bphi*P->d2bphiup_dxy;
    D2->ddG[0][0][1][3] =(P->dbphiup_dR*P->R_z+P->dbphiup_dth*P->th_z)*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x) + (P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z)*(P->dbphiup_dR*P->R_x+P->dbphiup_dth*P->th_x) - P->d2alpha02_dxz + P->bphiup*P->d2bphi_dxz + P->bphi*P->d2bphiup_dxz;

    D2->ddG[0][1][1][1] = (P->y1_2*(P->bphi*(-2*pow(2.0*P->x1_2,2) + P->rho2_2*2) + P->rho2_2*(2*2.0*P->x1_2*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x) - P->rho2_2*P->d2bphi_dx2)))/pow(P->rho2_2,3);
    D2->ddG[0][1][1][2] = (P->bphi*((P->rho2_2 - 2*P->y1_2*2.0*P->y1_2)*2.0*P->x1_2 + P->y1_2*P->rho2_2*0.0) + P->rho2_2*(P->y1_2*2.0*P->x1_2*(P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y) + (-P->rho2_2 + P->y1_2*2.0*P->y1_2)*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x) - P->y1_2*P->rho2_2*P->d2bphi_dxy))/pow(P->rho2_2,3);
    D2->ddG[0][1][1][3] = (P->y1_2*(2.0*P->x1_2*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z) - P->rho2_2*P->d2bphi_dxz))/pow(P->rho2_2,2);

    D2->ddG[0][1][2][1] =(P->bphi*((P->rho2_2 - 2*P->y1_2*2.0*P->y1_2)*2.0*P->x1_2 + P->y1_2*P->rho2_2*0.0) + P->rho2_2*(P->y1_2*2.0*P->x1_2*(P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y) + (-P->rho2_2 + P->y1_2*2.0*P->y1_2)*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x) - P->y1_2*P->rho2_2*P->d2bphi_dxy))/pow(P->rho2_2,3);
    D2->ddG[0][1][2][2] =(P->bphi*(-2*P->y1_2*pow(2.0*P->y1_2,2) + P->rho2_2*(2*2.0*P->y1_2 + P->y1_2*2)) + P->rho2_2*(-2*(P->rho2_2 - P->y1_2*2.0*P->y1_2)*(P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y) - P->y1_2*P->rho2_2*P->d2bphi_dy2))/pow(P->rho2_2,3);
    D2->ddG[0][1][2][3] = ((-P->rho2_2 + P->y1_2*2.0*P->y1_2)*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z) - P->y1_2*P->rho2_2*P->d2bphi_dyz)/pow(P->rho2_2,2);

    D2->ddG[0][1][3][1] =(P->y1_2*(2.0*P->x1_2*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z) - P->rho2_2*P->d2bphi_dxz))/pow(P->rho2_2,2);
    D2->ddG[0][1][3][2] =((-P->rho2_2 + P->y1_2*2.0*P->y1_2)*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z) - P->y1_2*P->rho2_2*P->d2bphi_dyz)/pow(P->rho2_2,2);
    D2->ddG[0][1][3][3] =-((P->y1_2*P->d2bphi_dz2)/P->rho2_2);

    D2->ddG[0][2][2][1] =(-(P->bphi*(2.0*P->y1_2*(P->rho2_2 - 2*P->x1_2*2.0*P->x1_2) + P->x1_2*P->rho2_2*0.0)) + P->rho2_2*(-(P->x1_2*(2.0*P->x1_2*(P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y) + 2.0*P->y1_2*(P->dbphi_dR*P->R_x+P->dbphi_dth*P->th_x))) + P->rho2_2*((P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y) + P->x1_2*P->d2bphi_dxy)))/pow(P->rho2_2,3);
    D2->ddG[0][2][2][2] =(P->x1_2*(P->bphi*(2*pow(2.0*P->y1_2,2) - P->rho2_2*2) + P->rho2_2*(-2*2.0*P->y1_2*(P->dbphi_dR*P->R_y+P->dbphi_dth*P->th_y) + P->rho2_2*P->d2bphi_dy2)))/pow(P->rho2_2,3);
    D2->ddG[0][2][2][3] =(P->x1_2*(-(2.0*P->y1_2*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z)) + P->rho2_2*P->d2bphi_dyz))/pow(P->rho2_2,2);

    D2->ddG[0][2][3][1] =((P->rho2_2 - P->x1_2*2.0*P->x1_2)*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z) + P->x1_2*P->rho2_2*P->d2bphi_dxz)/pow(P->rho2_2,2);
    D2->ddG[0][2][3][2] =(P->x1_2*(-(2.0*P->y1_2*(P->dbphi_dR*P->R_z+P->dbphi_dth*P->th_z)) + P->rho2_2*P->d2bphi_dyz))/pow(P->rho2_2,2);
    D2->ddG[0][2][3][3] =(P->x1_2*P->d2bphi_dz2)/P->rho2_2;

    D2->ddG[0][3][3][1] = 0.0;
    D2->ddG[0][3][3][2] = 0.0;
    D2->ddG[0][3][3][3] = 0.0;

    D2->ddG[1][1][1][1] = 2*(2*P->x1_2*P->fctGG + pow(P->x1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(2*P->fctGG + 4*P->x1_2*(P->dfctGG_dR*P->R_x) + pow(P->x1_2,2)*P->d2fctGG_dx2 + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->d2fctHH_dx2) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dx2;
    D2->ddG[1][1][1][2] = (P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y)*(2*P->x1_2*P->fctGG + pow(P->x1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x)) + (pow(P->x1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*P->y1_2*(2*P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(2*P->x1_2*(P->dfctGG_dR*P->R_y) + pow(P->x1_2,2)*P->d2fctGG_dxy + pow(P->bh_spin,2)*P->y1_2*(2*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x) + P->y1_2*P->d2fctHH_dxy)) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dxy;
    D2->ddG[1][1][1][3] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(2*P->x1_2*P->fctGG + pow(P->x1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x)) + (pow(P->x1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(2*P->x1_2*(P->dfctGG_dR*P->R_z) + pow(P->x1_2,2)*P->d2fctGG_dxz + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->d2fctHH_dxz) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dxz;

    D2->ddG[1][1][2][1] = (P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y)*(2*P->x1_2*P->fctGG + pow(P->x1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x)) + (pow(P->x1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*P->y1_2*(2*P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(2*P->x1_2*(P->dfctGG_dR*P->R_y) + pow(P->x1_2,2)*P->d2fctGG_dxy + pow(P->bh_spin,2)*P->y1_2*(2*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x) + P->y1_2*P->d2fctHH_dxy)) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dxy;
    D2->ddG[1][1][2][2] =2*(pow(P->x1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(pow(P->x1_2,2)*P->d2fctGG_dy2 + pow(P->bh_spin,2)*P->y1_2*(4*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y) + P->y1_2*P->d2fctHH_dy2)) + (1 + pow(P->x1_2,2)*P->fctGG)*P->d2psi4_2_dy2 + pow(P->bh_spin,2)*P->fctHH*(2*P->psi4_2 + P->y1_2*(4*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->y1_2*P->d2psi4_2_dy2));
    D2->ddG[1][1][2][3] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(pow(P->x1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*P->y1_2*(2*P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))) + (pow(P->x1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(pow(P->x1_2,2)*P->d2fctGG_dyz + pow(P->bh_spin,2)*P->y1_2*(2*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z) + P->y1_2*P->d2fctHH_dyz)) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dyz;

    D2->ddG[1][1][3][1] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(2*P->x1_2*P->fctGG + pow(P->x1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x)) + (pow(P->x1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(2*P->x1_2*(P->dfctGG_dR*P->R_z) + pow(P->x1_2,2)*P->d2fctGG_dxz + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->d2fctHH_dxz) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dxz;
    D2->ddG[1][1][3][2] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(pow(P->x1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*P->y1_2*(2*P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))) + (pow(P->x1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(pow(P->x1_2,2)*P->d2fctGG_dyz + pow(P->bh_spin,2)*P->y1_2*(2*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z) + P->y1_2*P->d2fctHH_dyz)) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dyz;
    D2->ddG[1][1][3][3] = 2*(pow(P->x1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->y1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->psi4_2*(pow(P->x1_2,2)*P->d2fctGG_dz2 + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->d2fctHH_dz2) + (1 + pow(P->x1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->y1_2,2)*P->fctHH)*P->d2psi4_2_dz2;

    D2->ddG[1][2][2][1] = P->psi4_2*(P->fctGG + P->y1_2*(P->dfctGG_dR*P->R_y) - pow(P->bh_spin,2)*(P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))) + P->y1_2*(P->fctGG - pow(P->bh_spin,2)*P->fctHH)*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->x1_2*P->y1_2*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y)*((P->dfctGG_dR*P->R_x) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x)) + P->x1_2*(P->fctGG + P->y1_2*(P->dfctGG_dR*P->R_y) - pow(P->bh_spin,2)*(P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->x1_2*P->psi4_2*((P->dfctGG_dR*P->R_x) + P->y1_2*P->d2fctGG_dxy - pow(P->bh_spin,2)*((P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x) + P->y1_2*P->d2fctHH_dxy)) + P->x1_2*P->y1_2*(P->fctGG - pow(P->bh_spin,2)*P->fctHH)*P->d2psi4_2_dxy;
    D2->ddG[1][2][2][2] = P->x1_2*(2*(P->fctGG + P->y1_2*(P->dfctGG_dR*P->R_y) - pow(P->bh_spin,2)*(P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(2*(P->dfctGG_dR*P->R_y) + P->y1_2*P->d2fctGG_dy2 - pow(P->bh_spin,2)*(2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y) + P->y1_2*P->d2fctHH_dy2)) + P->y1_2*(P->fctGG - pow(P->bh_spin,2)*P->fctHH)*P->d2psi4_2_dy2);
    D2->ddG[1][2][2][3] = P->x1_2*((P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->fctGG + P->y1_2*(P->dfctGG_dR*P->R_y) - pow(P->bh_spin,2)*(P->fctHH + P->y1_2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))) + P->y1_2*((P->dfctGG_dR*P->R_z) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*((P->dfctGG_dR*P->R_z) + P->y1_2*P->d2fctGG_dyz - pow(P->bh_spin,2)*((P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z) + P->y1_2*P->d2fctHH_dyz)) + P->y1_2*(P->fctGG - pow(P->bh_spin,2)*P->fctHH)*P->d2psi4_2_dyz);

    D2->ddG[1][2][3][1] = P->y1_2*(P->psi4_2*((P->dfctGG_dR*P->R_z) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z)) + (P->fctGG - pow(P->bh_spin,2)*P->fctHH)*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->x1_2*((P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*((P->dfctGG_dR*P->R_x) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x)) + ((P->dfctGG_dR*P->R_z) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(P->d2fctGG_dxz - pow(P->bh_spin,2)*P->d2fctHH_dxz) + (P->fctGG - pow(P->bh_spin,2)*P->fctHH)*P->d2psi4_2_dxz));
    D2->ddG[1][2][3][2] = P->x1_2*(P->psi4_2*((P->dfctGG_dR*P->R_z) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z)) + (P->fctGG - pow(P->bh_spin,2)*P->fctHH)*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->y1_2*((P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*((P->dfctGG_dR*P->R_y) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)) + ((P->dfctGG_dR*P->R_z) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(P->d2fctGG_dyz - pow(P->bh_spin,2)*P->d2fctHH_dyz) + (P->fctGG - pow(P->bh_spin,2)*P->fctHH)*P->d2psi4_2_dyz));
    D2->ddG[1][2][3][3] = P->x1_2*P->y1_2*(2*((P->dfctGG_dR*P->R_z) - pow(P->bh_spin,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->psi4_2*(P->d2fctGG_dz2 - pow(P->bh_spin,2)*P->d2fctHH_dz2) + (P->fctGG - pow(P->bh_spin,2)*P->fctHH)*P->d2psi4_2_dz2);

    D2->ddG[1][3][3][1] =P->x1_2*P->z1_2*((P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->dfctGG_dR*P->R_x) + (P->dfctGG_dR*P->R_z)*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x)) + P->psi4_2*(P->z1_2*(P->dfctGG_dR*P->R_z) + P->x1_2*(P->dfctGG_dR*P->R_x) + P->x1_2*P->z1_2*P->d2fctGG_dxz) + P->fctGG*(P->psi4_2 + P->z1_2*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->x1_2*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->x1_2*P->z1_2*P->d2psi4_2_dxz);
    D2->ddG[1][3][3][2] =P->x1_2*(P->z1_2*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->dfctGG_dR*P->R_y) + (P->fctGG + P->z1_2*(P->dfctGG_dR*P->R_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*((P->dfctGG_dR*P->R_y) + P->z1_2*P->d2fctGG_dyz) + P->z1_2*P->fctGG*P->d2psi4_2_dyz);
    D2->ddG[1][3][3][3] =P->x1_2*(2*(P->fctGG + P->z1_2*(P->dfctGG_dR*P->R_z))*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->psi4_2*(2*(P->dfctGG_dR*P->R_z) + P->z1_2*P->d2fctGG_dz2) + P->z1_2*P->fctGG*P->d2psi4_2_dz2);

    D2->ddG[2][2][2][1] = (P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y)*(pow(P->y1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*P->x1_2*(2*P->fctHH + P->x1_2*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x))) + (2*P->y1_2*P->fctGG + pow(P->y1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(P->y1_2*(2*(P->dfctGG_dR*P->R_x) + P->y1_2*P->d2fctGG_dxy) + pow(P->bh_spin,2)*P->x1_2*(2*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y) + P->x1_2*P->d2fctHH_dxy)) + (1 + pow(P->y1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->fctHH)*P->d2psi4_2_dxy;
    D2->ddG[2][2][2][2] = 2*(2*P->y1_2*P->fctGG + pow(P->y1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(2*P->fctGG + 4*P->y1_2*(P->dfctGG_dR*P->R_y) + pow(P->y1_2,2)*P->d2fctGG_dy2 + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->d2fctHH_dy2) + (1 + pow(P->y1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->fctHH)*P->d2psi4_2_dy2;
    D2->ddG[2][2][2][3] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(2*P->y1_2*P->fctGG + pow(P->y1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)) + (pow(P->y1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(2*P->y1_2*(P->dfctGG_dR*P->R_z) + pow(P->y1_2,2)*P->d2fctGG_dyz + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->d2fctHH_dyz) + (1 + pow(P->y1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->fctHH)*P->d2psi4_2_dyz;

    D2->ddG[2][2][3][1] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(pow(P->y1_2,2)*(P->dfctGG_dR*P->R_x) + pow(P->bh_spin,2)*P->x1_2*(2*P->fctHH + P->x1_2*(P->dfctHH_dR*P->R_x+P->dfctHH_dth*P->th_x))) + (pow(P->y1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*(pow(P->y1_2,2)*P->d2fctGG_dxz + pow(P->bh_spin,2)*P->x1_2*(2*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z) + P->x1_2*P->d2fctHH_dxz)) + (1 + pow(P->y1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->fctHH)*P->d2psi4_2_dxz;
    D2->ddG[2][2][3][2] = (P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(2*P->y1_2*P->fctGG + pow(P->y1_2,2)*(P->dfctGG_dR*P->R_y) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_y+P->dfctHH_dth*P->th_y)) + (pow(P->y1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->psi4_2*(2*P->y1_2*(P->dfctGG_dR*P->R_z) + pow(P->y1_2,2)*P->d2fctGG_dyz + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->d2fctHH_dyz) + (1 + pow(P->y1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->fctHH)*P->d2psi4_2_dyz;
    D2->ddG[2][2][3][3] = 2*(pow(P->y1_2,2)*(P->dfctGG_dR*P->R_z) + pow(P->bh_spin,2)*pow(P->x1_2,2)*(P->dfctHH_dR*P->R_z+P->dfctHH_dth*P->th_z))*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->psi4_2*(pow(P->y1_2,2)*P->d2fctGG_dz2 + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->d2fctHH_dz2) + (1 + pow(P->y1_2,2)*P->fctGG + pow(P->bh_spin,2)*pow(P->x1_2,2)*P->fctHH)*P->d2psi4_2_dz2;

    D2->ddG[2][3][3][1] = P->y1_2*(P->z1_2*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->dfctGG_dR*P->R_x) + (P->fctGG + P->z1_2*(P->dfctGG_dR*P->R_z))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->psi4_2*((P->dfctGG_dR*P->R_x) + P->z1_2*P->d2fctGG_dxz) + P->z1_2*P->fctGG*P->d2psi4_2_dxz);
    D2->ddG[2][3][3][2] = P->y1_2*P->z1_2*((P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->dfctGG_dR*P->R_y) + (P->dfctGG_dR*P->R_z)*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y)) + P->psi4_2*(P->z1_2*(P->dfctGG_dR*P->R_z) + P->y1_2*(P->dfctGG_dR*P->R_y) + P->y1_2*P->z1_2*P->d2fctGG_dyz) + P->fctGG*(P->psi4_2 + P->z1_2*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->y1_2*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->y1_2*P->z1_2*P->d2psi4_2_dyz);
    D2->ddG[2][3][3][3] =P->y1_2*(2*(P->fctGG + P->z1_2*(P->dfctGG_dR*P->R_z))*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->psi4_2*(2*(P->dfctGG_dR*P->R_z) + P->z1_2*P->d2fctGG_dz2) + P->z1_2*P->fctGG*P->d2psi4_2_dz2);

    D2->ddG[3][3][3][1] = pow(P->z1_2,2)*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->dfctGG_dR*P->R_x) + P->z1_2*(2*P->fctGG + P->z1_2*(P->dfctGG_dR*P->R_z))*(P->dpsi4_2_dR*P->R_x+P->dpsi4_2_dth*P->th_x) + P->z1_2*P->psi4_2*(2*(P->dfctGG_dR*P->R_x) + P->z1_2*P->d2fctGG_dxz) + (1 + pow(P->z1_2,2)*P->fctGG)*P->d2psi4_2_dxz;
    D2->ddG[3][3][3][2] = pow(P->z1_2,2)*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)*(P->dfctGG_dR*P->R_y) + P->z1_2*(2*P->fctGG + P->z1_2*(P->dfctGG_dR*P->R_z))*(P->dpsi4_2_dR*P->R_y+P->dpsi4_2_dth*P->th_y) + P->z1_2*P->psi4_2*(2*(P->dfctGG_dR*P->R_y) + P->z1_2*P->d2fctGG_dyz) + (1 + pow(P->z1_2,2)*P->fctGG)*P->d2psi4_2_dyz;
    D2->ddG[3][3][3][3] = P->z1_2*(2*(P->dfctGG_dR*P->R_z)*(2*P->psi4_2 + P->z1_2*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z)) + P->z1_2*P->psi4_2*P->d2fctGG_dz2) + P->d2psi4_2_dz2 + P->fctGG*(2*P->psi4_2 + P->z1_2*(4*(P->dpsi4_2_dR*P->R_z+P->dpsi4_2_dth*P->th_z) + P->z1_2*P->d2psi4_2_dz2));

  // symmetries of second derivatives to be implemented
    D2->ddG[1][0][1][1] =D2->ddG[0][1][1][1];
    D2->ddG[1][0][1][2] =D2->ddG[0][1][1][2];
    D2->ddG[1][0][1][3] =D2->ddG[0][1][1][3];
    D2->ddG[1][0][2][1] =D2->ddG[0][2][1][1];
    D2->ddG[1][0][2][2] =D2->ddG[0][2][1][2];
    D2->ddG[1][0][2][3] =D2->ddG[0][2][1][3];
    D2->ddG[1][0][3][1] =D2->ddG[0][3][1][1];
    D2->ddG[1][0][3][2] =D2->ddG[0][3][1][2];
    D2->ddG[1][0][3][3] =D2->ddG[0][3][1][3];
    D2->ddG[2][0][1][1] =D2->ddG[0][1][2][1];
    D2->ddG[2][0][1][2] =D2->ddG[0][1][2][2];
    D2->ddG[2][0][1][3] =D2->ddG[0][1][2][3];
    D2->ddG[2][0][2][1] =D2->ddG[0][2][2][1];
    D2->ddG[2][0][2][2] =D2->ddG[0][2][2][2];
    D2->ddG[2][0][2][3] =D2->ddG[0][2][2][3];
    D2->ddG[2][0][3][1] =D2->ddG[0][3][2][1];
    D2->ddG[2][0][3][2] =D2->ddG[0][3][2][2];
    D2->ddG[2][0][3][3] =D2->ddG[0][3][2][3];
    D2->ddG[3][0][1][1] =D2->ddG[0][1][3][1];
    D2->ddG[3][0][1][2] =D2->ddG[0][1][3][2];
    D2->ddG[3][0][1][3] =D2->ddG[0][1][3][3];
    D2->ddG[3][0][2][1] =D2->ddG[0][2][3][1];
    D2->ddG[3][0][2][2] =D2->ddG[0][2][3][2];
    D2->ddG[3][0][2][3] =D2->ddG[0][2][3][3];
    D2->ddG[3][0][3][1] =D2->ddG[0][3][3][1];
    D2->ddG[3][0][3][2] =D2->ddG[0][3][3][2];
    D2->ddG[3][0][3][3] =D2->ddG[0][3][3][3];
    D2->ddG[2][1][1][1] =D2->ddG[1][2][1][1];
    D2->ddG[2][1][1][2] =D2->ddG[1][2][1][2];
    D2->ddG[2][1][1][3] =D2->ddG[1][2][1][3];
    D2->ddG[2][1][2][1] =D2->ddG[1][2][2][1];
    D2->ddG[2][1][2][2] =D2->ddG[1][2][2][2];
    D2->ddG[2][1][2][3] =D2->ddG[1][2][2][3];
    D2->ddG[2][1][3][1] =D2->ddG[1][2][3][1];
    D2->ddG[2][1][3][2] =D2->ddG[1][2][3][2];
    D2->ddG[2][1][3][3] =D2->ddG[1][2][3][3];
    D2->ddG[3][1][1][1] =D2->ddG[1][3][1][1];
    D2->ddG[3][1][1][2] =D2->ddG[1][3][1][2];
    D2->ddG[3][1][1][3] =D2->ddG[1][3][1][3];
    D2->ddG[3][1][2][1] =D2->ddG[1][3][2][1];
    D2->ddG[3][1][2][2] =D2->ddG[1][3][2][2];
    D2->ddG[3][1][2][3] =D2->ddG[1][3][2][3];
    D2->ddG[3][1][3][1] =D2->ddG[1][3][3][1];
    D2->ddG[3][1][3][2] =D2->ddG[1][3][3][2];
    D2->ddG[3][1][3][3] =D2->ddG[1][3][3][3];
    D2->ddG[3][2][1][1] =D2->ddG[2][3][1][1];
    D2->ddG[3][2][1][2] =D2->ddG[2][3][1][2];
    D2->ddG[3][2][1][3] =D2->ddG[2][3][1][3];
    D2->ddG[3][2][2][1] =D2->ddG[2][3][2][1];
    D2->ddG[3][2][2][2] =D2->ddG[2][3][2][2];
    D2->ddG[3][2][2][3] =D2->ddG[2][3][2][3];
    D2->ddG[3][2][3][1] =D2->ddG[2][3][3][1];
    D2->ddG[3][2][3][2] =D2->ddG[2][3][3][2];
    D2->ddG[3][2][3][3] =D2->ddG[2][3][3][3];
}
