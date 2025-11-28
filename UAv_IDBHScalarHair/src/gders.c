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

  // second derivatives
    D2->ddG[0][0][1][1] =;
    D2->ddG[0][0][1][2] =;
    D2->ddG[0][0][1][3] =;

    D2->ddG[0][1][1][1] =;
    D2->ddG[0][1][1][2] =;
    D2->ddG[0][1][1][3] =;

    D2->ddG[0][1][2][1] =;
    D2->ddG[0][1][2][2] =;
    D2->ddG[0][1][2][3] =;

    D2->ddG[0][1][3][1] =;
    D2->ddG[0][1][3][2] =;
    D2->ddG[0][1][3][3] =;

    D2->ddG[0][2][2][1] =;
    D2->ddG[0][2][2][2] =;
    D2->ddG[0][2][2][3] =;

    D2->ddG[0][2][3][1] =;
    D2->ddG[0][2][3][2] =;
    D2->ddG[0][2][3][3] =;

    D2->ddG[0][3][3][1] =;
    D2->ddG[0][3][3][2] =;
    D2->ddG[0][3][3][3] =;

    D2->ddG[1][1][1][1] =;
    D2->ddG[1][1][1][2] =;
    D2->ddG[1][1][1][3] =;

    D2->ddG[1][1][2][1] =;
    D2->ddG[1][1][2][2] =;
    D2->ddG[1][1][2][3] =;

    D2->ddG[1][1][3][1] =;
    D2->ddG[1][1][3][2] =;
    D2->ddG[1][1][3][3] =;

    D2->ddG[1][2][2][1] =;
    D2->ddG[1][2][2][2] =;
    D2->ddG[1][2][2][3] =;

    D2->ddG[1][2][3][1] =;
    D2->ddG[1][2][3][2] =;
    D2->ddG[1][2][3][3] =;

    D2->ddG[1][3][3][1] =;
    D2->ddG[1][3][3][2] =;
    D2->ddG[1][3][3][3] =;

    D2->ddG[2][2][2][1] =;
    D2->ddG[2][2][2][2] =;
    D2->ddG[2][2][2][3] =;

    D2->ddG[2][2][3][1] =;
    D2->ddG[2][2][3][2] =;
    D2->ddG[2][2][3][3] =;

    D2->ddG[2][3][3][1] =;
    D2->ddG[2][3][3][2] =;
    D2->ddG[2][3][3][3] =;

    D2->ddG[3][3][3][1] =;
    D2->ddG[3][3][3][2] =;
    D2->ddG[3][3][3][3] =;

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
