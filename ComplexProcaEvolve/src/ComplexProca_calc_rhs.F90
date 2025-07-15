
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine ComplexProca_calc_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                hh(3,3), hu(3,3), trk, dethh, ch
  CCTK_REAL                lE1(3), lA1(3), lAphi1, lZeta1
  CCTK_REAL                lE2(3), lA2(3), lAphi2, lZeta2

  ! First derivatives
  CCTK_REAL                d1_alph(3), d1_beta(3,3)
  CCTK_REAL                d1_hh(3,3,3), d1_ch(3)
  CCTK_REAL                d1_lE1(3,3), d1_lA1(3,3), d1_lZeta1(3), d1_lAphi1(3)
  CCTK_REAL                d1_lE2(3,3), d1_lA2(3,3), d1_lZeta2(3), d1_lAphi2(3)

  ! Second derivatives
  CCTK_REAL                d2_lA1(3,3,3)
  CCTK_REAL                d2_lA2(3,3,3)

  ! Advection derivatives
  CCTK_REAL                ad1_lE1(3), ad1_lA1(3), ad1_lZeta1, ad1_lAphi1
  CCTK_REAL                ad1_lE2(3), ad1_lA2(3), ad1_lZeta2, ad1_lAphi2
  CCTK_REAL                d1_f1(3)   ! Place holder for the advection derivs
  CCTK_REAL                d1_f2(3)

  ! Auxiliary variables
  CCTK_REAL                cf1(3,3,3), cf2(3,3,3)

  ! Covaraint derivatives
  CCTK_REAL                cd_lA1(3,3), cd_dA1(3,3,3)
  CCTK_REAL                cd_lA2(3,3), cd_dA2(3,3,3)

  ! Right hand sides
  CCTK_REAL                rhs_lE1(3), rhs_lA1(3), rhs_lZeta1, rhs_lAphi1
  CCTK_REAL                rhs_lE2(3), rhs_lA2(3), rhs_lZeta2, rhs_lAphi2

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12, dxsq12, dysq12, dzsq12,         &
                           dxdy144, dxdz144, dydz144
  CCTK_REAL                odx60, ody60, odz60, odxsq180, odysq180, odzsq180,&
                           odxdy3600, odxdz3600, odydz3600
  CCTK_INT                 i, j, k
  CCTK_INT                 di, dj, dk
  CCTK_REAL, parameter ::  one = 1
  CCTK_INT                 a, b, c, m, n


  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  dxsq12 = 12*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
  dysq12 = 12*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2)
  dzsq12 = 12*CCTK_DELTA_SPACE(3)*CCTK_DELTA_SPACE(3)

  dxdy144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2)
  dxdz144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3)
  dydz144 = 144*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3)


  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  odxsq180 = 1 / (180*CCTK_DELTA_SPACE(1)**2)
  odysq180 = 1 / (180*CCTK_DELTA_SPACE(2)**2)
  odzsq180 = 1 / (180*CCTK_DELTA_SPACE(3)**2)

  odxdy3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2))
  odxdz3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3))
  odydz3600 = 1 / (3600*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3))


  ! convert ADM variables to BSSN-like ones
  call ComplexProca_adm2bssn(CCTK_PASS_FTOF)


  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(alph, beta, hh, hu, trk, dethh, ch,&
  !$OMP lE1, lA1, lAphi1, lZeta1,&
  !$OMP lE2, lA2, lAphi2, lZeta2,&
  !$OMP d1_alph, d1_beta, d1_hh, d1_ch,&
  !$OMP d1_lE1, d1_lA1, d1_lZeta1, d1_lAphi1,&
  !$OMP d2_lA1, ad1_lE1, ad1_lA1, ad1_lZeta1, ad1_lAphi1,&
  !$OMP d1_f1, cf1, cf2, cd_lA1, cd_dA1,&
  !$OMP rhs_lE1, rhs_lA1, rhs_lZeta1, rhs_lAphi1,&
  !$OMP d1_lE2, d1_lA2, d1_lZeta2, d1_lAphi2,&
  !$OMP d2_lA2, ad1_lE2, ad1_lA2, ad1_lZeta2, ad1_lAphi2,&
  !$OMP d1_f2, cd_lA2, cd_dA2,&
  !$OMP rhs_lE2, rhs_lA2, rhs_lZeta2, rhs_lAphi2,&
  !$OMP i, j, k,&
  !$OMP di, dj, dk,&
  !$OMP a, b, c, m, n)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------
    ch        = chi(i,j,k)

    hh(1,1)   = hxx(i,j,k)
    hh(1,2)   = hxy(i,j,k)
    hh(1,3)   = hxz(i,j,k)
    hh(2,2)   = hyy(i,j,k)
    hh(2,3)   = hyz(i,j,k)
    hh(3,3)   = hzz(i,j,k)
    hh(2,1)   = hh(1,2)
    hh(3,1)   = hh(1,3)
    hh(3,2)   = hh(2,3)

    trk       = tracek(i,j,k)

    alph      = alp(i,j,k)

    beta(1)   = betax(i,j,k)
    beta(2)   = betay(i,j,k)
    beta(3)   = betaz(i,j,k)

    lE1(1)     = E1x(i,j,k)
    lE1(2)     = E1y(i,j,k)
    lE1(3)     = E1z(i,j,k)

    lA1(1)     = A1x(i,j,k)
    lA1(2)     = A1y(i,j,k)
    lA1(3)     = A1z(i,j,k)

    lZeta1     = Zeta1(i,j,k)
    lAphi1     = Aphi1(i,j,k)
    
    lE2(1)     = E2x(i,j,k)
    lE2(2)     = E2y(i,j,k)
    lE2(3)     = E2z(i,j,k)

    lA2(1)     = A2x(i,j,k)
    lA2(2)     = A2y(i,j,k)
    lA2(3)     = A2z(i,j,k)

    lZeta2     = Zeta2(i,j,k)
    lAphi2     = Aphi2(i,j,k)
    !-------------------------------------------


    !------------ Invert 3-metric ----------------
    ! NOTE: deth = 1 by construction, but that is not satisfied numerically
    dethh =       hh(1,1) * hh(2,2) * hh(3,3)                              &
            + 2 * hh(1,2) * hh(1,3) * hh(2,3)                              &
            -     hh(1,1) * hh(2,3) ** 2                                   &
            -     hh(2,2) * hh(1,3) ** 2                                   &
            -     hh(3,3) * hh(1,2) ** 2
    hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     ) / dethh
    hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     ) / dethh
    hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     ) / dethh
    hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3)) / dethh
    hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2)) / dethh
    hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1)) / dethh
    hu(2,1) = hu(1,2)
    hu(3,1) = hu(1,3)
    hu(3,2) = hu(2,3)
    !-------------------------------------------

    if (derivs_order == 4) then

      !------------- Centered 1st derivatives -----------

      ! d1_ch(3)
      d1_ch(1) = (   -chi(i+2,j,k) + 8*chi(i+1,j,k)                          &
                  - 8*chi(i-1,j,k) +   chi(i-2,j,k) ) / dx12

      d1_ch(2) = (   -chi(i,j+2,k) + 8*chi(i,j+1,k)                          &
                  - 8*chi(i,j-1,k) +   chi(i,j-2,k) ) / dy12

      d1_ch(3) = (   -chi(i,j,k+2) + 8*chi(i,j,k+1)                          &
                  - 8*chi(i,j,k-1) +   chi(i,j,k-2) ) / dz12


      ! d1_hh(3,3,3)
      d1_hh(1,1,1) = (   -hxx(i+2,j,k) + 8*hxx(i+1,j,k)                      &
                      - 8*hxx(i-1,j,k) +   hxx(i-2,j,k) ) / dx12
      d1_hh(1,2,1) = (   -hxy(i+2,j,k) + 8*hxy(i+1,j,k)                      &
                      - 8*hxy(i-1,j,k) +   hxy(i-2,j,k) ) / dx12
      d1_hh(1,3,1) = (   -hxz(i+2,j,k) + 8*hxz(i+1,j,k)                      &
                      - 8*hxz(i-1,j,k) +   hxz(i-2,j,k) ) / dx12
      d1_hh(2,2,1) = (   -hyy(i+2,j,k) + 8*hyy(i+1,j,k)                      &
                      - 8*hyy(i-1,j,k) +   hyy(i-2,j,k) ) / dx12
      d1_hh(2,3,1) = (   -hyz(i+2,j,k) + 8*hyz(i+1,j,k)                      &
                      - 8*hyz(i-1,j,k) +   hyz(i-2,j,k) ) / dx12
      d1_hh(3,3,1) = (   -hzz(i+2,j,k) + 8*hzz(i+1,j,k)                      &
                      - 8*hzz(i-1,j,k) +   hzz(i-2,j,k) ) / dx12

      d1_hh(1,1,2) = (   -hxx(i,j+2,k) + 8*hxx(i,j+1,k)                      &
                      - 8*hxx(i,j-1,k) +   hxx(i,j-2,k) ) / dy12
      d1_hh(1,2,2) = (   -hxy(i,j+2,k) + 8*hxy(i,j+1,k)                      &
                      - 8*hxy(i,j-1,k) +   hxy(i,j-2,k) ) / dy12
      d1_hh(1,3,2) = (   -hxz(i,j+2,k) + 8*hxz(i,j+1,k)                      &
                      - 8*hxz(i,j-1,k) +   hxz(i,j-2,k) ) / dy12
      d1_hh(2,2,2) = (   -hyy(i,j+2,k) + 8*hyy(i,j+1,k)                      &
                      - 8*hyy(i,j-1,k) +   hyy(i,j-2,k) ) / dy12
      d1_hh(2,3,2) = (   -hyz(i,j+2,k) + 8*hyz(i,j+1,k)                      &
                      - 8*hyz(i,j-1,k) +   hyz(i,j-2,k) ) / dy12
      d1_hh(3,3,2) = (   -hzz(i,j+2,k) + 8*hzz(i,j+1,k)                      &
                      - 8*hzz(i,j-1,k) +   hzz(i,j-2,k) ) / dy12

      d1_hh(1,1,3) = (   -hxx(i,j,k+2) + 8*hxx(i,j,k+1)                      &
                      - 8*hxx(i,j,k-1) +   hxx(i,j,k-2) ) / dz12
      d1_hh(1,2,3) = (   -hxy(i,j,k+2) + 8*hxy(i,j,k+1)                      &
                      - 8*hxy(i,j,k-1) +   hxy(i,j,k-2) ) / dz12
      d1_hh(1,3,3) = (   -hxz(i,j,k+2) + 8*hxz(i,j,k+1)                      &
                      - 8*hxz(i,j,k-1) +   hxz(i,j,k-2) ) / dz12
      d1_hh(2,2,3) = (   -hyy(i,j,k+2) + 8*hyy(i,j,k+1)                      &
                      - 8*hyy(i,j,k-1) +   hyy(i,j,k-2) ) / dz12
      d1_hh(2,3,3) = (   -hyz(i,j,k+2) + 8*hyz(i,j,k+1)                      &
                      - 8*hyz(i,j,k-1) +   hyz(i,j,k-2) ) / dz12
      d1_hh(3,3,3) = (   -hzz(i,j,k+2) + 8*hzz(i,j,k+1)                      &
                      - 8*hzz(i,j,k-1) +   hzz(i,j,k-2) ) / dz12

      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)


      ! d1_alph(3)
      d1_alph(1) = (   -alp(i+2,j,k) + 8*alp(i+1,j,k)                        &
                    - 8*alp(i-1,j,k) +   alp(i-2,j,k) ) / dx12

      d1_alph(2) = (   -alp(i,j+2,k) + 8*alp(i,j+1,k)                        &
                    - 8*alp(i,j-1,k) +   alp(i,j-2,k) ) / dy12

      d1_alph(3) = (   -alp(i,j,k+2) + 8*alp(i,j,k+1)                        &
                    - 8*alp(i,j,k-1) +   alp(i,j,k-2) ) / dz12

      ! d1_beta (3,3)
      d1_beta(1,1)  = (   -betax(i+2,j,k) + 8*betax(i+1,j,k)                 &
                       - 8*betax(i-1,j,k) +   betax(i-2,j,k) ) / dx12
      d1_beta(2,1)  = (   -betay(i+2,j,k) + 8*betay(i+1,j,k)                 &
                       - 8*betay(i-1,j,k) +   betay(i-2,j,k) ) / dx12
      d1_beta(3,1)  = (   -betaz(i+2,j,k) + 8*betaz(i+1,j,k)                 &
                       - 8*betaz(i-1,j,k) +   betaz(i-2,j,k) ) / dx12

      d1_beta(1,2)  = (   -betax(i,j+2,k) + 8*betax(i,j+1,k)                 &
                       - 8*betax(i,j-1,k) +   betax(i,j-2,k) ) / dy12
      d1_beta(2,2)  = (   -betay(i,j+2,k) + 8*betay(i,j+1,k)                 &
                       - 8*betay(i,j-1,k) +   betay(i,j-2,k) ) / dy12
      d1_beta(3,2)  = (   -betaz(i,j+2,k) + 8*betaz(i,j+1,k)                 &
                       - 8*betaz(i,j-1,k) +   betaz(i,j-2,k) ) / dy12

      d1_beta(1,3)  = (   -betax(i,j,k+2) + 8*betax(i,j,k+1)                 &
                       - 8*betax(i,j,k-1) +   betax(i,j,k-2) ) / dz12
      d1_beta(2,3)  = (   -betay(i,j,k+2) + 8*betay(i,j,k+1)                 &
                       - 8*betay(i,j,k-1) +   betay(i,j,k-2) ) / dz12
      d1_beta(3,3)  = (   -betaz(i,j,k+2) + 8*betaz(i,j,k+1)                 &
                       - 8*betaz(i,j,k-1) +   betaz(i,j,k-2) ) / dz12

      ! d1_lE1(3,3)
      d1_lE1(1,1) = (   -E1x(i+2,j,k) + 8*E1x(i+1,j,k)               &
                    - 8*E1x(i-1,j,k) +   E1x(i-2,j,k) ) / dx12
      d1_lE1(2,1) = (   -E1y(i+2,j,k) + 8*E1y(i+1,j,k)               &
                    - 8*E1y(i-1,j,k) +   E1y(i-2,j,k) ) / dx12
      d1_lE1(3,1) = (   -E1z(i+2,j,k) + 8*E1z(i+1,j,k)               &
                    - 8*E1z(i-1,j,k) +   E1z(i-2,j,k) ) / dx12

      d1_lE1(1,2) = (   -E1x(i,j+2,k) + 8*E1x(i,j+1,k)               &
                    - 8*E1x(i,j-1,k) +   E1x(i,j-2,k) ) / dy12
      d1_lE1(2,2) = (   -E1y(i,j+2,k) + 8*E1y(i,j+1,k)               &
                    - 8*E1y(i,j-1,k) +   E1y(i,j-2,k) ) / dy12
      d1_lE1(3,2) = (   -E1z(i,j+2,k) + 8*E1z(i,j+1,k)               &
                    - 8*E1z(i,j-1,k) +   E1z(i,j-2,k) ) / dy12

      d1_lE1(1,3) = (   -E1x(i,j,k+2) + 8*E1x(i,j,k+1)               &
                    - 8*E1x(i,j,k-1) +   E1x(i,j,k-2) ) / dz12
      d1_lE1(2,3) = (   -E1y(i,j,k+2) + 8*E1y(i,j,k+1)               &
                    - 8*E1y(i,j,k-1) +   E1y(i,j,k-2) ) / dz12
      d1_lE1(3,3) = (   -E1z(i,j,k+2) + 8*E1z(i,j,k+1)               &
                    - 8*E1z(i,j,k-1) +   E1z(i,j,k-2) ) / dz12

      ! d1_lA1(3,3)
      d1_lA1(1,1) = (   -A1x(i+2,j,k) + 8*A1x(i+1,j,k)               &
                    - 8*A1x(i-1,j,k) +   A1x(i-2,j,k) ) / dx12
      d1_lA1(2,1) = (   -A1y(i+2,j,k) + 8*A1y(i+1,j,k)               &
                    - 8*A1y(i-1,j,k) +   A1y(i-2,j,k) ) / dx12
      d1_lA1(3,1) = (   -A1z(i+2,j,k) + 8*A1z(i+1,j,k)               &
                    - 8*A1z(i-1,j,k) +   A1z(i-2,j,k) ) / dx12

      d1_lA1(1,2) = (   -A1x(i,j+2,k) + 8*A1x(i,j+1,k)               &
                    - 8*A1x(i,j-1,k) +   A1x(i,j-2,k) ) / dy12
      d1_lA1(2,2) = (   -A1y(i,j+2,k) + 8*A1y(i,j+1,k)               &
                    - 8*A1y(i,j-1,k) +   A1y(i,j-2,k) ) / dy12
      d1_lA1(3,2) = (   -A1z(i,j+2,k) + 8*A1z(i,j+1,k)               &
                    - 8*A1z(i,j-1,k) +   A1z(i,j-2,k) ) / dy12

      d1_lA1(1,3) = (   -A1x(i,j,k+2) + 8*A1x(i,j,k+1)               &
                    - 8*A1x(i,j,k-1) +   A1x(i,j,k-2) ) / dz12
      d1_lA1(2,3) = (   -A1y(i,j,k+2) + 8*A1y(i,j,k+1)               &
                    - 8*A1y(i,j,k-1) +   A1y(i,j,k-2) ) / dz12
      d1_lA1(3,3) = (   -A1z(i,j,k+2) + 8*A1z(i,j,k+1)               &
                    - 8*A1z(i,j,k-1) +   A1z(i,j,k-2) ) / dz12

      ! d1_lZeta1(3)
      d1_lZeta1(1) = (   -Zeta1(i+2,j,k) + 8*Zeta1(i+1,j,k)                        &
                     - 8*Zeta1(i-1,j,k) +   Zeta1(i-2,j,k) ) / dx12

      d1_lZeta1(2) = (   -Zeta1(i,j+2,k) + 8*Zeta1(i,j+1,k)                        &
                     - 8*Zeta1(i,j-1,k) +   Zeta1(i,j-2,k) ) / dy12

      d1_lZeta1(3) = (   -Zeta1(i,j,k+2) + 8*Zeta1(i,j,k+1)                        &
                     - 8*Zeta1(i,j,k-1) +   Zeta1(i,j,k-2) ) / dz12

      ! d1_lAphi1(3)
      d1_lAphi1(1) = (   -Aphi1(i+2,j,k) + 8*Aphi1(i+1,j,k)                        &
                     - 8*Aphi1(i-1,j,k) +   Aphi1(i-2,j,k) ) / dx12

      d1_lAphi1(2) = (   -Aphi1(i,j+2,k) + 8*Aphi1(i,j+1,k)                        &
                     - 8*Aphi1(i,j-1,k) +   Aphi1(i,j-2,k) ) / dy12

      d1_lAphi1(3) = (   -Aphi1(i,j,k+2) + 8*Aphi1(i,j,k+1)                        &
                     - 8*Aphi1(i,j,k-1) +   Aphi1(i,j,k-2) ) / dz12

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_lA1(3,3,3)
      d2_lA1(1,1,1) = (  -A1x(i+2,j,k) + 16*A1x(i+1,j,k) - 30*A1x(i,j,k) &
                    + 16*A1x(i-1,j,k) -    A1x(i-2,j,k) ) / dxsq12
      d2_lA1(2,1,1) = (  -A1y(i+2,j,k) + 16*A1y(i+1,j,k) - 30*A1y(i,j,k) &
                    + 16*A1y(i-1,j,k) -    A1y(i-2,j,k) ) / dxsq12
      d2_lA1(3,1,1) = (  -A1z(i+2,j,k) + 16*A1z(i+1,j,k) - 30*A1z(i,j,k) &
                    + 16*A1z(i-1,j,k) -    A1z(i-2,j,k) ) / dxsq12

      d2_lA1(1,2,2) = (  -A1x(i,j+2,k) + 16*A1x(i,j+1,k) - 30*A1x(i,j,k) &
                    + 16*A1x(i,j-1,k) -    A1x(i,j-2,k) ) / dysq12
      d2_lA1(2,2,2) = (  -A1y(i,j+2,k) + 16*A1y(i,j+1,k) - 30*A1y(i,j,k) &
                    + 16*A1y(i,j-1,k) -    A1y(i,j-2,k) ) / dysq12
      d2_lA1(3,2,2) = (  -A1z(i,j+2,k) + 16*A1z(i,j+1,k) - 30*A1z(i,j,k) &
                    + 16*A1z(i,j-1,k) -    A1z(i,j-2,k) ) / dysq12

      d2_lA1(1,3,3) = (  -A1x(i,j,k+2) + 16*A1x(i,j,k+1) - 30*A1x(i,j,k) &
                    + 16*A1x(i,j,k-1) -    A1x(i,j,k-2) ) / dzsq12
      d2_lA1(2,3,3) = (  -A1y(i,j,k+2) + 16*A1y(i,j,k+1) - 30*A1y(i,j,k) &
                    + 16*A1y(i,j,k-1) -    A1y(i,j,k-2) ) / dzsq12
      d2_lA1(3,3,3) = (  -A1z(i,j,k+2) + 16*A1z(i,j,k+1) - 30*A1z(i,j,k) &
                    + 16*A1z(i,j,k-1) -    A1z(i,j,k-2) ) / dzsq12

      d2_lA1(1,1,2) = (  -A1x(i-2,j+2,k) +  8*A1x(i-1,j+2,k) -  8*A1x(i+1,j+2,k) +   A1x(i+2,j+2,k) &
                     + 8*A1x(i-2,j+1,k) - 64*A1x(i-1,j+1,k) + 64*A1x(i+1,j+1,k) - 8*A1x(i+2,j+1,k) &
                     - 8*A1x(i-2,j-1,k) + 64*A1x(i-1,j-1,k) - 64*A1x(i+1,j-1,k) + 8*A1x(i+2,j-1,k) &
                     +   A1x(i-2,j-2,k) -  8*A1x(i-1,j-2,k) +  8*A1x(i+1,j-2,k) -   A1x(i+2,j-2,k) ) / dxdy144
      d2_lA1(2,1,2) = (  -A1y(i-2,j+2,k) +  8*A1y(i-1,j+2,k) -  8*A1y(i+1,j+2,k) +   A1y(i+2,j+2,k) &
                     + 8*A1y(i-2,j+1,k) - 64*A1y(i-1,j+1,k) + 64*A1y(i+1,j+1,k) - 8*A1y(i+2,j+1,k) &
                     - 8*A1y(i-2,j-1,k) + 64*A1y(i-1,j-1,k) - 64*A1y(i+1,j-1,k) + 8*A1y(i+2,j-1,k) &
                     +   A1y(i-2,j-2,k) -  8*A1y(i-1,j-2,k) +  8*A1y(i+1,j-2,k) -   A1y(i+2,j-2,k) ) / dxdy144
      d2_lA1(3,1,2) = (  -A1z(i-2,j+2,k) +  8*A1z(i-1,j+2,k) -  8*A1z(i+1,j+2,k) +   A1z(i+2,j+2,k) &
                     + 8*A1z(i-2,j+1,k) - 64*A1z(i-1,j+1,k) + 64*A1z(i+1,j+1,k) - 8*A1z(i+2,j+1,k) &
                     - 8*A1z(i-2,j-1,k) + 64*A1z(i-1,j-1,k) - 64*A1z(i+1,j-1,k) + 8*A1z(i+2,j-1,k) &
                     +   A1z(i-2,j-2,k) -  8*A1z(i-1,j-2,k) +  8*A1z(i+1,j-2,k) -   A1z(i+2,j-2,k) ) / dxdy144

      d2_lA1(1,1,3) = (  -A1x(i-2,j,k+2) +  8*A1x(i-1,j,k+2) -  8*A1x(i+1,j,k+2) +   A1x(i+2,j,k+2) &
                     + 8*A1x(i-2,j,k+1) - 64*A1x(i-1,j,k+1) + 64*A1x(i+1,j,k+1) - 8*A1x(i+2,j,k+1) &
                     - 8*A1x(i-2,j,k-1) + 64*A1x(i-1,j,k-1) - 64*A1x(i+1,j,k-1) + 8*A1x(i+2,j,k-1) &
                     +   A1x(i-2,j,k-2) -  8*A1x(i-1,j,k-2) +  8*A1x(i+1,j,k-2) -   A1x(i+2,j,k-2) ) / dxdz144
      d2_lA1(2,1,3) = (  -A1y(i-2,j,k+2) +  8*A1y(i-1,j,k+2) -  8*A1y(i+1,j,k+2) +   A1y(i+2,j,k+2) &
                     + 8*A1y(i-2,j,k+1) - 64*A1y(i-1,j,k+1) + 64*A1y(i+1,j,k+1) - 8*A1y(i+2,j,k+1) &
                     - 8*A1y(i-2,j,k-1) + 64*A1y(i-1,j,k-1) - 64*A1y(i+1,j,k-1) + 8*A1y(i+2,j,k-1) &
                     +   A1y(i-2,j,k-2) -  8*A1y(i-1,j,k-2) +  8*A1y(i+1,j,k-2) -   A1y(i+2,j,k-2) ) / dxdz144
      d2_lA1(3,1,3) = (  -A1z(i-2,j,k+2) +  8*A1z(i-1,j,k+2) -  8*A1z(i+1,j,k+2) +   A1z(i+2,j,k+2) &
                     + 8*A1z(i-2,j,k+1) - 64*A1z(i-1,j,k+1) + 64*A1z(i+1,j,k+1) - 8*A1z(i+2,j,k+1) &
                     - 8*A1z(i-2,j,k-1) + 64*A1z(i-1,j,k-1) - 64*A1z(i+1,j,k-1) + 8*A1z(i+2,j,k-1) &
                     +   A1z(i-2,j,k-2) -  8*A1z(i-1,j,k-2) +  8*A1z(i+1,j,k-2) -   A1z(i+2,j,k-2) ) / dxdz144

      d2_lA1(1,2,3) = (  -A1x(i,j-2,k+2) +  8*A1x(i,j-1,k+2) -  8*A1x(i,j+1,k+2) +   A1x(i,j+2,k+2) &
                     + 8*A1x(i,j-2,k+1) - 64*A1x(i,j-1,k+1) + 64*A1x(i,j+1,k+1) - 8*A1x(i,j+2,k+1) &
                     - 8*A1x(i,j-2,k-1) + 64*A1x(i,j-1,k-1) - 64*A1x(i,j+1,k-1) + 8*A1x(i,j+2,k-1) &
                     +   A1x(i,j-2,k-2) -  8*A1x(i,j-1,k-2) +  8*A1x(i,j+1,k-2) -   A1x(i,j+2,k-2) ) / dydz144
      d2_lA1(2,2,3) = (  -A1y(i,j-2,k+2) +  8*A1y(i,j-1,k+2) -  8*A1y(i,j+1,k+2) +   A1y(i,j+2,k+2) &
                     + 8*A1y(i,j-2,k+1) - 64*A1y(i,j-1,k+1) + 64*A1y(i,j+1,k+1) - 8*A1y(i,j+2,k+1) &
                     - 8*A1y(i,j-2,k-1) + 64*A1y(i,j-1,k-1) - 64*A1y(i,j+1,k-1) + 8*A1y(i,j+2,k-1) &
                     +   A1y(i,j-2,k-2) -  8*A1y(i,j-1,k-2) +  8*A1y(i,j+1,k-2) -   A1y(i,j+2,k-2) ) / dydz144
      d2_lA1(3,2,3) = (  -A1z(i,j-2,k+2) +  8*A1z(i,j-1,k+2) -  8*A1z(i,j+1,k+2) +   A1z(i,j+2,k+2) &
                     + 8*A1z(i,j-2,k+1) - 64*A1z(i,j-1,k+1) + 64*A1z(i,j+1,k+1) - 8*A1z(i,j+2,k+1) &
                     - 8*A1z(i,j-2,k-1) + 64*A1z(i,j-1,k-1) - 64*A1z(i,j+1,k-1) + 8*A1z(i,j+2,k-1) &
                     +   A1z(i,j-2,k-2) -  8*A1z(i,j-1,k-2) +  8*A1z(i,j+1,k-2) -   A1z(i,j+2,k-2) ) / dydz144

      d2_lA1(:,2,1) = d2_lA1(:,1,2)
      d2_lA1(:,3,1) = d2_lA1(:,1,3)
      d2_lA1(:,3,2) = d2_lA1(:,2,3)
      
      ! d1_lE2(3,3)
      d1_lE2(1,1) = (   -E2x(i+2,j,k) + 8*E2x(i+1,j,k)               &
                    - 8*E2x(i-1,j,k) +   E2x(i-2,j,k) ) / dx12
      d1_lE2(2,1) = (   -E2y(i+2,j,k) + 8*E2y(i+1,j,k)               &
                    - 8*E2y(i-1,j,k) +   E2y(i-2,j,k) ) / dx12
      d1_lE2(3,1) = (   -E2z(i+2,j,k) + 8*E2z(i+1,j,k)               &
                    - 8*E2z(i-1,j,k) +   E2z(i-2,j,k) ) / dx12

      d1_lE2(1,2) = (   -E2x(i,j+2,k) + 8*E2x(i,j+1,k)               &
                    - 8*E2x(i,j-1,k) +   E2x(i,j-2,k) ) / dy12
      d1_lE2(2,2) = (   -E2y(i,j+2,k) + 8*E2y(i,j+1,k)               &
                    - 8*E2y(i,j-1,k) +   E2y(i,j-2,k) ) / dy12
      d1_lE2(3,2) = (   -E2z(i,j+2,k) + 8*E2z(i,j+1,k)               &
                    - 8*E2z(i,j-1,k) +   E2z(i,j-2,k) ) / dy12

      d1_lE2(1,3) = (   -E2x(i,j,k+2) + 8*E2x(i,j,k+1)               &
                    - 8*E2x(i,j,k-1) +   E2x(i,j,k-2) ) / dz12
      d1_lE2(2,3) = (   -E2y(i,j,k+2) + 8*E2y(i,j,k+1)               &
                    - 8*E2y(i,j,k-1) +   E2y(i,j,k-2) ) / dz12
      d1_lE2(3,3) = (   -E2z(i,j,k+2) + 8*E2z(i,j,k+1)               &
                    - 8*E2z(i,j,k-1) +   E2z(i,j,k-2) ) / dz12

      ! d1_lA2(3,3)
      d1_lA2(1,1) = (   -A2x(i+2,j,k) + 8*A2x(i+1,j,k)               &
                    - 8*A2x(i-1,j,k) +   A2x(i-2,j,k) ) / dx12
      d1_lA2(2,1) = (   -A2y(i+2,j,k) + 8*A2y(i+1,j,k)               &
                    - 8*A2y(i-1,j,k) +   A2y(i-2,j,k) ) / dx12
      d1_lA2(3,1) = (   -A2z(i+2,j,k) + 8*A2z(i+1,j,k)               &
                    - 8*A2z(i-1,j,k) +   A2z(i-2,j,k) ) / dx12

      d1_lA2(1,2) = (   -A2x(i,j+2,k) + 8*A2x(i,j+1,k)               &
                    - 8*A2x(i,j-1,k) +   A2x(i,j-2,k) ) / dy12
      d1_lA2(2,2) = (   -A2y(i,j+2,k) + 8*A2y(i,j+1,k)               &
                    - 8*A2y(i,j-1,k) +   A2y(i,j-2,k) ) / dy12
      d1_lA2(3,2) = (   -A2z(i,j+2,k) + 8*A2z(i,j+1,k)               &
                    - 8*A2z(i,j-1,k) +   A2z(i,j-2,k) ) / dy12

      d1_lA2(1,3) = (   -A2x(i,j,k+2) + 8*A2x(i,j,k+1)               &
                    - 8*A2x(i,j,k-1) +   A2x(i,j,k-2) ) / dz12
      d1_lA2(2,3) = (   -A2y(i,j,k+2) + 8*A2y(i,j,k+1)               &
                    - 8*A2y(i,j,k-1) +   A2y(i,j,k-2) ) / dz12
      d1_lA2(3,3) = (   -A2z(i,j,k+2) + 8*A2z(i,j,k+1)               &
                    - 8*A2z(i,j,k-1) +   A2z(i,j,k-2) ) / dz12

      ! d1_lZeta2(3)
      d1_lZeta2(1) = (   -Zeta2(i+2,j,k) + 8*Zeta2(i+1,j,k)                        &
                     - 8*Zeta2(i-1,j,k) +   Zeta2(i-2,j,k) ) / dx12

      d1_lZeta2(2) = (   -Zeta2(i,j+2,k) + 8*Zeta2(i,j+1,k)                        &
                     - 8*Zeta2(i,j-1,k) +   Zeta2(i,j-2,k) ) / dy12

      d1_lZeta2(3) = (   -Zeta2(i,j,k+2) + 8*Zeta2(i,j,k+1)                        &
                     - 8*Zeta2(i,j,k-1) +   Zeta2(i,j,k-2) ) / dz12

      ! d1_lAphi2(3)
      d1_lAphi2(1) = (   -Aphi2(i+2,j,k) + 8*Aphi2(i+1,j,k)                        &
                     - 8*Aphi2(i-1,j,k) +   Aphi2(i-2,j,k) ) / dx12

      d1_lAphi2(2) = (   -Aphi2(i,j+2,k) + 8*Aphi2(i,j+1,k)                        &
                     - 8*Aphi2(i,j-1,k) +   Aphi2(i,j-2,k) ) / dy12

      d1_lAphi2(3) = (   -Aphi2(i,j,k+2) + 8*Aphi2(i,j,k+1)                        &
                     - 8*Aphi2(i,j,k-1) +   Aphi2(i,j,k-2) ) / dz12

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_lA2(3,3,3)
      d2_lA2(1,1,1) = (  -A2x(i+2,j,k) + 16*A2x(i+1,j,k) - 30*A2x(i,j,k) &
                    + 16*A2x(i-1,j,k) -    A2x(i-2,j,k) ) / dxsq12
      d2_lA2(2,1,1) = (  -A2y(i+2,j,k) + 16*A2y(i+1,j,k) - 30*A2y(i,j,k) &
                    + 16*A2y(i-1,j,k) -    A2y(i-2,j,k) ) / dxsq12
      d2_lA2(3,1,1) = (  -A2z(i+2,j,k) + 16*A2z(i+1,j,k) - 30*A2z(i,j,k) &
                    + 16*A2z(i-1,j,k) -    A2z(i-2,j,k) ) / dxsq12

      d2_lA2(1,2,2) = (  -A2x(i,j+2,k) + 16*A2x(i,j+1,k) - 30*A2x(i,j,k) &
                    + 16*A2x(i,j-1,k) -    A2x(i,j-2,k) ) / dysq12
      d2_lA2(2,2,2) = (  -A2y(i,j+2,k) + 16*A2y(i,j+1,k) - 30*A2y(i,j,k) &
                    + 16*A2y(i,j-1,k) -    A2y(i,j-2,k) ) / dysq12
      d2_lA2(3,2,2) = (  -A2z(i,j+2,k) + 16*A2z(i,j+1,k) - 30*A2z(i,j,k) &
                    + 16*A2z(i,j-1,k) -    A2z(i,j-2,k) ) / dysq12

      d2_lA2(1,3,3) = (  -A2x(i,j,k+2) + 16*A2x(i,j,k+1) - 30*A2x(i,j,k) &
                    + 16*A2x(i,j,k-1) -    A2x(i,j,k-2) ) / dzsq12
      d2_lA2(2,3,3) = (  -A2y(i,j,k+2) + 16*A2y(i,j,k+1) - 30*A2y(i,j,k) &
                    + 16*A2y(i,j,k-1) -    A2y(i,j,k-2) ) / dzsq12
      d2_lA2(3,3,3) = (  -A2z(i,j,k+2) + 16*A2z(i,j,k+1) - 30*A2z(i,j,k) &
                    + 16*A2z(i,j,k-1) -    A2z(i,j,k-2) ) / dzsq12

      d2_lA2(1,1,2) = (  -A2x(i-2,j+2,k) +  8*A2x(i-1,j+2,k) -  8*A2x(i+1,j+2,k) +   A2x(i+2,j+2,k) &
                     + 8*A2x(i-2,j+1,k) - 64*A2x(i-1,j+1,k) + 64*A2x(i+1,j+1,k) - 8*A2x(i+2,j+1,k) &
                     - 8*A2x(i-2,j-1,k) + 64*A2x(i-1,j-1,k) - 64*A2x(i+1,j-1,k) + 8*A2x(i+2,j-1,k) &
                     +   A2x(i-2,j-2,k) -  8*A2x(i-1,j-2,k) +  8*A2x(i+1,j-2,k) -   A2x(i+2,j-2,k) ) / dxdy144
      d2_lA2(2,1,2) = (  -A2y(i-2,j+2,k) +  8*A2y(i-1,j+2,k) -  8*A2y(i+1,j+2,k) +   A2y(i+2,j+2,k) &
                     + 8*A2y(i-2,j+1,k) - 64*A2y(i-1,j+1,k) + 64*A2y(i+1,j+1,k) - 8*A2y(i+2,j+1,k) &
                     - 8*A2y(i-2,j-1,k) + 64*A2y(i-1,j-1,k) - 64*A2y(i+1,j-1,k) + 8*A2y(i+2,j-1,k) &
                     +   A2y(i-2,j-2,k) -  8*A2y(i-1,j-2,k) +  8*A2y(i+1,j-2,k) -   A2y(i+2,j-2,k) ) / dxdy144
      d2_lA2(3,1,2) = (  -A2z(i-2,j+2,k) +  8*A2z(i-1,j+2,k) -  8*A2z(i+1,j+2,k) +   A2z(i+2,j+2,k) &
                     + 8*A2z(i-2,j+1,k) - 64*A2z(i-1,j+1,k) + 64*A2z(i+1,j+1,k) - 8*A2z(i+2,j+1,k) &
                     - 8*A2z(i-2,j-1,k) + 64*A2z(i-1,j-1,k) - 64*A2z(i+1,j-1,k) + 8*A2z(i+2,j-1,k) &
                     +   A2z(i-2,j-2,k) -  8*A2z(i-1,j-2,k) +  8*A2z(i+1,j-2,k) -   A2z(i+2,j-2,k) ) / dxdy144

      d2_lA2(1,1,3) = (  -A2x(i-2,j,k+2) +  8*A2x(i-1,j,k+2) -  8*A2x(i+1,j,k+2) +   A2x(i+2,j,k+2) &
                     + 8*A2x(i-2,j,k+1) - 64*A2x(i-1,j,k+1) + 64*A2x(i+1,j,k+1) - 8*A2x(i+2,j,k+1) &
                     - 8*A2x(i-2,j,k-1) + 64*A2x(i-1,j,k-1) - 64*A2x(i+1,j,k-1) + 8*A2x(i+2,j,k-1) &
                     +   A2x(i-2,j,k-2) -  8*A2x(i-1,j,k-2) +  8*A2x(i+1,j,k-2) -   A2x(i+2,j,k-2) ) / dxdz144
      d2_lA2(2,1,3) = (  -A2y(i-2,j,k+2) +  8*A2y(i-1,j,k+2) -  8*A2y(i+1,j,k+2) +   A2y(i+2,j,k+2) &
                     + 8*A2y(i-2,j,k+1) - 64*A2y(i-1,j,k+1) + 64*A2y(i+1,j,k+1) - 8*A2y(i+2,j,k+1) &
                     - 8*A2y(i-2,j,k-1) + 64*A2y(i-1,j,k-1) - 64*A2y(i+1,j,k-1) + 8*A2y(i+2,j,k-1) &
                     +   A2y(i-2,j,k-2) -  8*A2y(i-1,j,k-2) +  8*A2y(i+1,j,k-2) -   A2y(i+2,j,k-2) ) / dxdz144
      d2_lA2(3,1,3) = (  -A2z(i-2,j,k+2) +  8*A2z(i-1,j,k+2) -  8*A2z(i+1,j,k+2) +   A2z(i+2,j,k+2) &
                     + 8*A2z(i-2,j,k+1) - 64*A2z(i-1,j,k+1) + 64*A2z(i+1,j,k+1) - 8*A2z(i+2,j,k+1) &
                     - 8*A2z(i-2,j,k-1) + 64*A2z(i-1,j,k-1) - 64*A2z(i+1,j,k-1) + 8*A2z(i+2,j,k-1) &
                     +   A2z(i-2,j,k-2) -  8*A2z(i-1,j,k-2) +  8*A2z(i+1,j,k-2) -   A2z(i+2,j,k-2) ) / dxdz144

      d2_lA2(1,2,3) = (  -A2x(i,j-2,k+2) +  8*A2x(i,j-1,k+2) -  8*A2x(i,j+1,k+2) +   A2x(i,j+2,k+2) &
                     + 8*A2x(i,j-2,k+1) - 64*A2x(i,j-1,k+1) + 64*A2x(i,j+1,k+1) - 8*A2x(i,j+2,k+1) &
                     - 8*A2x(i,j-2,k-1) + 64*A2x(i,j-1,k-1) - 64*A2x(i,j+1,k-1) + 8*A2x(i,j+2,k-1) &
                     +   A2x(i,j-2,k-2) -  8*A2x(i,j-1,k-2) +  8*A2x(i,j+1,k-2) -   A2x(i,j+2,k-2) ) / dydz144
      d2_lA2(2,2,3) = (  -A2y(i,j-2,k+2) +  8*A2y(i,j-1,k+2) -  8*A2y(i,j+1,k+2) +   A2y(i,j+2,k+2) &
                     + 8*A2y(i,j-2,k+1) - 64*A2y(i,j-1,k+1) + 64*A2y(i,j+1,k+1) - 8*A2y(i,j+2,k+1) &
                     - 8*A2y(i,j-2,k-1) + 64*A2y(i,j-1,k-1) - 64*A2y(i,j+1,k-1) + 8*A2y(i,j+2,k-1) &
                     +   A2y(i,j-2,k-2) -  8*A2y(i,j-1,k-2) +  8*A2y(i,j+1,k-2) -   A2y(i,j+2,k-2) ) / dydz144
      d2_lA2(3,2,3) = (  -A2z(i,j-2,k+2) +  8*A2z(i,j-1,k+2) -  8*A2z(i,j+1,k+2) +   A2z(i,j+2,k+2) &
                     + 8*A2z(i,j-2,k+1) - 64*A2z(i,j-1,k+1) + 64*A2z(i,j+1,k+1) - 8*A2z(i,j+2,k+1) &
                     - 8*A2z(i,j-2,k-1) + 64*A2z(i,j-1,k-1) - 64*A2z(i,j+1,k-1) + 8*A2z(i,j+2,k-1) &
                     +   A2z(i,j-2,k-2) -  8*A2z(i,j-1,k-2) +  8*A2z(i,j+1,k-2) -   A2z(i,j+2,k-2) ) / dydz144

      d2_lA2(:,2,1) = d2_lA2(:,1,2)
      d2_lA2(:,3,1) = d2_lA2(:,1,3)
      d2_lA2(:,3,2) = d2_lA2(:,2,3)


      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta(1) ) )
        dj = int( sign( one, beta(2) ) )
        dk = int( sign( one, beta(3) ) )

        ! ad1_lE1(3)
        d1_f1(1) = di * ( -3*E1x(i-di,j,k) - 10*E1x(i,j,k) + 18*E1x(i+di,j,k) &
                        - 6*E1x(i+2*di,j,k) + E1x(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*E1x(i,j-dj,k) - 10*E1x(i,j,k) + 18*E1x(i,j+dj,k) &
                        - 6*E1x(i,j+2*dj,k) + E1x(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*E1x(i,j,k-dk) - 10*E1x(i,j,k) + 18*E1x(i,j,k+dk) &
                        - 6*E1x(i,j,k+2*dk) + E1x(i,j,k+3*dk)) / dz12
        ad1_lE1(1) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * ( -3*E1y(i-di,j,k) - 10*E1y(i,j,k) + 18*E1y(i+di,j,k) &
                        - 6*E1y(i+2*di,j,k) + E1y(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*E1y(i,j-dj,k) - 10*E1y(i,j,k) + 18*E1y(i,j+dj,k) &
                        - 6*E1y(i,j+2*dj,k) + E1y(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*E1y(i,j,k-dk) - 10*E1y(i,j,k) + 18*E1y(i,j,k+dk) &
                        - 6*E1y(i,j,k+2*dk) + E1y(i,j,k+3*dk)) / dz12
        ad1_lE1(2) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * ( -3*E1z(i-di,j,k) - 10*E1z(i,j,k) + 18*E1z(i+di,j,k) &
                        - 6*E1z(i+2*di,j,k) + E1z(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*E1z(i,j-dj,k) - 10*E1z(i,j,k) + 18*E1z(i,j+dj,k) &
                        - 6*E1z(i,j+2*dj,k) + E1z(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*E1z(i,j,k-dk) - 10*E1z(i,j,k) + 18*E1z(i,j,k+dk) &
                        - 6*E1z(i,j,k+2*dk) + E1z(i,j,k+3*dk)) / dz12
        ad1_lE1(3) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        ! ad1_lA1(3)
        d1_f1(1) = di * ( -3*A1x(i-di,j,k) - 10*A1x(i,j,k) + 18*A1x(i+di,j,k) &
                        - 6*A1x(i+2*di,j,k) + A1x(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*A1x(i,j-dj,k) - 10*A1x(i,j,k) + 18*A1x(i,j+dj,k) &
                        - 6*A1x(i,j+2*dj,k) + A1x(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*A1x(i,j,k-dk) - 10*A1x(i,j,k) + 18*A1x(i,j,k+dk) &
                        - 6*A1x(i,j,k+2*dk) + A1x(i,j,k+3*dk)) / dz12
        ad1_lA1(1) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * ( -3*A1y(i-di,j,k) - 10*A1y(i,j,k) + 18*A1y(i+di,j,k) &
                        - 6*A1y(i+2*di,j,k) + A1y(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*A1y(i,j-dj,k) - 10*A1y(i,j,k) + 18*A1y(i,j+dj,k) &
                        - 6*A1y(i,j+2*dj,k) + A1y(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*A1y(i,j,k-dk) - 10*A1y(i,j,k) + 18*A1y(i,j,k+dk) &
                        - 6*A1y(i,j,k+2*dk) + A1y(i,j,k+3*dk)) / dz12
        ad1_lA1(2) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * ( -3*A1z(i-di,j,k) - 10*A1z(i,j,k) + 18*A1z(i+di,j,k) &
                        - 6*A1z(i+2*di,j,k) + A1z(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*A1z(i,j-dj,k) - 10*A1z(i,j,k) + 18*A1z(i,j+dj,k) &
                        - 6*A1z(i,j+2*dj,k) + A1z(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*A1z(i,j,k-dk) - 10*A1z(i,j,k) + 18*A1z(i,j,k+dk) &
                        - 6*A1z(i,j,k+2*dk) + A1z(i,j,k+3*dk)) / dz12
        ad1_lA1(3) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        ! ad1_lZeta1
        d1_f1(1) = di * ( -3*Zeta1(i-di,j,k) - 10*Zeta1(i,j,k) + 18*Zeta1(i+di,j,k)  &
                        - 6*Zeta1(i+2*di,j,k) + Zeta1(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*Zeta1(i,j-dj,k) - 10*Zeta1(i,j,k) + 18*Zeta1(i,j+dj,k)  &
                        - 6*Zeta1(i,j+2*dj,k) + Zeta1(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*Zeta1(i,j,k-dk) - 10*Zeta1(i,j,k) + 18*Zeta1(i,j,k+dk)  &
                        - 6*Zeta1(i,j,k+2*dk) + Zeta1(i,j,k+3*dk)) / dz12
        ad1_lZeta1 = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        ! ad1_lAphi1
        d1_f1(1) = di * ( -3*Aphi1(i-di,j,k) - 10*Aphi1(i,j,k) + 18*Aphi1(i+di,j,k)  &
                        - 6*Aphi1(i+2*di,j,k) + Aphi1(i+3*di,j,k)) / dx12
        d1_f1(2) = dj * ( -3*Aphi1(i,j-dj,k) - 10*Aphi1(i,j,k) + 18*Aphi1(i,j+dj,k)  &
                        - 6*Aphi1(i,j+2*dj,k) + Aphi1(i,j+3*dj,k)) / dy12
        d1_f1(3) = dk * ( -3*Aphi1(i,j,k-dk) - 10*Aphi1(i,j,k) + 18*Aphi1(i,j,k+dk)  &
                        - 6*Aphi1(i,j,k+2*dk) + Aphi1(i,j,k+3*dk)) / dz12
        ad1_lAphi1 = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)
        
        ! ad1_lE2(3)
        d1_f2(1) = di * ( -3*E2x(i-di,j,k) - 10*E2x(i,j,k) + 18*E2x(i+di,j,k) &
                        - 6*E2x(i+2*di,j,k) + E2x(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*E2x(i,j-dj,k) - 10*E2x(i,j,k) + 18*E2x(i,j+dj,k) &
                        - 6*E2x(i,j+2*dj,k) + E2x(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*E2x(i,j,k-dk) - 10*E2x(i,j,k) + 18*E2x(i,j,k+dk) &
                        - 6*E2x(i,j,k+2*dk) + E2x(i,j,k+3*dk)) / dz12
        ad1_lE2(1) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * ( -3*E2y(i-di,j,k) - 10*E2y(i,j,k) + 18*E2y(i+di,j,k) &
                        - 6*E2y(i+2*di,j,k) + E2y(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*E2y(i,j-dj,k) - 10*E2y(i,j,k) + 18*E2y(i,j+dj,k) &
                        - 6*E2y(i,j+2*dj,k) + E2y(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*E2y(i,j,k-dk) - 10*E2y(i,j,k) + 18*E2y(i,j,k+dk) &
                        - 6*E2y(i,j,k+2*dk) + E2y(i,j,k+3*dk)) / dz12
        ad1_lE2(2) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * ( -3*E2z(i-di,j,k) - 10*E2z(i,j,k) + 18*E2z(i+di,j,k) &
                        - 6*E2z(i+2*di,j,k) + E2z(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*E2z(i,j-dj,k) - 10*E2z(i,j,k) + 18*E2z(i,j+dj,k) &
                        - 6*E2z(i,j+2*dj,k) + E2z(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*E2z(i,j,k-dk) - 10*E2z(i,j,k) + 18*E2z(i,j,k+dk) &
                        - 6*E2z(i,j,k+2*dk) + E2z(i,j,k+3*dk)) / dz12
        ad1_lE2(3) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        ! ad1_lA2(3)
        d1_f2(1) = di * ( -3*A2x(i-di,j,k) - 10*A2x(i,j,k) + 18*A2x(i+di,j,k) &
                        - 6*A2x(i+2*di,j,k) + A2x(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*A2x(i,j-dj,k) - 10*A2x(i,j,k) + 18*A2x(i,j+dj,k) &
                        - 6*A2x(i,j+2*dj,k) + A2x(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*A2x(i,j,k-dk) - 10*A2x(i,j,k) + 18*A2x(i,j,k+dk) &
                        - 6*A2x(i,j,k+2*dk) + A2x(i,j,k+3*dk)) / dz12
        ad1_lA2(1) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * ( -3*A2y(i-di,j,k) - 10*A2y(i,j,k) + 18*A2y(i+di,j,k) &
                        - 6*A2y(i+2*di,j,k) + A2y(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*A2y(i,j-dj,k) - 10*A2y(i,j,k) + 18*A2y(i,j+dj,k) &
                        - 6*A2y(i,j+2*dj,k) + A2y(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*A2y(i,j,k-dk) - 10*A2y(i,j,k) + 18*A2y(i,j,k+dk) &
                        - 6*A2y(i,j,k+2*dk) + A2y(i,j,k+3*dk)) / dz12
        ad1_lA2(2) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * ( -3*A2z(i-di,j,k) - 10*A2z(i,j,k) + 18*A2z(i+di,j,k) &
                        - 6*A2z(i+2*di,j,k) + A2z(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*A2z(i,j-dj,k) - 10*A2z(i,j,k) + 18*A2z(i,j+dj,k) &
                        - 6*A2z(i,j+2*dj,k) + A2z(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*A2z(i,j,k-dk) - 10*A2z(i,j,k) + 18*A2z(i,j,k+dk) &
                        - 6*A2z(i,j,k+2*dk) + A2z(i,j,k+3*dk)) / dz12
        ad1_lA2(3) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        ! ad1_lZeta2
        d1_f2(1) = di * ( -3*Zeta2(i-di,j,k) - 10*Zeta2(i,j,k) + 18*Zeta2(i+di,j,k)  &
                        - 6*Zeta2(i+2*di,j,k) + Zeta2(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*Zeta2(i,j-dj,k) - 10*Zeta2(i,j,k) + 18*Zeta2(i,j+dj,k)  &
                        - 6*Zeta2(i,j+2*dj,k) + Zeta2(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*Zeta2(i,j,k-dk) - 10*Zeta2(i,j,k) + 18*Zeta2(i,j,k+dk)  &
                        - 6*Zeta2(i,j,k+2*dk) + Zeta2(i,j,k+3*dk)) / dz12
        ad1_lZeta2 = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        ! ad1_lAphi2
        d1_f2(1) = di * ( -3*Aphi2(i-di,j,k) - 10*Aphi2(i,j,k) + 18*Aphi2(i+di,j,k)  &
                        - 6*Aphi2(i+2*di,j,k) + Aphi2(i+3*di,j,k)) / dx12
        d1_f2(2) = dj * ( -3*Aphi2(i,j-dj,k) - 10*Aphi2(i,j,k) + 18*Aphi2(i,j+dj,k)  &
                        - 6*Aphi2(i,j+2*dj,k) + Aphi2(i,j+3*dj,k)) / dy12
        d1_f2(3) = dk * ( -3*Aphi2(i,j,k-dk) - 10*Aphi2(i,j,k) + 18*Aphi2(i,j,k+dk)  &
                        - 6*Aphi2(i,j,k+2*dk) + Aphi2(i,j,k+3*dk)) / dz12
        ad1_lAphi2 = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

      else

        ! ad1_lE1(3)
        ad1_lE1(1) = beta(1)*d1_lE1(1,1) + beta(2)*d1_lE1(1,2) + beta(3)*d1_lE1(1,3)
        ad1_lE1(2) = beta(1)*d1_lE1(2,1) + beta(2)*d1_lE1(2,2) + beta(3)*d1_lE1(2,3)
        ad1_lE1(3) = beta(1)*d1_lE1(3,1) + beta(2)*d1_lE1(3,2) + beta(3)*d1_lE1(3,3)

        ! ad1_lA1(3)
        ad1_lA1(1) = beta(1)*d1_lA1(1,1) + beta(2)*d1_lA1(1,2) + beta(3)*d1_lA1(1,3)
        ad1_lA1(2) = beta(1)*d1_lA1(2,1) + beta(2)*d1_lA1(2,2) + beta(3)*d1_lA1(2,3)
        ad1_lA1(3) = beta(1)*d1_lA1(3,1) + beta(2)*d1_lA1(3,2) + beta(3)*d1_lA1(3,3)

        ! ad1_lZeta1
        ad1_lZeta1 = beta(1)*d1_lZeta1(1) + beta(2)*d1_lZeta1(2) + beta(3)*d1_lZeta1(3)

        ! ad1_lAphi1
        ad1_lAphi1 = beta(1)*d1_lAphi1(1) + beta(2)*d1_lAphi1(2) + beta(3)*d1_lAphi1(3)
        
        ! ad1_lE2(3)
        ad1_lE2(1) = beta(1)*d1_lE2(1,1) + beta(2)*d1_lE2(1,2) + beta(3)*d1_lE2(1,3)
        ad1_lE2(2) = beta(1)*d1_lE2(2,1) + beta(2)*d1_lE2(2,2) + beta(3)*d1_lE2(2,3)
        ad1_lE2(3) = beta(1)*d1_lE2(3,1) + beta(2)*d1_lE2(3,2) + beta(3)*d1_lE2(3,3)

        ! ad1_lA2(3)
        ad1_lA2(1) = beta(1)*d1_lA2(1,1) + beta(2)*d1_lA2(1,2) + beta(3)*d1_lA2(1,3)
        ad1_lA2(2) = beta(1)*d1_lA2(2,1) + beta(2)*d1_lA2(2,2) + beta(3)*d1_lA2(2,3)
        ad1_lA2(3) = beta(1)*d1_lA2(3,1) + beta(2)*d1_lA2(3,2) + beta(3)*d1_lA2(3,3)

        ! ad1_lZeta2
        ad1_lZeta2 = beta(1)*d1_lZeta2(1) + beta(2)*d1_lZeta2(2) + beta(3)*d1_lZeta2(3)

        ! ad1_lAphi2
        ad1_lAphi2 = beta(1)*d1_lAphi2(1) + beta(2)*d1_lAphi2(2) + beta(3)*d1_lAphi2(3)

      end if
      !-------------------------------------------

    else if (derivs_order == 6) then

      !------------ Centered 1st derivatives -----
      ! d1_ch(3)
      d1_ch(1) = (  chi(i+3,j,k) - 9*chi(i+2,j,k) + 45*chi(i+1,j,k)          &
                  - chi(i-3,j,k) + 9*chi(i-2,j,k) - 45*chi(i-1,j,k) ) * odx60

      d1_ch(2) = (  chi(i,j+3,k) - 9*chi(i,j+2,k) + 45*chi(i,j+1,k)          &
                  - chi(i,j-3,k) + 9*chi(i,j-2,k) - 45*chi(i,j-1,k) ) * ody60

      d1_ch(3) = (  chi(i,j,k+3) - 9*chi(i,j,k+2) + 45*chi(i,j,k+1)          &
                  - chi(i,j,k-3) + 9*chi(i,j,k-2) - 45*chi(i,j,k-1) ) * odz60

      ! d1_hh(3,3,3)
      d1_hh(1,1,1) = (  hxx(i+3,j,k) - 9*hxx(i+2,j,k) + 45*hxx(i+1,j,k)      &
                      - hxx(i-3,j,k) + 9*hxx(i-2,j,k) - 45*hxx(i-1,j,k) ) * odx60
      d1_hh(1,2,1) = (  hxy(i+3,j,k) - 9*hxy(i+2,j,k) + 45*hxy(i+1,j,k)      &
                      - hxy(i-3,j,k) + 9*hxy(i-2,j,k) - 45*hxy(i-1,j,k) ) * odx60
      d1_hh(1,3,1) = (  hxz(i+3,j,k) - 9*hxz(i+2,j,k) + 45*hxz(i+1,j,k)      &
                      - hxz(i-3,j,k) + 9*hxz(i-2,j,k) - 45*hxz(i-1,j,k) ) * odx60
      d1_hh(2,2,1) = (  hyy(i+3,j,k) - 9*hyy(i+2,j,k) + 45*hyy(i+1,j,k)      &
                      - hyy(i-3,j,k) + 9*hyy(i-2,j,k) - 45*hyy(i-1,j,k) ) * odx60
      d1_hh(2,3,1) = (  hyz(i+3,j,k) - 9*hyz(i+2,j,k) + 45*hyz(i+1,j,k)      &
                      - hyz(i-3,j,k) + 9*hyz(i-2,j,k) - 45*hyz(i-1,j,k) ) * odx60
      d1_hh(3,3,1) = (  hzz(i+3,j,k) - 9*hzz(i+2,j,k) + 45*hzz(i+1,j,k)      &
                      - hzz(i-3,j,k) + 9*hzz(i-2,j,k) - 45*hzz(i-1,j,k) ) * odx60

      d1_hh(1,1,2) = (  hxx(i,j+3,k) - 9*hxx(i,j+2,k) + 45*hxx(i,j+1,k)      &
                      - hxx(i,j-3,k) + 9*hxx(i,j-2,k) - 45*hxx(i,j-1,k) ) * ody60
      d1_hh(1,2,2) = (  hxy(i,j+3,k) - 9*hxy(i,j+2,k) + 45*hxy(i,j+1,k)      &
                      - hxy(i,j-3,k) + 9*hxy(i,j-2,k) - 45*hxy(i,j-1,k) ) * ody60
      d1_hh(1,3,2) = (  hxz(i,j+3,k) - 9*hxz(i,j+2,k) + 45*hxz(i,j+1,k)      &
                      - hxz(i,j-3,k) + 9*hxz(i,j-2,k) - 45*hxz(i,j-1,k) ) * ody60
      d1_hh(2,2,2) = (  hyy(i,j+3,k) - 9*hyy(i,j+2,k) + 45*hyy(i,j+1,k)      &
                      - hyy(i,j-3,k) + 9*hyy(i,j-2,k) - 45*hyy(i,j-1,k) ) * ody60
      d1_hh(2,3,2) = (  hyz(i,j+3,k) - 9*hyz(i,j+2,k) + 45*hyz(i,j+1,k)      &
                      - hyz(i,j-3,k) + 9*hyz(i,j-2,k) - 45*hyz(i,j-1,k) ) * ody60
      d1_hh(3,3,2) = (  hzz(i,j+3,k) - 9*hzz(i,j+2,k) + 45*hzz(i,j+1,k)      &
                      - hzz(i,j-3,k) + 9*hzz(i,j-2,k) - 45*hzz(i,j-1,k) ) * ody60

      d1_hh(1,1,3) = (  hxx(i,j,k+3) - 9*hxx(i,j,k+2) + 45*hxx(i,j,k+1)      &
                      - hxx(i,j,k-3) + 9*hxx(i,j,k-2) - 45*hxx(i,j,k-1) ) * odz60
      d1_hh(1,2,3) = (  hxy(i,j,k+3) - 9*hxy(i,j,k+2) + 45*hxy(i,j,k+1)      &
                      - hxy(i,j,k-3) + 9*hxy(i,j,k-2) - 45*hxy(i,j,k-1) ) * odz60
      d1_hh(1,3,3) = (  hxz(i,j,k+3) - 9*hxz(i,j,k+2) + 45*hxz(i,j,k+1)      &
                      - hxz(i,j,k-3) + 9*hxz(i,j,k-2) - 45*hxz(i,j,k-1) ) * odz60
      d1_hh(2,2,3) = (  hyy(i,j,k+3) - 9*hyy(i,j,k+2) + 45*hyy(i,j,k+1)      &
                      - hyy(i,j,k-3) + 9*hyy(i,j,k-2) - 45*hyy(i,j,k-1) ) * odz60
      d1_hh(2,3,3) = (  hyz(i,j,k+3) - 9*hyz(i,j,k+2) + 45*hyz(i,j,k+1)      &
                      - hyz(i,j,k-3) + 9*hyz(i,j,k-2) - 45*hyz(i,j,k-1) ) * odz60
      d1_hh(3,3,3) = (  hzz(i,j,k+3) - 9*hzz(i,j,k+2) + 45*hzz(i,j,k+1)      &
                      - hzz(i,j,k-3) + 9*hzz(i,j,k-2) - 45*hzz(i,j,k-1) ) * odz60

      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)

      ! d1_alph(3)
      d1_alph(1) = (  alp(i+3,j,k) - 9*alp(i+2,j,k) + 45*alp(i+1,j,k) &
                    - alp(i-3,j,k) + 9*alp(i-2,j,k) - 45*alp(i-1,j,k) ) * odx60

      d1_alph(2) = (  alp(i,j+3,k) - 9*alp(i,j+2,k) + 45*alp(i,j+1,k) &
                    - alp(i,j-3,k) + 9*alp(i,j-2,k) - 45*alp(i,j-1,k) ) * ody60

      d1_alph(3) = (  alp(i,j,k+3) - 9*alp(i,j,k+2) + 45*alp(i,j,k+1) &
                    - alp(i,j,k-3) + 9*alp(i,j,k-2) - 45*alp(i,j,k-1) ) * odz60

      ! d1_beta(3,3)
      d1_beta(1,1) = (  betax(i+3,j,k) - 9*betax(i+2,j,k) + 45*betax(i+1,j,k) &
                      - betax(i-3,j,k) + 9*betax(i-2,j,k) - 45*betax(i-1,j,k) ) * odx60
      d1_beta(2,1) = (  betay(i+3,j,k) - 9*betay(i+2,j,k) + 45*betay(i+1,j,k) &
                      - betay(i-3,j,k) + 9*betay(i-2,j,k) - 45*betay(i-1,j,k) ) * odx60
      d1_beta(3,1) = (  betaz(i+3,j,k) - 9*betaz(i+2,j,k) + 45*betaz(i+1,j,k) &
                      - betaz(i-3,j,k) + 9*betaz(i-2,j,k) - 45*betaz(i-1,j,k) ) * odx60

      d1_beta(1,2) = (  betax(i,j+3,k) - 9*betax(i,j+2,k) + 45*betax(i,j+1,k) &
                      - betax(i,j-3,k) + 9*betax(i,j-2,k) - 45*betax(i,j-1,k) ) * ody60
      d1_beta(2,2) = (  betay(i,j+3,k) - 9*betay(i,j+2,k) + 45*betay(i,j+1,k) &
                      - betay(i,j-3,k) + 9*betay(i,j-2,k) - 45*betay(i,j-1,k) ) * ody60
      d1_beta(3,2) = (  betaz(i,j+3,k) - 9*betaz(i,j+2,k) + 45*betaz(i,j+1,k) &
                      - betaz(i,j-3,k) + 9*betaz(i,j-2,k) - 45*betaz(i,j-1,k) ) * ody60

      d1_beta(1,3) = (  betax(i,j,k+3) - 9*betax(i,j,k+2) + 45*betax(i,j,k+1) &
                      - betax(i,j,k-3) + 9*betax(i,j,k-2) - 45*betax(i,j,k-1) ) * odz60
      d1_beta(2,3) = (  betay(i,j,k+3) - 9*betay(i,j,k+2) + 45*betay(i,j,k+1) &
                      - betay(i,j,k-3) + 9*betay(i,j,k-2) - 45*betay(i,j,k-1) ) * odz60
      d1_beta(3,3) = (  betaz(i,j,k+3) - 9*betaz(i,j,k+2) + 45*betaz(i,j,k+1) &
                      - betaz(i,j,k-3) + 9*betaz(i,j,k-2) - 45*betaz(i,j,k-1) ) * odz60

      ! d1_lE1(3,3)
      d1_lE1(1,1) = (  E1x(i+3,j,k) - 9*E1x(i+2,j,k) + 45*E1x(i+1,j,k) &
                    - E1x(i-3,j,k) + 9*E1x(i-2,j,k) - 45*E1x(i-1,j,k) ) * odx60
      d1_lE1(2,1) = (  E1y(i+3,j,k) - 9*E1y(i+2,j,k) + 45*E1y(i+1,j,k) &
                    - E1y(i-3,j,k) + 9*E1y(i-2,j,k) - 45*E1y(i-1,j,k) ) * odx60
      d1_lE1(3,1) = (  E1z(i+3,j,k) - 9*E1z(i+2,j,k) + 45*E1z(i+1,j,k) &
                    - E1z(i-3,j,k) + 9*E1z(i-2,j,k) - 45*E1z(i-1,j,k) ) * odx60

      d1_lE1(1,2) = (  E1x(i,j+3,k) - 9*E1x(i,j+2,k) + 45*E1x(i,j+1,k) &
                    - E1x(i,j-3,k) + 9*E1x(i,j-2,k) - 45*E1x(i,j-1,k) ) * ody60
      d1_lE1(2,2) = (  E1y(i,j+3,k) - 9*E1y(i,j+2,k) + 45*E1y(i,j+1,k) &
                    - E1y(i,j-3,k) + 9*E1y(i,j-2,k) - 45*E1y(i,j-1,k) ) * ody60
      d1_lE1(3,2) = (  E1z(i,j+3,k) - 9*E1z(i,j+2,k) + 45*E1z(i,j+1,k) &
                    - E1z(i,j-3,k) + 9*E1z(i,j-2,k) - 45*E1z(i,j-1,k) ) * ody60

      d1_lE1(1,3) = (  E1x(i,j,k+3) - 9*E1x(i,j,k+2) + 45*E1x(i,j,k+1) &
                    - E1x(i,j,k-3) + 9*E1x(i,j,k-2) - 45*E1x(i,j,k-1) ) * odz60
      d1_lE1(2,3) = (  E1y(i,j,k+3) - 9*E1y(i,j,k+2) + 45*E1y(i,j,k+1) &
                    - E1y(i,j,k-3) + 9*E1y(i,j,k-2) - 45*E1y(i,j,k-1) ) * odz60
      d1_lE1(3,3) = (  E1z(i,j,k+3) - 9*E1z(i,j,k+2) + 45*E1z(i,j,k+1) &
                    - E1z(i,j,k-3) + 9*E1z(i,j,k-2) - 45*E1z(i,j,k-1) ) * odz60

      ! d1_lA1(3,3)
      d1_lA1(1,1) = (  A1x(i+3,j,k) - 9*A1x(i+2,j,k) + 45*A1x(i+1,j,k) &
                    - A1x(i-3,j,k) + 9*A1x(i-2,j,k) - 45*A1x(i-1,j,k) ) * odx60
      d1_lA1(2,1) = (  A1y(i+3,j,k) - 9*A1y(i+2,j,k) + 45*A1y(i+1,j,k) &
                    - A1y(i-3,j,k) + 9*A1y(i-2,j,k) - 45*A1y(i-1,j,k) ) * odx60
      d1_lA1(3,1) = (  A1z(i+3,j,k) - 9*A1z(i+2,j,k) + 45*A1z(i+1,j,k) &
                    - A1z(i-3,j,k) + 9*A1z(i-2,j,k) - 45*A1z(i-1,j,k) ) * odx60

      d1_lA1(1,2) = (  A1x(i,j+3,k) - 9*A1x(i,j+2,k) + 45*A1x(i,j+1,k) &
                    - A1x(i,j-3,k) + 9*A1x(i,j-2,k) - 45*A1x(i,j-1,k) ) * ody60
      d1_lA1(2,2) = (  A1y(i,j+3,k) - 9*A1y(i,j+2,k) + 45*A1y(i,j+1,k) &
                    - A1y(i,j-3,k) + 9*A1y(i,j-2,k) - 45*A1y(i,j-1,k) ) * ody60
      d1_lA1(3,2) = (  A1z(i,j+3,k) - 9*A1z(i,j+2,k) + 45*A1z(i,j+1,k) &
                    - A1z(i,j-3,k) + 9*A1z(i,j-2,k) - 45*A1z(i,j-1,k) ) * ody60

      d1_lA1(1,3) = (  A1x(i,j,k+3) - 9*A1x(i,j,k+2) + 45*A1x(i,j,k+1) &
                    - A1x(i,j,k-3) + 9*A1x(i,j,k-2) - 45*A1x(i,j,k-1) ) * odz60
      d1_lA1(2,3) = (  A1y(i,j,k+3) - 9*A1y(i,j,k+2) + 45*A1y(i,j,k+1) &
                    - A1y(i,j,k-3) + 9*A1y(i,j,k-2) - 45*A1y(i,j,k-1) ) * odz60
      d1_lA1(3,3) = (  A1z(i,j,k+3) - 9*A1z(i,j,k+2) + 45*A1z(i,j,k+1) &
                    - A1z(i,j,k-3) + 9*A1z(i,j,k-2) - 45*A1z(i,j,k-1) ) * odz60

      ! d1_lZeta1(3)
      d1_lZeta1(1) = (  Zeta1(i+3,j,k) - 9*Zeta1(i+2,j,k) + 45*Zeta1(i+1,j,k) &
                     - Zeta1(i-3,j,k) + 9*Zeta1(i-2,j,k) - 45*Zeta1(i-1,j,k) ) * odx60

      d1_lZeta1(2) = (  Zeta1(i,j+3,k) - 9*Zeta1(i,j+2,k) + 45*Zeta1(i,j+1,k) &
                     - Zeta1(i,j-3,k) + 9*Zeta1(i,j-2,k) - 45*Zeta1(i,j-1,k) ) * ody60

      d1_lZeta1(3) = (  Zeta1(i,j,k+3) - 9*Zeta1(i,j,k+2) + 45*Zeta1(i,j,k+1) &
                     - Zeta1(i,j,k-3) + 9*Zeta1(i,j,k-2) - 45*Zeta1(i,j,k-1) ) * odz60

      ! d1_lAphi1(3)
      d1_lAphi1(1) = (  Aphi1(i+3,j,k) - 9*Aphi1(i+2,j,k) + 45*Aphi1(i+1,j,k) &
                     - Aphi1(i-3,j,k) + 9*Aphi1(i-2,j,k) - 45*Aphi1(i-1,j,k) ) * odx60

      d1_lAphi1(2) = (  Aphi1(i,j+3,k) - 9*Aphi1(i,j+2,k) + 45*Aphi1(i,j+1,k) &
                     - Aphi1(i,j-3,k) + 9*Aphi1(i,j-2,k) - 45*Aphi1(i,j-1,k) ) * ody60

      d1_lAphi1(3) = (  Aphi1(i,j,k+3) - 9*Aphi1(i,j,k+2) + 45*Aphi1(i,j,k+1) &
                     - Aphi1(i,j,k-3) + 9*Aphi1(i,j,k-2) - 45*Aphi1(i,j,k-1) ) * odz60

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_A(3,3,3)
      d2_lA1(1,1,1) = (  2*A1x(i+3,j,k) - 27*A1x(i+2,j,k) + 270*A1x(i+1,j,k) - 490*A1x(i,j,k)&
                      + 2*A1x(i-3,j,k) - 27*A1x(i-2,j,k) + 270*A1x(i-1,j,k) ) * odxsq180
      d2_lA1(2,1,1) = (  2*A1y(i+3,j,k) - 27*A1y(i+2,j,k) + 270*A1y(i+1,j,k) - 490*A1y(i,j,k)&
                      + 2*A1y(i-3,j,k) - 27*A1y(i-2,j,k) + 270*A1y(i-1,j,k) ) * odxsq180
      d2_lA1(3,1,1) = (  2*A1z(i+3,j,k) - 27*A1z(i+2,j,k) + 270*A1z(i+1,j,k) - 490*A1z(i,j,k)&
                      + 2*A1z(i-3,j,k) - 27*A1z(i-2,j,k) + 270*A1z(i-1,j,k) ) * odxsq180

      d2_lA1(1,2,2) = (  2*A1x(i,j+3,k) - 27*A1x(i,j+2,k) + 270*A1x(i,j+1,k) - 490*A1x(i,j,k)&
                      + 2*A1x(i,j-3,k) - 27*A1x(i,j-2,k) + 270*A1x(i,j-1,k) ) * odysq180
      d2_lA1(2,2,2) = (  2*A1y(i,j+3,k) - 27*A1y(i,j+2,k) + 270*A1y(i,j+1,k) - 490*A1y(i,j,k)&
                      + 2*A1y(i,j-3,k) - 27*A1y(i,j-2,k) + 270*A1y(i,j-1,k) ) * odysq180
      d2_lA1(3,2,2) = (  2*A1z(i,j+3,k) - 27*A1z(i,j+2,k) + 270*A1z(i,j+1,k) - 490*A1z(i,j,k)&
                      + 2*A1z(i,j-3,k) - 27*A1z(i,j-2,k) + 270*A1z(i,j-1,k) ) * odysq180

      d2_lA1(1,3,3) = (  2*A1x(i,j,k+3) - 27*A1x(i,j,k+2) + 270*A1x(i,j,k+1) - 490*A1x(i,j,k)&
                      + 2*A1x(i,j,k-3) - 27*A1x(i,j,k-2) + 270*A1x(i,j,k-1) ) * odzsq180
      d2_lA1(2,3,3) = (  2*A1y(i,j,k+3) - 27*A1y(i,j,k+2) + 270*A1y(i,j,k+1) - 490*A1y(i,j,k)&
                      + 2*A1y(i,j,k-3) - 27*A1y(i,j,k-2) + 270*A1y(i,j,k-1) ) * odzsq180
      d2_lA1(3,3,3) = (  2*A1z(i,j,k+3) - 27*A1z(i,j,k+2) + 270*A1z(i,j,k+1) - 490*A1z(i,j,k)&
                      + 2*A1z(i,j,k-3) - 27*A1z(i,j,k-2) + 270*A1z(i,j,k-1) ) * odzsq180

      d2_lA1(1,1,2) = (    -A1x(i-3,j+3,k) +   9*A1x(i-2,j+3,k) -   45*A1x(i-1,j+3,k) +   45*A1x(i+1,j+3,k) -   9*A1x(i+2,j+3,k) +    A1x(i+3,j+3,k) &
                      +  9*A1x(i-3,j+2,k) -  81*A1x(i-2,j+2,k) +  405*A1x(i-1,j+2,k) -  405*A1x(i+1,j+2,k) +  81*A1x(i+2,j+2,k) -  9*A1x(i+3,j+2,k) &
                      - 45*A1x(i-3,j+1,k) + 405*A1x(i-2,j+1,k) - 2025*A1x(i-1,j+1,k) + 2025*A1x(i+1,j+1,k) - 405*A1x(i+2,j+1,k) + 45*A1x(i+3,j+1,k) &
                      + 45*A1x(i-3,j-1,k) - 405*A1x(i-2,j-1,k) + 2025*A1x(i-1,j-1,k) - 2025*A1x(i+1,j-1,k) + 405*A1x(i+2,j-1,k) - 45*A1x(i+3,j-1,k) &
                      -  9*A1x(i-3,j-2,k) +  81*A1x(i-2,j-2,k) -  405*A1x(i-1,j-2,k) +  405*A1x(i+1,j-2,k) -  81*A1x(i+2,j-2,k) +  9*A1x(i+3,j-2,k) &
                      +    A1x(i-3,j-3,k) -   9*A1x(i-2,j-3,k) +   45*A1x(i-1,j-3,k) -   45*A1x(i+1,j-3,k) +   9*A1x(i+2,j-3,k) -    A1x(i+3,j-3,k) ) * odxdy3600
      d2_lA1(2,1,2) = (    -A1y(i-3,j+3,k) +   9*A1y(i-2,j+3,k) -   45*A1y(i-1,j+3,k) +   45*A1y(i+1,j+3,k) -   9*A1y(i+2,j+3,k) +    A1y(i+3,j+3,k) &
                      +  9*A1y(i-3,j+2,k) -  81*A1y(i-2,j+2,k) +  405*A1y(i-1,j+2,k) -  405*A1y(i+1,j+2,k) +  81*A1y(i+2,j+2,k) -  9*A1y(i+3,j+2,k) &
                      - 45*A1y(i-3,j+1,k) + 405*A1y(i-2,j+1,k) - 2025*A1y(i-1,j+1,k) + 2025*A1y(i+1,j+1,k) - 405*A1y(i+2,j+1,k) + 45*A1y(i+3,j+1,k) &
                      + 45*A1y(i-3,j-1,k) - 405*A1y(i-2,j-1,k) + 2025*A1y(i-1,j-1,k) - 2025*A1y(i+1,j-1,k) + 405*A1y(i+2,j-1,k) - 45*A1y(i+3,j-1,k) &
                      -  9*A1y(i-3,j-2,k) +  81*A1y(i-2,j-2,k) -  405*A1y(i-1,j-2,k) +  405*A1y(i+1,j-2,k) -  81*A1y(i+2,j-2,k) +  9*A1y(i+3,j-2,k) &
                      +    A1y(i-3,j-3,k) -   9*A1y(i-2,j-3,k) +   45*A1y(i-1,j-3,k) -   45*A1y(i+1,j-3,k) +   9*A1y(i+2,j-3,k) -    A1y(i+3,j-3,k) ) * odxdy3600
      d2_lA1(3,1,2) = (    -A1z(i-3,j+3,k) +   9*A1z(i-2,j+3,k) -   45*A1z(i-1,j+3,k) +   45*A1z(i+1,j+3,k) -   9*A1z(i+2,j+3,k) +    A1z(i+3,j+3,k) &
                      +  9*A1z(i-3,j+2,k) -  81*A1z(i-2,j+2,k) +  405*A1z(i-1,j+2,k) -  405*A1z(i+1,j+2,k) +  81*A1z(i+2,j+2,k) -  9*A1z(i+3,j+2,k) &
                      - 45*A1z(i-3,j+1,k) + 405*A1z(i-2,j+1,k) - 2025*A1z(i-1,j+1,k) + 2025*A1z(i+1,j+1,k) - 405*A1z(i+2,j+1,k) + 45*A1z(i+3,j+1,k) &
                      + 45*A1z(i-3,j-1,k) - 405*A1z(i-2,j-1,k) + 2025*A1z(i-1,j-1,k) - 2025*A1z(i+1,j-1,k) + 405*A1z(i+2,j-1,k) - 45*A1z(i+3,j-1,k) &
                      -  9*A1z(i-3,j-2,k) +  81*A1z(i-2,j-2,k) -  405*A1z(i-1,j-2,k) +  405*A1z(i+1,j-2,k) -  81*A1z(i+2,j-2,k) +  9*A1z(i+3,j-2,k) &
                      +    A1z(i-3,j-3,k) -   9*A1z(i-2,j-3,k) +   45*A1z(i-1,j-3,k) -   45*A1z(i+1,j-3,k) +   9*A1z(i+2,j-3,k) -    A1z(i+3,j-3,k) ) * odxdy3600


      d2_lA1(1,1,3) = (    -A1x(i-3,j,k+3) +   9*A1x(i-2,j,k+3) -   45*A1x(i-1,j,k+3) +   45*A1x(i+1,j,k+3) -   9*A1x(i+2,j,k+3) +    A1x(i+3,j,k+3) &
                      +  9*A1x(i-3,j,k+2) -  81*A1x(i-2,j,k+2) +  405*A1x(i-1,j,k+2) -  405*A1x(i+1,j,k+2) +  81*A1x(i+2,j,k+2) -  9*A1x(i+3,j,k+2) &
                      - 45*A1x(i-3,j,k+1) + 405*A1x(i-2,j,k+1) - 2025*A1x(i-1,j,k+1) + 2025*A1x(i+1,j,k+1) - 405*A1x(i+2,j,k+1) + 45*A1x(i+3,j,k+1) &
                      + 45*A1x(i-3,j,k-1) - 405*A1x(i-2,j,k-1) + 2025*A1x(i-1,j,k-1) - 2025*A1x(i+1,j,k-1) + 405*A1x(i+2,j,k-1) - 45*A1x(i+3,j,k-1) &
                      -  9*A1x(i-3,j,k-2) +  81*A1x(i-2,j,k-2) -  405*A1x(i-1,j,k-2) +  405*A1x(i+1,j,k-2) -  81*A1x(i+2,j,k-2) +  9*A1x(i+3,j,k-2) &
                      +    A1x(i-3,j,k-3) -   9*A1x(i-2,j,k-3) +   45*A1x(i-1,j,k-3) -   45*A1x(i+1,j,k-3) +   9*A1x(i+2,j,k-3) -    A1x(i+3,j,k-3) ) * odxdz3600
      d2_lA1(2,1,3) = (    -A1y(i-3,j,k+3) +   9*A1y(i-2,j,k+3) -   45*A1y(i-1,j,k+3) +   45*A1y(i+1,j,k+3) -   9*A1y(i+2,j,k+3) +    A1y(i+3,j,k+3) &
                      +  9*A1y(i-3,j,k+2) -  81*A1y(i-2,j,k+2) +  405*A1y(i-1,j,k+2) -  405*A1y(i+1,j,k+2) +  81*A1y(i+2,j,k+2) -  9*A1y(i+3,j,k+2) &
                      - 45*A1y(i-3,j,k+1) + 405*A1y(i-2,j,k+1) - 2025*A1y(i-1,j,k+1) + 2025*A1y(i+1,j,k+1) - 405*A1y(i+2,j,k+1) + 45*A1y(i+3,j,k+1) &
                      + 45*A1y(i-3,j,k-1) - 405*A1y(i-2,j,k-1) + 2025*A1y(i-1,j,k-1) - 2025*A1y(i+1,j,k-1) + 405*A1y(i+2,j,k-1) - 45*A1y(i+3,j,k-1) &
                      -  9*A1y(i-3,j,k-2) +  81*A1y(i-2,j,k-2) -  405*A1y(i-1,j,k-2) +  405*A1y(i+1,j,k-2) -  81*A1y(i+2,j,k-2) +  9*A1y(i+3,j,k-2) &
                      +    A1y(i-3,j,k-3) -   9*A1y(i-2,j,k-3) +   45*A1y(i-1,j,k-3) -   45*A1y(i+1,j,k-3) +   9*A1y(i+2,j,k-3) -    A1y(i+3,j,k-3) ) * odxdz3600
      d2_lA1(3,1,3) = (    -A1z(i-3,j,k+3) +   9*A1z(i-2,j,k+3) -   45*A1z(i-1,j,k+3) +   45*A1z(i+1,j,k+3) -   9*A1z(i+2,j,k+3) +    A1z(i+3,j,k+3) &
                      +  9*A1z(i-3,j,k+2) -  81*A1z(i-2,j,k+2) +  405*A1z(i-1,j,k+2) -  405*A1z(i+1,j,k+2) +  81*A1z(i+2,j,k+2) -  9*A1z(i+3,j,k+2) &
                      - 45*A1z(i-3,j,k+1) + 405*A1z(i-2,j,k+1) - 2025*A1z(i-1,j,k+1) + 2025*A1z(i+1,j,k+1) - 405*A1z(i+2,j,k+1) + 45*A1z(i+3,j,k+1) &
                      + 45*A1z(i-3,j,k-1) - 405*A1z(i-2,j,k-1) + 2025*A1z(i-1,j,k-1) - 2025*A1z(i+1,j,k-1) + 405*A1z(i+2,j,k-1) - 45*A1z(i+3,j,k-1) &
                      -  9*A1z(i-3,j,k-2) +  81*A1z(i-2,j,k-2) -  405*A1z(i-1,j,k-2) +  405*A1z(i+1,j,k-2) -  81*A1z(i+2,j,k-2) +  9*A1z(i+3,j,k-2) &
                      +    A1z(i-3,j,k-3) -   9*A1z(i-2,j,k-3) +   45*A1z(i-1,j,k-3) -   45*A1z(i+1,j,k-3) +   9*A1z(i+2,j,k-3) -    A1z(i+3,j,k-3) ) * odxdz3600

      d2_lA1(1,2,3) = (    -A1x(i,j-3,k+3) +   9*A1x(i,j-2,k+3) -   45*A1x(i,j-1,k+3) +   45*A1x(i,j+1,k+3) -   9*A1x(i,j+2,k+3) +    A1x(i,j+3,k+3) &
                      +  9*A1x(i,j-3,k+2) -  81*A1x(i,j-2,k+2) +  405*A1x(i,j-1,k+2) -  405*A1x(i,j+1,k+2) +  81*A1x(i,j+2,k+2) -  9*A1x(i,j+3,k+2) &
                      - 45*A1x(i,j-3,k+1) + 405*A1x(i,j-2,k+1) - 2025*A1x(i,j-1,k+1) + 2025*A1x(i,j+1,k+1) - 405*A1x(i,j+2,k+1) + 45*A1x(i,j+3,k+1) &
                      + 45*A1x(i,j-3,k-1) - 405*A1x(i,j-2,k-1) + 2025*A1x(i,j-1,k-1) - 2025*A1x(i,j+1,k-1) + 405*A1x(i,j+2,k-1) - 45*A1x(i,j+3,k-1) &
                      -  9*A1x(i,j-3,k-2) +  81*A1x(i,j-2,k-2) -  405*A1x(i,j-1,k-2) +  405*A1x(i,j+1,k-2) -  81*A1x(i,j+2,k-2) +  9*A1x(i,j+3,k-2) &
                      +    A1x(i,j-3,k-3) -   9*A1x(i,j-2,k-3) +   45*A1x(i,j-1,k-3) -   45*A1x(i,j+1,k-3) +   9*A1x(i,j+2,k-3) -    A1x(i,j+3,k-3) ) * odydz3600
      d2_lA1(2,2,3) = (    -A1y(i,j-3,k+3) +   9*A1y(i,j-2,k+3) -   45*A1y(i,j-1,k+3) +   45*A1y(i,j+1,k+3) -   9*A1y(i,j+2,k+3) +    A1y(i,j+3,k+3) &
                      +  9*A1y(i,j-3,k+2) -  81*A1y(i,j-2,k+2) +  405*A1y(i,j-1,k+2) -  405*A1y(i,j+1,k+2) +  81*A1y(i,j+2,k+2) -  9*A1y(i,j+3,k+2) &
                      - 45*A1y(i,j-3,k+1) + 405*A1y(i,j-2,k+1) - 2025*A1y(i,j-1,k+1) + 2025*A1y(i,j+1,k+1) - 405*A1y(i,j+2,k+1) + 45*A1y(i,j+3,k+1) &
                      + 45*A1y(i,j-3,k-1) - 405*A1y(i,j-2,k-1) + 2025*A1y(i,j-1,k-1) - 2025*A1y(i,j+1,k-1) + 405*A1y(i,j+2,k-1) - 45*A1y(i,j+3,k-1) &
                      -  9*A1y(i,j-3,k-2) +  81*A1y(i,j-2,k-2) -  405*A1y(i,j-1,k-2) +  405*A1y(i,j+1,k-2) -  81*A1y(i,j+2,k-2) +  9*A1y(i,j+3,k-2) &
                     +    A1y(i,j-3,k-3) -   9*A1y(i,j-2,k-3) +   45*A1y(i,j-1,k-3) -   45*A1y(i,j+1,k-3) +   9*A1y(i,j+2,k-3) -    A1y(i,j+3,k-3) ) * odydz3600
      d2_lA1(3,2,3) = (    -A1z(i,j-3,k+3) +   9*A1z(i,j-2,k+3) -   45*A1z(i,j-1,k+3) +   45*A1z(i,j+1,k+3) -   9*A1z(i,j+2,k+3) +    A1z(i,j+3,k+3) &
                      +  9*A1z(i,j-3,k+2) -  81*A1z(i,j-2,k+2) +  405*A1z(i,j-1,k+2) -  405*A1z(i,j+1,k+2) +  81*A1z(i,j+2,k+2) -  9*A1z(i,j+3,k+2) &
                      - 45*A1z(i,j-3,k+1) + 405*A1z(i,j-2,k+1) - 2025*A1z(i,j-1,k+1) + 2025*A1z(i,j+1,k+1) - 405*A1z(i,j+2,k+1) + 45*A1z(i,j+3,k+1) &
                      + 45*A1z(i,j-3,k-1) - 405*A1z(i,j-2,k-1) + 2025*A1z(i,j-1,k-1) - 2025*A1z(i,j+1,k-1) + 405*A1z(i,j+2,k-1) - 45*A1z(i,j+3,k-1) &
                      -  9*A1z(i,j-3,k-2) +  81*A1z(i,j-2,k-2) -  405*A1z(i,j-1,k-2) +  405*A1z(i,j+1,k-2) -  81*A1z(i,j+2,k-2) +  9*A1z(i,j+3,k-2) &
                      +    A1z(i,j-3,k-3) -   9*A1z(i,j-2,k-3) +   45*A1z(i,j-1,k-3) -   45*A1z(i,j+1,k-3) +   9*A1z(i,j+2,k-3) -    A1z(i,j+3,k-3) ) * odydz3600

      d2_lA1(:,2,1) = d2_lA1(:,1,2)
      d2_lA1(:,3,1) = d2_lA1(:,1,3)
      d2_lA1(:,3,2) = d2_lA1(:,2,3)
      
      
      
      ! d1_lE2(3,3)
      d1_lE2(1,1) = (  E2x(i+3,j,k) - 9*E2x(i+2,j,k) + 45*E2x(i+1,j,k) &
                    - E2x(i-3,j,k) + 9*E2x(i-2,j,k) - 45*E2x(i-1,j,k) ) * odx60
      d1_lE2(2,1) = (  E2y(i+3,j,k) - 9*E2y(i+2,j,k) + 45*E2y(i+1,j,k) &
                    - E2y(i-3,j,k) + 9*E2y(i-2,j,k) - 45*E2y(i-1,j,k) ) * odx60
      d1_lE2(3,1) = (  E2z(i+3,j,k) - 9*E2z(i+2,j,k) + 45*E2z(i+1,j,k) &
                    - E2z(i-3,j,k) + 9*E2z(i-2,j,k) - 45*E2z(i-1,j,k) ) * odx60

      d1_lE2(1,2) = (  E2x(i,j+3,k) - 9*E2x(i,j+2,k) + 45*E2x(i,j+1,k) &
                    - E2x(i,j-3,k) + 9*E2x(i,j-2,k) - 45*E2x(i,j-1,k) ) * ody60
      d1_lE2(2,2) = (  E2y(i,j+3,k) - 9*E2y(i,j+2,k) + 45*E2y(i,j+1,k) &
                    - E2y(i,j-3,k) + 9*E2y(i,j-2,k) - 45*E2y(i,j-1,k) ) * ody60
      d1_lE2(3,2) = (  E2z(i,j+3,k) - 9*E2z(i,j+2,k) + 45*E2z(i,j+1,k) &
                    - E2z(i,j-3,k) + 9*E2z(i,j-2,k) - 45*E2z(i,j-1,k) ) * ody60

      d1_lE2(1,3) = (  E2x(i,j,k+3) - 9*E2x(i,j,k+2) + 45*E2x(i,j,k+1) &
                    - E2x(i,j,k-3) + 9*E2x(i,j,k-2) - 45*E2x(i,j,k-1) ) * odz60
      d1_lE2(2,3) = (  E2y(i,j,k+3) - 9*E2y(i,j,k+2) + 45*E2y(i,j,k+1) &
                    - E2y(i,j,k-3) + 9*E2y(i,j,k-2) - 45*E2y(i,j,k-1) ) * odz60
      d1_lE2(3,3) = (  E2z(i,j,k+3) - 9*E2z(i,j,k+2) + 45*E2z(i,j,k+1) &
                    - E2z(i,j,k-3) + 9*E2z(i,j,k-2) - 45*E2z(i,j,k-1) ) * odz60

      ! d1_lA2(3,3)
      d1_lA2(1,1) = (  A2x(i+3,j,k) - 9*A2x(i+2,j,k) + 45*A2x(i+1,j,k) &
                    - A2x(i-3,j,k) + 9*A2x(i-2,j,k) - 45*A2x(i-1,j,k) ) * odx60
      d1_lA2(2,1) = (  A2y(i+3,j,k) - 9*A2y(i+2,j,k) + 45*A2y(i+1,j,k) &
                    - A2y(i-3,j,k) + 9*A2y(i-2,j,k) - 45*A2y(i-1,j,k) ) * odx60
      d1_lA2(3,1) = (  A2z(i+3,j,k) - 9*A2z(i+2,j,k) + 45*A2z(i+1,j,k) &
                    - A2z(i-3,j,k) + 9*A2z(i-2,j,k) - 45*A2z(i-1,j,k) ) * odx60

      d1_lA2(1,2) = (  A2x(i,j+3,k) - 9*A2x(i,j+2,k) + 45*A2x(i,j+1,k) &
                    - A2x(i,j-3,k) + 9*A2x(i,j-2,k) - 45*A2x(i,j-1,k) ) * ody60
      d1_lA2(2,2) = (  A2y(i,j+3,k) - 9*A2y(i,j+2,k) + 45*A2y(i,j+1,k) &
                    - A2y(i,j-3,k) + 9*A2y(i,j-2,k) - 45*A2y(i,j-1,k) ) * ody60
      d1_lA2(3,2) = (  A2z(i,j+3,k) - 9*A2z(i,j+2,k) + 45*A2z(i,j+1,k) &
                    - A2z(i,j-3,k) + 9*A2z(i,j-2,k) - 45*A2z(i,j-1,k) ) * ody60

      d1_lA2(1,3) = (  A2x(i,j,k+3) - 9*A2x(i,j,k+2) + 45*A2x(i,j,k+1) &
                    - A2x(i,j,k-3) + 9*A2x(i,j,k-2) - 45*A2x(i,j,k-1) ) * odz60
      d1_lA2(2,3) = (  A2y(i,j,k+3) - 9*A2y(i,j,k+2) + 45*A2y(i,j,k+1) &
                    - A2y(i,j,k-3) + 9*A2y(i,j,k-2) - 45*A2y(i,j,k-1) ) * odz60
      d1_lA2(3,3) = (  A2z(i,j,k+3) - 9*A2z(i,j,k+2) + 45*A2z(i,j,k+1) &
                    - A2z(i,j,k-3) + 9*A2z(i,j,k-2) - 45*A2z(i,j,k-1) ) * odz60

      ! d1_lZeta2(3)
      d1_lZeta2(1) = (  Zeta2(i+3,j,k) - 9*Zeta2(i+2,j,k) + 45*Zeta2(i+1,j,k) &
                     - Zeta2(i-3,j,k) + 9*Zeta2(i-2,j,k) - 45*Zeta2(i-1,j,k) ) * odx60

      d1_lZeta2(2) = (  Zeta2(i,j+3,k) - 9*Zeta2(i,j+2,k) + 45*Zeta2(i,j+1,k) &
                     - Zeta2(i,j-3,k) + 9*Zeta2(i,j-2,k) - 45*Zeta2(i,j-1,k) ) * ody60

      d1_lZeta2(3) = (  Zeta2(i,j,k+3) - 9*Zeta2(i,j,k+2) + 45*Zeta2(i,j,k+1) &
                     - Zeta2(i,j,k-3) + 9*Zeta2(i,j,k-2) - 45*Zeta2(i,j,k-1) ) * odz60

      ! d1_lAphi2(3)
      d1_lAphi2(1) = (  Aphi2(i+3,j,k) - 9*Aphi2(i+2,j,k) + 45*Aphi2(i+1,j,k) &
                     - Aphi2(i-3,j,k) + 9*Aphi2(i-2,j,k) - 45*Aphi2(i-1,j,k) ) * odx60

      d1_lAphi2(2) = (  Aphi2(i,j+3,k) - 9*Aphi2(i,j+2,k) + 45*Aphi2(i,j+1,k) &
                     - Aphi2(i,j-3,k) + 9*Aphi2(i,j-2,k) - 45*Aphi2(i,j-1,k) ) * ody60

      d1_lAphi2(3) = (  Aphi2(i,j,k+3) - 9*Aphi2(i,j,k+2) + 45*Aphi2(i,j,k+1) &
                     - Aphi2(i,j,k-3) + 9*Aphi2(i,j,k-2) - 45*Aphi2(i,j,k-1) ) * odz60

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_A(3,3,3)
      d2_lA2(1,1,1) = (  2*A2x(i+3,j,k) - 27*A2x(i+2,j,k) + 270*A2x(i+1,j,k) - 490*A2x(i,j,k)&
                      + 2*A2x(i-3,j,k) - 27*A2x(i-2,j,k) + 270*A2x(i-1,j,k) ) * odxsq180
      d2_lA2(2,1,1) = (  2*A2y(i+3,j,k) - 27*A2y(i+2,j,k) + 270*A2y(i+1,j,k) - 490*A2y(i,j,k)&
                      + 2*A2y(i-3,j,k) - 27*A2y(i-2,j,k) + 270*A2y(i-1,j,k) ) * odxsq180
      d2_lA2(3,1,1) = (  2*A2z(i+3,j,k) - 27*A2z(i+2,j,k) + 270*A2z(i+1,j,k) - 490*A2z(i,j,k)&
                      + 2*A2z(i-3,j,k) - 27*A2z(i-2,j,k) + 270*A2z(i-1,j,k) ) * odxsq180

      d2_lA2(1,2,2) = (  2*A2x(i,j+3,k) - 27*A2x(i,j+2,k) + 270*A2x(i,j+1,k) - 490*A2x(i,j,k)&
                      + 2*A2x(i,j-3,k) - 27*A2x(i,j-2,k) + 270*A2x(i,j-1,k) ) * odysq180
      d2_lA2(2,2,2) = (  2*A2y(i,j+3,k) - 27*A2y(i,j+2,k) + 270*A2y(i,j+1,k) - 490*A2y(i,j,k)&
                      + 2*A2y(i,j-3,k) - 27*A2y(i,j-2,k) + 270*A2y(i,j-1,k) ) * odysq180
      d2_lA2(3,2,2) = (  2*A2z(i,j+3,k) - 27*A2z(i,j+2,k) + 270*A2z(i,j+1,k) - 490*A2z(i,j,k)&
                      + 2*A2z(i,j-3,k) - 27*A2z(i,j-2,k) + 270*A2z(i,j-1,k) ) * odysq180

      d2_lA2(1,3,3) = (  2*A2x(i,j,k+3) - 27*A2x(i,j,k+2) + 270*A2x(i,j,k+1) - 490*A2x(i,j,k)&
                      + 2*A2x(i,j,k-3) - 27*A2x(i,j,k-2) + 270*A2x(i,j,k-1) ) * odzsq180
      d2_lA2(2,3,3) = (  2*A2y(i,j,k+3) - 27*A2y(i,j,k+2) + 270*A2y(i,j,k+1) - 490*A2y(i,j,k)&
                      + 2*A2y(i,j,k-3) - 27*A2y(i,j,k-2) + 270*A2y(i,j,k-1) ) * odzsq180
      d2_lA2(3,3,3) = (  2*A2z(i,j,k+3) - 27*A2z(i,j,k+2) + 270*A2z(i,j,k+1) - 490*A2z(i,j,k)&
                      + 2*A2z(i,j,k-3) - 27*A2z(i,j,k-2) + 270*A2z(i,j,k-1) ) * odzsq180

      d2_lA2(1,1,2) = (    -A2x(i-3,j+3,k) +   9*A2x(i-2,j+3,k) -   45*A2x(i-1,j+3,k) +   45*A2x(i+1,j+3,k) -   9*A2x(i+2,j+3,k) +    A2x(i+3,j+3,k) &
                      +  9*A2x(i-3,j+2,k) -  81*A2x(i-2,j+2,k) +  405*A2x(i-1,j+2,k) -  405*A2x(i+1,j+2,k) +  81*A2x(i+2,j+2,k) -  9*A2x(i+3,j+2,k) &
                      - 45*A2x(i-3,j+1,k) + 405*A2x(i-2,j+1,k) - 2025*A2x(i-1,j+1,k) + 2025*A2x(i+1,j+1,k) - 405*A2x(i+2,j+1,k) + 45*A2x(i+3,j+1,k) &
                      + 45*A2x(i-3,j-1,k) - 405*A2x(i-2,j-1,k) + 2025*A2x(i-1,j-1,k) - 2025*A2x(i+1,j-1,k) + 405*A2x(i+2,j-1,k) - 45*A2x(i+3,j-1,k) &
                      -  9*A2x(i-3,j-2,k) +  81*A2x(i-2,j-2,k) -  405*A2x(i-1,j-2,k) +  405*A2x(i+1,j-2,k) -  81*A2x(i+2,j-2,k) +  9*A2x(i+3,j-2,k) &
                      +    A2x(i-3,j-3,k) -   9*A2x(i-2,j-3,k) +   45*A2x(i-1,j-3,k) -   45*A2x(i+1,j-3,k) +   9*A2x(i+2,j-3,k) -    A2x(i+3,j-3,k) ) * odxdy3600
      d2_lA2(2,1,2) = (    -A2y(i-3,j+3,k) +   9*A2y(i-2,j+3,k) -   45*A2y(i-1,j+3,k) +   45*A2y(i+1,j+3,k) -   9*A2y(i+2,j+3,k) +    A2y(i+3,j+3,k) &
                      +  9*A2y(i-3,j+2,k) -  81*A2y(i-2,j+2,k) +  405*A2y(i-1,j+2,k) -  405*A2y(i+1,j+2,k) +  81*A2y(i+2,j+2,k) -  9*A2y(i+3,j+2,k) &
                      - 45*A2y(i-3,j+1,k) + 405*A2y(i-2,j+1,k) - 2025*A2y(i-1,j+1,k) + 2025*A2y(i+1,j+1,k) - 405*A2y(i+2,j+1,k) + 45*A2y(i+3,j+1,k) &
                      + 45*A2y(i-3,j-1,k) - 405*A2y(i-2,j-1,k) + 2025*A2y(i-1,j-1,k) - 2025*A2y(i+1,j-1,k) + 405*A2y(i+2,j-1,k) - 45*A2y(i+3,j-1,k) &
                      -  9*A2y(i-3,j-2,k) +  81*A2y(i-2,j-2,k) -  405*A2y(i-1,j-2,k) +  405*A2y(i+1,j-2,k) -  81*A2y(i+2,j-2,k) +  9*A2y(i+3,j-2,k) &
                      +    A2y(i-3,j-3,k) -   9*A2y(i-2,j-3,k) +   45*A2y(i-1,j-3,k) -   45*A2y(i+1,j-3,k) +   9*A2y(i+2,j-3,k) -    A2y(i+3,j-3,k) ) * odxdy3600
      d2_lA2(3,1,2) = (    -A2z(i-3,j+3,k) +   9*A2z(i-2,j+3,k) -   45*A2z(i-1,j+3,k) +   45*A2z(i+1,j+3,k) -   9*A2z(i+2,j+3,k) +    A2z(i+3,j+3,k) &
                      +  9*A2z(i-3,j+2,k) -  81*A2z(i-2,j+2,k) +  405*A2z(i-1,j+2,k) -  405*A2z(i+1,j+2,k) +  81*A2z(i+2,j+2,k) -  9*A2z(i+3,j+2,k) &
                      - 45*A2z(i-3,j+1,k) + 405*A2z(i-2,j+1,k) - 2025*A2z(i-1,j+1,k) + 2025*A2z(i+1,j+1,k) - 405*A2z(i+2,j+1,k) + 45*A2z(i+3,j+1,k) &
                      + 45*A2z(i-3,j-1,k) - 405*A2z(i-2,j-1,k) + 2025*A2z(i-1,j-1,k) - 2025*A2z(i+1,j-1,k) + 405*A2z(i+2,j-1,k) - 45*A2z(i+3,j-1,k) &
                      -  9*A2z(i-3,j-2,k) +  81*A2z(i-2,j-2,k) -  405*A2z(i-1,j-2,k) +  405*A2z(i+1,j-2,k) -  81*A2z(i+2,j-2,k) +  9*A2z(i+3,j-2,k) &
                      +    A2z(i-3,j-3,k) -   9*A2z(i-2,j-3,k) +   45*A2z(i-1,j-3,k) -   45*A2z(i+1,j-3,k) +   9*A2z(i+2,j-3,k) -    A2z(i+3,j-3,k) ) * odxdy3600


      d2_lA2(1,1,3) = (    -A2x(i-3,j,k+3) +   9*A2x(i-2,j,k+3) -   45*A2x(i-1,j,k+3) +   45*A2x(i+1,j,k+3) -   9*A2x(i+2,j,k+3) +    A2x(i+3,j,k+3) &
                      +  9*A2x(i-3,j,k+2) -  81*A2x(i-2,j,k+2) +  405*A2x(i-1,j,k+2) -  405*A2x(i+1,j,k+2) +  81*A2x(i+2,j,k+2) -  9*A2x(i+3,j,k+2) &
                      - 45*A2x(i-3,j,k+1) + 405*A2x(i-2,j,k+1) - 2025*A2x(i-1,j,k+1) + 2025*A2x(i+1,j,k+1) - 405*A2x(i+2,j,k+1) + 45*A2x(i+3,j,k+1) &
                      + 45*A2x(i-3,j,k-1) - 405*A2x(i-2,j,k-1) + 2025*A2x(i-1,j,k-1) - 2025*A2x(i+1,j,k-1) + 405*A2x(i+2,j,k-1) - 45*A2x(i+3,j,k-1) &
                      -  9*A2x(i-3,j,k-2) +  81*A2x(i-2,j,k-2) -  405*A2x(i-1,j,k-2) +  405*A2x(i+1,j,k-2) -  81*A2x(i+2,j,k-2) +  9*A2x(i+3,j,k-2) &
                      +    A2x(i-3,j,k-3) -   9*A2x(i-2,j,k-3) +   45*A2x(i-1,j,k-3) -   45*A2x(i+1,j,k-3) +   9*A2x(i+2,j,k-3) -    A2x(i+3,j,k-3) ) * odxdz3600
      d2_lA2(2,1,3) = (    -A2y(i-3,j,k+3) +   9*A2y(i-2,j,k+3) -   45*A2y(i-1,j,k+3) +   45*A2y(i+1,j,k+3) -   9*A2y(i+2,j,k+3) +    A2y(i+3,j,k+3) &
                      +  9*A2y(i-3,j,k+2) -  81*A2y(i-2,j,k+2) +  405*A2y(i-1,j,k+2) -  405*A2y(i+1,j,k+2) +  81*A2y(i+2,j,k+2) -  9*A2y(i+3,j,k+2) &
                      - 45*A2y(i-3,j,k+1) + 405*A2y(i-2,j,k+1) - 2025*A2y(i-1,j,k+1) + 2025*A2y(i+1,j,k+1) - 405*A2y(i+2,j,k+1) + 45*A2y(i+3,j,k+1) &
                      + 45*A2y(i-3,j,k-1) - 405*A2y(i-2,j,k-1) + 2025*A2y(i-1,j,k-1) - 2025*A2y(i+1,j,k-1) + 405*A2y(i+2,j,k-1) - 45*A2y(i+3,j,k-1) &
                      -  9*A2y(i-3,j,k-2) +  81*A2y(i-2,j,k-2) -  405*A2y(i-1,j,k-2) +  405*A2y(i+1,j,k-2) -  81*A2y(i+2,j,k-2) +  9*A2y(i+3,j,k-2) &
                      +    A2y(i-3,j,k-3) -   9*A2y(i-2,j,k-3) +   45*A2y(i-1,j,k-3) -   45*A2y(i+1,j,k-3) +   9*A2y(i+2,j,k-3) -    A2y(i+3,j,k-3) ) * odxdz3600
      d2_lA2(3,1,3) = (    -A2z(i-3,j,k+3) +   9*A2z(i-2,j,k+3) -   45*A2z(i-1,j,k+3) +   45*A2z(i+1,j,k+3) -   9*A2z(i+2,j,k+3) +    A2z(i+3,j,k+3) &
                      +  9*A2z(i-3,j,k+2) -  81*A2z(i-2,j,k+2) +  405*A2z(i-1,j,k+2) -  405*A2z(i+1,j,k+2) +  81*A2z(i+2,j,k+2) -  9*A2z(i+3,j,k+2) &
                      - 45*A2z(i-3,j,k+1) + 405*A2z(i-2,j,k+1) - 2025*A2z(i-1,j,k+1) + 2025*A2z(i+1,j,k+1) - 405*A2z(i+2,j,k+1) + 45*A2z(i+3,j,k+1) &
                      + 45*A2z(i-3,j,k-1) - 405*A2z(i-2,j,k-1) + 2025*A2z(i-1,j,k-1) - 2025*A2z(i+1,j,k-1) + 405*A2z(i+2,j,k-1) - 45*A2z(i+3,j,k-1) &
                      -  9*A2z(i-3,j,k-2) +  81*A2z(i-2,j,k-2) -  405*A2z(i-1,j,k-2) +  405*A2z(i+1,j,k-2) -  81*A2z(i+2,j,k-2) +  9*A2z(i+3,j,k-2) &
                      +    A2z(i-3,j,k-3) -   9*A2z(i-2,j,k-3) +   45*A2z(i-1,j,k-3) -   45*A2z(i+1,j,k-3) +   9*A2z(i+2,j,k-3) -    A2z(i+3,j,k-3) ) * odxdz3600

      d2_lA2(1,2,3) = (    -A2x(i,j-3,k+3) +   9*A2x(i,j-2,k+3) -   45*A2x(i,j-1,k+3) +   45*A2x(i,j+1,k+3) -   9*A2x(i,j+2,k+3) +    A2x(i,j+3,k+3) &
                      +  9*A2x(i,j-3,k+2) -  81*A2x(i,j-2,k+2) +  405*A2x(i,j-1,k+2) -  405*A2x(i,j+1,k+2) +  81*A2x(i,j+2,k+2) -  9*A2x(i,j+3,k+2) &
                      - 45*A2x(i,j-3,k+1) + 405*A2x(i,j-2,k+1) - 2025*A2x(i,j-1,k+1) + 2025*A2x(i,j+1,k+1) - 405*A2x(i,j+2,k+1) + 45*A2x(i,j+3,k+1) &
                      + 45*A2x(i,j-3,k-1) - 405*A2x(i,j-2,k-1) + 2025*A2x(i,j-1,k-1) - 2025*A2x(i,j+1,k-1) + 405*A2x(i,j+2,k-1) - 45*A2x(i,j+3,k-1) &
                      -  9*A2x(i,j-3,k-2) +  81*A2x(i,j-2,k-2) -  405*A2x(i,j-1,k-2) +  405*A2x(i,j+1,k-2) -  81*A2x(i,j+2,k-2) +  9*A2x(i,j+3,k-2) &
                      +    A2x(i,j-3,k-3) -   9*A2x(i,j-2,k-3) +   45*A2x(i,j-1,k-3) -   45*A2x(i,j+1,k-3) +   9*A2x(i,j+2,k-3) -    A2x(i,j+3,k-3) ) * odydz3600
      d2_lA2(2,2,3) = (    -A2y(i,j-3,k+3) +   9*A2y(i,j-2,k+3) -   45*A2y(i,j-1,k+3) +   45*A2y(i,j+1,k+3) -   9*A2y(i,j+2,k+3) +    A2y(i,j+3,k+3) &
                      +  9*A2y(i,j-3,k+2) -  81*A2y(i,j-2,k+2) +  405*A2y(i,j-1,k+2) -  405*A2y(i,j+1,k+2) +  81*A2y(i,j+2,k+2) -  9*A2y(i,j+3,k+2) &
                      - 45*A2y(i,j-3,k+1) + 405*A2y(i,j-2,k+1) - 2025*A2y(i,j-1,k+1) + 2025*A2y(i,j+1,k+1) - 405*A2y(i,j+2,k+1) + 45*A2y(i,j+3,k+1) &
                      + 45*A2y(i,j-3,k-1) - 405*A2y(i,j-2,k-1) + 2025*A2y(i,j-1,k-1) - 2025*A2y(i,j+1,k-1) + 405*A2y(i,j+2,k-1) - 45*A2y(i,j+3,k-1) &
                      -  9*A2y(i,j-3,k-2) +  81*A2y(i,j-2,k-2) -  405*A2y(i,j-1,k-2) +  405*A2y(i,j+1,k-2) -  81*A2y(i,j+2,k-2) +  9*A2y(i,j+3,k-2) &
                     +    A2y(i,j-3,k-3) -   9*A2y(i,j-2,k-3) +   45*A2y(i,j-1,k-3) -   45*A2y(i,j+1,k-3) +   9*A2y(i,j+2,k-3) -    A2y(i,j+3,k-3) ) * odydz3600
      d2_lA2(3,2,3) = (    -A2z(i,j-3,k+3) +   9*A2z(i,j-2,k+3) -   45*A2z(i,j-1,k+3) +   45*A2z(i,j+1,k+3) -   9*A2z(i,j+2,k+3) +    A2z(i,j+3,k+3) &
                      +  9*A2z(i,j-3,k+2) -  81*A2z(i,j-2,k+2) +  405*A2z(i,j-1,k+2) -  405*A2z(i,j+1,k+2) +  81*A2z(i,j+2,k+2) -  9*A2z(i,j+3,k+2) &
                      - 45*A2z(i,j-3,k+1) + 405*A2z(i,j-2,k+1) - 2025*A2z(i,j-1,k+1) + 2025*A2z(i,j+1,k+1) - 405*A2z(i,j+2,k+1) + 45*A2z(i,j+3,k+1) &
                      + 45*A2z(i,j-3,k-1) - 405*A2z(i,j-2,k-1) + 2025*A2z(i,j-1,k-1) - 2025*A2z(i,j+1,k-1) + 405*A2z(i,j+2,k-1) - 45*A2z(i,j+3,k-1) &
                      -  9*A2z(i,j-3,k-2) +  81*A2z(i,j-2,k-2) -  405*A2z(i,j-1,k-2) +  405*A2z(i,j+1,k-2) -  81*A2z(i,j+2,k-2) +  9*A2z(i,j+3,k-2) &
                      +    A2z(i,j-3,k-3) -   9*A2z(i,j-2,k-3) +   45*A2z(i,j-1,k-3) -   45*A2z(i,j+1,k-3) +   9*A2z(i,j+2,k-3) -    A2z(i,j+3,k-3) ) * odydz3600

      d2_lA2(:,2,1) = d2_lA2(:,1,2)
      d2_lA2(:,3,1) = d2_lA2(:,1,3)
      d2_lA2(:,3,2) = d2_lA2(:,2,3)


      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta(1) ) )
        dj = int( sign( one, beta(2) ) )
        dk = int( sign( one, beta(3) ) )

        ! ad1_lE1(3)
        d1_f1(1) = di * (   2*E1x(i-2*di,j,k) - 24*E1x(i-di,j,k) - 35*E1x(i,j,k) + 80*E1x(i+di,j,k) &
                        - 30*E1x(i+2*di,j,k) + 8*E1x(i+3*di,j,k) - E1x(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*E1x(i,j-2*dj,k) - 24*E1x(i,j-dj,k) - 35*E1x(i,j,k) + 80*E1x(i,j+dj,k) &
                        - 30*E1x(i,j+2*dj,k) + 8*E1x(i,j+3*dj,k) - E1x(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*E1x(i,j,k-2*dk) - 24*E1x(i,j,k-dk) - 35*E1x(i,j,k) + 80*E1x(i,j,k+dk) &
                        - 30*E1x(i,j,k+2*dk) + 8*E1x(i,j,k+3*dk) - E1x(i,j,k+4*dk) ) * odz60
        ad1_lE1(1) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * (   2*E1y(i-2*di,j,k) - 24*E1y(i-di,j,k) - 35*E1y(i,j,k) + 80*E1y(i+di,j,k) &
                        - 30*E1y(i+2*di,j,k) +  8*E1y(i+3*di,j,k) -  E1y(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*E1y(i,j-2*dj,k) - 24*E1y(i,j-dj,k) - 35*E1y(i,j,k) + 80*E1y(i,j+dj,k) &
                        - 30*E1y(i,j+2*dj,k) +  8*E1y(i,j+3*dj,k) -  E1y(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*E1y(i,j,k-2*dk) - 24*E1y(i,j,k-dk) - 35*E1y(i,j,k) + 80*E1y(i,j,k+dk) &
                        - 30*E1y(i,j,k+2*dk) +  8*E1y(i,j,k+3*dk) -  E1y(i,j,k+4*dk) ) * odz60
        ad1_lE1(2) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * (   2*E1z(i-2*di,j,k) - 24*E1z(i-di,j,k) - 35*E1z(i,j,k) + 80*E1z(i+di,j,k) &
                        - 30*E1z(i+2*di,j,k) +  8*E1z(i+3*di,j,k) -  E1z(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*E1z(i,j-2*dj,k) - 24*E1z(i,j-dj,k) - 35*E1z(i,j,k) + 80*E1z(i,j+dj,k) &
                        - 30*E1z(i,j+2*dj,k) +  8*E1z(i,j+3*dj,k) -  E1z(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*E1z(i,j,k-2*dk) - 24*E1z(i,j,k-dk) - 35*E1z(i,j,k) + 80*E1z(i,j,k+dk) &
                        - 30*E1z(i,j,k+2*dk) +  8*E1z(i,j,k+3*dk) -  E1z(i,j,k+4*dk) ) * odz60
        ad1_lE1(3) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        ! ad1_lA1(3)
        d1_f1(1) = di * (   2*A1x(i-2*di,j,k) - 24*A1x(i-di,j,k) - 35*A1x(i,j,k) + 80*A1x(i+di,j,k) &
                        - 30*A1x(i+2*di,j,k) + 8*A1x(i+3*di,j,k) - A1x(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*A1x(i,j-2*dj,k) - 24*A1x(i,j-dj,k) - 35*A1x(i,j,k) + 80*A1x(i,j+dj,k) &
                        - 30*A1x(i,j+2*dj,k) + 8*A1x(i,j+3*dj,k) - A1x(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*A1x(i,j,k-2*dk) - 24*A1x(i,j,k-dk) - 35*A1x(i,j,k) + 80*A1x(i,j,k+dk) &
                        - 30*A1x(i,j,k+2*dk) + 8*A1x(i,j,k+3*dk) - A1x(i,j,k+4*dk) ) * odz60
        ad1_lA1(1) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * (   2*A1y(i-2*di,j,k) - 24*A1y(i-di,j,k) - 35*A1y(i,j,k) + 80*A1y(i+di,j,k) &
                        - 30*A1y(i+2*di,j,k) +  8*A1y(i+3*di,j,k) -  A1y(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*A1y(i,j-2*dj,k) - 24*A1y(i,j-dj,k) - 35*A1y(i,j,k) + 80*A1y(i,j+dj,k) &
                        - 30*A1y(i,j+2*dj,k) +  8*A1y(i,j+3*dj,k) -  A1y(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*A1y(i,j,k-2*dk) - 24*A1y(i,j,k-dk) - 35*A1y(i,j,k) + 80*A1y(i,j,k+dk) &
                        - 30*A1y(i,j,k+2*dk) +  8*A1y(i,j,k+3*dk) -  A1y(i,j,k+4*dk) ) * odz60
        ad1_lA1(2) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        d1_f1(1) = di * (   2*A1z(i-2*di,j,k) - 24*A1z(i-di,j,k) - 35*A1z(i,j,k) + 80*A1z(i+di,j,k) &
                        - 30*A1z(i+2*di,j,k) +  8*A1z(i+3*di,j,k) -  A1z(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*A1z(i,j-2*dj,k) - 24*A1z(i,j-dj,k) - 35*A1z(i,j,k) + 80*A1z(i,j+dj,k) &
                        - 30*A1z(i,j+2*dj,k) +  8*A1z(i,j+3*dj,k) -  A1z(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*A1z(i,j,k-2*dk) - 24*A1z(i,j,k-dk) - 35*A1z(i,j,k) + 80*A1z(i,j,k+dk) &
                        - 30*A1z(i,j,k+2*dk) +  8*A1z(i,j,k+3*dk) -  A1z(i,j,k+4*dk) ) * odz60
        ad1_lA1(3) = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        ! ad1_lZeta1
        d1_f1(1) = di * (   2*Zeta1(i-2*di,j,k) - 24*Zeta1(i-di,j,k) - 35*Zeta1(i,j,k) + 80*Zeta1(i+di,j,k) &
                        - 30*Zeta1(i+2*di,j,k) +  8*Zeta1(i+3*di,j,k) -  Zeta1(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*Zeta1(i,j-2*dj,k) - 24*Zeta1(i,j-dj,k) - 35*Zeta1(i,j,k) + 80*Zeta1(i,j+dj,k) &
                        - 30*Zeta1(i,j+2*dj,k) +  8*Zeta1(i,j+3*dj,k) -  Zeta1(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*Zeta1(i,j,k-2*dk) - 24*Zeta1(i,j,k-dk) - 35*Zeta1(i,j,k) + 80*Zeta1(i,j,k+dk) &
                        - 30*Zeta1(i,j,k+2*dk) +  8*Zeta1(i,j,k+3*dk) -  Zeta1(i,j,k+4*dk) ) * odz60
        ad1_lZeta1 = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)

        ! ad1_lAphi1
        d1_f1(1) = di * (   2*Aphi1(i-2*di,j,k) - 24*Aphi1(i-di,j,k) - 35*Aphi1(i,j,k) + 80*Aphi1(i+di,j,k) &
                        - 30*Aphi1(i+2*di,j,k) +  8*Aphi1(i+3*di,j,k) -  Aphi1(i+4*di,j,k) ) * odx60
        d1_f1(2) = dj * (   2*Aphi1(i,j-2*dj,k) - 24*Aphi1(i,j-dj,k) - 35*Aphi1(i,j,k) + 80*Aphi1(i,j+dj,k) &
                        - 30*Aphi1(i,j+2*dj,k) +  8*Aphi1(i,j+3*dj,k) -  Aphi1(i,j+4*dj,k) ) * ody60
        d1_f1(3) = dk * (   2*Aphi1(i,j,k-2*dk) - 24*Aphi1(i,j,k-dk) - 35*Aphi1(i,j,k) + 80*Aphi1(i,j,k+dk) &
                        - 30*Aphi1(i,j,k+2*dk) +  8*Aphi1(i,j,k+3*dk) -  Aphi1(i,j,k+4*dk) ) * odz60
        ad1_lAphi1 = beta(1)*d1_f1(1) + beta(2)*d1_f1(2) + beta(3)*d1_f1(3)
        
        
        ! ad1_lE2(3)
        d1_f2(1) = di * (   2*E2x(i-2*di,j,k) - 24*E2x(i-di,j,k) - 35*E2x(i,j,k) + 80*E2x(i+di,j,k) &
                        - 30*E2x(i+2*di,j,k) + 8*E2x(i+3*di,j,k) - E2x(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*E2x(i,j-2*dj,k) - 24*E2x(i,j-dj,k) - 35*E2x(i,j,k) + 80*E2x(i,j+dj,k) &
                        - 30*E2x(i,j+2*dj,k) + 8*E2x(i,j+3*dj,k) - E2x(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*E2x(i,j,k-2*dk) - 24*E2x(i,j,k-dk) - 35*E2x(i,j,k) + 80*E2x(i,j,k+dk) &
                        - 30*E2x(i,j,k+2*dk) + 8*E2x(i,j,k+3*dk) - E2x(i,j,k+4*dk) ) * odz60
        ad1_lE2(1) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * (   2*E2y(i-2*di,j,k) - 24*E2y(i-di,j,k) - 35*E2y(i,j,k) + 80*E2y(i+di,j,k) &
                        - 30*E2y(i+2*di,j,k) +  8*E2y(i+3*di,j,k) -  E2y(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*E2y(i,j-2*dj,k) - 24*E2y(i,j-dj,k) - 35*E2y(i,j,k) + 80*E2y(i,j+dj,k) &
                        - 30*E2y(i,j+2*dj,k) +  8*E2y(i,j+3*dj,k) -  E2y(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*E2y(i,j,k-2*dk) - 24*E2y(i,j,k-dk) - 35*E2y(i,j,k) + 80*E2y(i,j,k+dk) &
                        - 30*E2y(i,j,k+2*dk) +  8*E2y(i,j,k+3*dk) -  E2y(i,j,k+4*dk) ) * odz60
        ad1_lE2(2) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * (   2*E2z(i-2*di,j,k) - 24*E2z(i-di,j,k) - 35*E2z(i,j,k) + 80*E2z(i+di,j,k) &
                        - 30*E2z(i+2*di,j,k) +  8*E2z(i+3*di,j,k) -  E2z(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*E2z(i,j-2*dj,k) - 24*E2z(i,j-dj,k) - 35*E2z(i,j,k) + 80*E2z(i,j+dj,k) &
                        - 30*E2z(i,j+2*dj,k) +  8*E2z(i,j+3*dj,k) -  E2z(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*E2z(i,j,k-2*dk) - 24*E2z(i,j,k-dk) - 35*E2z(i,j,k) + 80*E2z(i,j,k+dk) &
                        - 30*E2z(i,j,k+2*dk) +  8*E2z(i,j,k+3*dk) -  E2z(i,j,k+4*dk) ) * odz60
        ad1_lE2(3) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        ! ad1_lA2(3)
        d1_f2(1) = di * (   2*A2x(i-2*di,j,k) - 24*A2x(i-di,j,k) - 35*A2x(i,j,k) + 80*A2x(i+di,j,k) &
                        - 30*A2x(i+2*di,j,k) + 8*A2x(i+3*di,j,k) - A2x(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*A2x(i,j-2*dj,k) - 24*A2x(i,j-dj,k) - 35*A2x(i,j,k) + 80*A2x(i,j+dj,k) &
                        - 30*A2x(i,j+2*dj,k) + 8*A2x(i,j+3*dj,k) - A2x(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*A2x(i,j,k-2*dk) - 24*A2x(i,j,k-dk) - 35*A2x(i,j,k) + 80*A2x(i,j,k+dk) &
                        - 30*A2x(i,j,k+2*dk) + 8*A2x(i,j,k+3*dk) - A2x(i,j,k+4*dk) ) * odz60
        ad1_lA2(1) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * (   2*A2y(i-2*di,j,k) - 24*A2y(i-di,j,k) - 35*A2y(i,j,k) + 80*A2y(i+di,j,k) &
                        - 30*A2y(i+2*di,j,k) +  8*A2y(i+3*di,j,k) -  A2y(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*A2y(i,j-2*dj,k) - 24*A2y(i,j-dj,k) - 35*A2y(i,j,k) + 80*A2y(i,j+dj,k) &
                        - 30*A2y(i,j+2*dj,k) +  8*A2y(i,j+3*dj,k) -  A2y(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*A2y(i,j,k-2*dk) - 24*A2y(i,j,k-dk) - 35*A2y(i,j,k) + 80*A2y(i,j,k+dk) &
                        - 30*A2y(i,j,k+2*dk) +  8*A2y(i,j,k+3*dk) -  A2y(i,j,k+4*dk) ) * odz60
        ad1_lA2(2) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        d1_f2(1) = di * (   2*A2z(i-2*di,j,k) - 24*A2z(i-di,j,k) - 35*A2z(i,j,k) + 80*A2z(i+di,j,k) &
                        - 30*A2z(i+2*di,j,k) +  8*A2z(i+3*di,j,k) -  A2z(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*A2z(i,j-2*dj,k) - 24*A2z(i,j-dj,k) - 35*A2z(i,j,k) + 80*A2z(i,j+dj,k) &
                        - 30*A2z(i,j+2*dj,k) +  8*A2z(i,j+3*dj,k) -  A2z(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*A2z(i,j,k-2*dk) - 24*A2z(i,j,k-dk) - 35*A2z(i,j,k) + 80*A2z(i,j,k+dk) &
                        - 30*A2z(i,j,k+2*dk) +  8*A2z(i,j,k+3*dk) -  A2z(i,j,k+4*dk) ) * odz60
        ad1_lA2(3) = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        ! ad1_lZeta2
        d1_f2(1) = di * (   2*Zeta2(i-2*di,j,k) - 24*Zeta2(i-di,j,k) - 35*Zeta2(i,j,k) + 80*Zeta2(i+di,j,k) &
                        - 30*Zeta2(i+2*di,j,k) +  8*Zeta2(i+3*di,j,k) -  Zeta2(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*Zeta2(i,j-2*dj,k) - 24*Zeta2(i,j-dj,k) - 35*Zeta2(i,j,k) + 80*Zeta2(i,j+dj,k) &
                        - 30*Zeta2(i,j+2*dj,k) +  8*Zeta2(i,j+3*dj,k) -  Zeta2(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*Zeta2(i,j,k-2*dk) - 24*Zeta2(i,j,k-dk) - 35*Zeta2(i,j,k) + 80*Zeta2(i,j,k+dk) &
                        - 30*Zeta2(i,j,k+2*dk) +  8*Zeta2(i,j,k+3*dk) -  Zeta2(i,j,k+4*dk) ) * odz60
        ad1_lZeta2 = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

        ! ad1_lAphi2
        d1_f2(1) = di * (   2*Aphi2(i-2*di,j,k) - 24*Aphi2(i-di,j,k) - 35*Aphi2(i,j,k) + 80*Aphi2(i+di,j,k) &
                        - 30*Aphi2(i+2*di,j,k) +  8*Aphi2(i+3*di,j,k) -  Aphi2(i+4*di,j,k) ) * odx60
        d1_f2(2) = dj * (   2*Aphi2(i,j-2*dj,k) - 24*Aphi2(i,j-dj,k) - 35*Aphi2(i,j,k) + 80*Aphi2(i,j+dj,k) &
                        - 30*Aphi2(i,j+2*dj,k) +  8*Aphi2(i,j+3*dj,k) -  Aphi2(i,j+4*dj,k) ) * ody60
        d1_f2(3) = dk * (   2*Aphi2(i,j,k-2*dk) - 24*Aphi2(i,j,k-dk) - 35*Aphi2(i,j,k) + 80*Aphi2(i,j,k+dk) &
                        - 30*Aphi2(i,j,k+2*dk) +  8*Aphi2(i,j,k+3*dk) -  Aphi2(i,j,k+4*dk) ) * odz60
        ad1_lAphi2 = beta(1)*d1_f2(1) + beta(2)*d1_f2(2) + beta(3)*d1_f2(3)

      else

        ! ad1_lE1(3)
        ad1_lE1(1) = beta(1)*d1_lE1(1,1) + beta(2)*d1_lE1(1,2) + beta(3)*d1_lE1(1,3)
        ad1_lE1(2) = beta(1)*d1_lE1(2,1) + beta(2)*d1_lE1(2,2) + beta(3)*d1_lE1(2,3)
        ad1_lE1(3) = beta(1)*d1_lE1(3,1) + beta(2)*d1_lE1(3,2) + beta(3)*d1_lE1(3,3)

        ! ad1_lA1(3)
        ad1_lA1(1) = beta(1)*d1_lA1(1,1) + beta(2)*d1_lA1(1,2) + beta(3)*d1_lA1(1,3)
        ad1_lA1(2) = beta(1)*d1_lA1(2,1) + beta(2)*d1_lA1(2,2) + beta(3)*d1_lA1(2,3)
        ad1_lA1(3) = beta(1)*d1_lA1(3,1) + beta(2)*d1_lA1(3,2) + beta(3)*d1_lA1(3,3)

        ! ad1_lZeta1
        ad1_lZeta1 = beta(1)*d1_lZeta1(1) + beta(2)*d1_lZeta1(2) + beta(3)*d1_lZeta1(3)

        ! ad1_lAphi1
        ad1_lAphi1 = beta(1)*d1_lAphi1(1) + beta(2)*d1_lAphi1(2) + beta(3)*d1_lAphi1(3)
        
        
        ! ad1_lE2(3)
        ad1_lE2(1) = beta(1)*d1_lE2(1,1) + beta(2)*d1_lE2(1,2) + beta(3)*d1_lE2(1,3)
        ad1_lE2(2) = beta(1)*d1_lE2(2,1) + beta(2)*d1_lE2(2,2) + beta(3)*d1_lE2(2,3)
        ad1_lE2(3) = beta(1)*d1_lE2(3,1) + beta(2)*d1_lE2(3,2) + beta(3)*d1_lE2(3,3)

        ! ad1_lA2(3)
        ad1_lA2(1) = beta(1)*d1_lA2(1,1) + beta(2)*d1_lA2(1,2) + beta(3)*d1_lA2(1,3)
        ad1_lA2(2) = beta(1)*d1_lA2(2,1) + beta(2)*d1_lA2(2,2) + beta(3)*d1_lA2(2,3)
        ad1_lA2(3) = beta(1)*d1_lA2(3,1) + beta(2)*d1_lA2(3,2) + beta(3)*d1_lA2(3,3)

        ! ad1_lZeta2
        ad1_lZeta2 = beta(1)*d1_lZeta2(1) + beta(2)*d1_lZeta2(2) + beta(3)*d1_lZeta2(3)

        ! ad1_lAphi2
        ad1_lAphi2 = beta(1)*d1_lAphi2(1) + beta(2)*d1_lAphi2(2) + beta(3)*d1_lAphi2(3)

      end if
      !-------------------------------------------

    else
      call CCTK_WARN(0, "derivs_order not yet implemented.")
    end if

    !------------ Christoffel symbols ----------
    cf1 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do
    cf1(:,2,1) = cf1(:,1,2)
    cf1(:,3,1) = cf1(:,1,3)
    cf1(:,3,2) = cf1(:,2,3)

    cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    cf2(:,2,1) = cf2(:,1,2)
    cf2(:,3,1) = cf2(:,1,3)
    cf2(:,3,2) = cf2(:,2,3)
    !-------------------------------------------


    !------------ Covariant derivatives --------
    cd_lA1 = d1_lA1
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          cd_lA1(a,b) = cd_lA1(a,b) - cf2(m,a,b) * lA1(m)
       end do
      end do
    end do

    ! note that this is *not* D_i D_j A_k (ie, the covariant derivative of D_j A_k)!
    ! but since this variable only appears in the combination (D_k D^i A^k - D_k D^k A^i),
    ! the terms that are computed herein are enough.

    cd_dA1 = d2_lA1
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            cd_dA1(a,b,c) = cd_dA1(a,b,c) - cf2(m,a,c) * d1_lA1(m,b)         &
                                        - cf2(m,b,c) * d1_lA1(a,m)
           end do
         end do
       end do
     end do
     
     
     cd_lA2 = d1_lA2
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          cd_lA2(a,b) = cd_lA2(a,b) - cf2(m,a,b) * lA2(m)
       end do
      end do
    end do

    cd_dA2 = d2_lA2
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            cd_dA2(a,b,c) = cd_dA2(a,b,c) - cf2(m,a,c) * d1_lA2(m,b)         &
                                        - cf2(m,b,c) * d1_lA2(a,m)
           end do
         end do
       end do
     end do
    !-------------------------------------------


    !--------- Evolution of E, A, Aphi, Zeta ----------

    ! rhs_lE1
    rhs_lE1  = ad1_lE1

    do a = 1, 3
      do m = 1, 3
        rhs_lE1(a) = rhs_lE1(a) - lE1(m) * d1_beta(a,m)
      end do
    end do

    rhs_lE1  = rhs_lE1 + alph * trk * lE1

    do a = 1, 3
       do b = 1, 3
          rhs_lE1(a) = rhs_lE1(a) + alph * ch**conf_fac_exponent * hu(a,b) * d1_lZeta1(b)      &
                                + alph * mu*mu * ch**conf_fac_exponent * hu(a,b) * lA1(b)
          do c = 1, 3
             do m = 1, 3
                rhs_lE1(a) = rhs_lE1(a)                                                       &
                     + ch**(2*conf_fac_exponent) * hu(a,b) * hu(c,m) * d1_alph(m)           &
                                                 * ( d1_lA1(c,b) - d1_lA1(b,c) )              &
                     + alph * ch**(2*conf_fac_exponent) * hu(a,b) * hu(c,m)                 &
                                                 * ( cd_dA1(m,b,c) - cd_dA1(b,m,c) )          &
                     + 0.50d0*conf_fac_exponent * alph * ch**(2*conf_fac_exponent - 1)      &
                        * hu(a,b) * hu(c,m) * ( d1_lA1(m,b)*d1_ch(c) - d1_lA1(b,c)*d1_ch(m) )
             end do
          end do
       end do
    end do


    ! rhs_lA1
    rhs_lA1  = ad1_lA1

    do a = 1, 3
      do m = 1, 3
        rhs_lA1(a) = rhs_lA1(a) + lA1(m) * d1_beta(m,a)
      end do
    end do

    rhs_lA1  = rhs_lA1 - alph * d1_lAphi1 - lAphi1 * d1_alph

    do a = 1,3
       do b = 1,3
          rhs_lA1(a) = rhs_lA1(a) - alph * hh(a,b) * lE1(b) * ch**(-conf_fac_exponent)
       end do
    end do


    ! rhs_lAphi1
    rhs_lAphi1  = ad1_lAphi1  + alph * trk * lAphi1 - alph * lZeta1

    do a = 1, 3
        do b = 1, 3
           rhs_lAphi1 = rhs_lAphi1                                                     &
                     + 0.5d0 * alph * conf_fac_exponent * ch**(conf_fac_exponent-1)  &
                             * hu(a,b) * lA1(a) * d1_ch(b)                            &
                     - alph * ch**conf_fac_exponent * hu(a,b) * cd_lA1(a,b)           &
                     - ch**conf_fac_exponent * hu(a,b) * lA1(a) * d1_alph(b)
        end do
    end do


    ! rhs_lZeta1
    rhs_lZeta1 = ad1_lZeta1 - alph * kappa * lZeta1 + alph * mu*mu * lAphi1

    do a = 1, 3
       rhs_lZeta1 = rhs_lZeta1 + alph * d1_lE1(a,a)                         &
                 - 1.5d0 * conf_fac_exponent * alph * lE1(a) * d1_ch(a) / ch
    end do

    rhs_lZeta1 = Zeta_Omega_fac * rhs_lZeta1
    
    ! rhs_lE2
    rhs_lE2  = ad1_lE2

    do a = 1, 3
      do m = 1, 3
        rhs_lE2(a) = rhs_lE2(a) - lE2(m) * d1_beta(a,m)
      end do
    end do

    rhs_lE2  = rhs_lE2 + alph * trk * lE2

    do a = 1, 3
       do b = 1, 3
          rhs_lE2(a) = rhs_lE2(a) + alph * ch**conf_fac_exponent * hu(a,b) * d1_lZeta2(b)      &
                                + alph * mu*mu * ch**conf_fac_exponent * hu(a,b) * lA2(b)
          do c = 1, 3
             do m = 1, 3
                rhs_lE2(a) = rhs_lE2(a)                                                       &
                     + ch**(2*conf_fac_exponent) * hu(a,b) * hu(c,m) * d1_alph(m)           &
                                                 * ( d1_lA2(c,b) - d1_lA2(b,c) )              &
                     + alph * ch**(2*conf_fac_exponent) * hu(a,b) * hu(c,m)                 &
                                                 * ( cd_dA2(m,b,c) - cd_dA2(b,m,c) )          &
                     + 0.50d0*conf_fac_exponent * alph * ch**(2*conf_fac_exponent - 1)      &
                        * hu(a,b) * hu(c,m) * ( d1_lA2(m,b)*d1_ch(c) - d1_lA2(b,c)*d1_ch(m) )
             end do
          end do
       end do
    end do


    ! rhs_lA2
    rhs_lA2  = ad1_lA2

    do a = 1, 3
      do m = 1, 3
        rhs_lA2(a) = rhs_lA2(a) + lA2(m) * d1_beta(m,a)
      end do
    end do

    rhs_lA2  = rhs_lA2 - alph * d1_lAphi2 - lAphi2 * d1_alph

    do a = 1,3
       do b = 1,3
          rhs_lA2(a) = rhs_lA2(a) - alph * hh(a,b) * lE2(b) * ch**(-conf_fac_exponent)
       end do
    end do


    ! rhs_lAphi2
    rhs_lAphi2  = ad1_lAphi2  + alph * trk * lAphi2 - alph * lZeta2

    do a = 1, 3
        do b = 1, 3
           rhs_lAphi2 = rhs_lAphi2                                                     &
                     + 0.5d0 * alph * conf_fac_exponent * ch**(conf_fac_exponent-1)  &
                             * hu(a,b) * lA2(a) * d1_ch(b)                            &
                     - alph * ch**conf_fac_exponent * hu(a,b) * cd_lA2(a,b)           &
                     - ch**conf_fac_exponent * hu(a,b) * lA2(a) * d1_alph(b)
        end do
    end do


    ! rhs_lZeta2
    rhs_lZeta2 = ad1_lZeta2 - alph * kappa * lZeta2 + alph * mu*mu * lAphi2

    do a = 1, 3
       rhs_lZeta2 = rhs_lZeta2 + alph * d1_lE2(a,a)                         &
                 - 1.5d0 * conf_fac_exponent * alph * lE2(a) * d1_ch(a) / ch
    end do

    rhs_lZeta2 = Zeta_Omega_fac * rhs_lZeta2


    !-------------------------------------------

    ! if( abs(y(i,j,k)) < 1.0d-05 .and. abs(z(i,j,k)) < 1.0d-05 ) then
    !    write(*,*) 'i, j, k    = ', i, j, k
    !    write(*,*) 'x          = ', x(i,j,k)
    !    write(*,*) 'Ex         = ', lE(1)
    !    write(*,*) 'rhs_lE     = ', rhs_lE
    !    write(*,*) 'rhs_lA     = ', rhs_lA
    !    write(*,*) 'rhs_lAphi  = ', rhs_lAphi
    !    write(*,*) 'rhs_lZeta  = ', rhs_lZeta
    !    call flush(6)
    ! end if


    !-------- Write to the gridfunctions -------

    rhs_E1x(i,j,k) = rhs_lE1(1)
    rhs_E1y(i,j,k) = rhs_lE1(2)
    rhs_E1z(i,j,k) = rhs_lE1(3)

    rhs_A1x(i,j,k) = rhs_lA1(1)
    rhs_A1y(i,j,k) = rhs_lA1(2)
    rhs_A1z(i,j,k) = rhs_lA1(3)

    rhs_Zeta1(i,j,k) = rhs_lZeta1

    rhs_Aphi1(i,j,k) = rhs_lAphi1
    
    
    rhs_E2x(i,j,k) = rhs_lE2(1)
    rhs_E2y(i,j,k) = rhs_lE2(2)
    rhs_E2z(i,j,k) = rhs_lE2(3)

    rhs_A2x(i,j,k) = rhs_lA2(1)
    rhs_A2y(i,j,k) = rhs_lA2(2)
    rhs_A2z(i,j,k) = rhs_lA2(3)

    rhs_Zeta2(i,j,k) = rhs_lZeta2

    rhs_Aphi2(i,j,k) = rhs_lAphi2

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine ComplexProca_calc_rhs
!
!=============================================================================
!

subroutine ComplexProca_calc_rhs_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, parameter :: one  = 1.0
  CCTK_REAL, parameter :: zero = 0.0
  CCTK_INT ierr

  ierr = NewRad_Apply(cctkGH, E1x, rhs_E1x, E01(1), one, n_E(1))
  ierr = NewRad_Apply(cctkGH, E1y, rhs_E1y, E01(2), one, n_E(2))
  ierr = NewRad_Apply(cctkGH, E1z, rhs_E1z, E01(3), one, n_E(3))

  ierr = NewRad_Apply(cctkGH, A1x, rhs_A1x, A01(1), one, n_A(1))
  ierr = NewRad_Apply(cctkGH, A1y, rhs_A1y, A01(2), one, n_A(2))
  ierr = NewRad_Apply(cctkGH, A1z, rhs_A1z, A01(3), one, n_A(3))

  ierr = NewRad_Apply(cctkGH, Aphi1, rhs_Aphi1, Aphi01, one, n_Aphi)
  ierr = NewRad_Apply(cctkGH, Zeta1, rhs_Zeta1, zero, one, n_Zeta)
  
  ierr = NewRad_Apply(cctkGH, E2x, rhs_E2x, E02(1), one, n_E(1))
  ierr = NewRad_Apply(cctkGH, E2y, rhs_E2y, E02(2), one, n_E(2))
  ierr = NewRad_Apply(cctkGH, E2z, rhs_E2z, E02(3), one, n_E(3))

  ierr = NewRad_Apply(cctkGH, A2x, rhs_A2x, A02(1), one, n_A(1))
  ierr = NewRad_Apply(cctkGH, A2y, rhs_A2y, A02(2), one, n_A(2))
  ierr = NewRad_Apply(cctkGH, A2z, rhs_A2z, A02(3), one, n_A(3))

  ierr = NewRad_Apply(cctkGH, Aphi2, rhs_Aphi2, Aphi02, one, n_Aphi)
  ierr = NewRad_Apply(cctkGH, Zeta2, rhs_Zeta2, zero, one, n_Zeta)

end subroutine ComplexProca_calc_rhs_bdry
