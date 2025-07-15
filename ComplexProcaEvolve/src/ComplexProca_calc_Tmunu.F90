#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine ComplexProca_calc_Tmunu( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                gg(3,3), gu(3,3), detgg
  CCTK_REAL                lE1(3), lB1(3), lA1(3), lAphi1
  CCTK_REAL                Ed1(3), Bd1(3)
  CCTK_REAL                lE2(3), lB2(3), lA2(3), lAphi2
  CCTK_REAL                Ed2(3), Bd2(3)
  CCTK_REAL                Tab(4,4)

  ! First derivatives
  CCTK_REAL                d1_lA1(3,3)
  CCTK_REAL                d1_lA2(3,3)

  ! Auxiliary variables
  CCTK_REAL                eps_lc_d(3,3,3), eps_lc_u(3,3,3)

  ! Matter variables
  CCTK_REAL                srcE, srcjdi(3), srcSij(3,3)

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12
  CCTK_REAL                odx60, ody60, odz60
  CCTK_REAL                aux

  CCTK_REAL, parameter ::  one  = 1
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi4  = 4*pi
  CCTK_REAL, parameter ::  pi8  = 8*pi
  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c, m

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i,a,b,c,m,aux,&
  !$OMP                                 alph,beta,&
  !$OMP                                 gg,gu,detgg,&
  !$OMP                                 lE1,lB1,lA1,lAphi1,&
  !$OMP                                 lE2,lB2,lA2,lAphi2,&
  !$OMP                                 Ed1,Bd1,Ed2,Bd2,Tab,&
  !$OMP                                 eps_lc_d,eps_lc_u,&
  !$OMP                                 d1_lA1,d1_lA2, &
  !$OMP                                 srcE, srcjdi, srcSij)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           !------------ Get local variables ----------

           alph      = alp(i,j,k)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           gg(1,1)   = gxx(i,j,k)
           gg(1,2)   = gxy(i,j,k)
           gg(1,3)   = gxz(i,j,k)
           gg(2,2)   = gyy(i,j,k)
           gg(2,3)   = gyz(i,j,k)
           gg(3,3)   = gzz(i,j,k)
           gg(2,1)   = gg(1,2)
           gg(3,1)   = gg(1,3)
           gg(3,2)   = gg(2,3)

           lE1(1)     = E1x(i,j,k)
           lE1(2)     = E1y(i,j,k)
           lE1(3)     = E1z(i,j,k)

           lA1(1)     = A1x(i,j,k)
           lA1(2)     = A1y(i,j,k)
           lA1(3)     = A1z(i,j,k)

           lAphi1     = Aphi1(i,j,k)

          lE2(1)     = E2x(i,j,k)
          lE2(2)     = E2y(i,j,k)
          lE2(3)     = E2z(i,j,k)

          lA2(1)     = A2x(i,j,k)
          lA2(2)     = A2y(i,j,k)
          lA2(3)     = A2z(i,j,k)

          lAphi2     = Aphi2(i,j,k)

           Ed1(1)     = gg(1,1) * lE1(1) + gg(1,2) * lE1(2) + gg(1,3) * lE1(3)
           Ed1(2)     = gg(2,1) * lE1(1) + gg(2,2) * lE1(2) + gg(2,3) * lE1(3)
           Ed1(3)     = gg(3,1) * lE1(1) + gg(3,2) * lE1(2) + gg(3,3) * lE1(3)

           Ed2(1)     = gg(1,1) * lE2(1) + gg(1,2) * lE2(2) + gg(1,3) * lE2(3)
           Ed2(2)     = gg(2,1) * lE2(1) + gg(2,2) * lE2(2) + gg(2,3) * lE2(3)
           Ed2(3)     = gg(3,1) * lE2(1) + gg(3,2) * lE2(2) + gg(3,3) * lE2(3)


           !------------ Invert 3-metric ----------------
           detgg   =     gg(1,1) * gg(2,2) * gg(3,3)                              &
                   + 2 * gg(1,2) * gg(1,3) * gg(2,3)                              &
                   -     gg(1,1) * gg(2,3) ** 2                                   &
                   -     gg(2,2) * gg(1,3) ** 2                                   &
                   -     gg(3,3) * gg(1,2) ** 2

           gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detgg
           gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detgg
           gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detgg
           gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detgg
           gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detgg
           gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detgg
           gu(2,1) = gu(1,2)
           gu(3,1) = gu(1,3)
           gu(3,2) = gu(2,3)
           !-------------------------------------------------


           !------------- Centered 1st derivatives ----------

           if (derivs_order == 4) then

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

           else if (derivs_order == 6) then

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


           else
             call CCTK_WARN(0, "derivs_order not yet implemented.")
           end if


           !------------ Levi-Civita tensor ----------
           eps_lc_u        =  0
           eps_lc_u(1,2,3) =  1
           eps_lc_u(2,3,1) =  1
           eps_lc_u(3,1,2) =  1
           eps_lc_u(3,2,1) = -1
           eps_lc_u(2,1,3) = -1
           eps_lc_u(1,3,2) = -1
           eps_lc_u = eps_lc_u / sqrt(detgg)
           eps_lc_d = eps_lc_u * detgg
           !------------------------------------------


           ! magnetic field B (here used as an auxiliary variable)
           lB1 = 0
           do a = 1, 3
             do b = 1, 3
               do m = 1, 3
                 lB1(a) = lB1(a) + eps_lc_u(a,b,m) * d1_lA1(m,b)
               end do
             end do
           end do
           
           lB2 = 0
           do a = 1, 3
             do b = 1, 3
               do m = 1, 3
                 lB2(a) = lB2(a) + eps_lc_u(a,b,m) * d1_lA2(m,b)
               end do
             end do
           end do

           Bd1(1) = gg(1,1) * lB1(1) + gg(1,2) * lB1(2) + gg(1,3) * lB1(3)
           Bd1(2) = gg(2,1) * lB1(1) + gg(2,2) * lB1(2) + gg(2,3) * lB1(3)
           Bd1(3) = gg(3,1) * lB1(1) + gg(3,2) * lB1(2) + gg(3,3) * lB1(3)
           
           Bd2(1) = gg(1,1) * lB2(1) + gg(1,2) * lB2(2) + gg(1,3) * lB2(3)
           Bd2(2) = gg(2,1) * lB2(1) + gg(2,2) * lB2(2) + gg(2,3) * lB2(3)
           Bd2(3) = gg(3,1) * lB2(1) + gg(3,2) * lB2(2) + gg(3,3) * lB2(3)
           !-------------------------------------------

           !------------ Matter terms -----------------
           !
           ! mu = 0, 1, 2, 3; i,a = 1,2,3
           !
           ! n_mu = (-alph, 0, 0, 0)
           ! n^mu = (1, -betax, -betay, -betaz)/alph
           !
           ! rho = n^mu n^nu T_{mu nu}
           !     = (T_{00} - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alph^2
           !
           ! j_a = -h_a^mu n^nu T_{mu nu}
           !     = -(T_{a 0} - beta^j T_{a j})/alph
           !
           ! S_{a b} = h_{a mu} h_{b nu} T^{mu nu} = T_{a b}


           ! srcE = rho = n^mu n^nu T_{mu nu}
           srcE = mu*mu * lAphi1*lAphi1 + mu*mu * lAphi2*lAphi2
           do a = 1, 3
              do b = 1, 3
                 srcE = srcE + ( lE1(a) * lE1(b) + lB1(a) * lB1(b) ) * gg(a,b)    &
                             + mu*mu * lA1(a) * lA1(b) * gu(a,b) &
                             + ( lE2(a) * lE2(b) + lB2(a) * lB2(b) ) * gg(a,b)    &
                             + mu*mu * lA2(a) * lA2(b) * gu(a,b)
              end do
           end do
           srcE = srcE / pi8

           !srcjdi = j_a
           srcjdi = mu*mu * lAphi1 * lA1 + mu*mu * lAphi2 * lA2
           do a = 1, 3
              do b = 1, 3
                 do m = 1, 3
                    srcjdi(a) = srcjdi(a) + eps_lc_d(a,b,m) * lE1(b) * lB1(m) &
                    		 + eps_lc_d(a,b,m) * lE2(b) * lB2(m)
                 end do
              end do
           end do
           srcjdi = srcjdi / pi4


           ! srcSij = S_{a b}

           aux = 0
           do a = 1, 3
              do b = 1, 3
                 aux = aux + ( lE1(a) * lE1(b) + lB1(a) * lB1(b) ) * gg(a,b)    &
                           - mu*mu * lA1(a) * lA1(b) * gu(a,b) &
                           + ( lE2(a) * lE2(b) + lB2(a) * lB2(b) ) * gg(a,b)    &
                           - mu*mu * lA2(a) * lA2(b) * gu(a,b)
              end do
           end do

           srcSij = 0.5 * (aux + mu*mu * lAphi1*lAphi1) * gg &
                        + 0.5 * (mu*mu * lAphi2*lAphi2) * gg
           do a = 1, 3
              do b = 1, 3
                 srcSij(a,b) = srcSij(a,b) - Ed1(a) * Ed1(b) - Bd1(a) * Bd1(b)  &
                             + mu*mu * lA1(a) * lA1(b) &
                             - Ed2(a) * Ed2(b) - Bd2(a) * Bd2(b)  &
                             + mu*mu * lA2(a) * lA2(b)
              end do
           end do
           srcSij = srcSij / pi4

           !------------------------------------------


           ! now to fill in the stress-energy tensor. note that we use Tab(4,4)
           ! for T_{0 0}
           !
           ! T_{a b} = S_{a b}
           !
           ! T_{0 0} = alph^2 rho - 2 alph beta^a j_a + beta^a beta^b S_{a b}
           !
           ! T_{0 a} = -alph j_a + beta^b S_{a b}

           Tab(1:3,1:3) = srcSij(1:3,1:3)

           Tab(1:3,4) = -alph * srcjdi(1:3)
           do b = 1, 3
              Tab(1:3,4) = Tab(1:3,4) + beta(b) * srcSij(1:3,b)
           end do
           Tab(4,1:3) = Tab(1:3,4)

           Tab(4,4) = alph**2 * srcE
           do a = 1, 3
              Tab(4,4) = Tab(4,4) - 2 * alph * beta(a) * srcjdi(a)
              do b = 1, 3
                 Tab(4,4) = Tab(4,4) + beta(a) * beta(b) * srcSij(a,b)
              end do
           end do

           ! and finally store it in the Tmunu variables
           eTtt(i,j,k) = eTtt(i,j,k) + Tab(4,4)
           eTtx(i,j,k) = eTtx(i,j,k) + Tab(4,1)
           eTty(i,j,k) = eTty(i,j,k) + Tab(4,2)
           eTtz(i,j,k) = eTtz(i,j,k) + Tab(4,3)
           eTxx(i,j,k) = eTxx(i,j,k) + Tab(1,1)
           eTxy(i,j,k) = eTxy(i,j,k) + Tab(1,2)
           eTxz(i,j,k) = eTxz(i,j,k) + Tab(1,3)
           eTyy(i,j,k) = eTyy(i,j,k) + Tab(2,2)
           eTyz(i,j,k) = eTyz(i,j,k) + Tab(2,3)
           eTzz(i,j,k) = eTzz(i,j,k) + Tab(3,3)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine ComplexProca_calc_Tmunu
