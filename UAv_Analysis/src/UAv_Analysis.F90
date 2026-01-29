
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

subroutine UAv_Analysis_gfs( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL alph, beta(3), Tab(4,4)
  CCTK_REAL gd(3,3), gu(3,3), detgd

  CCTK_REAL aux, S, rho
  CCTK_REAL mom(3)
  
  ! names x0, y0, z0 used as members of the thorn
  ! They are set in dedicated functions, to be called before this routine
  CCTK_REAL x1, y1, z1

  CCTK_INT  i, j, k, m, n

  CCTK_INT type_bits, state_outside

  logical docalc

  type_bits     = -1
  state_outside = -1
  
  if (do_analysis_every .le. 0) then
     return
  end if
  
  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
     return
  endif

  if (excise_horizon /= 0) then

     call SpaceMask_GetTypeBits(type_bits, "mask")
     call SpaceMask_GetStateBits(state_outside, "mask", "outside")

     if (type_bits < 0) then
        call CCTK_WARN(0, "Thorn AHFinderDirect not activated, but excise_horizon requires it.")
     end if

     if (state_outside < 0) then
        call CCTK_WARN(0, "Error in obtaining AHFinderDirect GetStateBits")
     end if

  end if

!   write(*,*) 'Checking origin coordinates for the analysis in UAv_Analysis'
!   write(*,*) 'x0 = ', x0
!   write(*,*) 'y0 = ', y0
!   write(*,*) 'z0 = ', z0

  dE_gf_volume   = 0
  dJx_gf_volume  = 0
  dJy_gf_volume  = 0
  dJz_gf_volume  = 0
  dIxx_gf_volume = 0
  dIxy_gf_volume = 0
  dIxz_gf_volume = 0
  dIyy_gf_volume = 0
  dIyz_gf_volume = 0
  dIzz_gf_volume = 0
  if (compute_density_rho == 1) then
     density_rho    = 0
  end if
  if (compute_density_p == 1) then
     density_px     = 0
     density_py     = 0
     density_pz     = 0
  end if

  ! Note that these loops will also exclude at least one layer of points on the physical boundary,
  ! (see https://lists.einsteintoolkit.org/pipermail/users/2024-September/009465.html )
  ! but as of now we choose to live with that
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    ! checking if outside the horizon, if asking for it to be excised
    docalc = .true.
    if (excise_horizon /= 0) then
       if (.not. SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, state_outside)) then
          docalc = .false.
       end if
    end if

    ! if inside the horizon, no need to compute the rest (continue with the next
    ! iteration of the do loop)
    if (.not. docalc) cycle

    !--------------Get local variables ----------
    gd(1,1) = gxx(i,j,k)
    gd(1,2) = gxy(i,j,k)
    gd(1,3) = gxz(i,j,k)
    gd(2,2) = gyy(i,j,k)
    gd(2,3) = gyz(i,j,k)
    gd(3,3) = gzz(i,j,k)
    gd(2,1) = gd(1,2)
    gd(3,1) = gd(1,3)
    gd(3,2) = gd(2,3)

    alph    = alp(i,j,k)

    beta(1) = betax(i,j,k)
    beta(2) = betay(i,j,k)
    beta(3) = betaz(i,j,k)

    x1      = x(i,j,k) - x0
    y1      = y(i,j,k) - y0
    z1      = z(i,j,k) - z0

    ! stress-energy tensor variables
    Tab = 0
    if (stress_energy_state /= 0) then
       Tab(4,4) = eTtt(i,j,k)
       Tab(4,1) = eTtx(i,j,k)
       Tab(4,2) = eTty(i,j,k)
       Tab(4,3) = eTtz(i,j,k)
       Tab(1,1) = eTxx(i,j,k)
       Tab(1,2) = eTxy(i,j,k)
       Tab(1,3) = eTxz(i,j,k)
       Tab(2,2) = eTyy(i,j,k)
       Tab(2,3) = eTyz(i,j,k)
       Tab(3,3) = eTzz(i,j,k)
       Tab(1,4) = Tab(4,1)
       Tab(2,4) = Tab(4,2)
       Tab(3,4) = Tab(4,3)
       Tab(2,1) = Tab(1,2)
       Tab(3,1) = Tab(1,3)
       Tab(3,2) = Tab(2,3)
    end if
    !--------------------------------------------


    !-------------- Invert metric ---------------
    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2 * gd(1,2) * gd(1,3) * gd(2,3)                                &
            -     gd(1,1) * gd(2,3) ** 2                                     &
            -     gd(2,2) * gd(1,3) ** 2                                     &
            -     gd(3,3) * gd(1,2) ** 2
    gu(1,1) = (gd(2,2) * gd(3,3) - gd(2,3) ** 2     ) / detgd
    gu(2,2) = (gd(1,1) * gd(3,3) - gd(1,3) ** 2     ) / detgd
    gu(3,3) = (gd(1,1) * gd(2,2) - gd(1,2) ** 2     ) / detgd
    gu(1,2) = (gd(1,3) * gd(2,3) - gd(1,2) * gd(3,3)) / detgd
    gu(1,3) = (gd(1,2) * gd(2,3) - gd(1,3) * gd(2,2)) / detgd
    gu(2,3) = (gd(1,3) * gd(1,2) - gd(2,3) * gd(1,1)) / detgd
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
    !--------------------------------------------

    ! Eulerian energy density
    rho = Tab(4,4)
    do m = 1, 3
       rho = rho - 2 * beta(m) * Tab(m,4)
       do n = 1, 3
          rho = rho + beta(m) * beta(n) * Tab(m,n)
       end do
    end do
    rho = rho / ( alph * alph )

    if (compute_density_rho == 1) then
      density_rho(i,j,k) = rho
    end if
    
    ! momentum density
    do n = 1, 3
      mom(n) = Tab(4,n)
      do m = 1, 3
         mom(n) = mom(n) - beta(m) * Tab(m,n)
      end do
      mom(n) = - mom(n) / alph
    end do
   
    if (compute_density_p == 1) then
      density_px(i,j,k) = mom(1)
      density_py(i,j,k) = mom(2)
      density_pz(i,j,k) = mom(3)
    end if
   

    aux = 0
    do m = 1, 3
      do n = 1, 3
        aux = aux + beta(m) * beta(n) * Tab(m,n)
      end do
    end do
    aux = (Tab(4,4) - aux) / alph

    S = 0
    do m = 1, 3
       do n = 1, 3
          S = S + gu(m,n) * Tab(m,n)
       end do
    end do


    ! dE_gf_volume = (alpha h^ij T_ij + T_tt / alpha - beta^i beta^j T_ij / alpha) sqrt(detgd)

    dE_gf_volume(i,j,k)   = (alph * S + aux) * sqrt(detgd)

    ! dJz = (-y p_x + x p_y) sqrt(detgd)        + permutations
    dJz_gf_volume(i,j,k)  = (-y1 * mom(1) + x1 * mom(2)) * sqrt(detgd)
    dJx_gf_volume(i,j,k)  = (-z1 * mom(2) + y1 * mom(3)) * sqrt(detgd)
    dJy_gf_volume(i,j,k)  = (-x1 * mom(3) + z1 * mom(1)) * sqrt(detgd)
    
    ! dI_ij = rho * x^i x^j * alpha * sqrt(detgd)
    dIxx_gf_volume(i,j,k) = alph * rho * x1 * x1 * sqrt(detgd)
    dIxy_gf_volume(i,j,k) = alph * rho * x1 * y1 * sqrt(detgd)
    dIxz_gf_volume(i,j,k) = alph * rho * x1 * z1 * sqrt(detgd)
    dIyy_gf_volume(i,j,k) = alph * rho * y1 * y1 * sqrt(detgd)
    dIyz_gf_volume(i,j,k) = alph * rho * y1 * z1 * sqrt(detgd)
    dIzz_gf_volume(i,j,k) = alph * rho * z1 * z1 * sqrt(detgd)

  end do
  end do
  end do
end subroutine UAv_Analysis_gfs


subroutine UAv_Analysis_IntegrateVol( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  ! num_out_vals: number of output values for a given reduction
  CCTK_INT, PARAMETER :: num_in_fields = 10, num_out_vals = 1  
  CCTK_REAL out_vals(num_in_fields*num_out_vals)
  
  CCTK_INT ierr
  CCTK_INT reduction_handle, varid(num_in_fields)

  CCTK_INT i
  CCTK_REAL dV

  character(len=*), PARAMETER :: thorn_str = "UAv_Analysis::"
  CCTK_INT, PARAMETER :: thorn_strlen = LEN(thorn_str), var_strlen = 14 ! 14 for dIxy_gf_volume (largest so far)
  CCTK_INT, PARAMETER :: full_strlen = thorn_strlen + var_strlen 
  character(len=full_strlen), dimension(num_in_fields) :: varnames

  varnames = [character(len=full_strlen) :: &
               thorn_str//"dE_gf_volume  ", &
               thorn_str//"dJx_gf_volume ", &
               thorn_str//"dJy_gf_volume ", &
               thorn_str//"dJz_gf_volume ", &
               thorn_str//"dIxx_gf_volume", &
               thorn_str//"dIxy_gf_volume", &
               thorn_str//"dIxz_gf_volume", &
               thorn_str//"dIyy_gf_volume", &
               thorn_str//"dIyz_gf_volume", &
               thorn_str//"dIzz_gf_volume"]


  if (do_analysis_every .le. 0) then
     return
  end if

  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
     return
  end if
  
  call CCTK_ReductionHandle(reduction_handle, 'sum')
  if (reduction_handle < 0) then
     call CCTK_WARN(0, 'Could not obtain a handle for sum reduction')
  end if


  ! Get var IDs
  do i = 1, num_in_fields
      call CCTK_VarIndex(varid(i), varnames(i))
      if (varid(i) < 0) then
         call CCTK_WARN(0, 'Could not get index to grid array '//varnames(i))
      end if
  end do


  ! Call reduction
  ! do a sum over all processors
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, num_out_vals, CCTK_VARIABLE_REAL, &
       out_vals, num_in_fields, &
       varid(1), & ! E
       varid(2), varid(3), varid(4), & ! J_i
       varid(5), varid(6), varid(7), varid(8), varid(9), varid(10)) ! I_ij
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing the auxiliary XX_gf_volume grid functions.')
  end if

  ! the multiplication with the volume element needs to be done here
  dV = cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  do i = 1,num_in_fields
      out_vals(i) = out_vals(i) * dV
  end do

  total_energy = out_vals(1)

  total_angular_momentum_x = out_vals(2)
  total_angular_momentum_y = out_vals(3)
  total_angular_momentum_z = out_vals(4)

  Ixx = out_vals(5)
  Ixy = out_vals(6)
  Ixz = out_vals(7)
  Iyy = out_vals(8)
  Iyz = out_vals(9)
  Izz = out_vals(10)

  ! write(*,*) 'total_energy = ', total_energy

  ! note that the integrated values just obtained do *not* take into account
  ! any grid symmetries that may be present. so one needs to remember to
  ! multiply by the appropriate factor (ie, 8 in the case of using just an
  ! octant, etc) a posteriori.

end subroutine UAv_Analysis_IntegrateVol
