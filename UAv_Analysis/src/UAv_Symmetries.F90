#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine UAv_Analysis_Symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_Analysis::density_rho" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "UAv_Analysis::density_px" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "UAv_Analysis::density_py" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "UAv_Analysis::density_pz" )

  ! Initialize the symmetry factors for the volume integrals
  call UAv_Analysis_InitSymmetryFactors( CCTK_PASS_FTOF )

end subroutine UAv_Analysis_Symmetries


! Wrapper to initialize the symmetry factors for the volume integrals
subroutine UAv_Analysis_InitSymmetryFactors( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  character(len=100) :: symmetry_fstr
  CCTK_INT :: symmetry_fstr_len

  if (CCTK_Equals(symmetry, "none") /= 0) then
    call UAv_Analysis_InitSymmetryFactors_None( CCTK_PASS_FTOF )
  
  else if (CCTK_Equals(symmetry, "bitant_x") /= 0) then
    call UAv_Analysis_InitSymmetryFactors_BitantX( CCTK_PASS_FTOF )
  
  else if (CCTK_Equals(symmetry, "bitant_y") /= 0) then
    call UAv_Analysis_InitSymmetryFactors_BitantY( CCTK_PASS_FTOF )
  
  else if (CCTK_Equals(symmetry, "bitant_z") /= 0) then
    call UAv_Analysis_InitSymmetryFactors_BitantZ( CCTK_PASS_FTOF )
  
  ! Default to none, but issue warning
  else 
    call CCTK_FortranString(symmetry_fstr_len, symmetry, symmetry_fstr)
    call CCTK_WARN (CCTK_WARN_ALERT, "Unknown symmetry factors for: '"//symmetry_fstr//"'. Defaulting to 'symmetry=none'." )
    call UAv_Analysis_InitSymmetryFactors_None( CCTK_PASS_FTOF )
  end if

end subroutine UAv_Analysis_InitSymmetryFactors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines to initialize the symmetry factors for the volume integrals for different symmetries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! None / Default
subroutine UAv_Analysis_InitSymmetryFactors_None( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  sym_factor_dE       = 1

  sym_factor_dJz      = 1
  sym_factor_dJx      = 1
  sym_factor_dJy      = 1

  sym_factor_drho     = 1
  sym_factor_dCoMx    = 1
  sym_factor_dCoMy    = 1
  sym_factor_dCoMz    = 1

  sym_factor_dpx      = 1
  sym_factor_dpy      = 1
  sym_factor_dpz      = 1

  sym_factor_dIxx     = 1
  sym_factor_dIxy     = 1
  sym_factor_dIxz     = 1
  sym_factor_dIyy     = 1
  sym_factor_dIyz     = 1
  sym_factor_dIzz     = 1

end subroutine UAv_Analysis_InitSymmetryFactors_None

! Bitant even symmetry across x=0 plane
subroutine UAv_Analysis_InitSymmetryFactors_BitantX( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  sym_factor_dE       = 2

  sym_factor_dJx      = 2
  sym_factor_dJy      = 0
  sym_factor_dJz      = 0

  sym_factor_drho     = 2
  sym_factor_dCoMx    = 0
  sym_factor_dCoMy    = 2
  sym_factor_dCoMz    = 2

  sym_factor_dpx      = 0
  sym_factor_dpy      = 2
  sym_factor_dpz      = 2

  sym_factor_dIxx     = 2
  sym_factor_dIxy     = 0
  sym_factor_dIxz     = 0
  sym_factor_dIyy     = 2
  sym_factor_dIyz     = 2
  sym_factor_dIzz     = 2

end subroutine UAv_Analysis_InitSymmetryFactors_BitantX

! Bitant even symmetry across y=0 plane
subroutine UAv_Analysis_InitSymmetryFactors_BitantY( CCTK_ARGUMENTS )
  
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  sym_factor_dE       = 2

  sym_factor_dJx      = 0
  sym_factor_dJy      = 2
  sym_factor_dJz      = 0

  sym_factor_drho     = 2
  sym_factor_dCoMx    = 2
  sym_factor_dCoMy    = 0
  sym_factor_dCoMz    = 2

  sym_factor_dpx      = 2
  sym_factor_dpy      = 0
  sym_factor_dpz      = 2

  sym_factor_dIxx     = 2
  sym_factor_dIxy     = 0
  sym_factor_dIxz     = 2
  sym_factor_dIyy     = 2
  sym_factor_dIyz     = 0
  sym_factor_dIzz     = 2

end subroutine UAv_Analysis_InitSymmetryFactors_BitantY

! Bitant even symmetry across z=0 plane
subroutine UAv_Analysis_InitSymmetryFactors_BitantZ( CCTK_ARGUMENTS )
  
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  sym_factor_dE       = 2

  sym_factor_dJx      = 0
  sym_factor_dJy      = 0
  sym_factor_dJz      = 2

  sym_factor_drho     = 2
  sym_factor_dCoMx    = 2
  sym_factor_dCoMy    = 2
  sym_factor_dCoMz    = 0

  sym_factor_dpx      = 2
  sym_factor_dpy      = 2
  sym_factor_dpz      = 0

  sym_factor_dIxx     = 2
  sym_factor_dIxy     = 2
  sym_factor_dIxz     = 0
  sym_factor_dIyy     = 2
  sym_factor_dIyz     = 0
  sym_factor_dIzz     = 2

end subroutine UAv_Analysis_InitSymmetryFactors_BitantZ