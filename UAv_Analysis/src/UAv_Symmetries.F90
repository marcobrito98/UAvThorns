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

end subroutine UAv_Analysis_Symmetries
