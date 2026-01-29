#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine UAv_Analysis_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr
  CCTK_INT, parameter :: one = 1
  CCTK_INT, parameter :: bndsize = 3

  if (do_analysis_every .le. 0) then
     return
  end if

  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
     return
  endif

  if (compute_density_rho == 1) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
          "UAv_Analysis::density_rho", "flat")
     if (ierr < 0)                                                           &
          call CCTK_ERROR("Failed to register BC for UAv_Analysis::density_rho!")
  end if

  if (compute_density_p == 1) then 
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
          "UAv_Analysis::density_p", "flat")
     if (ierr < 0)                                                           &
          call CCTK_ERROR("Failed to register BC for UAv_Analysis::density_p!")
  end if

end subroutine UAv_Analysis_Boundaries
