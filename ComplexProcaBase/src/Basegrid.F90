! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine complexproca_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ComplexProcaBase::E1x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ComplexProcaBase::E1y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ComplexProcaBase::E1z" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ComplexProcaBase::A1x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ComplexProcaBase::A1y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ComplexProcaBase::A1z" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ComplexProcaBase::Zeta1" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ComplexProcaBase::Aphi1" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ComplexProcaBase::E2x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ComplexProcaBase::E2y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ComplexProcaBase::E2z" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ComplexProcaBase::A2x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ComplexProcaBase::A2y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ComplexProcaBase::A2z" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ComplexProcaBase::Zeta2" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ComplexProcaBase::Aphi2" )

end subroutine complexproca_symmetries
