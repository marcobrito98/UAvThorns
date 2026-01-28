#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine ComplexProca_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1

  ! The outgoing (radiative) boundary conditions are being handled from calls to
  ! the NewRad infraestructure. Here we register all BCs as 'none', which
  ! enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ComplexProcaBase::E1i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::E1i!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ComplexProcaBase::A1i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::A1i!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ComplexProcaBase::Aphi1", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::Aphi1!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ComplexProcaBase::Zeta1", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::Zeta1!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ComplexProcaBase::E2i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::E2i!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ComplexProcaBase::A2i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::A2i!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ComplexProcaBase::Aphi2", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::Aphi2!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ComplexProcaBase::Zeta2", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ComplexProcaBase::Zeta2!")

end subroutine ComplexProca_Boundaries
