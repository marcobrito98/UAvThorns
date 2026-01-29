#include "cctk_Arguments.h"
#include "cctk.h"

subroutine ComplexProca_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_E1x    = 0
  rhs_E1y    = 0
  rhs_E1z    = 0

  rhs_A1x    = 0
  rhs_A1y    = 0
  rhs_A1z    = 0

  rhs_Aphi1  = 0

  rhs_Zeta1  = 0
  
  rhs_E2x    = 0
  rhs_E2y    = 0
  rhs_E2z    = 0

  rhs_A2x    = 0
  rhs_A2y    = 0
  rhs_A2z    = 0

  rhs_Aphi2  = 0

  rhs_Zeta2  = 0

end subroutine ComplexProca_zero_rhs
