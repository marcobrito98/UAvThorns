/* ParamCheck.c : Check that the parameters provided make sense                  */
/* ============================================================================= */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void UAv_IDBHScalarHair_ParamCheck(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
  We check if the parameters are compatible.
  Namely, we check that "HairyBH" and "ScalarBS" are used consistently, 
  along with the appropriate angular frequency and horizon radius parameters.
  */

  // Hairy BH simulation: use OmegaH and rH, not omega_BS.
  if (CCTK_Equals(initial_data, "HairyBH")) {
    // Consistent keywords
    if (CCTK_Equals(initial_lapse, "ScalarBS")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = HairyBH' and 'initial_lapse = ScalarBS' is not allowed.");
    }
    if (CCTK_Equals(initial_shift, "ScalarBS")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = HairyBH' and 'initial_shift = ScalarBS' is not allowed.");
    }

    // Angular frequency
    if (omega_BS > 1e-16) {
      CCTK_PARAMWARN("Using 'initial_data = HairyBH' with a non-zero 'omega_BS' is not allowed. "
                   "Unset 'omega_BS' and check that you set 'OmegaH' and 'rH' properly.");
    }
  }

  // Scalar BS simulation: use omega_BS, not OmegaH nor rH.
  if (CCTK_Equals(initial_data, "ScalarBS")) {
    // Consistent keywords
    if (CCTK_Equals(initial_lapse, "HairyBH")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = ScalarBS' and 'initial_lapse = HairyBH' is not allowed.");
    }
    if (CCTK_Equals(initial_shift, "HairyBH")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = ScalarBS' and 'initial_shift = HairyBH' is not allowed.");
    }

    // Angular frequency
    if (OmegaH > 1e-16) {
      CCTK_PARAMWARN("Using 'initial_data = ScalarBS' with a non-zero 'OmegaH' is not allowed. "
                   "Unset 'OmegaH' and check that you set 'omega_BS' properly.");
    }
    // Horizon radius
    if (rH > 1e-16) {
      CCTK_PARAMWARN("Using 'initial_data = ScalarBS' with a non-zero 'rH' is not allowed. "
                   "Unset 'rH'.");
    }
  }
}
