/* ParamCheck.c : Check that the parameters provided make sense                  */
/* ============================================================================= */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void UAv_IDBHProcaHair_ParamCheck(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
  We check if the parameters are compatible.
  Namely, we check that "ProcaHairyBH" and "ProcaBS" are used consistently, 
  along with the appropriate angular frequency and horizon radius parameters.
  */

  // Hairy BH simulation: use OmegaH and rH, not omega_BS.
  if (CCTK_Equals(initial_data, "ProcaHairyBH")) {
    // Consistent keywords
    if (CCTK_Equals(initial_lapse, "ProcaBS")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = ProcaHairyBH' and 'initial_lapse = ProcaBS' is not allowed.");
    }
    if (CCTK_Equals(initial_shift, "ProcaBS")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = ProcaHairyBH' and 'initial_shift = ProcaBS' is not allowed.");
    }

    // Angular frequency
    if (omega_BS > 1e-16) {
      CCTK_PARAMWARN("Using 'initial_data = ProcaHairyBH' with a non-zero 'omega_BS' is not allowed. "
                   "Unset 'omega_BS' and check that you set 'OmegaH' and 'rH' properly.");
    }
  }

  // Proca BS simulation: use omega_BS, not OmegaH nor rH.
  if (CCTK_Equals(initial_data, "ProcaBS")) {
    // Consistent keywords
    if (CCTK_Equals(initial_lapse, "ProcaHairyBH")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = ProcaBS' and 'initial_lapse = ProcaHairyBH' is not allowed.");
    }
    if (CCTK_Equals(initial_shift, "ProcaHairyBH")) {
      CCTK_PARAMWARN("Using parameter 'initial_data = ProcaBS' and 'initial_shift = ProcaHairyBH' is not allowed.");
    }

    // Angular frequency
    if (OmegaH > 1e-16) {
      CCTK_PARAMWARN("Using 'initial_data = ProcaBS' with a non-zero 'OmegaH' is not allowed. "
                   "Unset 'OmegaH' and check that you set 'omega_BS' properly.");
    }
    // Horizon radius
    if (rH > 1e-16) {
      CCTK_PARAMWARN("Using 'initial_data = ProcaBS' with a non-zero 'rH' is not allowed. "
                   "Unset 'rH'.");
    }
  }
}
