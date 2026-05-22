#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void TargetTracker_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        
        if (which_surface_to_store_info[itarget] != -1 && which_surface_to_store_info[itarget] >= nsurfaces) {
            CCTK_VPARAMWARN("For target %d, surface index is greater than the number of spherical surfaces.", itarget);
        }
        
        
        // --------------------------------
        // CHECK TARGET VARIABLES
        // --------------------------------
        CCTK_INT index;
        // x
        index = CCTK_VarIndex (target_x[itarget]);
        if (!CCTK_Equals(target_x[itarget], "fixed")) {
            if (index < 0) {
                CCTK_VPARAMWARN("Could not get index of target %d x variable: %s.", itarget, target_x[itarget]);
            }
            else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
                CCTK_VPARAMWARN("Target %d x variable: %s is not a grid scalar.", itarget, target_x[itarget]);
            }
        }
        // y
        index = CCTK_VarIndex (target_y[itarget]);
        if (!CCTK_Equals(target_y[itarget], "fixed")) {
            if (index < 0) {
                CCTK_VPARAMWARN("Could not get index of target %d y variable: %s.", itarget, target_y[itarget]);
            }
            else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
                CCTK_VPARAMWARN("Target %d y variable: %s is not a grid scalar.", itarget, target_y[itarget]);
            }
        }
        // z
        index = CCTK_VarIndex (target_z[itarget]);
        if (!CCTK_Equals(target_z[itarget], "fixed")) {
            if (index < 0) {
                CCTK_VPARAMWARN("Could not get index of target %d z variable: %s.", itarget, target_z[itarget]);
            }
            else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
                CCTK_VPARAMWARN("Target %d z variable: %s is not a grid scalar.", itarget, target_z[itarget]);
            }
        }
    }
}
