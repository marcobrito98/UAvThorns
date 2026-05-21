#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "TargetTracker.h"

///////////////////////////////////////////////////////////////////////////////
// First initialization of one target
///////////////////////////////////////////////////////////////////////////////
void InitializeOneTarget (CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Get the index of variables. It's not supposed to change during the simulation (I think).
    // Validity of parameters should have been checked in ParamCheck.
    // Since validity has been checked, we can use the negative return value for "fixed" targets.
    // WARNING: Since the parameters are steerable, make sure the rest is robust
    
    
    // x source
    target_id_x[itarget] = CCTK_VarIndex (target_x[itarget]);
    // y source
    target_id_y[itarget] = CCTK_VarIndex (target_y[itarget]);
    // z source
    target_id_z[itarget] = CCTK_VarIndex (target_z[itarget]);

    // initial position
    target_loc_x[itarget] = initial_x[itarget];
    target_loc_y[itarget] = initial_y[itarget];
    target_loc_z[itarget] = initial_z[itarget];

    // Set initial active status of the target based on the current parameters.
    TargetActivationCondition(CCTK_PASS_CTOC, itarget);
    
    if (is_tracked[itarget]) {
        
        CCTK_VINFO("Initialized %s target %d with sources:", is_active[itarget] ? "active" : "inactive", itarget);
        
        // Check if fixed target in each dimension

        if (target_id_x[itarget] >= 0) {
            CCTK_VINFO("x: %s", target_x[itarget]);
        } else {
            CCTK_VINFO("x = %g (fixed)", initial_x[itarget]);
        }

        if (target_id_y[itarget] >= 0) {
            CCTK_VINFO("y: %s", target_y[itarget]);
        } else {
            CCTK_VINFO("y = %g (fixed)", initial_y[itarget]);
        }

        if (target_id_z[itarget] >= 0) {
            CCTK_VINFO("z: %s", target_z[itarget]);
        } else {
            CCTK_VINFO("z = %g (fixed)", initial_z[itarget]);
        }

    } // end if track
}

///////////////////////////////////////////////////////////////////////////////
// Helper function to recover the name of the target variable for one dimension (x, y, z) from its checkpointed ID.
///////////////////////////////////////////////////////////////////////////////
void RecoverOneTargetNameOneDim (CCTK_ARGUMENTS, const struct TargetInfoBundleOneDim bundle) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Unpack bundle
    const CCTK_INT itarget          = bundle.itarget;
    const CCTK_INT* ptr_current_id  = bundle.ptr_current_id;
    const char* tgt_name            = bundle.tgt_name;
    const char* dim_name            = bundle.dim_name;

    const CCTK_INT recovery_ID      = CCTK_VarIndex (tgt_name);
    char checkpointed_Name[150];
    if (*ptr_current_id >= 0) {
        sprintf(checkpointed_Name, "%s::%s", CCTK_ImpFromVarI(*ptr_current_id), CCTK_VarName(*ptr_current_id));
    }
    else {
        sprintf(checkpointed_Name, "fixed");
    }

    // Output info
    if (recovery_ID != *ptr_current_id || !CCTK_Equals(tgt_name, checkpointed_Name)) {
        char message[1000];
        sprintf(message, "At recovery of target %d, parameter 'target_%s' is '%s' (ID: %d),\n"
                         "    but the variable with chekcpointed index 'target_id_%s[%d]' = %d is '%s'.\n"
                         "    Setting target_%s[%d] = %s.",
                         itarget, dim_name, tgt_name, recovery_ID,
                         dim_name, itarget, *ptr_current_id, checkpointed_Name,
                         dim_name, itarget, checkpointed_Name);
        if (verbose) {
            CCTK_VINFO("%s", message);
        }
        else {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "%s", message);
        }
    } //enf if mismatch

    // Set parameter to checkpointed value
    // A recovered negative index should mean a fixed target
    // The ParameterSet will take effect at the next DECLARE_CCTK_PARAMETERS
    char param_name[100];
    sprintf(param_name, "target_%s[%d]", dim_name, itarget);
    CCTK_ParameterSet(param_name, "TargetTracker", (*ptr_current_id >= 0) ? checkpointed_Name : "fixed");
}

///////////////////////////////////////////////////////////////////////////////
// (Partial) recovery of parameters for one target (after restart from checkpoint)
// This is not completely robust for all parameters and all scenarios.
///////////////////////////////////////////////////////////////////////////////
void RecoverOneTarget (CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // The ParameterSet will take effect at the next DECLARE_CCTK_PARAMETERS

    /* Parameter: track 
     * Most important to recover since it controls the tracking and is always steerable.
     */
    // Info
    if (track[itarget] != is_tracked[itarget]) {
        char message[1000];
        sprintf(message, "At recovery of target %d, parameter 'track' is '%s', but the internal flag 'is_tracked' is '%s'.\n"
                         "    Setting track[%d] = %s.",
                         itarget, track[itarget] ? "yes" : "no", is_tracked[itarget] ? "yes" : "no",
                         itarget, is_tracked[itarget] ? "yes" : "no");
        if (verbose) {
            CCTK_VINFO("%s", message);
        }
        else {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "%s", message);
        }
    }
    // Set
    char param_name[100];
    sprintf(param_name, "track[%d]", itarget);
    CCTK_ParameterSet(param_name, "TargetTracker", is_tracked[itarget] ? "yes" : "no");

    //////////////////////////////////////////////

    /* Parameters: target variables
     * We recover the names from their checkpointed IDs.
     * NOTE: If these parameters become not always steerable, this step is probably unncessary.
     */
    const struct TargetInfoBundleOneDim bundle_x = {itarget, &target_id_x[itarget], target_x[itarget], "x"};
    RecoverOneTargetNameOneDim(CCTK_PASS_CTOC, bundle_x);
    const struct TargetInfoBundleOneDim bundle_y = {itarget, &target_id_y[itarget], target_y[itarget], "y"};
    RecoverOneTargetNameOneDim(CCTK_PASS_CTOC, bundle_y);
    const struct TargetInfoBundleOneDim bundle_z = {itarget, &target_id_z[itarget], target_z[itarget], "z"};
    RecoverOneTargetNameOneDim(CCTK_PASS_CTOC, bundle_z);

    //////////////////////////////////////////////

    /* Parameters that we don't try to recover:
     * - force_params_at_recovery: STEERABLE=RECOVER
     * - track_every: STEERABLE=RECOVER
     * - start/stop_tracking_after_time: STEERABLE=RECOVER
     * - which_surface_to_store_info: STEERABLE=NEVER
     * - initial_x/y/z: STEERABLE=RECOVER, and they are ignored anyway if force_params_at_recovery = no
     * - verbose: STEERABLE=ALWAYS, but it's not critical (and the user should know if they try to steer it mid-run)
     */



}

///////////////////////////////////////////////////////////////////////////////
// Routine called at the first initialization of the simulation (i.e. not recovery)
///////////////////////////////////////////////////////////////////////////////
void TargetTracker_Initialization(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        InitializeOneTarget(CCTK_PASS_CTOC, itarget);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Routine called at recovery
///////////////////////////////////////////////////////////////////////////////
void TargetTracker_Recovery(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    /* Fundamental ambiguity: have parameters changed during the run and the parameter is reused identically,
     * or has the user purposefully changed the parameter file?
     * Having a fully consistent and robust architecture seems involved, so we have to make assumptions and concessions.
     *
     * By default, we want to keep the recovered parameters, but without having to fully duplicate parameters into the interface.
     * (I don't fully understand if we can get the checkpointed values of (steerable) parameters and store them somehow.)
     * The user can force a re-read of parameters at recovery (for those that we try to infer).
     */
    
    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        if (force_params_at_recovery[itarget]) {
            if (verbose) {
                CCTK_VINFO("Forcing re-initialization of target %d from parameter file at recovery.", itarget);
            }
            InitializeOneTarget(CCTK_PASS_CTOC, itarget);
        }
        else {
            RecoverOneTarget(CCTK_PASS_CTOC, itarget);
        }
    }
}
