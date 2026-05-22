#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

/*
 * Contains helper functions for the TargetTracker thorn.
 */

/////////////////////////////////////////////////////////////////////////////////
// Helper struc to pass arguments to functions
//////////////////////////////////////////////////////////////////////////////////
struct TargetInfoBundleOneDim {
    CCTK_INT itarget;
    CCTK_INT* ptr_current_id;
    const char* tgt_name;
    const char* dim_name;
};

////////////////////////////////////////////////////////////////////////////////
// Helper function for factorization and consistency.
// Gives the condition for a target to be active based on the current parameters
// and updates the is_active status of the target accordingly.
////////////////////////////////////////////////////////////////////////////////

inline void TargetActivationCondition(CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Update value of internal is_tracked flag based on track parameter (which is STEERABLE=ALWAYS)
    is_tracked[itarget] = track[itarget];

    is_active[itarget] = ( 
        is_tracked[itarget]
        // &&  track_every[itarget] > 0     // Should be satisfied by construction
        &&  (cctk_time >= start_tracking_after_time[itarget])
        &&  (cctk_time <= stop_tracking_after_time[itarget])
    );
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Helper function to print an error message and trigger termination after the time step.
// Deactivates the concerned target.
/////////////////////////////////////////////////////////////////////////////////////////
inline void TerminateHelper(CCTK_ARGUMENTS, const char* message, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;

    CCTK_VWARN(1, "%s", message);
    CCTK_VWARN(1, "Deactivating target %d and triggering termination at iteration %d (simulation time %g).", itarget, cctk_iteration, cctk_time);
    CCTK_TerminateNext(cctkGH);

    is_active[itarget] = 0;
    return;
}