#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// Initialize auxiliary members
void UAv_Initialization (CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  // -------------------------------
  // MULTIPATCH
  // -------------------------------

  // For information, see WARNING in param.ccl

  // Setting the flags has to be done here, cannot be done in PARAMCHECK

  // Trick to check multipatch usage without having to inherit the thorn
  *is_multipatch = CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification");
  
  // default (non Llama grid): don't use volume form
  *use_volume_form = 0; 

  // TODO: If there's ever a flag in Coordinates that indicates that the volume_form is trivial (e.g. is_cartesian),
  // we could check it and set *use_volume_form = 0 in that case.

  if (*is_multipatch > 0) {
    // If the simulation uses multipatch, make sure that the volume form variable is stored.
    // Supposedly, this should be handled correctly in Coordinates, 
    // so this should be redundant with the ParamCheck on store_volume_form.
    CCTK_INT* volume_form_state_ptr = CCTK_VarDataPtr(cctkGH, 0, "Coordinates::volume_form_state");
    // CCTK_VINFO("volume_form_state = %d", *volume_form_state_ptr);
    
    if (volume_form_state_ptr != NULL) {
      // No volume form actually stored (should not happen if ParamCheck passed)
      if (*volume_form_state_ptr != 1) {
        CCTK_ERROR("The patch system that you are using does not store the volume form, although you set Coordinates::store_volume_form = yes. "     
                   "UAv_Analysis thorn will not be able to use the volume form for integrations. " 
                   "There is likely a problem in thorn Coordinates.");
      }
      // else, we're good to go
      else {
        CCTK_INFO("Using volume form from multipatch system for integrations in UAv_Analysis.");
        *use_volume_form = 1;
      }
    }
    else {
      CCTK_ERROR("Problem acquiring pointer to Coordinates::volume_form_state variable.");
    }
  }


  // -------------------------------
  // ORIGIN TRACKING
  // -------------------------------

  // We can already initialize the coordinates
  *x0 = origin_x;
  *y0 = origin_y;
  *z0 = origin_z;

  // Initialize origin tracking if needed
  // origin_from_grid_scalar_index not allocated in schedule.ccl if track_origin_from_grid_scalar = no
  if (track_origin_from_grid_scalar) {
    // Get the index of variables. It's not supposed to change during the simulation (I think).
    // Validity of parameters should have been checked in ParamCheck.
    // Since validity has been checked, a negative value here should mean "fixed".

    // x source
    *origin_from_grid_scalar_index_x = CCTK_VarIndex (track_origin_source_x);
    // y source
    *origin_from_grid_scalar_index_y = CCTK_VarIndex (track_origin_source_y);
    // z source
    *origin_from_grid_scalar_index_z = CCTK_VarIndex (track_origin_source_z);
  }
  
  // Info
  CCTK_VINFO("Origin used in the analysis:");
  // x
  if (track_origin_from_grid_scalar && *origin_from_grid_scalar_index_x >= 0) {
    CCTK_VINFO("x0: %s", track_origin_source_x);
  }
  else {
    CCTK_VINFO("x0 = %g (fixed)", *x0);
  }
  // y
  if (track_origin_from_grid_scalar && *origin_from_grid_scalar_index_y >= 0) {
    CCTK_VINFO("y0: %s", track_origin_source_y);
  }
  else {
    CCTK_VINFO("y0 = %g (fixed)", *y0);
  }
  // z
  if (track_origin_from_grid_scalar && *origin_from_grid_scalar_index_z >= 0) {
    CCTK_VINFO("z0: %s", track_origin_source_z);
  }
  else {
    CCTK_VINFO("z0 = %g (fixed)", *z0);
  }
}