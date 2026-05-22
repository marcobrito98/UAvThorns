/* ParamCheck.c : Check that the parameters provided make sense                  */
/* ============================================================================= */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void UAv_Analysis_ParamCheck(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // -------------------------------
  // MULTIPATCH
  // -------------------------------
  
  // For information, see WARNING in param.ccl

  // For multipatch / Llama grids, we need the volume form computed in Coordinates
  // This is activated by the parameter Coordinates::store_volume_form = yes

  // Trick to check multipatch usage without having to inherit the thorn
  if( CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification") ) {
    
    // -1 is the "error" return value of ParameterGetType() (see repos/flesh/src/main/Parameters.c and repos/flesh/src/include/cctk_Parameter.h)
    CCTK_INT type = -1; 
    const CCTK_INT* store_volume_form_ptr = CCTK_ParameterGet("store_volume_form", "Coordinates", &type);

    if (store_volume_form_ptr == NULL || type != PARAMETER_BOOLEAN) {
      CCTK_ERROR("Problem acquiring pointer to Coordinates::store_volume_form parameter.");
    }
    else if (!(*store_volume_form_ptr)) {
      CCTK_VPARAMWARN("You are using a multipatch system, but you set Coordinates::store_volume_form = no. "   
                      "UAv_Analysis thorn needs to use the volume form for integrations. " 
                      "Please set Coordinates::store_volume_form = yes.");
    }
    // else, we're good to go
  }


  // -------------------------------
  // ORIGIN TRACKING
  // -------------------------------

  // Check validity of grid scalars provided
  if (track_origin_from_grid_scalar) {
    CCTK_INT index;

    // x
    index = CCTK_VarIndex (track_origin_source_x);
    if (!CCTK_Equals(track_origin_source_x, "fixed")) {
      if (index < 0) {
        CCTK_VPARAMWARN("Could not get index of chosen track_origin_source_x: %s.", track_origin_source_x);
      }
      else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
        CCTK_VPARAMWARN("Chosen track_origin_source_x: %s is not a grid scalar.", track_origin_source_x);
      }
    }
    
    // y
    index = CCTK_VarIndex (track_origin_source_y);
    if (!CCTK_Equals(track_origin_source_y, "fixed")) {
      if (index < 0) {
        CCTK_VPARAMWARN("Could not get index of chosen track_origin_source_y: %s.", track_origin_source_y);
      }
      else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
        CCTK_VPARAMWARN("Chosen track_origin_source_y: %s is not a grid scalar.", track_origin_source_y);
      }
    }
    
    // z
    index = CCTK_VarIndex (track_origin_source_z);
    if (!CCTK_Equals(track_origin_source_z, "fixed")) {
      if (index < 0) {
        CCTK_VPARAMWARN("Could not get index of chosen track_origin_source_z: %s.", track_origin_source_z);
      }
      else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
        CCTK_VPARAMWARN("Chosen track_origin_source_z: %s is not a grid scalar.", track_origin_source_z);
      }
    }
  }

}