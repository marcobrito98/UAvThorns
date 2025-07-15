/* ParamCheck.c : Check that the parameters provided make sense                  */
/* ============================================================================= */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void UAv_Analysis_ParamCheck(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Check validity of grid scalars provided
  if (track_origin_from_grid_scalar) {
    CCTK_INT index;
    // x
    index = CCTK_VarIndex (track_origin_source_x);
    if (index < 0) {
      CCTK_VPARAMWARN("Could not get index of chosen track_origin_source_x: %s.", track_origin_source_x);
    }
    if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
      CCTK_VPARAMWARN("Chosen track_origin_source_x: %s is not a grid scalar.", track_origin_source_x);
    }
    // y
    index = CCTK_VarIndex (track_origin_source_y);
    if (index < 0) {
      CCTK_VPARAMWARN("Could not get index of chosen track_origin_source_y: %s.", track_origin_source_y);
    }
    if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
      CCTK_VPARAMWARN("Chosen track_origin_source_y: %s is not a grid scalar.", track_origin_source_y);
    }
    // z
    index = CCTK_VarIndex (track_origin_source_z);
    if (index < 0) {
      CCTK_VPARAMWARN("Could not get index of chosen track_origin_source_z: %s.", track_origin_source_z);
    }
    if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
      CCTK_VPARAMWARN("Chosen track_origin_source_z: %s is not a grid scalar.", track_origin_source_z);
    }
  }

}