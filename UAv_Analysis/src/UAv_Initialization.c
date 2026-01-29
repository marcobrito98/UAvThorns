#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// Initialize auxiliary members
void UAv_Initialization (CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (track_origin_from_grid_scalar) {
    // Get the index of variables. It's not supposed to change during the simulation (I think).
    // Validity of parameters should have been checked in ParamCheck.
    // Seems a bit redundant to do this affectation here and not in ParamCheck, but more in the logic.
    
    // x source
    *origin_from_grid_scalar_index_x = CCTK_VarIndex (track_origin_source_x);
    // y source
    *origin_from_grid_scalar_index_y = CCTK_VarIndex (track_origin_source_y);
    // z source
    *origin_from_grid_scalar_index_z = CCTK_VarIndex (track_origin_source_z);

    CCTK_VINFO("Tracking origin used in the analysis with grid scalars.");
    CCTK_VINFO("x0 = %s", track_origin_source_x);
    CCTK_VINFO("y0 = %s", track_origin_source_y);
    CCTK_VINFO("z0 = %s", track_origin_source_z);
  }
  else { // no tracking from grid scalar
    // origin_from_grid_scalar_index not allocated in schedule.ccl in that case

    // We can already initialize the coordinates
    *x0 = origin_x;
    *y0 = origin_y;
    *z0 = origin_z;

    CCTK_VINFO("Using fixed origin in the analysis.");
    CCTK_VINFO("x0 = %g", *x0);
    CCTK_VINFO("y0 = %g", *y0);
    CCTK_VINFO("z0 = %g", *z0);
  }
}