#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


/*
  NOTES:
  - If we want to go more general than fixed coordinates or tracking from a grid scalar,
  maybe we could have a UAv_Set_origin() function, that would call or incorporate this one.
  - It could be tempting to just set the x0, y0, z0 pointers to the corresponding variables,
  at Initialization [though we would need to be really sure that the pointers don't change
  during a simulation, in any case]. However, these are const pointers, so actually they're
  read-only, we can't set them. This function takes a negligible time anyway.
*/

// if (track_origin_from_grid_scalar) in schedule 
void UAv_Track_origin (CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (do_analysis_every <= 0) return;
  if (cctk_iteration % do_analysis_every != 0) return;

  // Track the coordinates of the origin from the chosen grid scalars
  *x0 = * (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, *origin_from_grid_scalar_index_x);
  *y0 = * (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, *origin_from_grid_scalar_index_y);
  *z0 = * (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, *origin_from_grid_scalar_index_z);


  // // TODO: Could be added to some verbose parameter...
  // CCTK_VINFO ("Tracking origin in UAv_Analysis:");
  // CCTK_VINFO ("x0 = %g", *x0);
  // CCTK_VINFO ("y0 = %g", *y0);
  // CCTK_VINFO ("z0 = %g", *z0);
}