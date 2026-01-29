
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

// Register a mask with values "inside", "buffer", "outside", to be used with
// AHFinderDirect
int UAv_Analysis_RegisterMask(void)
{
  int ierr;

  const char *state_list[3] = {"inside", "buffer", "outside"};

  ierr = SpaceMask_RegisterType("mask", 3, state_list);

  if (ierr)
    CCTK_WARN(0, "Failed to register the mask!");

  return 0;
}
