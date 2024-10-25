
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

#define SMALL (1.e-9)

void UAv_ID_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL []);


void UAv_IDProcaBS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* 
    TODO: Implement Proca star.

    Proca Boson Star limit (rH=0) not implemented yet.
    To do it, copy and adapt from UAv_IDBHScalarHair/src/ScalarBS.c
    and UAv_IDBHProcaHair/src/BHProcaHair.c.
    In the scalar case, just taking the BH routine with rH == 0 was not
    correct in a few places, although on the few numerical tests performed 
    it didn't seem too critical. Thus the dedicated separate routine.
    This is probably the case here too, it needs to be checked and done properly.
  */

  CCTK_ERROR("Function UAv_IDProcaBS not implemented yet.");

  return;
}
