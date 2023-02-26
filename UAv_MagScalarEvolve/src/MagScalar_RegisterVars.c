
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void MagScalar_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs, var;

  // register evolution and rhs gridfunction groups with MoL

  /* metric and extrinsic curvature */
  group = CCTK_GroupIndex("ADMBase::lapse");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::shift");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  /* Ei and rhs_Ei */
  group = CCTK_GroupIndex("MagScalarBase::Ei");
  rhs   = CCTK_GroupIndex("MagScalarEvolve::rhs_Ei");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ai and rhs_Ai */
  group = CCTK_GroupIndex("MagScalarBase::Ai");
  rhs   = CCTK_GroupIndex("MagScalarEvolve::rhs_Ai");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Aphi and rhs_Aphi */
  var   = CCTK_VarIndex("MagScalarBase::Aphi");
  rhs   = CCTK_VarIndex("MagScalarEvolve::rhs_Aphi");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Zeta and rhs_Zeta */
  var   = CCTK_VarIndex("MagScalarBase::Zeta");
  rhs   = CCTK_VarIndex("MagScalarEvolve::rhs_Zeta");
  ierr += MoLRegisterEvolved(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
