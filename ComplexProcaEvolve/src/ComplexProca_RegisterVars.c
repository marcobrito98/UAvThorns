
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ComplexProca_RegisterVars(CCTK_ARGUMENTS)
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
  group = CCTK_GroupIndex("ComplexProcaBase::E1i");
  rhs   = CCTK_GroupIndex("ComplexProcaEvolve::rhs_E1i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ai and rhs_Ai */
  group = CCTK_GroupIndex("ComplexProcaBase::A1i");
  rhs   = CCTK_GroupIndex("ComplexProcaEvolve::rhs_A1i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Aphi and rhs_Aphi */
  var   = CCTK_VarIndex("ComplexProcaBase::Aphi1");
  rhs   = CCTK_VarIndex("ComplexProcaEvolve::rhs_Aphi1");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Zeta and rhs_Zeta */
  var   = CCTK_VarIndex("ComplexProcaBase::Zeta1");
  rhs   = CCTK_VarIndex("ComplexProcaEvolve::rhs_Zeta1");
  ierr += MoLRegisterEvolved(var, rhs);
  
  /* Ei and rhs_Ei */
  group = CCTK_GroupIndex("ComplexProcaBase::E2i");
  rhs   = CCTK_GroupIndex("ComplexProcaEvolve::rhs_E2i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ai and rhs_Ai */
  group = CCTK_GroupIndex("ComplexProcaBase::A2i");
  rhs   = CCTK_GroupIndex("ComplexProcaEvolve::rhs_A2i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Aphi and rhs_Aphi */
  var   = CCTK_VarIndex("ComplexProcaBase::Aphi2");
  rhs   = CCTK_VarIndex("ComplexProcaEvolve::rhs_Aphi2");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Zeta and rhs_Zeta */
  var   = CCTK_VarIndex("ComplexProcaBase::Zeta2");
  rhs   = CCTK_VarIndex("ComplexProcaEvolve::rhs_Zeta2");
  ierr += MoLRegisterEvolved(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
