#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void UAv_TPDelta_CaptureBase(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        gxx_ref[ind] = gxx[ind];
        gxy_ref[ind] = gxy[ind];
        gxz_ref[ind] = gxz[ind];
        gyy_ref[ind] = gyy[ind];
        gyz_ref[ind] = gyz[ind];
        gzz_ref[ind] = gzz[ind];

        kxx_ref[ind] = kxx[ind];
        kxy_ref[ind] = kxy[ind];
        kxz_ref[ind] = kxz[ind];
        kyy_ref[ind] = kyy[ind];
        kyz_ref[ind] = kyz[ind];
        kzz_ref[ind] = kzz[ind];
      }
    }
  }
}

void UAv_TPDelta_Assemble(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        if (apply_to_metric) {
          gxx[ind] = gxx_ref[ind] + delta_scale * (gxx[ind] - gxx_ref[ind]);
          gxy[ind] = gxy_ref[ind] + delta_scale * (gxy[ind] - gxy_ref[ind]);
          gxz[ind] = gxz_ref[ind] + delta_scale * (gxz[ind] - gxz_ref[ind]);
          gyy[ind] = gyy_ref[ind] + delta_scale * (gyy[ind] - gyy_ref[ind]);
          gyz[ind] = gyz_ref[ind] + delta_scale * (gyz[ind] - gyz_ref[ind]);
          gzz[ind] = gzz_ref[ind] + delta_scale * (gzz[ind] - gzz_ref[ind]);
        }

        if (apply_to_curv) {
          kxx[ind] = kxx_ref[ind] + delta_scale * (kxx[ind] - kxx_ref[ind]);
          kxy[ind] = kxy_ref[ind] + delta_scale * (kxy[ind] - kxy_ref[ind]);
          kxz[ind] = kxz_ref[ind] + delta_scale * (kxz[ind] - kxz_ref[ind]);
          kyy[ind] = kyy_ref[ind] + delta_scale * (kyy[ind] - kyy_ref[ind]);
          kyz[ind] = kyz_ref[ind] + delta_scale * (kyz[ind] - kyz_ref[ind]);
          kzz[ind] = kzz_ref[ind] + delta_scale * (kzz[ind] - kzz_ref[ind]);
        }
      }
    }
  }
}
