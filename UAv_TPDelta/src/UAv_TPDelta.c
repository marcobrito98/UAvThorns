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

          CCTK_REAL tmp[4][4];
          for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < 4; ++b) {
              tmp[a][b] = 0.0;
            }
          }

          //isto vem do novo two punctures
          tmp[1][1] = gxx[ind];
          tmp[1][2] = gxy[ind];
          tmp[1][3] = gxz[ind];
          tmp[2][1] = gxy[ind];
          tmp[2][2] = gyy[ind];
          tmp[2][3] = gyz[ind];
          tmp[3][1] = gxz[ind];
          tmp[3][2] = gyz[ind];
          tmp[3][3] = gzz[ind];
 
          

          // voltar a por a metrica de uma so puncture
          gxx[ind] = gxx_ref[ind]; //+ delta_scale * (gxx[ind] - gxx_ref[ind]);
          gxy[ind] = gxy_ref[ind]; //+ delta_scale * (gxy[ind] - gxy_ref[ind]);
          gxz[ind] = gxz_ref[ind]; //+ delta_scale * (gxz[ind] - gxz_ref[ind]);
          gyy[ind] = gyy_ref[ind]; //+ delta_scale * (gyy[ind] - gyy_ref[ind]);
          gyz[ind] = gyz_ref[ind]; //+ delta_scale * (gyz[ind] - gyz_ref[ind]);
          gzz[ind] = gzz_ref[ind]; //+ delta_scale * (gzz[ind] - gzz_ref[ind]);

          //agora volto a gravar o novo two punctures nestas variaveis
          gxx_ref[ind] = tmp[1][1];
          gxy_ref[ind] = tmp[1][2];
          gxz_ref[ind] = tmp[1][3];
          gyy_ref[ind] = tmp[2][2];
          gyz_ref[ind] = tmp[2][3];
          gzz_ref[ind] = tmp[3][3];
        }

        if (apply_to_curv) {

          CCTK_REAL tmp[4][4];
          for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < 4; ++b) {
              tmp[a][b] = 0.0;
            }
          }

          //isto vem do novo two punctures
          tmp[1][1] = kxx[ind];
          tmp[1][2] = kxy[ind];
          tmp[1][3] = kxz[ind];
          tmp[2][1] = kxy[ind];
          tmp[2][2] = kyy[ind];
          tmp[2][3] = kyz[ind];
          tmp[3][1] = kxz[ind];
          tmp[3][2] = kyz[ind];
          tmp[3][3] = kzz[ind];

          kxx[ind] = kxx_ref[ind]; //+ delta_scale * (kxx[ind] - kxx_ref[ind]);
          kxy[ind] = kxy_ref[ind]; //+ delta_scale * (kxy[ind] - kxy_ref[ind]);
          kxz[ind] = kxz_ref[ind]; //+ delta_scale * (kxz[ind] - kxz_ref[ind]);
          kyy[ind] = kyy_ref[ind]; //+ delta_scale * (kyy[ind] - kyy_ref[ind]);
          kyz[ind] = kyz_ref[ind]; //+ delta_scale * (kyz[ind] - kyz_ref[ind]);
          kzz[ind] = kzz_ref[ind]; //+ delta_scale * (kzz[ind] - kzz_ref[ind]);

          //agora volto a gravar o novo two punctures nestas variaveis
          kxx_ref[ind] = tmp[1][1];
          kxy_ref[ind] = tmp[1][2];
          kxz_ref[ind] = tmp[1][3];
          kyy_ref[ind] = tmp[2][2];
          kyz_ref[ind] = tmp[2][3];
          kzz_ref[ind] = tmp[3][3];

        }
      }
    }
  }
}
