
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#define MAXBUF 1000  // max number of characters per line in input file. this
                     // should be safe

/* utility routine to check whether given array of chars is an empty line (or
   a line with only whitespace) */
static bool is_empty(const char *s)
{
  while (*s != '\0') {
    if (!isspace(*s))
      return false;
    s++;
  }
  return true;
}

void UAv_IDBHProcaHair_read_data(CCTK_INT *NF_p, CCTK_INT *NX_p, CCTK_REAL Xtmp[], CCTK_REAL thtmp[],
               CCTK_REAL F1[], CCTK_REAL F2[], CCTK_REAL F0[], CCTK_REAL W[],
               CCTK_REAL H1[], CCTK_REAL H2[], CCTK_REAL H3[], CCTK_REAL V[])
{
  DECLARE_CCTK_PARAMETERS;

  FILE *infile;
  /* open input file */
  infile = fopen(infilename, "rt");
  if (infile == NULL) {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
    "Unable to open file %s\n", infilename);
  } else {
    CCTK_VInfo(CCTK_THORNSTRING, "Reading data file %s", infilename);
  }

  /* read data from input file */
  char buf[MAXBUF];
  CCTK_INT NF = 0;     // NF will be the size of the full array
  CCTK_INT NX = 0;     // NX will be the number of X points
  bool first_block = true;
  while (fgets(buf, MAXBUF, infile) > 0) {

    // skip comments
    if (buf[0] == '#')
      continue;

    // empty lines mark a new theta coordinate. use that to count total number
    // of X points.
    if (is_empty(buf)) {
      if (first_block) {
        NX = NF;
        first_block = false;
      }
      continue;
    }

    /* printf("%s\n", buf); */
    sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &Xtmp[NF], &thtmp[NF], &F1[NF], &F2[NF], &F0[NF], &W[NF], &H1[NF], &H2[NF], &H3[NF], &V[NF]);

    // take into account different normalization used for the stress-energy
    // tensor in the input files, which may assume G = 1,
    // whereas for ComplexProcaEvolve thorns it is assumed that 4 pi G = 1.
    if (normalization_Tmunu == 0) { // 4 pi G = 1 in input file
      H1[NF] /= (Xtmp[NF]/(1-Xtmp[NF]));
      H2[NF] *= 1;
      H3[NF] *= 1;
       V[NF] *= 1;
    } else if (normalization_Tmunu == 1) { // G = 1 in input file
      H1[NF] *= 2.0*sqrt(M_PI);
      H2[NF] *= 2.0*sqrt(M_PI);
      H3[NF] *= 2.0*sqrt(M_PI);
       V[NF] *= 2.0*sqrt(M_PI);
    } else if (normalization_Tmunu == 2) { // G = 1 and Tmunu has factor of 0.5 in input file
      H1[NF] *= sqrt(2.0*M_PI);
      H2[NF] *= sqrt(2.0*M_PI);
      H3[NF] *= sqrt(2.0*M_PI);
       V[NF] *= sqrt(2.0*M_PI);
    }

    NF++;
  }

  // the following is a hack to set correctly the number of X points when we
  // only have one block (with the logic above and with just one block of data,
  // NX may never be assigned to NF after NF is properly incremented).
  if (NX == 0)
    NX = NF;

  *NF_p = NF;
  *NX_p = NX;

  fclose(infile);

  CCTK_INFO("Data file read.");

  return;
}
