/* TwoPunctures:  File  "utilities.h"*/

#include <math.h>

#include "cctk.h"

#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164	/* Pi/2*/
#define Piq 0.78539816339744830961566084582	/* Pi/4*/

#define TINY 1.0e-20
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define nrerror UAV2_nrerror
#define ivector UAV2_ivector
#define dvector UAV2_dvector
#define imatrix UAV2_imatrix
#define dmatrix UAV2_dmatrix
#define d3tensor UAV2_d3tensor
#define free_ivector UAV2_free_ivector
#define free_dvector UAV2_free_dvector
#define free_imatrix UAV2_free_imatrix
#define free_dmatrix UAV2_free_dmatrix
#define free_d3tensor UAV2_free_d3tensor

#define minimum2 UAV2_minimum2
#define minimum3 UAV2_minimum3
#define maximum2 UAV2_maximum2
#define maximum3 UAV2_maximum3
#define pow_int UAV2_pow_int

#define chebft_Zeros UAV2_chebft_Zeros
#define chebft_Extremes UAV2_chebft_Extremes
#define chder UAV2_chder
#define chebev UAV2_chebev
#define fourft UAV2_fourft
#define fourder UAV2_fourder
#define fourder2 UAV2_fourder2
#define fourev UAV2_fourev

#define norm1 UAV2_norm1
#define norm2 UAV2_norm2
#define scalarproduct UAV2_scalarproduct

void nrerror (char error_text[]);
int *ivector (long nl, long nh);
CCTK_REAL *dvector (long nl, long nh);
int **imatrix (long nrl, long nrh, long ncl, long nch);
CCTK_REAL **dmatrix (long nrl, long nrh, long ncl, long nch);
CCTK_REAL ***d3tensor (long nrl, long nrh, long ncl, long nch, long ndl,
		    long ndh);
void free_ivector (int *v, long nl, long nh);
void free_dvector (CCTK_REAL *v, long nl, long nh);
void free_imatrix (int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix (CCTK_REAL **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor (CCTK_REAL ***t, long nrl, long nrh, long ncl, long nch,
		    long ndl, long ndh);

int minimum2 (int i, int j);
int minimum3 (int i, int j, int k);
int maximum2 (int i, int j);
int maximum3 (int i, int j, int k);
int pow_int (int mantisse, int exponent);

void chebft_Zeros (CCTK_REAL u[], int n, int inv);
void chebft_Extremes (CCTK_REAL u[], int n, int inv);
void chder (CCTK_REAL *c, CCTK_REAL *cder, int n);
CCTK_REAL chebev (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c[], int m, CCTK_REAL x);
void fourft (CCTK_REAL *u, int N, int inv);
void fourder (CCTK_REAL u[], CCTK_REAL du[], int N);
void fourder2 (CCTK_REAL u[], CCTK_REAL d2u[], int N);
CCTK_REAL fourev (CCTK_REAL *u, int N, CCTK_REAL x);


CCTK_REAL norm1 (CCTK_REAL *v, int n);
CCTK_REAL norm2 (CCTK_REAL *v, int n);
CCTK_REAL scalarproduct (CCTK_REAL *v, CCTK_REAL *w, int n);

void UAV2_nrerror (char error_text[]);
int *UAV2_ivector (long nl, long nh);
CCTK_REAL *UAV2_dvector (long nl, long nh);
int **UAV2_imatrix (long nrl, long nrh, long ncl, long nch);
CCTK_REAL **UAV2_dmatrix (long nrl, long nrh, long ncl, long nch);
CCTK_REAL ***UAV2_d3tensor (long nrl, long nrh, long ncl, long nch, long ndl,
			    long ndh);
void UAV2_free_ivector (int *v, long nl, long nh);
void UAV2_free_dvector (CCTK_REAL *v, long nl, long nh);
void UAV2_free_imatrix (int **m, long nrl, long nrh, long ncl, long nch);
void UAV2_free_dmatrix (CCTK_REAL **m, long nrl, long nrh, long ncl, long nch);
void UAV2_free_d3tensor (CCTK_REAL ***t, long nrl, long nrh, long ncl, long nch,
			      long ndl, long ndh);

int UAV2_minimum2 (int i, int j);
int UAV2_minimum3 (int i, int j, int k);
int UAV2_maximum2 (int i, int j);
int UAV2_maximum3 (int i, int j, int k);
int UAV2_pow_int (int mantisse, int exponent);

void UAV2_chebft_Zeros (CCTK_REAL u[], int n, int inv);
void UAV2_chebft_Extremes (CCTK_REAL u[], int n, int inv);
void UAV2_chder (CCTK_REAL *c, CCTK_REAL *cder, int n);
CCTK_REAL UAV2_chebev (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c[], int m, CCTK_REAL x);
void UAV2_fourft (CCTK_REAL *u, int N, int inv);
void UAV2_fourder (CCTK_REAL u[], CCTK_REAL du[], int N);
void UAV2_fourder2 (CCTK_REAL u[], CCTK_REAL d2u[], int N);
CCTK_REAL UAV2_fourev (CCTK_REAL *u, int N, CCTK_REAL x);

CCTK_REAL UAV2_norm1 (CCTK_REAL *v, int n);
CCTK_REAL UAV2_norm2 (CCTK_REAL *v, int n);
CCTK_REAL UAV2_scalarproduct (CCTK_REAL *v, CCTK_REAL *w, int n);
