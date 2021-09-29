#include "mex.h"
#include <math.h>
#include <string.h>


void egrss_err_param(const char * func, int info) {
  fprintf(stderr,"** On entry to %s, parameter number %d had an illegal value\n", func, -info);
  return;
}

int egrss_dpotrf(
  const int p,
  const int n,
  const double *restrict U,
  const int ldu,
  double *restrict V,
  const int ldv,
  double *restrict d,
  const int incd,
  double *restrict workspace
) {

  int info=0;
  double tmp;

  /* Test input parameters */
  if (p < 0) info = -1;
  if (n < 0) info = -2;
  if (U==NULL) info = -3;
  if (ldu < p) info = -4;
  if (V==NULL) info = -5;
  if (ldv < p) info = -6;
  if (d == NULL && incd != 0) info = -8;
  if (d != NULL && incd == 0) info = -8;
  if (workspace==NULL) info = -9;
  if (info != 0) {
    const char *func = __func__;
    egrss_err_param(func, info);
    return info;
  }

  /* Initialize workspace */
  for (int k=0;k<p*p;k++) workspace[k] = 0.0;

  /* Cholesky factorization */
  if (d==NULL) {d = &tmp;}
  for (int i=0; i<n; i++) {

    /* v_i = v_i - P*u_i */
    for (int j=0;j<p;j++) {
      for (int k=0;k<p;k++) {
        V[k] -= workspace[p*j+k]*U[j];
      }
    }

    /* Scale v_i by 1/(u_i'*v_i) or 1/(u_i'*v_i + d_i) */
    tmp = 0.0;
    for (int k=0;k<p;k++) *d += U[k]*V[k];
    if (*d<=0.0) return i+1;
    *d = sqrt(*d);
    tmp = 1.0/(*d);
    for (int k=0;k<p;k++) V[k] *= tmp;

    /* P = P + v_i*v_i' */
    for (int j=0;j<p;j++) {
      for (int k=0;k<p;k++) {
        workspace[p*j+k] += V[j]*V[k];
      }
    }
    U += ldu;
    V += ldv;
    d += incd;
  }
  return info;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* Check for proper number of arguments. */
  if ( nrhs > 3 ) {
    mexErrMsgIdAndTxt("egrss:potrf:invalidNumInputs","Too many input arguments.");
	}
	else if (nrhs < 2) {
		mexErrMsgIdAndTxt("egrss:potrf:invalidNumInputs","Not enough input arguments.");
	}
	if ((nrhs == 3 && nlhs != 2) || nlhs > 2) {
			mexErrMsgIdAndTxt("egrss:potrf:invalidNumOutputs","Invalid number of output arguments.");
	}

  /* Check input types */
	if ( !mxIsDouble(prhs[0]) ||
    	 !mxIsDouble(prhs[1]) ||
			 (nrhs == 3 && !mxIsDouble(prhs[2])) ) {
		mexErrMsgIdAndTxt("egrss:potrf:invalidInput","Inputs must be real-valued numeric arrays.");
	}

	mwSize p = mxGetM(prhs[0]);
	mwSize n = mxGetN(prhs[0]);

  /* Validate input arguments */
  if (p != mxGetM(prhs[1]) || n != mxGetN(prhs[1])) {
    mexErrMsgIdAndTxt("egrss:potrf:matchdims","Dimensions of Ut and Vt do not match.");
  }
	if (nrhs == 3 && !(mxGetM(prhs[2]) == 1 || mxGetN(prhs[2]) == 1)) {
    mexErrMsgIdAndTxt("egrss:potrf:matchdims","Size of d is invalid.");
	}
  if (nrhs == 3 && (mxGetNumberOfElements(prhs[2]) != n && mxGetNumberOfElements(prhs[2]) != 1)) {
    mexErrMsgIdAndTxt("egrss:potrf:matchdims","Length of d is invalid.");
  }

	/* Allocate/initialize output array W */
  plhs[0] = mxDuplicateArray(prhs[1]);
	double * c = NULL;
	int inc = 0;
  if (nlhs == 2) {
		/* Allocate/initialize output array c */
		if (mxGetNumberOfElements(prhs[2]) == n) {
			plhs[1] = mxDuplicateArray(prhs[2]);
			c = mxGetPr(plhs[1]);
		}
		else {
			/* Input is a scalar */
			plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
			c = mxGetPr(plhs[1]);
			double val = mxGetScalar(prhs[2]);
			for (int i=0;i<n;i++) c[i] = val;
		}
		inc = 1;
	}

	/* Allocate workspace */
	mxArray * work = mxCreateDoubleMatrix(p,p,mxREAL);
  double * ws = mxGetPr(work);

	/* Compute factorization and clean up */
	int info = egrss_dpotrf(p,n,mxGetPr(prhs[0]),p,mxGetPr(plhs[0]),p,c,inc,ws);
	mxDestroyArray(work);
	if (info) {
		mexErrMsgIdAndTxt("egrss:potrf:failure","Factorization failed.");
	}
}