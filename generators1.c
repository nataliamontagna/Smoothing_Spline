#include "mex.h"
#include <math.h>
#include <string.h>


void egrss_err_param(const char * func, int info) {
  fprintf(stderr,"** On entry to %s, parameter number %d had an illegal value\n", func, -info);
  return;
}

int egrss_dsplkgr(
    const int p,
    const int n,
    double *restrict U,
    const int ldu,
    double *restrict V,
    const int ldv,
    const double *restrict t,
    const int inct
) {

    int info = 0;

    /* Test input parameters */
    if (p < 0) info = -1;
    if (n < 0) info = -2;
    if (U==NULL) info = -3;
    if (ldu < p) info = -4;
    if (V==NULL) info = -5;
    if (ldv < p) info = -6;
    if (t==NULL) info = -7;
    if (inct == 0) info = -8;
    if (info != 0) {
      const char *func = __func__;
      egrss_err_param(func, info);
      return info;
    }

    // Check monotonicity
    int increasing = t[0]<=t[inct*(n-1)];

    // Swap U and V is sequence is decreasing
    if (!increasing) {
        double * tmp = U;
        U = V;
        V = tmp;
    }

    /* Compute U and V*/
    double thold = *t;
    for (int i=0;i<n;i++) {
        if ((increasing && *t>=thold) || (!increasing && *t<=thold)) {
            thold = *t;
        }
        else {
            /* Not a monotonic sequence */
            info = i+1;
            break;
        }
        double ti = (thold >= 0) ? thold : 0.0;;

        U[p-1] = 1.0;
        for (int k=p-2;k>=0;k--) {
            U[k] = U[k+1]*ti/(p-1-k);
        }
        V[0] = U[0]*ti/p;
        for (int k=1;k<p;k++) {
            V[k] = -V[k-1]*ti/(p+k);
        }

        U += ldu;
        V += ldv;
        t += inct;
    }

    return info;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

   /* Check for proper number of arguments. */
	if ( nrhs > 2 ) {
 		mexErrMsgIdAndTxt("egrss:generators:invalidNumInputs","Too many input arguments.");
	}
	else if (nrhs < 2) {
		mexErrMsgIdAndTxt("egrss:generators:invalidNumInputs","Not enough input arguments.");
	}
	if (nlhs != 2) {
		mexErrMsgIdAndTxt("egrss:potrf:invalidNumOutputs","Invalid number of output arguments.");
	}

	/* Validate input arguments */
   if (!mxIsDouble(prhs[0]) || !(mxGetM(prhs[0]) == 1 || mxGetN(prhs[0]) == 1)) {
	 	mexErrMsgIdAndTxt("egrss:generators:matchdims","Input t must be a vector.");
	}
   if (!mxIsScalar(prhs[1]) || mxGetScalar(prhs[1]) < 1 || (mxGetScalar(prhs[1]) != (int)mxGetScalar(prhs[1]))) {
		mexErrMsgIdAndTxt("egrss:generators:invalidInput","Input p must be a positive integer.");
	}

	mwSize n = mxGetM(prhs[0]) > mxGetN(prhs[0]) ? mxGetM(prhs[0]) : mxGetN(prhs[0]);
	mwSize p = (mwSize) mxGetScalar(prhs[1]);
  double * t = mxGetPr(prhs[0]);

  /* Allocate output arrays Ut and Vt */
 	plhs[0] = mxCreateDoubleMatrix(p,n,mxREAL);
 	plhs[1] = mxCreateDoubleMatrix(p,n,mxREAL);

	/* Compute SS generators */
	int info;
	info = egrss_dsplkgr(p,n,mxGetPr(plhs[0]),p,mxGetPr(plhs[1]),p,t,1);
	if (info > 0) {
		mexErrMsgIdAndTxt("egrss:generators:failure","Input t must be monotonic.");
	}

}