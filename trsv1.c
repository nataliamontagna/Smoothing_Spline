#include "mex.h"
#include <math.h>
#include <string.h>


void egrss_err_param(const char * func, int info) {
  fprintf(stderr,"** On entry to %s, parameter number %d had an illegal value\n", func, -info);
  return;
}

int egrss_dtrsv(
  const char trans,
  const int p,
  const int n,
  const double *restrict U,
  const int ldu,
  const double *restrict W,
  const int ldw,
  const double *restrict c,
  const int incc,
  double *restrict b,
  const int incb,
  double *restrict workspace
) {

    int info=0;

    /* Test input parameters */
    if (trans != 'N' && trans != 'T') info = -1;
    if (p < 0) info = -2;
    if (n < 0) info = -3;
    if (U==NULL) info = -4;
    if (ldu < p) info = -5;
    if (W==NULL) info = -6;
    if (ldw < p) info = -7;
    if (c == NULL && incc != 0) info = -9;
    if (b == NULL) info = -10;
    if (incb == 0) info = -11;
    if (workspace==NULL) info = -12;
    if (info != 0) {
      const char *func = __func__;
      egrss_err_param(func, info);
      return info;
    }

    /* Initialize workspace */
    for (int k=0;k<p;k++) workspace[k] = 0.0;
    double * z = workspace;

    /* Solve L*x = b or L'*x = b */
    if (c != NULL) {
      if (trans == 'N') {
        for (int i=0; i<n; i++) {
          /* Compute u_i'*z */
          double dot=0.0;
          for (int k=0;k<p;k++) {
            dot += U[k]*z[k];
          }

          /* Update b_i := (b_i-u_i'*z)/c_i */
          *b = ((*b)-dot)/(*c);

          /* Update z */
          for (int k=0;k<p;k++) {
            z[k] += W[k]*(*b);
          }

          U += ldu; W += ldw;
          b += incb; c += incc;
        }
      } else {
        /* trans == 'T' */
        U += (n-1)*ldu; W += (n-1)*ldw;
        b += (n-1)*incb; c += (n-1)*incc;
        for (int i=0; i<n; i++) {
          /* Compute v_i'*z */
          double dot=0.0;
          for (int k=0;k<p;k++) {
            dot += W[k]*z[k];
          }

          /* Update b_i */
          *b = ((*b) - dot)/(*c);

          /* Update z */
          for (int k=0;k<p;k++) {
            z[k] += U[k]*(*b);
          }

          U -= ldu; W -= ldw;
          b -= incb; c -= incc;
        }
      }
    } else {
      /* d == NULL */
      if (trans == 'N') {
        for (int i=0;i<n;i++) {

          /* Compute b_i := (b_i - u_i'*z)/(v_i'*u_i) */
          double dot = 0.0;
          for (int k=0;k<p;k++) {
            *b -= U[k]*z[k];
            dot += U[k]*W[k];
          }
          *b /= dot;

          /* Update z */
          for (int k=0;k<p;k++) {
            z[k] += W[k]*(*b);
          }

          U += ldu; W += ldw;
          b += incb;
        }
      } else {
        /* trans == 'T' */
        U += (n-1)*ldu; W += (n-1)*ldw;
        b += (n-1)*incb;
        for (int i=0;i<n;i++) {

          /* Compute b_i := (b_i - u_i'*z)/(v_i'*u_i) */
          double dot = 0.0;
          for (int k=0;k<p;k++) {
            *b -= W[k]*z[k];
            dot += U[k]*W[k];
          }
          *b /= dot;

          /* Update z */
          for (int k=0;k<p;k++) {
            z[k] += U[k]*(*b);
          }

          U -= ldu; W -= ldw;
          b -= incb;
        }
      }
    }

    return info;
  }



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

   /* Check for proper number of arguments. */
   if ( nrhs > 5 ) {
      mexErrMsgIdAndTxt("egrss:trsv:invalidNumInputs","Too many input arguments.");
  	}
  	else if (nrhs < 3) {
  		mexErrMsgIdAndTxt("egrss:trsv:invalidNumInputs","Not enough input arguments.");
  	}
  	if (nlhs > 1) {
  		mexErrMsgIdAndTxt("egrss:trsv:invalidNumOutputs","Invalid number of output arguments.");
  	}

    /* Check input types */
  	if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) ) {
  		mexErrMsgIdAndTxt("egrss:trsv:invalidInput","Inputs must be real-valued numeric arrays.");
  	}

  	mwSize p = mxGetM(prhs[0]);
  	mwSize n = mxGetN(prhs[0]);

    /* Validate input arguments */
    if (p != mxGetM(prhs[1]) || n != mxGetN(prhs[1])) {
      mexErrMsgIdAndTxt("egrss:trsv:matchdims","Dimensions of Ut and Wt do not match.");
    }
    char trans = 'N';
    double * c = NULL;
    int inc = 0;
    if (nrhs == 3) {
      if (!(mxGetM(prhs[2]) == n && mxGetN(prhs[2]) == 1)) {
        mexErrMsgIdAndTxt("egrss:trsv:matchdims","Dimensions of b are invalid.");
    	 }
    }
    else if (nrhs == 4) {
      if (mxIsChar(prhs[3])){
        trans = *mxGetChars(prhs[3]);
        if (mxGetNumberOfElements(prhs[3]) != 1 || (trans != 'N' && trans != 'T')) {
          mexErrMsgIdAndTxt("egrss:trsv:invalidInput","Invalid input: expected 'N' or 'T'.");
        }
      }
      else if (mxIsDouble(prhs[3])) {
        if (!(mxGetM(prhs[3]) == n && mxGetN(prhs[3]) == 1)) {
          mexErrMsgIdAndTxt("egrss:trsv:matchdims","Dimensions of b are invalid.");
        }
        c = mxGetPr(prhs[2]);
        inc = 1;
      }
      else {
        mexErrMsgIdAndTxt("egrss:trsv:invalidInput","Invalid input: expected real-valued numeric array or a charater.");
      }
    }
    else if (nrhs == 5) {
      if (mxIsChar(prhs[4])){
         trans = *mxGetChars(prhs[4]);
         if (mxGetNumberOfElements(prhs[4]) != 1 || (trans != 'N' && trans != 'T')) {
           mexErrMsgIdAndTxt("egrss:trsv:invalidInput","Invalid input: expected 'N' or 'T'.");
         }
      }
      else {
        mexErrMsgIdAndTxt("egrss:trsv:invalidInput","Invalid input: expected 'N' or 'T'.");
      }
      if (mxIsDouble(prhs[3])) {
        if (!(mxGetM(prhs[3]) == n && mxGetN(prhs[3]) == 1)) {
          mexErrMsgIdAndTxt("egrss:trsv:matchdims","Dimensions of b are invalid.");
        }
        c = mxGetPr(prhs[2]);
        inc = 1;
      }
      else {
        mexErrMsgIdAndTxt("egrss:trsv:invalidInput","Invalid input: expected real-valued numeric array.");
      }
    }

  	/* Allocate/initialize output array b */
   plhs[0] = mxDuplicateArray(c == NULL ? prhs[2] : prhs[3]);
   double * b = mxGetPr(plhs[0]);

  	/* Allocate workspace */
  	mxArray * work = mxCreateDoubleMatrix(p,1,mxREAL);
   double * ws = mxGetPr(work);

  	/* Compute matrix-vector product and clean up */
  	int info = egrss_dtrsv(trans,p,n,mxGetPr(prhs[0]),p,mxGetPr(prhs[1]),p,c,inc,b,1,ws);
  	mxDestroyArray(work);
  	if (info) {
  		mexErrMsgIdAndTxt("egrss:trsv:failure","Factorization failed.");
  	}
}