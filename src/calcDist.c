/*Beta-binomial distribution in C */
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
/* n1 & n2 the length of the array x1,y1 and x2,y2 */

SEXP calcDist(SEXP x1, SEXP y1,SEXP x2, SEXP y2){
	
  int n0 = length(x1);
  int n1 = length(x2);
  double *X1 = REAL(x1);
  double *X2 = REAL(x2);
  double *Y1 = REAL(y1);
  double *Y2 = REAL(y2);
  SEXP r, ans;

  PROTECT(r = allocMatrix(REALSXP, n0, n1));
  PROTECT(ans = allocVector(REALSXP, 1));
  double *Ans = REAL(ans);
  double *R = REAL(r);
    for (int j = 0; j < n1; j++) {
      for(int i = 0; i < n0; i++){
 	Ans[0] = sqrt( pow(X1[i]-X2[j],2.0) + pow(Y1[i]-Y2[j],2.0) );
	R[i+n0*j] = Ans[0];
	/*SET_VECTOR_ELT(r, i+n0*j,ans);*/ 
      }
    }
    
    UNPROTECT(2);

    return r;
}



