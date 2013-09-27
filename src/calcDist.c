/*Beta-binomial distribution in C */
#include <R.h>
#include <Rmath.h>
/* n1 & n2 the length of the array x1,y1 and x2,y2 */

void calcDist(double *x1, double *y1,double *x2, double *y2, double *r, 
	      int *n1){
  int n = n1[0];

    for (int i = 0; i < n; i++) {
      r[i] = sqrt( pow(x1[i]-x2[0],2.0) + pow(y1[i]-y2[0],2.0));
    }

}



