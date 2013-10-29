/*Beta-binomial distribution in C */
#include <R.h>
#include <Rmath.h>

/* l: The length of the array k */

void dbbf(int *lin, int *k, int *nin,  double *ain, double *bin, double *x)
{
  int l = lin[0];
  int n = nin[0];
  double a = ain[0];
  double b = bin[0];
  int i;
                           

  for (i = 0; i < l; i++) {
    x[i] = exp(lgamma(n+1) - lgamma(k[i]+1) - gamma((n-k[i])+1) + lbeta(k[i]+a, n-k[i] + b) - lbeta(a,b));
  }

}



