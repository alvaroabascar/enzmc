#define VERBOSE 0
#ifndef __NLR__
  #define __NLR__

double chisquare(int n, double yi[n], double yfit[n], double sig[n]);

double lvmrq(
           int n,                       /* number of points */
           int m,                       /* number of parameters */
           int mfit,                    /* number of parameters to adjust */
           int nvars,                   /* number of independent variables */
           double xi[n][nvars],         /* data points */
           double yi[n],
           double a[m],                 /* parameters (guess) */
           int *fit[m],                 /* fit[i]=1 --> adjust a[i].
                                         * fit[i]=0 --> keep a[i] fixed. */
           double f(double x[], double params[]), /* model function */
           double covar[][mfit],
           double *sigp[n],             /* pointer to array of deviations, set
                                           to NULL if they are not known */
           double results[2]);          /* results[0] = final chi2.
                                         * results[1] = final increment in chi2
                                         */


void buildAlphaBeta(int n, int mfit, int lambda, double dyda[][mfit],
                    double alpha[mfit][mfit], double beta[mfit][1],
                    double sig[n], double yi[n], double yfit[n]);
#endif
