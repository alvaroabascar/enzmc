#include "lvmrq.h"
#include <gaussjbs.h>
#include <matrix.h>
#include <mathlib.h>
#include <stdarg.h>

#define CONVERGENCE 1e-20
#define ITERLIM 500
#define LAMBDA_FACTOR 10
#define LAMBDA_START 1e-3
#define CHISQLIM 1

/* lvmrq: non-linear regression algorithm, using the Levenberg-Marquardt
 * method.
 *
 * Input parameters:
 *
 * "n"                     -> number of data points
 * "m"                     -> number of parameters of the model
 * "nvars"                 -> number of independent variables in the model
 * "xi[n][nvars]","yi[n]"  -> data points
 * "a[m]" -> is the set of "m" parameters, from which the first "mfit" will
 *           be adjusted, and the next (m - mfit) will be kept fixed.
 * "*fit[m]" -> pointer to array indicating what parameters to adjust. If fit[i]
 * is set to 0, the parameter is kept fixed. If it is set to 1, it is adjusted.
 * "f(double xi[], double params[])" -> the model function, which accepts an
 *                                      array of independent variables "xi[]"
 *                                      and an array of parameters "params[]"
 * "covar[mfit][mfit]" -> a matrix into which the covariances will be stored
 * "knownsig" -> set it to 1 if an array of standard deviations is passed. If
 *               not, set it to 0.
 * "..." -> pass an array of "n" doubles containing the standard deviation of
 *          each point, in case that it is known.
 *
 * Return values:
 * Once the function has been called, "a[m]" contains the adjusted parameters,
 * and covar[mfit][mfit] contains the matrix of covariances.
 * the return value is the the last value of chi square
 * results[2] holds the value of chi2 and the convergence value
 */
double lvmrq(
           int n,                     /* number of points */
           int m,                     /* number of parameters */
           int mfit,
           int nvars,                 /* number of independent variables */
           double xi[n][nvars],       /* data points */
           double yi[n],
           double a0[m],               /* parameters (guess) */
           int *fitp[m],               /* ptr to array indicating what parameters
                                       * to fix (0) and what ones to adjust (1)
                                       */
           double f(double x[], double params[]), /* model function */
           double covar[][mfit],
           double *sigp[n],           /* pointer to array of deviations, set
                                         it to NULL if they are not known */
           double results[2])
{
    /* definitions here */
    int i, j, k, iters = 0;
    int fit[m],       /* indicates what parameters to fit */
        a_order[m];   /* indicates the original order of the params */

    double chisq,
           a[m],               /* reordered parameters, adjustable first */
           alpha[mfit][mfit],  /* (1/2) * hessian matrix */
           beta[mfit][1],      /* (-1/2)*gradient vector */
           yfit[n],            /* will hold the fitted values for each xi */
           dyda[n][mfit],      /* derivatives. Each row contains the derivatives
                                * at one point with respect to each parameter */
           lambda, ainc[mfit], /* ainc = correction over the parameters */
           anew[m],            /* new values of the parameters (a + ainc) */
           sig[n],             /* array of standard deviations */
           tmp, conv;
    /* Build the array of deviations */
    if (sigp != NULL) { /* if deviations are known */
        for (i = 0; i < n; i++) {
            sig[i] = (*sigp)[i];
        }
    } else { /* deviations are not known, set them to 1 */
        for (i = 0; i < n; i++) {
            sig[i] = 1;
        }
    }
    /* Build the array of parameters to fit (1) or fix (0) */
    if (fitp != NULL) {
        for (i = 0; i < m; i++) {
            fit[i] = (*fitp)[i];
        }
    } else {
        for (i = 0; i < m; i++) {
            fit[i] = 1;
        }
    }
    /* Build an array with the adjustable parameters first */
    for (i = j = 0; i < m; i++) {
        if (fit[i] == 1) {
            a[j] = a0[i];
            a_order[j] = i;
            j++;
        }
    }
    for (i = 0; i < m; i++) {
        if (fit[i] == 0) {
            a[j] = a0[i];
            a_order[j] = i;
            j++;
        }
    }

    /* this function is used as an interface to the provided functions,
     * to calculate the fitted values and the derivatives
     */
    void func(int n, int m, int mfit, int nvars, double xi[n][nvars],
              double yi[n], double a[], double dyda[n][mfit], double yfit[n])
    {
        int i, k;
        double a_local[m];
        for (i = 0; i < m; i++) {
            a_local[i] = a[a_order[i]];
        }
        for (i = 0; i < n; i++) {
            /* obtain fitted ys */
            yfit[i] = f(xi[i], a_local);
            for (k = 0; k < mfit; k++) {
            /* obtain derivatives with respect to each parameter a[k] at each
             * abscise xi. dyda has dimensions n by mfit. Each row contains
             * the derivatives at one point with respect to each parameter */
                dyda[i][k] = dfda(xi[i], a_local, f, a_order[k]);
            }
        }
    }

    /* call func to fill yfit and dyda */
    func(n, m, mfit, nvars, xi, yi, a, dyda, yfit);
    /* calculate the initial value of chi square */
    chisq = chisquare(n, yi, yfit, sig);
    /* set lambda to a low value (ex. 0.001) */
    lambda = LAMBDA_START;
    conv = CONVERGENCE + 1; /* just to make true the following conditional */
    while (iters < ITERLIM && (conv > CONVERGENCE || conv < 0) && chisq) {
        iters++;
        /* calculate alpha and beta, note that alpha is symmetric */
        buildAlphaBeta(n, mfit, lambda, dyda, alpha, beta, sig, yi, yfit);
        /* solve alpha*ainc = beta to get ainc
         * beta contains the correction over a */
        gaussj(mfit, 1, alpha, beta); 
        /* update anew */
        for (i = 0; i < m; i++) {
            anew[i] = i < mfit ? (a[i] + beta[i][0]) : a[i];
        }
        /* calculate chisquare(a + ainc) */
        func(n, m, mfit, nvars, xi, yi, anew, dyda, yfit);
        tmp = chisquare(n, yi, yfit, sig);
        conv = chisq - tmp;
        /* worse fit */
        if (conv <= 0 && chisq != 0) {
            lambda *= LAMBDA_FACTOR;
            func(n, m, mfit, nvars, xi, yi, a, dyda, yfit);
            /* new iteration with the same parameter set "a", but now
             * yfit and dyda correspond to "anew". Logic says I should
             * recover the old yfit and dyda, but it only seems to harm
             * the performance */
        /* better fit */
        } else if (conv > CONVERGENCE) {
            vcopy(mfit, a, anew); /* copy anew to a */
            chisq = tmp;          /* update chisq */
            lambda /= LAMBDA_FACTOR;
        }
    }

    /* if the deviations were not known, estimate them as the variance of the
     * residues */
    if (!sigp) {
        tmp = 0;
        for (i = 0; i < n; i++) {
            sig[0] = (yi[i] - f(xi[i], a));
            sig[0] *= sig[0];
            tmp += sig[0];
        }
        tmp /= n;
        tmp = sqrt(tmp);
        for (i = 0; i < n; i++) {
            sig[i] = tmp;
        }
    }

    /* build alpha with lambda = 0 and calculate its inverse (the matrix
     * of covariances) */
    lambda = 0;
    func(n, m, mfit, nvars, xi, yi, anew, dyda, yfit); /* update dyda, yfit */
    buildAlphaBeta(n, mfit, lambda, dyda, alpha, beta, sig, yi, yfit);
    for (i = 0; i < mfit; i++) { /* build identity matrix */
        for (j = 0; j < mfit; j++) {
            covar[i][j] = (i != j ? 0 : 1);
        }
    }
    gaussj(mfit, mfit, alpha, covar); /* builds inverse of alpha, saves it in covar */

    for (i = 0; i < m; i++) {
        a0[a_order[i]] = a[i];
    }
    if (VERBOSE) {
        printf("\n****[ lvmrq configuration values ]****\n");
        printf("\n-Maximum iterations:         %d\n", ITERLIM);
        printf("\n-Convergence limit :         %.4e\n", CONVERGENCE);
        printf("\n-Initial value of lambda:      %.4e\n", LAMBDA_START);
        printf("\n-Lambda modificator factor:  %d\n", LAMBDA_FACTOR);
        printf("****************************************\n");
        printf("\n****[ lvmrq result values ]****\n");
        printf("\n-Number of iterations: %d\n", iters);
        printf("\n-Final value of chi2: %.4e\n", chisq);
        printf("\n-Convergence value (chi2(iters) - chi2(iters-1)): %.10e\n", conv);
        printf("\n-Final set of parameters:\n");
        vprint(m, a0);
        printf("\n-Covariance matrix:\n");
        mprint(mfit, mfit, covar);
        printf("****************************************\n");
    }
    results[0] = chisq;
    results[1] = -conv;
    return iters;
}
           
/* Given n values (yi[0...n]) with standard deviations sig[0...n], and
 * n fitted values (y[0...n]), results the chi square of this data set:
 *
 * chi square = sum([(yi-yfit)^2]/sig)
 */
double chisquare(int n, double yi[n], double yfit[n], double sig[n])
{
    int i;
    double chisq, tmp;
    for (i = chisq = 0; i < n; i++) {
        tmp = (yi[i] - yfit[i]) / sig[i];
        chisq += tmp * tmp;
    }
    return chisq;
}

/* buildAlphaBeta: builds alpha (1/2*hessian) and beta (-1/2*gradient)
 */
void buildAlphaBeta(int n, int mfit, int lambda, double dyda[][mfit],
                    double alpha[mfit][mfit], double beta[mfit][1],
                    double sig[n], double yi[n], double yfit[n])
{
    double tmp;
    int i, j, k;
    for (i = 0; i < mfit; i++) {
        for (j = 0; j <= i; j++) {
            alpha[i][j] =  0;
            for (k = 0; k < n; k++) {
                tmp = 1/sig[k];
                tmp *= tmp;
                alpha[i][j] += dyda[k][i]*dyda[k][j] * tmp;
            }
            if (i != j) {
                alpha[j][i] = alpha[i][j];
            } else {
                alpha[i][j] *= (1+lambda);
            }
        }
        beta[i][0] = 0;
        for (k = 0; k < n; k++) {
            beta[i][0] += ((yi[k] - yfit[k]) * dyda[k][i]) / (sig[k]*sig[k]);
        }
    }
}
