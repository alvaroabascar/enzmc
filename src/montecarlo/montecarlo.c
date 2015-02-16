#include <lvmrq.h>
#include <random.h>
#include <stdio.h>
#include <matrix.h>
#include <time.h>

/* Input:
 *
 *      modelo
 *      desviacion
 *      valor real parametros
 *      ¿guess parametros?
 *      numero de puntos
 *      puntos
 *
 * Output:
 *      media de los parámetros
 *      varianza parámetros
 *      numero de convergencias
 *
 * Proceso:
 *      1. Obtener curva perfecta
 *      2. Añadir error aleatorio
 *      3. Regresión no lineal sobre curva pseudo-experimental
 *      4. Guardar parámetros
 *      5. Repetir pasos 2-4 N veces
 *      6. Calcular media y varianza
 *
 */

#define abs(x) (x >= 0 ? x : -1*x)

int montecarlo(double model(double X[], double p[]), /* model function */
                double params[], /* real value of the parameters */
                double guess[], /* guess of the parameters */
                double dev, /* an estimation of the deviation */
                int nsims, /* number of simulations */
                int n, /* number of points */
                int nvars, /* number of independent variables of the model */
                int m, /* number of parameters of the model */
                int mfit, /* number of adjustable parameters */
                int fit[m], /* ptr to array which indicates what parameters
                              * will be adjusted (fit[i] = 1) and what ones
                              * will be fixed (fit[i] = 0) */
                double *xi[nvars], /* data points (indep. variables) */
                double *params_mean,
                double *variances,
                FILE *fp)
{
    int i, j, nsuccess, niters, skip;
    long seed;
    double y[n]; /* dependent variable values at the n given points */
    double yi[n]; /* same, with error added */
    double params_guess[m];
    double sig[n];
    double *sigp[n];
    double covar[mfit][mfit];
    double results[2];
    /* build array of y values and array of deviations */
    for (i = 0; i < n; i++) {
        y[i] = model(xi[i], params);
        sig[i] = dev;
    }
    *sigp = sig;
    for (i = 0; i < m; i++) {
        params_mean[i] = variances[i] = 0;
    }
    if (fp != NULL) {
        fprintf(fp, "y = ");
        vector_printf(fp, n, y);
    }
    seed = time(NULL); /* set seed */
    for (i = nsuccess = 0; i < nsims; i++) {
        printf("\rRunning simulations: %d%%",  i*100 / nsims);
        if (fp != NULL)
            fprintf(fp, "\n- Sim. num. %d\n", i);
        /* add error */
        for (j = 0; j < n; j++) {
            yi[j] = y[j] + boxmuller(&seed) * dev;
        }
        if (fp != NULL) {
            fprintf(fp, "yi = ");
            vector_printf(fp, n, yi);
        }
        vcopy(m, params_guess, guess); /* original guess array */
        niters = lvmrq(n, m, mfit, nvars, xi, yi, params_guess, fit, model, covar, sigp, results);
        /* Check whether the result is valid. A large variance (in the
         * covariance matrix) indicates that something went wrong */
        skip = 0;
        for (j = 0; j < m; j++)
            if (covar[j][j] > 10*params[j] || abs(params_guess[j]) > 100*abs(params[j]))
                skip = 1;
        if (skip)
            continue;
        nsuccess++;
        for (j = 0; j < m; j++) {
            params_mean[j] += params_guess[j]; /* add new value */
            variances[j] += (params_guess[j] - params[j])*
                            (params_guess[j] - params[j]);
        }
        if (fp != NULL) {
            fprintf(fp, "- Number of iterations:  %d\n", niters);
            fprintf(fp, "- Parameters:\n");
            vector_printf(fp, m, params_guess);
            fprintf(fp, "- Chi2: %.4e\n", results[0]);
            fprintf(fp, "- Chi2(niters) - Chi2(niters-1):  %.4e\n", results[1]);
            fprintf(fp, "- Matrix of covariances:\n");
            mfprint(fp, mfit, mfit, covar);
        }
    }
    printf("\rRunning simulations: 100%%\n");
    for (i = 0; i < m; i++) {
        params_mean[i] /= nsuccess;
        variances[i] /= nsuccess;
    }
    return nsuccess;
}
