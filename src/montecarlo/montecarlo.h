int montecarlo(double model(double X[], double p[]), /* model function */
                double params[], /* real value of the parameters */
                double guess[], /* guess of the parameters */
                double dev, /* an estimation of the deviation */
                int nsims, /* number of simulations */
                int n, /* number of points */
                int nvars, /* number of independent variables of the model */
                int m, /* number of parameters of the model */
                int mfit, /* number of adjustable parameters */
                int fit[m],  /* array which indicates what parameters
                              * will be adjusted (fit[i] = 1) and what ones
                              * will be fixed (fit[i] = 0) */
                double xi[n][nvars], /* data points (indep. variables) */
                double params_mean[m], 
                double variances[m],
                FILE *fp);
