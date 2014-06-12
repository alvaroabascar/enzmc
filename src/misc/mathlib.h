#ifndef __MATH_INC__
#define __MATH_INC__
#include <math.h>
#define abs(x) ((x > 0) ? x : -1*x)

double sum(int n, double data[n]);
double mean(int n, double data[n]);
double std(int n, double data[n]);
double max(int n, double data[n]);
double dfda(double xi[],
            double params[],
            double fun(double x[], double params[]),
            int k);

#endif
