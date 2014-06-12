#include "math.h"

double mean(int n, double data[n])
{
    int i;
    double sum = 0;
    for (i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / n;
}

double sum(int n, double data[n])
{
    double sum;
    for (sum = 0; n > 0; n--) {
        sum += data[n];
    }
    return sum;
}

double std(int n, double data[n])
{
    int i;
    double std, tmp;
    double meanv = mean(n, data);
    for (i = std = 0; i < n; i++) {
        tmp = data[i] - meanv;
        std += tmp * tmp;
    }
    return sqrt(std/(n-1));
}

double max(int n, double data[n])
{
    double max;
    for (max = 0; n > 0; n--) {
        if (data[n-1] > max) {
            max = data[n-1];
        }
    }
    return max;
}

/* diff:
 * calculates the derivate of a function with respect to its kth parameter, at
 * the point defined by the array of independent variables xi[]
 */
double dfda(double xi[],
            double params[],
            double fun(double x[], double params[]), /* func to derivate */
            int k) /* parameter with respect to which derivate */
{
    double diff;
    double backupk = params[k]; /* backup kth parameter */
    double h = 1e-4;
    /* calculate function at current point */
    params[k] += h; /* function at ak = ak+h */
    diff = fun(xi, params);
    params[k] = backupk - h; /* function at ak = ak-h */
    diff = (diff - fun(xi, params))/(2*h);
    params[k] = backupk;
    return diff;
}


