#include <stdio.h>
#include "../models/michaelis.c"
#include "math.c"
#include "matrix.c"
#include <math.h>

double fun(double x[], double a[])
{
    return exp(a[0]);
}

double dfun(double x[], double a[], int k)
{
    return exp(a[0]);
}

#define N 100
int main()
{
    int i;
    double a[2] = {5, 10}; /* Vm, Km */
    double S[] =  {10};
    double ana, num, mic;
    FILE *fana, *fnum, *fmic;
    fana = fopen("ana", "w");
    fnum = fopen("num", "w");
    fmic = fopen("mic", "w");
    for (i = 0; i < N; i++) {
        a[0] = (double) i;
        ana = michaelis_dyda(S, a, 0);
        num = dfda(S, a, michaelis, 0);
        mic = michaelis(S, a);
        fprintf(fana, "%.8e\n", ana);
        fprintf(fnum, "%.8e\n", num);
        fprintf(fmic, "%.8e\n", mic);
    }
}
