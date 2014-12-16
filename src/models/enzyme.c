#include "enzyme.h"
#include <stdio.h>
#include <math.h>

/* For an explanation of the models, please refer to the header file */

double michaelis(double S[1], double p[2])
{
    /* S[0] = [S]
     * p[0] = Vmax, p[1] = Km */
    return p[0] * S[0] / (p[1] + S[0]);
}

double alberty(double S[2], double p[3])
{
    /* S[0] = [AX], S[1] = [B]
     * p[0] = Vmax, p[1] = KmAX, p[2] = KmB, p[3] = KsA */
    return p[0] * S[0] * S[1] / \
           (p[1]*S[1] + p[2]*S[0] + S[0]*S[1] + p[3]*p[2]);
}

double pingpong(double S[2], double p[3])
{
    /* S[0] = [AX], S[1] = [B]
     * p[0] = Vmax, p[1] = KmAX, p[2] = KmB */
    return p[0] * S[0] * S[1] /
       (p[1]*S[1] + p[2]*S[0] + S[0]*S[1]);
}

double mixed(double X[2], double p[4])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIa, p[3] = KIb
     */
    return p[0]*X[0] /
           (p[1]*(1+X[1]/p[2]) + X[0]*(1+X[1]/p[3]));
}

double competitive(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIa
     */
    return p[0]*X[0] /
           (p[1]*(1+X[1]/p[2]) + X[0]);
}

double uncompetitive(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIb
     */
    return p[0]*X[0] /
           (p[1] + X[0]*(1+X[1]/p[2]));
}

double noncompetitive(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIb
     */
    return p[0]*X[0] /
           ((p[1] + X[0])*(1+X[1]/p[2]));
}

double ph(double X[2], double p[6])
{
    /* X[0] = [S], X[1] = [H+]
     * p[0] = Vmax, p[1] = Km
     * p[2] = Ka1, p[3] = Ka2, p[4] = Ka3, p[5] = Ka4
     */
     return p[0]*X[0] /
            (p[1]*(1+X[1]/p[2] + p[4]/X[1])+X[0]*(1+X[1]/p[3]+p[5]/X[1]));
}

double michaelistemp(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = T2
     * p[0] = Vmax, p[1] = Km, p[2] = Ea, p[3] = T1 
     */
     return X[0]*p[0]*exp((-p[2]/8.3144621)*(1/p[3] - 1/X[1])) /
                                                            (p[2] + X[0]);
}

double michaelis_inactiv(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = t
     * p[0] = Vmax, p[1] = Km, p[2] = kt
     */
    return X[0]*p[0]*exp(-p[2]*X[1])/(p[1]+X[0]);
}
