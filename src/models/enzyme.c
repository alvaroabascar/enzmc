#include "enzyme.h"
#include <stdio.h>
#include <math.h>

/* Michaelis-Menten model:
 * Parameters:
 *              S[1]      -> array of one element (substrate concentration)
 *              p[2] -> array of two elements: { Vmax, Km }
 * Output:
 *              The initial reaction rate for the given [S], Vmax and Km.
 */
double michaelis(double S[1], double p[2])
{
    return p[0] * S[0] / (p[1] + S[0]);
}

/* Michaelis-Menten derivative:
 * Parameters:
 *              S[1]      -> array of one element (substrate concentration)
 *              p[2] -> array of two elements: { Vmax, Km }
 *              k         -> 0 = derivative with respect to Vmax
 *                           1 = derivative with respect to Km
 * Output:
 *              if k = 0 : dv/dVmax
 *              if k = 1 : dv/dKm
 */
double michaelis_dyda(double S[1], double p[2], int k)
{
    double tmp;
    /* p[0] = vmax */
    if (k == 0) {
        return S[0] / (p[1] + S[0]);
    /* p[1] = Km */
    } else if (k == 1) {
        return - p[0] * S[0] / ((p[1] + S[0]) * (p[1] + S[0]));
    }
}

/********************* Multi Substrate ************************/

/* Alberty equation (multi-substrate kinetic)
 *
 * E = Enzyme
 * B, AX = Substrates
 *
 *    +B       +AX                   -A        -BX    --->
 * E <---> EB <---> EAXB <---> EABX <---> EBX <---> E
 *    -B       -AX                   +A        +BX    <---
 *
 * Equation:
 * v0 = Vmax[AX][B] / (KmAX[B] + KmB[AX] + [AX][B] + KsAX*KmB)
 *
 * Parameters:
 *              S[2]      -> Array of independent variables {[AX], [B]}
 *              p[3] -> Array of parameters {Vmax, KmAX, KmB}
 * Output:
 *              Initial reaction rate for the given [AX],[B],Vmax,KmAX,KmB
 */
double alberty(double S[2], double p[4])
{
    /* S[0] = [AX], S[1] = [B]
     * p[0] = Vmax, p[1] = KmAX, p[2] = KmB, p[3] = KsA */
    return p[0] * S[0] * S[1] / \
           (p[1]*S[1] + p[2]*S[0] + S[0]*S[1] + p[3]*p[2]);
}

/* Double displacement equation
 * (Analogous to Alberty, without the last term in the denominator
 * 
 * Equation:
 * v0 = Vmax[AX][B] / (KmAX[B] + KmB[AX] + [AX][B])
 *
 * Parameters:
 *              S[2]      -> Array of independent variables {[AX], [B]}
 *              p[3] -> Array of parameters {Vmax, KmAX, KmB}
 *
 * Output:
 *              Initial reaction rate for the given [AX],[B],Vmax,KmAX,KmB
 */
double pingpong(double S[2], double p[3])
{
    /* S[0] = [AX], S[1] = [B]
     * p[0] = Vmax, p[1] = KmAX, p[2] = KmB */
    return p[0] * S[0] * S[1] /
       (p[1]*S[1] + p[2]*S[0] + S[0]*S[1]);
}

/*********************** Inhibition ***********************/

/* Mixed inhibition
 *
 * Model:
 *          Ks         k2
 * E + S <------> ES -----> E + P
 *
 *          KIa
 * E + I <------> EI
 *
 *          Ks'           KIb
 * EI + S <------> ESI <------> ES + I
 *
 * Equation:
 *
 * v0 = Vmax[S] / (Km*(1+[I]/KIa) + [S]*(1+[I]/KIb))
 *
 * Parameters:
 *
 *              X[2] -> {[S], [I]}
 *              p[4] -> {Vmax, Km, KIa, KIb}
 *
 * Output:
 *              The initial reaction rate
 */
double mixed(double X[2], double p[4])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIa, p[3] = KIb
     */
    return p[0]*X[0] /
           (p[1]*(1+X[1]/p[2]) + X[0]*(1+X[1]/p[3]));
}

/* Competitive inhibition
 *
 * Analogous to mixed, with KIb = infinite
 *
 * Equation:
 *
 * v0 = Vmax[S] / (Km*(1+[I]/KIa) + [S])
 *
 * Parameters:
 *              X[2] -> {[S], [I]}
 *              p[4] -> {Vmax, Km, KIa}
 *
 * Output:
 *              The initial reaction rate
 */
double competitive(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIa
     */
    return p[0]*X[0] /
           (p[1]*(1+X[1]/p[2]) + X[0]);
}

/* Uncompetitive inhibition
 *
 * Analogous to mixed, with KIa = infinite
 *
 * Equation:
 *
 * v0 = Vmax[S] / (Km + [S]*(1+[I]/KIb))
 *
 * Parameters:
 *              X[2] -> {[S], [I]}
 *              p[4] -> {Vmax, Km, KIb}
 *
 * Output:
 *              The initial reaction rate
 */
double uncompetitive(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIb
     */
    return p[0]*X[0] /
           (p[1] + X[0]*(1+X[1]/p[2]));
}

/* Non-competitive inhibition
 *
 * Analogous to mixed, with KIa = KIb
 *
 * Equation:
 *
 * v0 = Vmax[S] / ((Km + [S])*(1+[I]/KIb)))
 *
 * Parameters:
 *              X[2] -> {[S], [I]}
 *              p[4] -> {Vmax, Km, KIb}
 *
 * Output:
 *              The initial reaction rate
 */
double noncompetitive(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = [I]
     * p[0] = Vmax, p[1] = Km, p[2] = KIb
     */
    return p[0]*X[0] /
           ((p[1] + X[0])*(1+X[1]/p[2]));
}

/* Effect of pH
 *
 *
 * Model:
 *         Ks       k2
 * E + S <----> ES ----> E + P
 *
 *           Ka1                      Ka2                     Ks'
 * E + H+  <-----> EH  ,  ES + H+  <-----> HES ,    HE + S <-----> HES
 *
 *           Ka3                     Ka4                      K's
 * E + OH- <-----> EOH ,  ES + OH- <-----> ESOH ,   EOH + S <-----> ESOH
 *
 * Equation:
 *            v0 = Vm*[S] /
 *                 (Km*(1+[H+]/Ka1 + Ka3/[H+]) + [S]*(1 + [H+]/Ka2 + Ka4/[H+]))
 *
 * Parameters:
 *              X[2] -> Array of concentrations {[S], [H+]}
 *              p[6] -> Array of parameters {Vmax, Km, Ka1, Ka2, Ka3, Ka4}
 * Output:
 *              Initial reaction rate.
 */
double ph(double X[2], double p[6])
{
    /* X[0] = [S], X[1] = [H+]
     * p[0] = Vmax, p[1] = Km
     * p[2] = Ka1, p[3] = Ka2, p[4] = Ka3, p[5] = Ka4
     */
     return p[0]*X[0] /
            (p[1]*(1+X[1]/p[2] + p[4]/X[1])+X[0]*(1+X[1]/p[3]+p[5]/X[1]));
}

/* Effect of the temperature on the rate of a michaelis-menten reaction
 *
 *                    kcat
 * E + S  <-----> ES -----> E + P
 *
 * Equation:
 *           v0 = [S]*vmax(T1)*exp(-Ea/R*(1/T1 - 1/T2))/ (Km + [S])
 * 
 *
 * Parameters:
 *           X[2] -> Array of {[S], T2}
 *           p[2] -> Array of parameters {AE, Ea, Km, T1}
 *
 * Note: T1 is the temperature at which Vmax is refferred, and must be
 *       fixed
 */
double michaelistemp(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = T2
     * p[0] = Vmax, p[1] = Km, p[2] = Ea, p[3] = T1 
     */
     return X[0]*p[0]*exp((-p[2]/8.3144621)*(1/p[3] - 1/X[1])) /
                                                            (p[2] + X[0]);
}

/* Michaelis-Menten kinetics accounting for thermal inactivation of the enzyme
 *
 * Equation:
 *           v0 = S*vmax*exp(-kt*t) / (Km + S)
 *
 * Parameters:
 *           X[2] -> {[S], [t]} (substrate concentrations and times)
 *           p[3] -> {vmax, Km, kt} (kt = first order inactivation rate constant
 */
double michaelisinactiv(double X[2], double p[3])
{
    /* X[0] = [S], X[1] = t
     * p[0] = Vmax, p[1] = Km, p[2] = kt
     */
    return X[0]*p[0]*exp(-p[2]*X[1])/(p[1]+X[0]);
}
