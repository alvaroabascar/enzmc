#ifndef _ENZYME_

double michaelis(double S[1], double p[2]);
/* Michaelis-Menten
 * Parameters:
 *              S[1]      -> array of one element (substrate concentration)
 *              p[2] -> array of two elements: { Vmax, Km }
 * Output:
 *              The initial reaction rate for the given [S], Vmax and Km.
 */

/********************* Multi Substrate models ************************/

double alberty(double S[2], double p[3]);
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

double pingpong(double S[2], double p[3]);
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

/*********************** Inhibition ***********************/

double mixed(double X[2], double p[4]);
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

double competitive(double X[2], double p[3]);
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

double uncompetitive(double X[2], double p[3]);
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

double noncompetitive(double X[2], double p[3]);
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

double ph(double X[2], double p[6]);
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

double michaelistemp(double X[2], double p[3]);
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

double michaelis_inactiv(double X[2], double p[3]);
/* Michaelis-Menten kinetics accounting for thermal inactivation of the enzyme
 *
 * Equation:
 *           v0 = S*vmax*exp(-kt*t) / (Km + S)
 *
 * Parameters:
 *           X[2] -> {[S], [t]} (substrate concentrations and times)
 *           p[3] -> {vmax, Km, kt} (kt = first order inactivation rate constant
 */

#define _ENZYME_
#endif
