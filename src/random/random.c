#include <stdio.h>
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long *idum) {
    long k;
    float ans;

    *idum ^= MASK;
    k = (*idum)/IQ;
    *idum = IA * (*idum - k*IQ) - IR*k;
    if (*idum < 0) {
        * idum += IM;
    }
    ans = AM*(*idum);
    *idum ^= MASK;
    return ans;
}

float boxmuller(int seed)
{
    long dum = 0;
    float x1, x2, y1, r, fac;
    static float y2;
    static int flag = 0;

    // Init seed
    if (seed < 0) {
        dum = -1 * seed;
    }
    if (flag) {
        flag = 0;
        return y2;
    }
    do {
        x1 = ran0(&dum)*2 - 1;
        x2 = ran0(&dum)*2 - 1;
        r = x1*x1 + x2*x2;
    } while (r >= 1);
    fac = sqrt((-2*log(r)/r));
    y1 = x1*fac;
    y2 = x2*fac;
    flag = 1;
    return y1;
}
