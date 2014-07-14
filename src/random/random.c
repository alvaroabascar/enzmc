#include <stdio.h>
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long *seed) {
    long k;
    float ans;

    *seed ^= MASK;
    k = (*seed)/IQ;
    *seed = IA * (*seed - k*IQ) - IR*k;
    if (*seed < 0) {
        * seed += IM;
    }
    ans = AM*(*seed);
    *seed ^= MASK;
    return ans;
}

float boxmuller(long *seed)
{
    float x1, x2, y1, r, fac;
    static float y2;
    static int flag = 0;

    if (flag) {
        flag = 0;
        return y2;
    }
    do {
        x1 = ran0(seed)*2 - 1;
        x2 = ran0(seed)*2 - 1;
        r = x1*x1 + x2*x2;
    } while (r >= 1);
    fac = sqrt((-2*log(r)/r));
    y1 = x1*fac;
    y2 = x2*fac;
    flag = 1;
    return y1;
}
