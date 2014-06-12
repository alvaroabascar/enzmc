#include "random.h"

float boxmuller(int seed)
{
    float x1, x2, y1, r, fac;
    static float y2;
    static int flag = 0;
    if (seed < 0) {
        srand(-1*seed);
    }
    if (flag) {
        flag = 0;
        return y2;
    }
    do {
        x1 = (float)rand()/((float)RAND_MAX/2) - 1;
        x2 = (float)rand()/((float)RAND_MAX/2) - 1;
        r = x1*x1 + x2*x2;
    } while (r >= 1);
    fac = sqrt((-2*log(r)/r));
    y1 = x1*fac;
    y2 = x2*fac;
    flag = 1;
    return y1;
}
