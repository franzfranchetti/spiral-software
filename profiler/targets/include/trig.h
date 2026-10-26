
#ifndef TRIG_H_INCLUDED
#define TRIG_64_H_INCLUDED
#include <math.h>
#include <stdlib.h>

#define PI    3.14159265358979323846

/* More accurate cosine */
static double cospi(double a) { /*return cos(PI*a); */
    double sign = +1.0;
    while (a >= 2.0) a -= 2.0;
    while (a < 0.0) a += 2.0;
    if (a >= 1.0) { a -= 1.0; sign = -1.0; }

    if (a == 0.0) return sign;
    else if (a == 0.25) return sign * 0.70710678118654752440;
    else if (a == 0.5) return 0.0;

    if (a > 0.5) { a -= 0.5; return -sign * sin(PI * a); }
    else { return sign * cos(PI * a); }
}

static double sinpi(double a) {/*return sin(PI*a); */
    return cospi(0.5 - a);
}

static _Dcomplex omega(int N, int k) {
    _Dcomplex w;

    w = _DCOMPLEX_(cos(2 * PI * k / N), sin(2 * PI * k / N));
    return w;
}

#endif

