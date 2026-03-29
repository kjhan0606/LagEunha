#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double postype;

// Yoshida 4th order symplectic integrator coefficients
const postype w1 = 0.675603595979828817023843904485730;
const postype w2 = -0.175603595979828817023843904485730;
const postype w3 = 0.675603595979828817023843904485730 / 2.0;


void acceleration(postype x, postype y, postype eps2, postype *ax, postype *ay) {
    postype r2 = x * x + y * y;
    postype r3 = pow(r2 + eps2, -1.5L);
    *ax = -x * r3;
    *ay = -y * r3;
}

int main(int argc, char **argv) {
    postype eps2 = 0.1 * 0.1;
    postype x = 1.5, y = 0, vx, vy, ax, ay;
    postype time = 0, Dtime;
    postype r = sqrt(x * x + y * y);
    postype vamp = r * pow(r * r + eps2, -0.75L);
    vx = -y / r * vamp;
    vy = x / r * vamp;

    do {
        r = sqrt(x * x + y * y);
        vamp = sqrt(vx * vx + vy * vy);
        Dtime = 2 * M_PI * r / vamp / 100.0;
        
        // Yoshida 4th order symplectic integration steps
		x += w3 * vx * Dtime;
        y += w3 * vy * Dtime;
//      i = 1
        acceleration(x, y, eps2, &ax, &ay);
        vx += w1 * ax * Dtime;
        vy += w1 * ay * Dtime;
        x += w2 * vx * Dtime;
        y += w2 * vy * Dtime;


//      i = 2
        acceleration(x, y, eps2, &ax, &ay);
        vx += w2 * ax * Dtime;
        vy += w2 * ay * Dtime;

        x += w1 * vx * Dtime;
        y += w1 * vy * Dtime;
//      i = 3
//
        acceleration(x, y, eps2, &ax, &ay);
        vx += w3 * ax * Dtime;
        vy += w3 * ay * Dtime;

        x += w3 * vx * Dtime;
        y += w3 * vy * Dtime;

        time += Dtime;
        printf("%g %g\n", x, y);
    } while (time < 120);

    return 0;
}
