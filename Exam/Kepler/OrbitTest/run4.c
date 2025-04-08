// Leap-frog 2nd sympletic form
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double postype;

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

    acceleration(x, y, eps2, &ax, &ay);
    vx += 0.5 * ax * Dtime;
    vy += 0.5 * ay * Dtime;

    do {
        r = sqrt(x * x + y * y);
        vamp = sqrt(vx * vx + vy * vy);
        Dtime = 2 * M_PI * r / vamp / 100.0;
        
        // Leapfrog integration
        x += vx * Dtime;
        y += vy * Dtime;
        
        acceleration(x, y, eps2, &ax, &ay);
        vx += ax * Dtime;
        vy += ay * Dtime;
        
        time += Dtime;
        printf("%g %g\n", x, y);
    } while (time < 120);

    return 0;
}
