// RK4 by Open AI
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
    postype x = 0.2, y = 0, vx, vy, ax, ay;
    postype time = 0, Dtime;
    postype r = sqrt(x * x + y * y);
    postype vamp = r * pow(r * r + eps2, -0.75L);
    vx = -y / r * vamp;
    vy = x / r * vamp;

    do {
        r = sqrt(x * x + y * y);
        vamp = sqrt(vx * vx + vy * vy);
        Dtime = 2 * M_PI * r / vamp / 100.0;
        
        // RK4 coefficients
        postype k1x, k1y, k1vx, k1vy;
        postype k2x, k2y, k2vx, k2vy;
        postype k3x, k3y, k3vx, k3vy;
        postype k4x, k4y, k4vx, k4vy;
        
        acceleration(x, y, eps2, &ax, &ay);
        k1x = vx * Dtime;
        k1y = vy * Dtime;
        k1vx = ax * Dtime;
        k1vy = ay * Dtime;
        
        acceleration(x + k1x / 2, y + k1y / 2, eps2, &ax, &ay);
        k2x = (vx + k1vx / 2) * Dtime;
        k2y = (vy + k1vy / 2) * Dtime;
        k2vx = ax * Dtime;
        k2vy = ay * Dtime;
        
        acceleration(x + k2x / 2, y + k2y / 2, eps2, &ax, &ay);
        k3x = (vx + k2vx / 2) * Dtime;
        k3y = (vy + k2vy / 2) * Dtime;
        k3vx = ax * Dtime;
        k3vy = ay * Dtime;
        
        acceleration(x + k3x, y + k3y, eps2, &ax, &ay);
        k4x = (vx + k3vx) * Dtime;
        k4y = (vy + k3vy) * Dtime;
        k4vx = ax * Dtime;
        k4vy = ay * Dtime;
        
        // Update position and velocity using RK4
        x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        vx += (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6;
        vy += (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6;
        
        time += Dtime;
        printf("%g %g\n", x, y);
    } while (time < 120);

    return 0;
}
