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
    postype r3 = pow(r2 + eps2, -1.5);
    *ax = -x * r3;
    *ay = -y * r3;
}

void integrate(postype *x, postype *y, postype *vx, postype *vy, postype eps2, postype dt) {
    postype ax, ay;
    
    // First step
    *x += w3 * (*vx) * dt;
    *y += w3 * (*vy) * dt;
    acceleration(*x, *y, eps2, &ax, &ay);
    *vx += w3 * ax * dt;
    *vy += w3 * ay * dt;
    
    // Second step
    *x += w2 * (*vx) * dt;
    *y += w2 * (*vy) * dt;
    acceleration(*x, *y, eps2, &ax, &ay);
    *vx += w2 * ax * dt;
    *vy += w2 * ay * dt;
    
    // Third step
    *x += w1 * (*vx) * dt;
    *y += w1 * (*vy) * dt;
    acceleration(*x, *y, eps2, &ax, &ay);
    *vx += w1 * ax * dt;
    *vy += w1 * ay * dt;
    
    // Fourth step
    *x += w3 * (*vx) * dt;
    *y += w3 * (*vy) * dt;
    acceleration(*x, *y, eps2, &ax, &ay);
    *vx += w3 * ax * dt;
    *vy += w3 * ay * dt;
}

int main() {
    postype eps2 = 0.01;
    postype x = 1.5, y = 0.0;
    postype r = sqrt(x * x + y * y);
    postype vamp = r * pow(r * r + eps2, -0.75L);
    postype vx = -y / r * vamp;
    postype vy = x / r * vamp;
    postype time = 0.0, Dtime;
    
    while (time < 120.0) {
        r = sqrt(x * x + y * y);
        vamp = sqrt(vx * vx + vy * vy);
        Dtime = 2 * M_PI * r / vamp / 100.0;
        integrate(&x, &y, &vx, &vy, eps2, Dtime);
        time += Dtime;
        printf("%g %g\n", x, y);
    }
    
    return 0;
}
