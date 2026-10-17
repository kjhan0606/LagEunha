#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
/*
#include "eunha.h"
*/



float midinf(float (*)(SimParameters *, float), SimParameters *, float , float , int );
float midsql(float (*)(SimParameters *, float), SimParameters *, float , float , int );
float midsqu(float (*)(SimParameters *, float), SimParameters *, float , float , int );
float midexp(float (*)(SimParameters *, float), SimParameters *, float , float , int );

/* Quintessence model is based on Sefusatti & Vernizzi 2011, arXiv:1101.1026v3 */
/* Note: read just below equation(29) . The Growthfactor for quintessence model is give in 
 * integral form (Eq. 29)*/

float simparqsimp(float (*func)(SimParameters *, float), float, float, SimParameters *);
float simparqromo(float (*)(SimParameters *, float), float , float ,
		    float (*)(float(*)(SimParameters *, float), SimParameters *, float, float, int), SimParameters *);

float wlambda(SimParameters *, float );
float DEfunc(SimParameters *, float );
float gzDE(SimParameters *, float );
float HofEz(SimParameters *, float );
float ap(SimParameters *, float );
float apold(SimParameters *, float );
float Omega_matter(SimParameters *, float );
float Omega_lambda(SimParameters *, float );
float Cw(SimParameters *, float );
float CwOveraH3(SimParameters *, float );
float Dall(SimParameters *, float );
float Dmslope(SimParameters *, float );
float growthall(SimParameters *, float );
float coscon(SimParameters *, float );
float coscon2(SimParameters *, float );
float growth(SimParameters *, float );
float growth2(SimParameters *, float );
float growthall(SimParameters *, float );
float growthgen(SimParameters *, float ); 
float XgrowthM(SimParameters *, float ); 
float DplusAllQuint0(SimParameters *, float ); 
