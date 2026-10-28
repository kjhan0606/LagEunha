#include "reversible_memory.h"
int rm_probe(int n,const double *dx,const double *dy,int owner,double coverage,
             double target,double *row,double *work,double *meta) {
    return rm_quadratic_row(n,dx,dy,owner,coverage,target,row,work,meta);
}
