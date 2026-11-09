#ifndef GFS_MATERIAL_DAMPING_BOUND_H
#define GFS_MATERIAL_DAMPING_BOUND_H
#include <math.h>
/* Dimension-independent block Gershgorin accumulation for one rank-one
 * damping row. ids must be UNIQUE after gathering duplicate face entries.
 * vector is a dim-strided per-cell array; mass and row are per-cell arrays.
 * Overwrites vector[dim*i] with its mass-normalized norm as workspace.
 * Caller owns allocation, validation of cell indices, and parallel reduction.
 */
static inline int gfs_material_damping_row(int dim,int used,const int*ids,
 double*vector,const double*mass,double coefficient,double*row){
 if(dim<1||dim>3||used<0||!ids||!vector||!mass||!row||
    !(coefficient>=0)||!isfinite(coefficient))return 0;
 double sum=0;
 for(int j=0;j<used;j++){
  int i=ids[j];if(!(mass[i]>0)||!isfinite(mass[i]))return 0;
  /* hypot(0,x) is |x|; retain hypot for the remaining components
   * to preserve overflow/underflow protection in every dimension. */
  double length=fabs(vector[dim*i]);for(int a=1;a<dim;a++)length=hypot(length,vector[dim*i+a]);
  double norm=length/sqrt(mass[i]);if(!isfinite(norm))return 0;
  vector[dim*i]=norm;sum+=norm;
 }
 if(!isfinite(sum))return 0;
 for(int j=0;j<used;j++){int i=ids[j];row[i]+=coefficient*vector[dim*i]*sum;if(!isfinite(row[i]))return 0;}
 return 1;
}
#endif
