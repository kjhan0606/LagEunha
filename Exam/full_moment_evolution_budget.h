#ifndef EXAM_FULL_MOMENT_EVOLUTION_BUDGET_H
#define EXAM_FULL_MOMENT_EVOLUTION_BUDGET_H
#include <math.h>
#include <float.h>
/* final_step is the global endpoint index. Limit new work per batch, rather
 * than forbidding any continuation beyond global step 200 (t=1 for h=.005).
 * Caller output is unchanged on failure; no timestep adaptation is permitted. */
static inline int full_moment_evolution_budget(double time,double h,int final_step,int *start_step){
 if(!start_step||!isfinite(time)||time<0||!isfinite(h)||h<=0||h>.005||final_step<1||final_step>1000000)return 1;
 double index=time/h;
 if(!isfinite(index)||index<0||index>final_step)return 1;
 int first=(int)llround(index);
 if(first<0||first>=final_step||final_step-first>200
  ||fabs(time-first*h)>32*DBL_EPSILON*fmax(1.,time))return 1;
 *start_step=first;return 0;
}
#endif
