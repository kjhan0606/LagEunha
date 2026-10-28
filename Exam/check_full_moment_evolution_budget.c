#include "full_moment_evolution_budget.h"
#include <stdio.h>
int main(void){
 struct {double time,h;int end,first,valid;} cases[]={
  {0,.005,200,0,1},{.05,.005,200,10,1},{1.,.005,400,200,1},
  {5.,.005,1200,1000,1},{9.,.005,2000,1800,1},
  {1.,.005,401,0,0},{1.,.005,200,0,0},{1.,.005,199,0,0},
  {.0011,.001,2,0,0},{NAN,.005,200,0,0},{INFINITY,.005,200,0,0},
  {1.,0,200,0,0},{1.,.006,200,0,0},{1.,DBL_MIN,200,0,0},
  {-1.,.005,200,0,0},{0,.005,1000001,0,0}
 };
 size_t count=sizeof cases/sizeof cases[0];
 for(size_t i=0;i<count;i++){
  int first=123456;int status=full_moment_evolution_budget(cases[i].time,cases[i].h,cases[i].end,&first);
  if(cases[i].valid?(status||first!=cases[i].first):(!status||first!=123456))return (int)i+1;
 }
 puts("{\"status\":\"passed_native_continuation_batch_budget\",\"cases\":16,\"maximum_additional_steps\":200,\"maximum_actual_RHS\":400,\"t1_and_t10_global_indices_allowed\":true,\"physical_providers\":0,\"accepted_time_steps\":0}");
 return 0;
}
