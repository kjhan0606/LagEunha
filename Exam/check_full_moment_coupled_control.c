/* Native fixtures: no geometry, fluid provider, GPU, MPI or accepted step. */
#include "full_moment_coupled_control.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
typedef struct {
 int calls,mode,finalized;
 double elapsed,rows[10][24],trials[10][24],mass,other_mass,factors[4],seed;
 int64_t id;long double Psi;
} Context;
static int count(void *p){return ((Context*)p)->calls;}
static double elapsed(void *p){return ((Context*)p)->elapsed;}
static int evaluate(void *p,const FullMomentState *base,const FullMomentState *trial,int iteration,
 FullMomentState *mapped,FullMomentCoupledResidual *residual,long double *solver,size_t solver_count){
 Context *c=p;c->calls+=c->mode==3?2:1;
 if(c->mode==2)return 17;
 assert(solver_count==21);memset(solver,0,solver_count*sizeof(*solver));
 memcpy(c->trials[iteration],trial->state,24*sizeof(double));
 memcpy(c->rows[iteration],trial->state,24*sizeof(double));
 if(c->mode==9)c->rows[iteration][0]=iteration==0?1.:iteration==1?.9:trial->state[0];
 else if(c->mode==12&&iteration==4)c->rows[iteration][0]=trial->state[0];
 else c->rows[iteration][0]=c->mode==1||c->mode==6?trial->state[0]:.01*trial->state[0]+.99;
 if(c->mode==10)c->rows[iteration][5]=iteration==0?1.:iteration==1?2.:trial->state[5];
 *mapped=*base;mapped->state=c->rows[iteration];mapped->header.time=2.;mapped->header.epoch=10+iteration;
 if(c->mode==7)mapped->mass=&c->other_mass;
 memset(residual,0,sizeof *residual);
 residual->maximum[FM_RES_X]=fabs(mapped->state[0]-trial->state[0]);
 solver[0]=(long double)mapped->state[0]-trial->state[0];
 if(c->mode==9){
  if(iteration==0)solver[0]=1.;
  if(iteration==1)solver[0]=-.1L;
  if(iteration==2)solver[0]=0.;
  residual->maximum[FM_RES_X]=fabsl(solver[0]);
 }
 if(c->mode==10){
  solver[0]=iteration==0?1.L:iteration==1?8.L/7.:0.;
  residual->maximum[FM_RES_X]=iteration==2?0.:1.;
 }
 if(c->mode==11){
  residual->metric_epoch=iteration?1:0;
  residual->metric_reset_cells=iteration==1?1:0;
 }
 if(c->mode==12){
  const long double sequence[5]={1.,.5,.6,.7,0.};
  solver[0]=sequence[iteration<5?iteration:4];
  residual->maximum[FM_RES_X]=fabsl(solver[0]);
 }
 for(size_t k=0;k<solver_count;k++)residual->vector_norm2+=solver[k]*solver[k];
 if(c->mode==1)residual->maximum[FM_RES_THERMAL]=1e-3;
 if(c->mode==6)c->elapsed=100.;
 return 0;
}
static int finalize(void *p,const FullMomentState *base,const FullMomentState *trial,
 const FullMomentState *mapped,const FullMomentCoupledResidual *residual,int iteration){
 Context *c=p;(void)base;(void)iteration;c->finalized++;
 assert(full_moment_coupled_converged(residual));
 assert(mapped->state!=trial->state);
 if(c->mode==9){assert(fabs(c->trials[2][0]-10./11.)<1e-12);return 0;}
 if(c->mode==10){assert(fabs(c->trials[2][5]-1.)<1e-12);return 0;}
 if(c->mode==12){assert(fabs(mapped->state[0]-trial->state[0])<1e-14);return 0;}
 assert(fabs(mapped->state[0]-(.01*trial->state[0]+.99))<1e-12);
 if(c->mode==4)return 9;
 if(c->mode==5)c->calls++;
 return 0;
}
static FullMomentState make_base(Context *c,double *state){
 FullMomentState s={0};s.header.cells=1;s.header.state_elements=24;s.header.time=1.;s.header.epoch=1;s.header.run_token=42;
 c->mass=1.;c->other_mass=2.;c->factors[0]=c->factors[3]=1.;
 s.ids=&c->id;s.mass=&c->mass;s.state=state;s.factors=c->factors;s.Psi=&c->Psi;s.weight_rate_seed=&c->seed;return s;
}
static void make_physical_raw(double state[24]){
 memset(state,0,24*sizeof(*state));
 state[5]=state[8]=1.;state[15]=state[18]=1.;state[19]=2.;
}
static void check_type2_fit(void){
 double gamma[3];int depth,rank;long double predicted,condition;
 long double current[2]={1.,0.},h0[2]={0.,1.},h1[2]={-1.,-1.};
 const long double *history[3]={h0,h1,NULL};
 assert(full_moment_coupled_type2_fit(2,2,current,history,gamma,&depth,&rank,&predicted,&condition));
 assert(depth==2&&rank==2&&fabs(gamma[0]-2./3.)<1e-12&&fabs(gamma[1]-1./3.)<1e-12);
 assert(predicted<1e-24L&&isfinite(condition));

 long double dependent_h1[2]={-1.,2.};history[1]=dependent_h1;
 assert(full_moment_coupled_type2_fit(2,2,current,history,gamma,&depth,&rank,&predicted,&condition));
 assert(depth==1&&rank==1&&fabs(gamma[0]-.5)<1e-12);

 long double large_current[2]={1.e6L,0.},near_history[2]={1.e6L-1.e-3L,0.};
 history[0]=near_history;history[1]=NULL;
 assert(!full_moment_coupled_type2_fit(2,1,large_current,history,gamma,&depth,&rank,&predicted,&condition));
 long double invalid[2]={NAN,0.};history[0]=invalid;
 assert(!full_moment_coupled_type2_fit(2,1,current,history,gamma,&depth,&rank,&predicted,&condition));
}
int main(void){
 check_type2_fit();
 for(int mode=0;mode<9;mode++){
  Context c={0};c.mode=mode;double raw[24];make_physical_raw(raw);
  FullMomentState base=make_base(&c,raw),out={0};out.header.magic=93;
  FullMomentState expected=out;FullMomentCoupledControlReport report={0};
  FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=evaluate;callbacks.finalize=finalize;
  callbacks.provider_count=count;callbacks.elapsed_seconds=elapsed;
  callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=150;callbacks.wall_limit_seconds=10.;
  if(mode==8)c.elapsed=100.;
  int status=full_moment_coupled_control(&base,&callbacks,&out,&report);
  if(mode==0){
   assert(status==FM_COUPLED_OK&&report.converged&&report.final_checks_passed);
   assert(c.calls==3&&c.finalized==1&&report.iterations==3);
   assert(fabs(out.state[0]-1.)<1e-10);full_moment_state_close(&out);
  }else{
   assert(status!=FM_COUPLED_OK&&memcmp(&out,&expected,sizeof out)==0);
   if(mode==1)assert(status==FM_COUPLED_NOT_CONVERGED&&c.calls==10&&c.finalized==0);
   if(mode==2)assert(status==FM_COUPLED_EVALUATION_FAILED&&c.calls==1);
   if(mode==3)assert(status==FM_COUPLED_COUNTER_VIOLATION&&c.calls==2);
   if(mode==4)assert(status==FM_COUPLED_FINAL_CHECK_FAILED&&c.finalized==1);
   if(mode==5)assert(status==FM_COUPLED_COUNTER_VIOLATION&&c.finalized==1);
   if(mode==6)assert(status==FM_COUPLED_WALL_BUDGET&&c.finalized==0);
   if(mode==7)assert(status==FM_COUPLED_EVALUATION_FAILED&&report.callback_status==-1);
   if(mode==8)assert(status==FM_COUPLED_WALL_BUDGET&&c.calls==0);
  }
  assert(raw[0]==0.&&c.mass==1.);
 }
 {
  Context c={0};c.mode=9;double raw[24];make_physical_raw(raw);
  FullMomentState base=make_base(&c,raw),out={0};FullMomentCoupledControlReport report={0};
  FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=evaluate;callbacks.finalize=finalize;
  callbacks.provider_count=count;callbacks.elapsed_seconds=elapsed;
  callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=150;callbacks.wall_limit_seconds=10.;
  assert(full_moment_coupled_control(&base,&callbacks,&out,&report)==FM_COUPLED_OK);
  assert(c.calls==3&&report.iterations==3&&c.finalized==1);
  assert(report.anderson_depth[1]==1&&report.anderson_qr_rank[1]==1);
  assert(fabs(report.anderson_gamma[1][0]-1./11.)<1e-12&&report.anderson_gamma[1][1]==0.);
  assert(report.anderson_lambda[1]==1.&&fabs(c.trials[2][0]-10./11.)<1e-12);
  assert(out.state[20]==raw[20]&&out.state[21]==raw[21]&&out.state[22]==raw[22]&&out.state[23]==raw[23]);
  assert(out.weight_rate_seed[0]==c.seed&&out.mass[0]==c.mass);
  full_moment_state_close(&out);
 }
 {
  Context c={0};c.mode=10;double raw[24];make_physical_raw(raw);
  FullMomentState base=make_base(&c,raw),out={0};FullMomentCoupledControlReport report={0};
  FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=evaluate;callbacks.finalize=finalize;
  callbacks.provider_count=count;callbacks.elapsed_seconds=elapsed;
  callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=150;callbacks.wall_limit_seconds=10.;
  assert(full_moment_coupled_control(&base,&callbacks,&out,&report)==FM_COUPLED_OK);
  assert(report.anderson_depth[1]==1&&fabs(report.anderson_gamma[1][0]-8.)<1e-12);
  assert(report.anderson_lambda[1]==.125&&report.anderson_fallback[1]==FM_AA_FALLBACK_POSITIVITY_BACKTRACK);
  full_moment_state_close(&out);
 }
 {
  Context c={0};c.mode=11;double raw[24];make_physical_raw(raw);
  FullMomentState base=make_base(&c,raw),out={0};FullMomentCoupledControlReport report={0};
  FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=evaluate;callbacks.finalize=finalize;
  callbacks.provider_count=count;callbacks.elapsed_seconds=elapsed;
  callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=150;callbacks.wall_limit_seconds=10.;
  assert(full_moment_coupled_control(&base,&callbacks,&out,&report)==FM_COUPLED_OK);
  assert(report.anderson_depth[1]==0&&report.anderson_fallback[1]==FM_AA_FALLBACK_METRIC_RESET);
  full_moment_state_close(&out);
 }
 {
  Context c={0};c.mode=12;double raw[24];make_physical_raw(raw);
  FullMomentState base=make_base(&c,raw),out={0};FullMomentCoupledControlReport report={0};
  FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=evaluate;callbacks.finalize=finalize;
  callbacks.provider_count=count;callbacks.elapsed_seconds=elapsed;
  callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=150;callbacks.wall_limit_seconds=10.;
  assert(full_moment_coupled_control(&base,&callbacks,&out,&report)==FM_COUPLED_OK);
  assert(report.prediction_miss_streak[2]==1);
  assert(report.prediction_miss_streak[3]==0&&report.anderson_fallback[3]==FM_AA_FALLBACK_PREDICTION_RESET);
  full_moment_state_close(&out);
 }
 Context c={0};c.mode=1;c.calls=149;double raw[24];make_physical_raw(raw);
 FullMomentState base=make_base(&c,raw),out={0};FullMomentCoupledControlReport report={0};
 FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=evaluate;callbacks.finalize=finalize;
 callbacks.provider_count=count;callbacks.elapsed_seconds=elapsed;
 callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=150;callbacks.wall_limit_seconds=10.;
 assert(full_moment_coupled_control(&base,&callbacks,&out,&report)==FM_COUPLED_PROVIDER_BUDGET);
 assert(c.calls==150&&report.iterations==1&&out.state==NULL);
 printf("{\"status\":\"passed_native_coupled_control_fixtures\",\"cases\":18,\"type2_depth_two_qr\":true,\"rank_deficiency_drops_oldest\":true,\"coefficient_limit_fallback\":true,\"controller_type2_integration\":true,\"positivity_backtracking\":true,\"metric_epoch_reset\":true,\"prediction_miss_history_reset\":true,\"fixed_mass_ids_and_ledger_preserved\":true,\"cholesky_rebuilt\":true,\"step_provider_limit\":10,\"total_provider_limit\":150,\"accepted_time_steps\":0,\"scope\":\"Scalar fake-callback controller fixtures only; MPI adapter and physical gate remain to be validated\"}\n");
 return 0;
}
