#include "full_moment_coupled_control.h"
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

static void coupled_heap_state_close(FullMomentState *s){
 if(!s)return;
 free(s->state);free(s->factors);free(s->Psi);free(s->weight_rate_seed);
 memset(s,0,sizeof *s);
}
static void coupled_history_clear(FullMomentState maps[3],long double *residuals[3],int *count){
 for(int i=0;i<3;i++){full_moment_state_close(&maps[i]);free(residuals[i]);residuals[i]=NULL;}
 *count=0;
}
static int coupled_dot(size_t n,const long double *a,const long double *b,long double *out){
 long double value=0.;
 for(size_t i=0;i<n;i++)value+=a[i]*b[i];
 if(!isfinite(value))return 1;
 *out=value;return 0;
}
static int coupled_norm2(size_t n,const long double *a,long double *out){
 long double value=0.;
 for(size_t i=0;i<n;i++)value+=a[i]*a[i];
 if(!isfinite(value)||value<0.)return 1;
 *out=value;return 0;
}
static int coupled_collective_status(const FullMomentCoupledCallbacks *c,int local_status){
 return c->collective_status?c->collective_status(c->context,local_status):local_status;
}
static long double coupled_predicted_norm2(size_t n,const long double *current,
 const long double *const history[3],int depth,const double gamma[3],double lambda){
 long double norm=0.;
 for(size_t k=0;k<n;k++){
  long double value=current[k];
  for(int j=0;j<depth;j++){
   const long double *newer=j?history[j-1]:current;
   value-=(long double)lambda*gamma[j]*(newer[k]-history[j][k]);
  }
  norm+=value*value;
 }
 return norm;
}

int full_moment_coupled_type2_fit(size_t count,int maximum_depth,
 const long double *current,const long double *const history[3],
 double gamma[3],int *used_depth,int *qr_rank,
 long double *predicted_norm2,long double *condition){
 if(gamma)memset(gamma,0,3*sizeof(*gamma));
 if(used_depth)*used_depth=0;
 if(qr_rank)*qr_rank=0;
 if(predicted_norm2)*predicted_norm2=INFINITY;
 if(condition)*condition=INFINITY;
 if(!count||!current||!history||!gamma||!used_depth||!qr_rank||!predicted_norm2||!condition)return 0;
 if(maximum_depth>3)maximum_depth=3;
 if(maximum_depth<1)return 0;
 for(int j=0;j<maximum_depth;j++)if(!history[j])return 0;
 long double current_norm2=0.;if(coupled_norm2(count,current,&current_norm2)||!(current_norm2>0.))return 0;
 const long double rank_tolerance=64.*sqrtl(LDBL_EPSILON);
 for(int depth=maximum_depth;depth>0;){
  if(count>SIZE_MAX/(size_t)depth/sizeof(long double))return 0;
  long double *work=malloc((size_t)depth*count*sizeof(*work));
  if(!work)return 0;
  long double R[3][3]={{0.}};int rank=0,dependent=0;
  for(int j=0;j<depth;j++){
   const long double *newer=j?history[j-1]:current;
   const long double *older=history[j];
   long double *column=work+(size_t)j*count,original_norm2=0.;
   for(size_t k=0;k<count;k++){
    column[k]=newer[k]-older[k];original_norm2+=column[k]*column[k];
   }
   long double original_norm=sqrtl(original_norm2);
   if(!isfinite(original_norm2)||!isfinite(original_norm)||!(original_norm>0.)){
    dependent=1;break;
   }
   for(int pass=0;pass<2;pass++)for(int i=0;i<j;i++){
    long double coefficient=0.;
    if(coupled_dot(count,work+(size_t)i*count,column,&coefficient)){dependent=1;break;}
    R[i][j]+=coefficient;
    for(size_t k=0;k<count;k++)column[k]-=coefficient*work[(size_t)i*count+k];
   }
   if(dependent)break;
   long double residual_norm2=0.;
   for(size_t k=0;k<count;k++)residual_norm2+=column[k]*column[k];
   long double residual_norm=sqrtl(residual_norm2);
   if(!isfinite(residual_norm2)||!isfinite(residual_norm)||
      residual_norm<=rank_tolerance*original_norm){dependent=1;break;}
   R[j][j]=residual_norm;
   for(size_t k=0;k<count;k++)column[k]/=residual_norm;
   rank++;
  }
  if(dependent||rank<depth){
   free(work);depth=rank;continue; /* drop the oldest dependent columns */
  }
  double trial_gamma[3]={0.};int valid=1;
  for(int i=depth-1;i>=0;i--){
   long double rhs=0.;
   if(coupled_dot(count,work+(size_t)i*count,current,&rhs)){valid=0;break;}
   for(int j=i+1;j<depth;j++)rhs-=(long double)R[i][j]*trial_gamma[j];
   if(!(R[i][i]>0.)||!isfinite(rhs)){valid=0;break;}
   trial_gamma[i]=(double)(rhs/R[i][i]);
   if(!isfinite(trial_gamma[i])||fabs(trial_gamma[i])>10.){valid=0;break;}
  }
  long double condition_estimate=INFINITY;
  if(valid){
   long double inverse[3][3]={{0.}},norm_R=0.,norm_inverse=0.;
   for(int i=0;i<depth;i++){
    long double row=0.;for(int j=i;j<depth;j++)row+=fabsl(R[i][j]);
    if(row>norm_R)norm_R=row;
   }
   for(int column=0;column<depth;column++)for(int i=depth-1;i>=0;i--){
    long double rhs=i==column?1.:0.;
    for(int j=i+1;j<depth;j++)rhs-=R[i][j]*inverse[j][column];
    inverse[i][column]=rhs/R[i][i];
   }
   for(int i=0;i<depth;i++){
    long double row=0.;for(int j=0;j<depth;j++)row+=fabsl(inverse[i][j]);
    if(row>norm_inverse)norm_inverse=row;
   }
   condition_estimate=norm_R*norm_inverse;
   if(!isfinite(condition_estimate))valid=0;
  }
  long double prediction=INFINITY;
  if(valid){
   const long double *hist[3]={history[0],history[1],history[2]};
   prediction=coupled_predicted_norm2(count,current,hist,depth,trial_gamma,1.);
   if(!isfinite(prediction)||prediction>current_norm2*(1.+1.e-10L))valid=0;
  }
  if(!valid){free(work);depth--;continue;}
  for(int j=0;j<depth;j++)gamma[j]=trial_gamma[j];
  *used_depth=depth;*qr_rank=rank;*predicted_norm2=prediction;
  *condition=condition_estimate;
  free(work);return isfinite(*condition)?1:0;
 }
 return 0;
}

/* Type-II candidate, damped from the current map output. Only per-cell raw
 * U and Psi are accelerated. IDs/mass stay fixed, the four global ledger
 * entries and weight-rate seed come from the current map, and factors are
 * reconstructed. W is projected to the accepted base's common gauge. */
static int coupled_type2_candidate(const FullMomentState *base,
 const FullMomentState *mapped,const FullMomentState maps[3],int depth,
 const double gamma[3],double lambda,FullMomentState *out){
 if(!base||!mapped||!maps||!gamma||!out||out->state||!isfinite(lambda)||lambda<0.||lambda>1.)return 1;
 size_t n=(size_t)base->header.cells,raw_count=20*n,state_count=raw_count+4;
 FullMomentState s={0};s.header=mapped->header;s.ids=base->ids;s.mass=base->mass;
 s.state=malloc(state_count*sizeof(*s.state));s.factors=malloc(4*n*sizeof(*s.factors));
 s.Psi=malloc(n*sizeof(*s.Psi));s.weight_rate_seed=malloc(n*sizeof(*s.weight_rate_seed));
 if(!s.state||!s.factors||!s.Psi||!s.weight_rate_seed){coupled_heap_state_close(&s);return 2;}
 memcpy(s.state,mapped->state,state_count*sizeof(*s.state));
 memcpy(s.weight_rate_seed,mapped->weight_rate_seed,n*sizeof(*s.weight_rate_seed));
 for(size_t k=0;k<raw_count;k++){
  long double aa=mapped->state[k];
  for(int j=0;j<depth;j++){
   const double *newer=j?maps[j-1].state:mapped->state;
   const double *older=maps[j].state;
   aa-=(long double)gamma[j]*((long double)newer[k]-older[k]);
  }
  long double value=(long double)mapped->state[k]+lambda*(aa-mapped->state[k]);
  s.state[k]=(double)value;
  if(!isfinite(s.state[k])){coupled_heap_state_close(&s);return 3;}
 }
 for(size_t i=0;i<n;i++){
  long double aa=mapped->Psi[i];
  for(int j=0;j<depth;j++){
   const long double *newer=j?maps[j-1].Psi:mapped->Psi;
   const long double *older=maps[j].Psi;
   aa-=(long double)gamma[j]*(newer[i]-older[i]);
  }
  s.Psi[i]=(long double)mapped->Psi[i]+lambda*(aa-mapped->Psi[i]);
  if(!isfinite(s.Psi[i])||!isfinite(s.weight_rate_seed[i])){coupled_heap_state_close(&s);return 3;}
 }
 long double gauge=0.;
 for(size_t i=0;i<n;i++)gauge+=(long double)s.state[20*i+2]-base->state[20*i+2];
 gauge/=n;
 if(!isfinite(gauge)){coupled_heap_state_close(&s);return 3;}
 for(size_t i=0;i<n;i++)s.state[20*i+2]=(double)((long double)s.state[20*i+2]-gauge);
 for(size_t i=0;i<n;i++){
  FullMomentCenteredState decoded;
  if(full_moment_centered_decode(s.state+20*i,s.mass[i],&decoded)){
   coupled_heap_state_close(&s);return 4;
  }
  long double l00=sqrtl(decoded.R[0]),l10=decoded.R[2]/l00;
  long double rest=decoded.R[3]-l10*l10;
  if(!(rest>0.)||!isfinite(rest)){coupled_heap_state_close(&s);return 4;}
  long double l11=sqrtl(rest);
  s.factors[4*i]=(double)l00;s.factors[4*i+1]=0.;
  s.factors[4*i+2]=(double)l10;s.factors[4*i+3]=(double)l11;
  if(!isfinite(s.factors[4*i])||!isfinite(s.factors[4*i+2])||!isfinite(s.factors[4*i+3])||
     s.factors[4*i]<=0.||s.factors[4*i+3]<=0.){coupled_heap_state_close(&s);return 4;}
 }
 *out=s;return 0;
}

int full_moment_coupled_control(const FullMomentState *base,
 const FullMomentCoupledCallbacks *c,FullMomentState *out,FullMomentCoupledControlReport *report){
 if(!base||!c||!out||out==base||!report||!base->header.cells||base->header.cells>INT_MAX/21||
    !base->state||!base->ids||!base->mass||!base->factors||!base->Psi||!base->weight_rate_seed||out->state||
    !c->evaluate||!c->finalize||!c->provider_count||!c->elapsed_seconds||
    c->maximum_step_providers<1||c->maximum_step_providers>10||
    c->maximum_total_providers<1||c->maximum_total_providers>150||
    !isfinite(c->wall_limit_seconds)||c->wall_limit_seconds<=0.)return FM_COUPLED_INVALID;
 size_t n=(size_t)base->header.cells;
 if(base->header.state_elements!=20*n+4||n>SIZE_MAX/21)return FM_COUPLED_INVALID;
 size_t residual_count=21*n;
 FullMomentCoupledControlReport r={0};
 r.providers_before=c->provider_count(c->context);r.providers_after=r.providers_before;
 if(r.providers_before<0){*report=r;return FM_COUPLED_INVALID;}
 FullMomentState history_maps[3]={0};long double *history_residuals[3]={NULL,NULL,NULL};
 FullMomentState mapped={0},trial_owned={0};const FullMomentState *trial=base;
 int history_count=0,active_metric_epoch=-1,prediction_miss_streak=0,prediction_valid=0;
 int prediction_metric_epoch=-1;long double previous_prediction=0.;
 int status=FM_COUPLED_NOT_CONVERGED;
 for(int iteration=0;iteration<c->maximum_step_providers;iteration++){
  double elapsed=c->elapsed_seconds(c->context);
  int wall_expired=!isfinite(elapsed)||elapsed>=c->wall_limit_seconds;
  if(coupled_collective_status(c,wall_expired)){status=FM_COUPLED_WALL_BUDGET;break;}
  int before=c->provider_count(c->context);
  if(coupled_collective_status(c,before>=c->maximum_total_providers)){
   status=FM_COUPLED_PROVIDER_BUDGET;break;
  }
  long double *current_residual=malloc(residual_count*sizeof(*current_residual));
  if(coupled_collective_status(c,!current_residual)){
   free(current_residual);status=FM_COUPLED_EVALUATION_FAILED;r.callback_status=-3;break;
  }
  memset(&r.last_residual,0,sizeof r.last_residual);
  r.callback_status=c->evaluate(c->context,base,trial,iteration,&mapped,&r.last_residual,
                                current_residual,residual_count);
  r.iterations=iteration+1;r.providers_after=c->provider_count(c->context);
  if(coupled_collective_status(c,r.callback_status!=0)){
   free(current_residual);status=FM_COUPLED_EVALUATION_FAILED;break;
  }
  if(coupled_collective_status(c,r.providers_after!=before+1)){
   free(current_residual);status=FM_COUPLED_COUNTER_VIOLATION;break;
  }
  int map_invalid=mapped.header.cells!=base->header.cells||mapped.header.state_elements!=20*n+4||
     mapped.header.run_token!=base->header.run_token||mapped.header.epoch<=base->header.epoch||
     !isfinite(mapped.header.time)||mapped.header.time<=base->header.time||
     !mapped.state||!mapped.ids||!mapped.mass||!mapped.factors||!mapped.Psi||!mapped.weight_rate_seed||
     memcmp(mapped.ids,base->ids,n*sizeof(*base->ids))||memcmp(mapped.mass,base->mass,n*sizeof(*base->mass));
  if(coupled_collective_status(c,map_invalid)){
   free(current_residual);status=FM_COUPLED_EVALUATION_FAILED;r.callback_status=-1;break;
  }
  long double vector_norm=0.;
  if(coupled_collective_status(c,coupled_norm2(residual_count,current_residual,&vector_norm))){
   free(current_residual);status=FM_COUPLED_EVALUATION_FAILED;r.callback_status=-4;break;
  }
  r.last_residual.vector_norm2=vector_norm;
  r.iteration_residual[iteration]=r.last_residual;
  int reset_history=0;
  if(r.last_residual.metric_reset_cells>0||
     (active_metric_epoch>=0&&r.last_residual.metric_epoch!=active_metric_epoch)){
   reset_history=1;prediction_valid=0;prediction_miss_streak=0;
   r.anderson_fallback[iteration]=FM_AA_FALLBACK_METRIC_RESET;
  }
  active_metric_epoch=r.last_residual.metric_epoch;
  if(!reset_history&&prediction_valid&&prediction_metric_epoch==active_metric_epoch){
   int miss=previous_prediction>0.?
    vector_norm>4.*previous_prediction:
    vector_norm>4096.*LDBL_EPSILON;
   prediction_miss_streak=miss?prediction_miss_streak+1:0;
   if(prediction_miss_streak>=2){
    reset_history=1;prediction_valid=0;prediction_miss_streak=0;
    r.anderson_fallback[iteration]=FM_AA_FALLBACK_PREDICTION_RESET;
   }
  }
  if(reset_history)coupled_history_clear(history_maps,history_residuals,&history_count);
  r.prediction_miss_streak[iteration]=prediction_miss_streak;
  r.anderson_predicted_norm2[iteration]=vector_norm;
  if(full_moment_coupled_converged(&r.last_residual)){
   r.converged=1;
   elapsed=c->elapsed_seconds(c->context);
   wall_expired=!isfinite(elapsed)||elapsed>=c->wall_limit_seconds;
   if(coupled_collective_status(c,wall_expired)){free(current_residual);status=FM_COUPLED_WALL_BUDGET;break;}
   r.callback_status=c->finalize(c->context,base,trial,&mapped,&r.last_residual,iteration);
   int after=c->provider_count(c->context);
   if(coupled_collective_status(c,after!=r.providers_after)){r.providers_after=after;free(current_residual);status=FM_COUPLED_COUNTER_VIOLATION;break;}
   if(coupled_collective_status(c,r.callback_status!=0)){free(current_residual);status=FM_COUPLED_FINAL_CHECK_FAILED;break;}
   elapsed=c->elapsed_seconds(c->context);
   wall_expired=!isfinite(elapsed)||elapsed>=c->wall_limit_seconds;
   if(coupled_collective_status(c,wall_expired)){free(current_residual);status=FM_COUPLED_WALL_BUDGET;break;}
   r.final_checks_passed=1;*out=mapped;memset(&mapped,0,sizeof mapped);
   free(current_residual);status=FM_COUPLED_OK;break;
  }
  int maximum_depth=history_count<3?history_count:3;
  const long double *older_residuals[3]={history_residuals[0],history_residuals[1],history_residuals[2]};
  double gamma[3]={0.};int used_depth=0,qr_rank=0,fit=0;
  long double predicted=vector_norm,condition=INFINITY;
  if(maximum_depth>0&&!reset_history){
   if(c->type2_fit)fit=c->type2_fit(c->context,residual_count,maximum_depth,
     current_residual,older_residuals,gamma,&used_depth,&qr_rank,&predicted,&condition);
   else fit=full_moment_coupled_type2_fit(residual_count,maximum_depth,current_residual,
     older_residuals,gamma,&used_depth,&qr_rank,&predicted,&condition);
   if(fit<0){free(current_residual);status=FM_COUPLED_EVALUATION_FAILED;r.callback_status=-5;break;}
  }
  if(fit>0&&(used_depth<1||used_depth>maximum_depth||qr_rank!=used_depth||
     !isfinite(predicted)||!isfinite(condition))){fit=0;used_depth=qr_rank=0;}
  r.anderson_depth[iteration]=fit?used_depth:0;
  r.anderson_qr_rank[iteration]=fit?qr_rank:0;
  r.anderson_condition[iteration]=fit?(double)condition:INFINITY;
  for(int j=0;j<3;j++)r.anderson_gamma[iteration][j]=fit?gamma[j]:0.;
  FullMomentState next_trial={0};int candidate_status=1;double lambda=0.;
  int fallback=r.anderson_fallback[iteration];
  if(fit){
   int candidate_found=0;
   for(int backtrack=0;backtrack<=32;backtrack++){
    lambda=ldexp(1.,-backtrack);
    candidate_status=coupled_type2_candidate(base,&mapped,history_maps,used_depth,gamma,lambda,&next_trial);
    int any_candidate_failure=coupled_collective_status(c,candidate_status!=0);
    if(!any_candidate_failure){candidate_found=1;break;}
    coupled_heap_state_close(&next_trial);
   }
   if(!candidate_found){
    lambda=0.;candidate_status=coupled_type2_candidate(base,&mapped,history_maps,used_depth,gamma,lambda,&next_trial);
    int any_candidate_failure=coupled_collective_status(c,candidate_status!=0);
    if(!any_candidate_failure)candidate_status=0;
    else candidate_status=1;
    fallback=FM_AA_FALLBACK_PURE_MAP;
   }else if(lambda<1.)fallback=FM_AA_FALLBACK_POSITIVITY_BACKTRACK;
  }else{
   lambda=0.;candidate_status=coupled_type2_candidate(base,&mapped,history_maps,0,gamma,lambda,&next_trial);
   int any_candidate_failure=coupled_collective_status(c,candidate_status!=0);
   if(!any_candidate_failure)candidate_status=0;
   if(!fallback)fallback=maximum_depth?FM_AA_FALLBACK_QR_OR_COEFFICIENT:FM_AA_FALLBACK_NO_HISTORY;
  }
  if(coupled_collective_status(c,candidate_status!=0)){
   free(current_residual);coupled_heap_state_close(&next_trial);
   status=FM_COUPLED_EVALUATION_FAILED;r.callback_status=-6;break;
  }
  r.anderson_fallback[iteration]=fallback;r.anderson_lambda[iteration]=lambda;
  if(fit){
   const long double *hist[3]={history_residuals[0],history_residuals[1],history_residuals[2]};
   predicted=coupled_predicted_norm2(residual_count,current_residual,hist,used_depth,gamma,lambda);
   if(!isfinite(predicted)){
    free(current_residual);coupled_heap_state_close(&next_trial);
    status=FM_COUPLED_EVALUATION_FAILED;r.callback_status=-7;break;
   }
   r.anderson_predicted_norm2[iteration]=predicted;
   prediction_valid=lambda>0.;prediction_metric_epoch=active_metric_epoch;previous_prediction=predicted;
  }else {prediction_valid=0;prediction_miss_streak=0;}
  coupled_heap_state_close(&trial_owned);trial_owned=next_trial;memset(&next_trial,0,sizeof next_trial);
  for(int j=2;j>0;j--){
   full_moment_state_close(&history_maps[j]);free(history_residuals[j]);
   history_maps[j]=history_maps[j-1];memset(&history_maps[j-1],0,sizeof history_maps[j-1]);
   history_residuals[j]=history_residuals[j-1];history_residuals[j-1]=NULL;
  }
  history_maps[0]=mapped;memset(&mapped,0,sizeof mapped);
  history_residuals[0]=current_residual;current_residual=NULL;
  if(history_count<3)history_count++;
  trial=&trial_owned;
 }
 full_moment_state_close(&mapped);coupled_heap_state_close(&trial_owned);
 coupled_history_clear(history_maps,history_residuals,&history_count);
 *report=r;return status;
}
