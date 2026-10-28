/* RK4 bookkeeping, rejection and fieldwise checkpoint regression.
 * Synthetic nonzero hydro force/heat plus auxiliary oscillator. This is not
 * the production integrator, and does not claim full snapshot I/O coverage.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "reversible_memory_state.h"
#include "reversible_memory_checkpoint.h"
#define REQUIRE(test) do { if(!(test)){fprintf(stderr,"failed line %d\n",__LINE__);exit(2);} } while(0)
typedef struct {double x,v,thermal;rm_owner_state memory;} state;
/* m=1, B=2, a_hydro=.3, thermal rate=.07. Ehydro_dot=.3*v+.07.
 * Auxiliary acceleration MUST NOT enter the hydro kinetic subtraction.
 */
static double thermal_rate(double v) {return (.3*v+.07)-v*.3;}
static void trial(state *s,double dt) {
    state initial=*s;double kx[4],kv[4],ku[4];
    REQUIRE(rms_begin(&s->memory,dt)==0);
    for(int stage=1;stage<=4;stage++) {
        double a=.3-2*s->memory.z[0];
        kx[stage-1]=dt*s->v;kv[stage-1]=dt*a;ku[stage-1]=dt*thermal_rate(s->v);
        REQUIRE(rms_stage(&s->memory,stage,2*s->v,0)==0);
        if(stage<4){double f=stage==3?1.:.5;
            s->x=initial.x+f*kx[stage-1];s->v=initial.v+f*kv[stage-1];s->thermal=initial.thermal+f*ku[stage-1];}
    }
    s->x=initial.x+(kx[0]+2*kx[1]+2*kx[2]+kx[3])/6;
    s->v=initial.v+(kv[0]+2*kv[1]+2*kv[2]+kv[3])/6;
    s->thermal=initial.thermal+(ku[0]+2*ku[1]+2*ku[2]+ku[3])/6;
}
static void accepted(state *s,double dt){trial(s,dt);REQUIRE(rms_accept(&s->memory)==0);}
static void compare(const state *a,const state *b){
    REQUIRE(a->x==b->x && a->v==b->v && a->thermal==b->thermal);
    REQUIRE(a->memory.id==b->memory.id && a->memory.phase==b->memory.phase);
    for(int j=0;j<2;j++)REQUIRE(a->memory.z[j]==b->memory.z[j]);
}
int main(void){
    state initial={0};initial.x=.4;initial.v=.7;initial.thermal=.2;
    REQUIRE(rms_init(&initial.memory,539,.03,-.04)==0);
    state direct=initial,retry=initial;accepted(&direct,.005);trial(&retry,.01);
    uint64_t id;double values[2];REQUIRE(rms_pack(&retry.memory,&id,values)==1);
    REQUIRE(rms_reject(&retry.memory)==0);
    retry.x=initial.x;retry.v=initial.v;retry.thermal=initial.thermal;
    compare(&retry,&initial);accepted(&retry,.005);compare(&retry,&direct);
    /* Round-trip through a real sidecar; hydro pairing is represented by a
     * supplied fingerprint here, not a production snapshot writer. */
    state uninterrupted=direct,restarted=direct;
    REQUIRE(rms_pack(&direct.memory,&id,values)==0);
    memset(&restarted.memory,0,sizeof(restarted.memory));
    REQUIRE(rms_init(&restarted.memory,id,values[0],values[1])==0);
    REQUIRE(rms_checkpoint_write("memory_state_test.bin",&direct.memory,1,12345,.005)==0);
    restarted.memory.z[0]=987.;restarted.memory.z[1]=-876.;
    REQUIRE(rms_checkpoint_read("memory_state_test.bin",&restarted.memory,1,12345,.005)==0);
    compare(&direct,&restarted);
    rm_owner_state checkpoint_saved=restarted.memory;
    REQUIRE(rms_checkpoint_read("memory_state_test.bin",&restarted.memory,1,12346,.005)==1);
    REQUIRE(memcmp(&checkpoint_saved,&restarted.memory,sizeof(checkpoint_saved))==0);
    REQUIRE(rms_checkpoint_read("memory_state_test.bin",&restarted.memory,1,12345,.006)==1);
    REQUIRE(memcmp(&checkpoint_saved,&restarted.memory,sizeof(checkpoint_saved))==0);
    FILE *source=fopen("memory_state_test.bin","rb"),*truncated=fopen("memory_state_truncated.bin","wx");
    REQUIRE(source && truncated);
    for(int j=0;j<60;j++){int byte=fgetc(source);REQUIRE(byte!=EOF);REQUIRE(fputc(byte,truncated)!=EOF);}
    REQUIRE(fclose(source)==0 && fclose(truncated)==0);
    REQUIRE(rms_checkpoint_read("memory_state_truncated.bin",&restarted.memory,1,12345,.005)==1);
    REQUIRE(memcmp(&checkpoint_saved,&restarted.memory,sizeof(checkpoint_saved))==0);
    for(int k=0;k<20;k++){accepted(&uninterrupted,.005);accepted(&restarted,.005);}
    compare(&uninterrupted,&restarted);
    REQUIRE(fabs(uninterrupted.thermal-(initial.thermal+.07*.105))<1e-14);
    rm_owner_state probe=initial.memory,saved;
    REQUIRE(rms_stage(&probe,1,0,0)==1);REQUIRE(rms_begin(&probe,.01)==0);
    saved=probe;REQUIRE(rms_stage(&probe,1,NAN,0)==1);
    REQUIRE(memcmp(&probe,&saved,sizeof(probe))==0);
    REQUIRE(rms_stage(&probe,2,0,0)==1);REQUIRE(rms_accept(&probe)==1);
    REQUIRE(rms_reject(&probe)==0);
    /* Permuting records preserves persistent ID/state associations. */
    rm_owner_state records[3];
    for(int j=0;j<3;j++)REQUIRE(rms_init(records+j,100+j,.1*j,-.2*j)==0);
    rm_owner_state tmp=records[0];records[0]=records[2];records[2]=tmp;
    for(int j=0;j<3;j++)REQUIRE(records[j].z[0]==.1*(records[j].id-100));
    printf("{\"scope\":\"C auxiliary state bookkeeping, actual sidecar codec and synthetic coupled RK4 only; production hooks not connected\",\"rejection_restores_state\":true,\"checkpoint_roundtrip_bitwise_trajectory\":true,\"hash_time_truncation_rejected_without_mutation\":true,\"thermal_energy_isolated\":true,\"invalid_stage_and_nan_rejected\":true,\"persistent_id_permutation\":true,\"final_thermal\":%.17g,\"final_auxiliary_energy\":%.17g}\n",uninterrupted.thermal,rms_energy(&uninterrupted.memory));
    return 0;
}
