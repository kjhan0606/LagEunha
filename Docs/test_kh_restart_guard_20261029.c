#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "eunha.h"
#include "Exam/gfs_moment.h"
#include "Exam/gfs_thermal64.h"
#include "Exam/KH/kh_restart_guard.h"

int main(void){
    treevorostressrk4particletype p;
    memset(&p,0,sizeof(p));
    if(!kh_restart_stress_stride_valid(0,sizeof(treevorork4particletype)) ||
       !kh_restart_stress_stride_valid(1,sizeof(p)) ||
       kh_restart_stress_stride_valid(1,sizeof(p)-1) ||
       kh_restart_stress_stride_valid(1,sizeof(p)+1) ||
       kh_restart_has_extended_state(&p)) return 1;
    p.stress.E_inv_yy=GFS_T64_TAG;
    if(!kh_restart_has_extended_state(&p)) return 2;
    memset(&p,0,sizeof(p));
    if(sizeof(PosType)!=sizeof(double) || offsetof(RK4,w2backup)!=24*sizeof(double)) return 3;
    GfsMomentPayload m;
    memset(&m,0,sizeof(m));
    m.tag=GFS_MOMENT_TAG;
    memcpy(&p.rk4,&m,sizeof(m));
    if(!kh_restart_has_extended_state(&p)) return 4;
    puts("restart guard: legacy/current accepted; unknown stress and tagged state detected");
    return 0;
}
