#ifndef KH_RESTART_GUARD_H
#define KH_RESTART_GUARD_H

/* Call with a complete current-ABI stress record, never an unknown stride. */
static inline int kh_restart_has_extended_state(const treevorostressrk4particletype *p){
    int moment_tagged=0;
    if(sizeof(PosType)==sizeof(double) && offsetof(RK4,w2backup)==24*sizeof(double)){
        GfsMomentPayload moment;
        memcpy(&moment,&p->rk4,sizeof(moment));
        moment_tagged=(moment.tag==GFS_MOMENT_TAG);
    }
    return GFS_T64_TAGGED(*p) || moment_tagged;
}

static inline int kh_restart_stress_stride_valid(int has_stress, size_t stride){
    return !has_stress || stride==sizeof(treevorostressrk4particletype);
}
#endif
