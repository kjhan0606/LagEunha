#ifndef LAGFORCE_MOMENT_GRAVITY_STEP_H
#define LAGFORCE_MOMENT_GRAVITY_STEP_H
#include "moment_quadrupole_gravity.h"
#include "moment_affine_kick.h"
#include "moment_material_drift.h"
#ifdef __CUDACC__
#define MS_INLINE static __host__ __device__ __forceinline__
#else
#define MS_INLINE static inline
#endif
/* Diagnostic KDK of the quadrupole moment Hamiltonian only.
 * Input raw17,rc,dt; output raw17 with Psi recomputed from that Hamiltonian.
 * No pressure, mesh, remap, weight law, or external-potential ledger reuse.
 * Caller must explicitly initialize this approximate Psi before comparison.
 */
MS_INLINE int ms_kick(const double *q,double rc,double dt,double *out) {
    double gravity_in[18],source[12],in[25];
    for(int k=0;k<17;k++){gravity_in[k]=q[k];in[k]=q[k];}
    gravity_in[17]=rc;if(mq_source(gravity_in,source))return 1;
    for(int k=0;k<6;k++)in[17+k]=source[k];in[23]=dt;in[24]=0;
    if(mk_kick(in,out))return 1;
    out[15]=q[15]+0.5*((out[12]-q[12])+(out[14]-q[14]));
    return 0;
}
MS_INLINE int ms_step(const double *in,double *out) {
    double first[17],drift[17],middle[17],ginput[18],source[12],rc=in[17],dt=in[18];
    if(ms_kick(in,rc,0.5*dt,first))return 1;
    for(int k=0;k<16;k++)drift[k]=first[k];drift[16]=dt;
    if(md_drift(drift,middle))return 1;middle[16]=first[16];
    if(ms_kick(middle,rc,0.5*dt,out))return 1;
    for(int k=0;k<17;k++)ginput[k]=out[k];ginput[17]=rc;
    if(mq_source(ginput,source))return 1;out[16]=source[7];return 0;
}
#undef MS_INLINE
#endif
