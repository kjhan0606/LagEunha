#ifndef MATERIAL_FULL_INITIAL_BATCH_H
#define MATERIAL_FULL_INITIAL_BATCH_H
#include <stdint.h>
#include "material_shared_trace_batch.h"
#include "gfs_enriched_euler.h"
struct MFDoubleDouble {double hi,lo;};
enum {MF_INITIAL_METRICS=16,MF_BATCH_TRIANGLES=4096};
struct MFInitialBlock {
    MFDoubleDouble mass[ST_MAX*ST_MAX],force[2][ST_MAX],metric[MF_INITIAL_METRICS];
    double minJ;uint64_t points;int bad_point,bad_kind;
};
enum {MF_CURV_METRICS=12};
struct MFCurvatureBlock {
    MFDoubleDouble metric[MF_CURV_METRICS];
    double max_chi,minJ,max_velocity_error;
    uint64_t points;int bad_point,max_chi_point;
};
extern "C" int mf_initial_create(const STTriangle*,size_t,const STCoefficient*,size_t,
    const GfsEnrichQuadrature*,const GfsEnrichConfig*,void**);
extern "C" int mf_initial_batch(void*,size_t,size_t,MFInitialBlock*,double*);
extern "C" void mf_initial_destroy(void*);
extern "C" int mf_curvature_upload(void*,const STCoefficient*,size_t);
extern "C" int mf_curvature_batch(void*,size_t,size_t,MFCurvatureBlock*,double*);
// Three passes at an actual changed RK state. Generating coefficients/delta0
// remain immutable; source delta/delta_dot/mean pressure are stage-local.
enum {MF_EVOL_METRICS=19};
struct MFEvolvedCell {double delta,delta_dot,mean_deltaP,generating_delta;};
struct MFEvolvedBlock {
    MFDoubleDouble force[2][ST_MAX],metric[MF_EVOL_METRICS];
    double minJ,min_j,min_U,max_grad_v,min_length,max_sound,max_omega,max_J;
    uint64_t points;int bad_point,bad_kind;
};
extern "C" int mf_evolved_bind(void*,const STCoefficient*,size_t,size_t);
extern "C" int mf_evolved_stage(void*,const STCoefficient*,size_t,const MFEvolvedCell*,size_t);
extern "C" int mf_evolved_context_range(void*,const MFEvolvedCell*,size_t,size_t);
extern "C" int mf_evolved_batch(void*,int,size_t,size_t,MFEvolvedBlock*,double*);
// Optional bounded point tile: cache fixed shapes and actual current fields
// for one batch of complete source cells across its three dependent passes.
extern "C" int mf_evolved_cache_reserve(void*,size_t);
extern "C" int mf_evolved_cache_build(void*,size_t,size_t,const uint32_t*,MFEvolvedBlock*,double*);
extern "C" int mf_evolved_cache_batch(void*,int,size_t,size_t,MFEvolvedBlock*,double*);
#endif
