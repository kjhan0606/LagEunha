#ifndef LAGFORCE_MOMENT_STATE_LAYOUT_H
#define LAGFORCE_MOMENT_STATE_LAYOUT_H
#include <math.h>

/* Canonical planar raw volume moments, also used for signed moment rates.
 * D_ab = integral rho x_a v_b; Dxy and Dyx are independent.
 * U=Ehyd-(Pixx+Piyy)/2: unresolved velocity variance is not heat.
 * Psi is separate and explicitly tagged. A diagnostic reconstruction is not
 * a transported state. Legacy wm17 bytes and particle structs are unchanged.
 * Generator x[2], additive power W[length^2], factor L[4], stable ID, frozen
 * flag and stage/geometry epochs belong outside this raw array. No geometry
 * weight conversion, physical-state validation or serialization is done here.
 */
enum LagMomentRaw18Index {
    LMR_M=0, LMR_QX=1, LMR_QY=2,
    LMR_QXX=3, LMR_QXY=4, LMR_QYX=5, LMR_QYY=6,
    LMR_PX=7, LMR_PY=8,
    LMR_DXX=9, LMR_DXY=10, LMR_DYX=11, LMR_DYY=12,
    LMR_PIXX=13, LMR_PIXY=14, LMR_PIYX=15, LMR_PIYY=16,
    LMR_EHYD=17, LMR_COUNT=18
};
enum LagMomentPsiMode {
    LMP_ABSENT=0, LMP_DIAGNOSTIC_RECONSTRUCTION=1, LMP_TRANSPORTED_MOMENT=2
};
typedef struct LagMomentPsi {
    int mode;
    double value;
} LagMomentPsi;
enum LagMomentAdapterStatus {
    LMA_OK=0, LMA_NULL_POINTER=1, LMA_NONFINITE=2,
    LMA_PSI_ABSENT=3, LMA_INVALID_PSI_MODE=4, LMA_ASYMMETRIC=5
};

#ifdef __CUDACC__
#define LML_INLINE static __host__ __device__ __forceinline__
#else
#define LML_INLINE static inline
#endif

LML_INLINE int lml_check_psi_mode(int mode) {
    if(mode==LMP_ABSENT)return LMA_PSI_ABSENT;
    if(mode!=LMP_DIAGNOSTIC_RECONSTRUCTION && mode!=LMP_TRANSPORTED_MOMENT)
        return LMA_INVALID_PSI_MODE;
    return LMA_OK;
}

/* Exact, transactional conversion. All validation precedes output writes.
 * Q and Pi must be exactly symmetric; no averaging or implicit projection.
 * No positivity checks: signed RHS rates use this same layout.
 * The result metadata retains the input Psi mode; the adapter never upgrades
 * a reconstructed diagnostic to a transported moment. Temporary arrays also
 * permit input/output array overlap without corrupting unread input. Psi
 * input/output may alias each other, but not the double output array.
 */
LML_INLINE int lml_raw18_to_wm17(const double *raw,const LagMomentPsi *psi,
                                double *out,LagMomentPsi *out_psi) {
    if(!raw || !psi || !out || !out_psi)return LMA_NULL_POINTER;
    const int status=lml_check_psi_mode(psi->mode);
    if(status)return status;
    if(!isfinite(psi->value))return LMA_NONFINITE;
    for(int k=0;k<LMR_COUNT;k++)if(!isfinite(raw[k]))return LMA_NONFINITE;
    if(raw[LMR_QXY]!=raw[LMR_QYX] || raw[LMR_PIXY]!=raw[LMR_PIYX] ||
       (raw[LMR_QXY]==0. && signbit(raw[LMR_QXY])!=signbit(raw[LMR_QYX])) ||
       (raw[LMR_PIXY]==0. && signbit(raw[LMR_PIXY])!=signbit(raw[LMR_PIYX])))
        return LMA_ASYMMETRIC;
    const double tmp[17]={raw[0],raw[1],raw[2],raw[3],raw[4],raw[6],
        raw[7],raw[8],raw[9],raw[10],raw[11],raw[12],
        raw[13],raw[14],raw[16],raw[17],psi->value};
    const int mode=psi->mode;const double value=psi->value;
    for(int k=0;k<17;k++)out[k]=tmp[k];
    out_psi->mode=mode;out_psi->value=value;
    return LMA_OK;
}

/* Legacy wm17 has a Psi value but no semantic tag in its array bytes.
 * The caller must supply the source's documented mode explicitly.
 */
LML_INLINE int lml_wm17_to_raw18(const double *wm,int source_psi_mode,
                                double *out,LagMomentPsi *out_psi) {
    if(!wm || !out || !out_psi)return LMA_NULL_POINTER;
    const int status=lml_check_psi_mode(source_psi_mode);
    if(status)return status;
    for(int k=0;k<17;k++)if(!isfinite(wm[k]))return LMA_NONFINITE;
    const double tmp[18]={wm[0],wm[1],wm[2],wm[3],wm[4],wm[4],wm[5],
        wm[6],wm[7],wm[8],wm[9],wm[10],wm[11],wm[12],wm[13],wm[13],wm[14],wm[15]};
    const double value=wm[16];
    for(int k=0;k<LMR_COUNT;k++)out[k]=tmp[k];
    out_psi->mode=source_psi_mode;out_psi->value=value;
    return LMA_OK;
}
#undef LML_INLINE
#endif
