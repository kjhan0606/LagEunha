#ifndef LAGFORCE_TCF_FACE_TERMS_H
#define LAGFORCE_TCF_FACE_TERMS_H
/* FP64 default TCF face algebra. Caller supplies the actual reconstructed Q
 * and common face geometry. Reconstruction, ownership, boundary policy and
 * accumulation remain caller responsibilities. Keep normal/tangent separate
 * so production accumulation order and the historical NOTAN option survive. */
#ifdef __CUDACC__
#define TCF_TERMS_INLINE static __host__ __device__ __forceinline__
#else
#define TCF_TERMS_INLINE static inline
#endif
typedef struct {
    double normal_x,normal_y,tangent_x,tangent_y,energy,volume,heat;
} TcfFaceTerms;
TCF_TERMS_INLINE TcfFaceTerms tcf_face_terms(double A,double d,double ex,double ey,
        double tx,double ty,double Pi,double Pj,double vix,double viy,
        double vjx,double vjy,double ui,double uj,double Q){
    TcfFaceTerms o;
    double Pb=.5*(Pi+Pj);
    o.normal_x=-A*(Pb+Q)*ex;o.normal_y=-A*(Pb+Q)*ey;
    o.tangent_x=-A*(Pj-Pi)*tx/d;o.tangent_y=-A*(Pj-Pi)*ty/d;
    o.energy=-A*(.5*(Pj*ui+Pi*uj)+(Pj*(vix*tx+viy*ty)-Pi*(vjx*tx+vjy*ty))/d+.5*Q*(ui+uj));
    o.volume=A*(.5*(ui+uj)+((vix-vjx)*tx+(viy-vjy)*ty)/d);
    o.heat=.5*A*Q*(ui-uj);
    return o;
}
#undef TCF_TERMS_INLINE
#endif
