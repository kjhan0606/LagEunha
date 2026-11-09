#ifndef GFS_MATERIAL_FACE_RESPONSE_H
#define GFS_MATERIAL_FACE_RESPONSE_H
#include "gfs_exact_riemann_ale.h"
#include "gfs_exact_wave_impedance.h"
#ifdef __CUDACC__
#define GFS_MATERIAL_HD static __host__ __device__ inline
#else
#define GFS_MATERIAL_HD static inline
#endif
/* Slip-wall velocity is an affine map, u_g = R u + c, where
 * R = I - 2 n n^T and c = 2 (n dot wall_velocity) n.
 * Use the SAME R in the force pullback. The returned offset power is the
 * term u_g dot F_g - u dot R^T F_g; it belongs in the boundary ledger.
 * These helpers do not prescribe pressure ghosts or replace wall geometry.
 */
GFS_MATERIAL_HD int gfs_material_reflect_vector(int dim,const double n[3],
 const double v[3],double out[3]){
 if(dim<1||dim>3||!n||!v||!out)return 0;
 double norm=0,dot=0,q[3]={0};
 for(int a=0;a<dim;a++){
  if(!isfinite(n[a])||!isfinite(v[a]))return 0;
  norm+=n[a]*n[a];dot+=n[a]*v[a];
 }
 if(fabs(norm-1)>128*DBL_EPSILON||!isfinite(dot))return 0;
 for(int a=0;a<dim;a++){q[a]=v[a]-2*dot*n[a];if(!isfinite(q[a]))return 0;}
 for(int a=0;a<3;a++)out[a]=q[a];return 1;
}
GFS_MATERIAL_HD int gfs_material_slip_velocity(int dim,const double n[3],
 const double wall[3],const double u[3],double ghost[3]){
 if(!wall||!ghost)return 0;double q[3],wn=0;
 if(!gfs_material_reflect_vector(dim,n,u,q))return 0;
 for(int a=0;a<dim;a++){if(!isfinite(wall[a]))return 0;wn+=n[a]*wall[a];}
 for(int a=0;a<dim;a++){q[a]+=2*wn*n[a];if(!isfinite(q[a]))return 0;}
 for(int a=0;a<3;a++)ghost[a]=q[a];return 1;
}
GFS_MATERIAL_HD int gfs_material_slip_force_pullback(int dim,const double n[3],
 const double wall[3],const double ghost_force[3],double real_force[3],double *offset_power){
 if(!wall||!real_force||!offset_power)return 0;double q[3],wn=0,fn=0;
 if(!gfs_material_reflect_vector(dim,n,ghost_force,q))return 0;
 for(int a=0;a<dim;a++){if(!isfinite(wall[a]))return 0;wn+=n[a]*wall[a];fn+=n[a]*ghost_force[a];}
 double power=2*wn*fn;if(!isfinite(power))return 0;
 for(int a=0;a<3;a++)real_force[a]=q[a];*offset_power=power;return 1;
}
/* Primitive-specific reconstruction weights, shared by host/device gather
 * and matching adjoint scatter. Caller validates finite alpha in [0,1]. */
typedef struct {double rho,pressure,velocity,reversible;} GfsMaterialTraceWeights;
GFS_MATERIAL_HD GfsMaterialTraceWeights gfs_material_trace_weights(
 double base,int owner_entry,const double *alpha3){
 GfsMaterialTraceWeights w={base,base,base,base};
 if(alpha3){w.rho=alpha3[0]*base+(1-alpha3[0])*owner_entry;
  w.pressure=alpha3[1]*base+(1-alpha3[1])*owner_entry;
  w.velocity=alpha3[2]*base+(1-alpha3[2])*owner_entry;}
#ifndef GFS_MATERIAL_SPLIT_REVERSIBLE
 w.reversible=w.velocity;
#elif GFS_MATERIAL_SPLIT_REVERSIBLE == 2
 /* Reuse the thermodynamic reconstruction in the reversible pair.
  * Velocity limiting affects only the dissipative jump. */
 w.reversible=w.pressure;
#elif GFS_MATERIAL_SPLIT_REVERSIBLE != 1
#error "Unsupported experimental reversible trace policy"
#endif
 return w;
}
/* Experimental operator coefficients, not a production integrator.
 * Return 1: material contact; 2: vacuum fan (no shared contact); 0: invalid.
 * Heat coefficients must be assembled with the matching adjoint force and
 * volume maps. They alone do not specify mesh motion or prove entropy accuracy.
 */
typedef struct {
    double pstar,ustar,zleft,zright,alpha_left,alpha_right;
    double mean_speed,pressure_speed,stress,heat_left,heat_right;
} GfsMaterialFaceResponse;
GFS_MATERIAL_HD int gfs_material_face_response(int dim,GfsFvPrimitive l,
        GfsFvPrimitive r,const double normal[3],double gamma,
        GfsMaterialFaceResponse *out){
    if(!out)return 0;
    GfsMaterialFaceResponse q={0};*out=q;
    GfsExactFan f;
    if(!gfs_exact_riemann_prepare(dim,l,r,normal,gamma,&f))return 0;
    if(f.star.vacuum)return 2;
    if(!gfs_exact_wave_impedance(l.rho,l.p,f.star.p,gamma,&q.zleft)||
       !gfs_exact_wave_impedance(r.rho,r.p,f.star.p,gamma,&q.zright))return 0;
    double z=q.zleft+q.zright;
    if(!isfinite(z)||!(z>0))return 0;
    q.alpha_left=q.zleft/z;q.alpha_right=q.zright/z;
    q.pstar=f.star.p;q.ustar=f.star.u;
    q.mean_speed=q.alpha_left*f.ul+q.alpha_right*f.ur;
    q.pressure_speed=(l.p-r.p)/z;
    double jump=f.ur-f.ul;
    q.stress=q.zleft*q.alpha_right*jump;
    double dl=q.pressure_speed+q.alpha_right*jump;
    double dr=q.pressure_speed-q.alpha_left*jump;
    q.heat_left=q.zleft*dl*dl;q.heat_right=q.zright*dr*dr;
    if(!isfinite(q.mean_speed)||!isfinite(q.pressure_speed)||
       !isfinite(q.stress)||!isfinite(q.heat_left)||!isfinite(q.heat_right))return 0;
    *out=q;return 1;
}
/* Integrated coefficients for the two reconstruction transposes and owners.
 * Scatter side_force using the EXACT velocity reconstruction coefficients,
 * and side_volume using the EXACT pressure reconstruction coefficients. These
 * maps may differ; density reconstruction enters only the face response.
 * Scatter owner terms only to the two face owners. No mass-rate output exists.
 * The pressure argument here is the owner CELL pressure, not reconstructed
 * face pressure: this is essential for matching -P_cell*dV_cell in Udot.
 */
typedef struct {
    double side_force[2][3],side_volume[2];
    double pressure_force[2][3],dissipative_force[2][3];
    double owner_volume[2],owner_heat[2];
} GfsMaterialFaceExchange;
GFS_MATERIAL_HD int gfs_material_face_exchange_drive(int dim,
        const double normal[3],double area,double dp,
        const GfsMaterialFaceResponse *r,
        GfsMaterialFaceExchange *out){
    if(!r||!out||dim<1||dim>3||!(area>0)||!isfinite(area)||
       !isfinite(dp))return 0;
    GfsMaterialFaceExchange q={0};double norm=0;
    for(int a=0;a<dim;a++)norm+=normal[a]*normal[a];
    if(!isfinite(norm)||fabs(norm-1)>128*DBL_EPSILON)return 0;
    for(int side=0;side<2;side++){
        double sign=side?-1:1,alpha=side?r->alpha_right:r->alpha_left;
        for(int a=0;a<dim;a++){
            q.side_force[side][a]=area*(alpha*dp+sign*r->stress)*normal[a];
            q.pressure_force[side][a]=area*alpha*dp*normal[a];
            q.dissipative_force[side][a]=area*sign*r->stress*normal[a];
            if(!isfinite(q.side_force[side][a]))return 0;
        }
        q.side_volume[side]=sign*area*r->pressure_speed;
        q.owner_volume[side]=sign*area*r->mean_speed;
        q.owner_heat[side]=area*(side?r->heat_right:r->heat_left);
        if(!isfinite(q.side_volume[side])||!isfinite(q.owner_volume[side])||
           !isfinite(q.owner_heat[side])||q.owner_heat[side]<0)return 0;
    }
    *out=q;return 1;
}
/* Physical-pressure wrapper retains the original admissibility contract.
 * A generalized potential drive can be gauge shifted, so the low-level
 * exchange accepts its difference only; it is never used as an EOS pressure. */
GFS_MATERIAL_HD int gfs_material_face_exchange(int dim,
        const double normal[3],double area,double owner_pressure_left,
        double owner_pressure_right,const GfsMaterialFaceResponse *r,
        GfsMaterialFaceExchange *out){
    if(!(owner_pressure_left>0)||!(owner_pressure_right>0)||
       !isfinite(owner_pressure_left)||!isfinite(owner_pressure_right))return 0;
    return gfs_material_face_exchange_drive(dim,normal,area,
        owner_pressure_left-owner_pressure_right,r,out);
}
/* Pressure force and pressure work use one reversible trace; damping and
 * its positive heat use the limited jump trace. */
GFS_MATERIAL_HD double gfs_material_scatter_force(GfsMaterialTraceWeights w,
 const GfsMaterialFaceExchange*x,int side,int axis){
#ifdef GFS_MATERIAL_SPLIT_REVERSIBLE
 return w.reversible*x->pressure_force[side][axis]+w.velocity*x->dissipative_force[side][axis];
#else
 return w.velocity*x->side_force[side][axis];
#endif
}
#undef GFS_MATERIAL_HD
#endif
