#include "material_compact_thermo_cuda.h"
#include <cfloat>
#include <climits>
namespace {
// Same compensated summation algebra as the user's preserved native CUDA TU.
using DD=MFDoubleDouble;
__device__ DD dd(double a){return {a,0};}
__device__ DD add(DD a,DD b){double s=a.hi+b.hi,z=s-a.hi,e=(a.hi-(s-z))+(b.hi-z);e+=a.lo+b.lo;double h=s+e;return {h,e-(h-s)};}
__global__ void evolved(const compact_geometry::Patch*patch,size_t np,const STCoefficient*c,
    const STCoefficient*generating,size_t nc,const compact_p3::ThermalDensity*field,
    const MFEvolvedCell*cell,size_t cells,const uint64_t*offset,
    const compact_thermo::QuadraturePoint*quadrature,size_t points,
    size_t start,size_t count,int pass,compact_thermo::Config cfg,MFEvolvedBlock*out){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;int tid=threadIdx.x;
    __shared__ double force[2][64][12],metric[64][MF_EVOL_METRICS];
    __shared__ double jac[64],jhat[64],heat[64],grad[64],length[64],sound[64],omega[64];
    __shared__ uint64_t first,nq,cell_id;
    __shared__ int n,bad[64],invalid;
    if(tid==0){
        out[at]={};out[at].minJ=out[at].min_j=out[at].min_U=out[at].min_length=DBL_MAX;out[at].bad_point=-1;
        invalid=s>=np;first=0;nq=0;cell_id=0;n=0;
        if(!invalid){
            first=offset[s];uint64_t end=offset[s+1];cell_id=patch[s].polynomial.cell;n=patch[s].curved?12:10;
            invalid=first>=end||end>points||end-first>INT_MAX||cell_id>=cells||patch[s].curved>1;
            if(!invalid)nq=end-first;
        }
        if(invalid){out[at].bad_point=0;out[at].bad_kind=1;}
    }
    __syncthreads();if(invalid)return;const auto&T=patch[s];
    DD sum={0,0},fx={0,0},fy={0,0};
    for(uint64_t base=0;base<nq;base+=64){
        if(tid<64){
            for(int k=0;k<MF_EVOL_METRICS;k++)metric[tid][k]=0;
            for(int d=0;d<2;d++)for(int a=0;a<12;a++)force[d][tid][a]=0;
            jac[tid]=jhat[tid]=heat[tid]=length[tid]=DBL_MAX;grad[tid]=sound[tid]=omega[tid]=0;bad[tid]=0;
            uint64_t qindex=base+tid;
            if(qindex<nq){
                const auto&q=quadrature[first+qindex];compact_geometry::Basis basis;GfsMaterialPoint12 p,p0;
                double total=0;bool valid=q.weight>0&&isfinite(q.weight);
                for(double l:q.l){valid=valid&&l>=0&&isfinite(l);total+=l;}
                valid=valid&&fabs(total-1)<=32*DBL_EPSILON&&compact_geometry::basis(T,q.l,&basis);
                if(valid)for(int a=0;a<basis.n;a++)valid=valid&&basis.id[a]<nc;
                if(!valid)bad[tid]=2;
                else if(!compact_geometry::point(basis,c,&p))bad[tid]=3;
                else{
                    double B[10];compact_p3::bernstein(q.l,B);
                    double rho=compact_p3::field(field[s].mass,B),e0=compact_p3::field(field[s].internal_energy,B);
                    if(!(rho>0)||!isfinite(rho)||!(e0>0)||!isfinite(e0))bad[tid]=4;
                    else{
                        double jd=p.J*(p.grad_v[0]+p.grad_v[3]),dm=q.weight*rho;
                        jac[tid]=p.J;metric[tid][0]=q.weight;metric[tid][1]=q.weight*p.J;metric[tid][2]=q.weight*jd;metric[tid][3]=dm;
                        double x=p.position[0],y=p.position[1],vx=p.velocity[0],vy=p.velocity[1];
                        double z=x*x+y*y+cfg.eps*cfg.eps,g=pow(z,-1.5),gx=-x*g,gy=-y*g;
                        grad[tid]=sqrt(p.grad_v[0]*p.grad_v[0]+p.grad_v[1]*p.grad_v[1]+p.grad_v[2]*p.grad_v[2]+p.grad_v[3]*p.grad_v[3]);omega[tid]=sqrt(g);
                        double refh=DBL_MAX;for(int a=0;a<3;a++)refh=fmin(refh,1/hypot(T.polynomial.grad[a][0],T.polynomial.grad[a][1]));
                        double tr=p.F[0]*p.F[0]+p.F[1]*p.F[1]+p.F[2]*p.F[2]+p.F[3]*p.F[3];
                        double maxsv=sqrt(.5*(tr+sqrt(fmax(0.,tr*tr-4*p.J*p.J))));length[tid]=refh*p.J/maxsv;
                        if(pass>0){
                            const auto&C=cell[cell_id];
                            if(!compact_geometry::point(basis,generating,&p0))bad[tid]=5;
                            else{
                                auto thermal=compact_thermo::evaluate(cfg,e0,rho,p0.J+C.generating_delta,p.J+C.delta,jd+C.delta_dot);
                                if(!thermal.valid)bad[tid]=6;
                                else{
                                    jhat[tid]=thermal.j;heat[tid]=q.weight*thermal.energy_reference;sound[tid]=sqrt(thermal.sound_squared);
                                    metric[tid][4]=q.weight*thermal.delta_pressure;metric[tid][5]=heat[tid];metric[tid][6]=q.weight*thermal.energy_dot_reference;metric[tid][7]=q.weight*thermal.j;
                                    metric[tid][8]=.5*dm*(vx*vx+vy*vy);metric[tid][9]=dm*(x*vy-y*vx);metric[tid][10]=-dm/sqrt(z);metric[tid][11]=dm*(vx*gx+vy*gy);metric[tid][12]=dm*(x*gy-y*gx);
                                    metric[tid][13]=dm*x;metric[tid][14]=dm*y;metric[tid][15]=dm*vx;metric[tid][16]=dm*vy;metric[tid][17]=dm*(x*x+y*y);metric[tid][18]=2*dm*(x*vx+y*vy);
                                    if(pass==2){
                                        double residual=thermal.delta_pressure-C.mean_deltaP,sigma[3]={-residual,0,-residual};
                                        double A[12],Ax[12],Ay[12],F[2][12]={},u=0,gwork=0;
                                        for(int a=0;a<n;a++){A[a]=basis.N[a].v;Ax[a]=basis.N[a].x;Ay[a]=basis.N[a].y;}
                                        if(!gfs_material12_stress_at_point(n,&p,A,Ax,Ay,rho,q.weight,sigma,gx,gy,F,&u,&gwork))bad[tid]=7;
                                        for(int d=0;d<2;d++)for(int a=0;a<n;a++)force[d][tid][a]=F[d][a];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        __syncthreads();int count_here=int(nq-base<64?nq-base:64);
        if(tid<MF_EVOL_METRICS)for(int k=0;k<count_here;k++)sum=add(sum,dd(metric[k][tid]));
        if(tid<n&&pass==2)for(int k=0;k<count_here;k++){fx=add(fx,dd(force[0][k][tid]));fy=add(fy,dd(force[1][k][tid]));}
        if(tid==0)for(int k=0;k<count_here;k++){
            if(bad[k]&&out[at].bad_point<0){out[at].bad_point=int(base)+k;out[at].bad_kind=bad[k];}if(!bad[k])out[at].points++;
            out[at].minJ=fmin(out[at].minJ,jac[k]);out[at].max_J=fmax(out[at].max_J,jac[k]);out[at].min_j=fmin(out[at].min_j,jhat[k]);out[at].min_U=fmin(out[at].min_U,heat[k]);
            out[at].max_grad_v=fmax(out[at].max_grad_v,grad[k]);out[at].min_length=fmin(out[at].min_length,length[k]);out[at].max_sound=fmax(out[at].max_sound,sound[k]);out[at].max_omega=fmax(out[at].max_omega,omega[k]);
        }
        __syncthreads();
    }
    if(tid<MF_EVOL_METRICS)out[at].metric[tid]=sum;
    if(tid<n&&pass==2){out[at].force[0][tid]=fx;out[at].force[1][tid]=fy;}
}
}
extern "C" int compact_material_batch_device(const compact_geometry::Patch*p,size_t np,
    const STCoefficient*c,const STCoefficient*g,size_t nc,const compact_p3::ThermalDensity*f,
    const MFEvolvedCell*cell,size_t cells,const uint64_t*offset,
    const compact_thermo::QuadraturePoint*q,size_t points,size_t start,size_t count,int pass,
    compact_thermo::Config cfg,MFEvolvedBlock*out,cudaStream_t stream){
    if(!count)return int(cudaSuccess);
    if(!p||!c||!g||!nc||!f||!cell||!cells||!offset||!q||!points||!out||start>np||count>np-start||count>MF_BATCH_TRIANGLES||pass<0||pass>2)return int(cudaErrorInvalidValue);
    if(!(cfg.gamma>1)||!(cfg.background_pressure>0)||!(cfg.eps>0)||!isfinite(cfg.gamma)||!isfinite(cfg.background_pressure)||!isfinite(cfg.eps))return int(cudaErrorInvalidValue);
    evolved<<<unsigned(count),128,0,stream>>>(p,np,c,g,nc,f,cell,cells,offset,q,points,start,count,pass,cfg,out);
    return int(cudaPeekAtLastError());
}
