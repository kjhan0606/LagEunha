// Small-remainder and exact-predicate fixture. No accepted source/run is read.
#include "check_compact_relative_probe.h"
#include "material_compact_predicate.h"
#include <quadmath.h>
#include <iostream>
#include <iomanip>
namespace cg=compact_geometry;using Q=__float128;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void cuda_need(cudaError_t e){need(e==cudaSuccess,cudaGetErrorString(e));}
Q L(int node,Q u){Q out=1;for(int j=0;j<4;j++)if(j!=node)out*= (3*u-j)/Q(node-j);return out;}
Q dL(int node,Q u){Q out=0;for(int j=0;j<4;j++)if(j!=node){Q term=Q(3)/Q(node-j);for(int k=0;k<4;k++)if(k!=node&&k!=j)term*=(3*u-k)/Q(node-k);out+=term;}return out;}
std::array<Q,2>reference(const cg::FrozenTrace&f,const RelativeInput&q,int k){
    Q u=Q(q.lo)+Q(q.t)*(Q(q.hi)-q.lo),value=0,derivative=0;
    for(int j=0;j<4;j++){value+=Q(f.coefficient[j].x[k])*L(j,u);derivative+=Q(f.coefficient[j].x[k])*dL(j,u)*(Q(q.hi)-q.lo);}
    value-=(1-Q(q.t))*Q(q.B[k])+Q(q.t)*Q(q.C[k]);derivative-=Q(q.C[k])-q.B[k];return {value,derivative};
}
Q shape(int a,const double*l){int index[3];gfs_material_p3_index(a,index);Q out=1;for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)out*=(3*Q(l[k])-j)/Q(index[k]-j);return out;}
template<class T>struct Device {T*ptr=nullptr;size_t n;Device(size_t count):n(count){cuda_need(cudaMalloc(reinterpret_cast<void**>(&ptr),n*sizeof(T)));}~Device(){cudaFree(ptr);}void put(const T*x){cuda_need(cudaMemcpy(ptr,x,n*sizeof(T),cudaMemcpyHostToDevice));}std::vector<T>get(){std::vector<T>out(n);cuda_need(cudaMemcpy(out.data(),ptr,n*sizeof(T),cudaMemcpyDeviceToHost));return out;}};
int main(){try{
    cg::Patch p;p.curved=1;p.modes.component=0;p.modes.id[0]=10;p.modes.id[1]=11;p.modes.scale[0]=ldexp(1.,-36);p.modes.scale[1]=ldexp(1.,-37);
    double vertices[3][2]={{12,-11.90625},{12,-12},{12.09375,-12}};for(int a=0;a<3;a++)for(int k=0;k<2;k++)p.polynomial.vertex[a][k]=vertices[a][k];compact_p3::geometry(p.polynomial);
    for(int a=0;a<10;a++)p.polynomial.id[a]=a;
    STEdge e={};e.rank=0;e.A[0]=12;e.A[1]=-12;e.D[0]=.09375;e.eps=.047;std::vector<STCoefficient>old(4);
    for(int j=0;j<4;j++){e.poly[j]=j;old[j].x[0]=12+j*.03125;old[j].x[1]=-12;}
    old[1].x[0]+=p.modes.scale[0];old[2].x[0]-=.5*p.modes.scale[0];old[1].x[1]-=p.modes.scale[0];old[2].x[1]+=.25*p.modes.scale[0];
    p.trace=cg::freeze_trace(e,old,0,1,5);cg::interpolate_chart(p);
    std::vector<RelativeInput>input;
    for(auto interval:{std::array<double,2>{0,1},std::array<double,2>{.13,.83},std::array<double,2>{.75,.25}})for(int j=0;j<=128;j++){
        RelativeInput q={};q.t=j/128.;q.lo=interval[0];q.hi=interval[1];STJet b[2],c[2];need(cg::trace(p.trace,STJet(q.lo),b)&&cg::trace(p.trace,STJet(q.hi),c),"frozen endpoint values");
        for(int k=0;k<2;k++){q.B[k]=b[k].v;q.C[k]=c[k].v;}double s=j%2?.75:.5;q.l[0]=1-s;q.l[1]=s*(1-q.t);q.l[2]=s*q.t;input.push_back(q);
    }
    Q new_error=0,old_error=0,mode_error=0,endpoint_residual=0;std::vector<RelativeOutput>cpu(input.size());
    for(size_t j=0;j<input.size();j++){
        const auto&q=input[j];auto&v=cpu[j];cg::Basis basis;
        v.valid=cg::relative_trace(p.trace,STJet(q.t,1,0),q.lo,q.hi,q.B,q.C,v.relative)&&cg::basis(p,q.l,&basis);need(v.valid,"stable CPU relative geometry");
        STJet gamma[2];need(cg::trace(p.trace,STJet(q.lo)+STJet(q.t,1,0)*STJet(q.hi-q.lo),gamma),"old absolute trace reference");
        for(int k=0;k<2;k++){
            v.mode[k]=basis.N[10+k];auto exact=reference(p.trace,q,k);STJet chord=(STJet(1)-STJet(q.t,1,0))*STJet(q.B[k])+STJet(q.t,1,0)*STJet(q.C[k]);STJet unstable=gamma[k]-chord;
            new_error=fmaxq(new_error,fabsq(Q(v.relative[k].v)-exact[0])/Q(p.modes.scale[k]));new_error=fmaxq(new_error,fabsq(Q(v.relative[k].x)-exact[1])/Q(p.modes.scale[k]));
            old_error=fmaxq(old_error,fabsq(Q(unstable.v)-exact[0])/Q(p.modes.scale[k]));old_error=fmaxq(old_error,fabsq(Q(unstable.x)-exact[1])/Q(p.modes.scale[k]));
            if(q.t==0||q.t==1)endpoint_residual=fmaxq(endpoint_residual,fabsq(Q(v.relative[k].v)));
            RelativeInput full={};full.lo=0;full.hi=1;full.t=double(Q(q.l[2])/(Q(q.l[1])+q.l[2]));for(int d=0;d<2;d++){full.B[d]=vertices[1][d];full.C[d]=vertices[2][d];}
            Q s=Q(q.l[1])+q.l[2],h=s*s*s*s*reference(p.trace,full,k)[0];for(int a=0;a<10;a++)h-=shape(a,q.l)*Q(p.interpolation[a][k]);h/=Q(p.modes.scale[k]);
            mode_error=fmaxq(mode_error,fabsq(Q(v.mode[k].v)-h));
        }
    }
    need(new_error<1e-12Q&&mode_error<1e-12Q&&old_error>1000*new_error&&endpoint_residual>0,"independent binary128 geometry difference/mode accuracy and retained endpoint residual");
    uint64_t predicates=0;
    for(int j=0;j<128;j++)for(int delta=-2;delta<=2;delta++){
        double A[2]={12-j*ldexp(1.,-20),-12+j*ldexp(1.,-21)},B[2]={A[0]+.09375,A[1]+.0625},P[2]={A[0]+.046875,A[1]+.03125+delta*ldexp(1.,-49)};
        Q det=(Q(B[0])-A[0])*(Q(P[1])-A[1])-(Q(B[1])-A[1])*(Q(P[0])-A[0]);int sign=det>0?1:det<0?-1:0;
        need(compact_predicate::orientation(A,B,P)==sign&&compact_predicate::orientation(B,A,P)==-sign,"exact orientation and zero classification");predicates+=2;
    }
    cuda_need(cudaSetDevice(0));Device<cg::Patch>dp(1);dp.put(&p);Device<RelativeInput>dq(input.size());dq.put(input.data());Device<RelativeOutput>dv(input.size());
    cuda_need(cudaError_t(compact_relative_probe(dp.ptr,dq.ptr,input.size(),dv.ptr)));cuda_need(cudaDeviceSynchronize());auto gpu=dv.get();double gpu_error=0;
    for(size_t j=0;j<gpu.size();j++){
        need(gpu[j].valid,"stable GPU geometry");for(int k=0;k<2;k++){
            gpu_error=std::max(gpu_error,fabs(gpu[j].relative[k].v-cpu[j].relative[k].v)/p.modes.scale[k]);gpu_error=std::max(gpu_error,fabs(gpu[j].relative[k].x-cpu[j].relative[k].x)/p.modes.scale[k]);
            for(auto pair:{std::array<double,2>{gpu[j].mode[k].v,cpu[j].mode[k].v},std::array<double,2>{gpu[j].mode[k].x,cpu[j].mode[k].x},std::array<double,2>{gpu[j].mode[k].y,cpu[j].mode[k].y}})gpu_error=std::max(gpu_error,fabs(pair[0]-pair[1]));
        }
    }
    need(gpu_error<2e-11,"CPU/Ada tiny geometry and derivative agreement");
    std::cout<<std::setprecision(17)<<"{\"status\":\"RELATIVE_TINY_GEOMETRY_BINARY128_ADA_PREDICATE_PASS_NO_PHYSICAL_RUN\",\"samples\":"<<input.size()<<",\"predicate_checks\":"<<predicates<<",\"coordinate_scale\":12,\"mode_scale\":"<<p.modes.scale[0]<<",\"old_normalized_max_error\":"<<double(old_error)<<",\"new_normalized_max_error\":"<<double(new_error)<<",\"mode_value_error\":"<<double(mode_error)<<",\"CPU_GPU_max_normalized_difference\":"<<gpu_error<<",\"retained_endpoint_residual\":"<<double(endpoint_residual)<<",\"actual_source_executed\":false,\"new_physical_step\":false}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
