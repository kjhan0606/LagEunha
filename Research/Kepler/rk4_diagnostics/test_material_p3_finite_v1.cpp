/* Finite-geometry RHS control, using preserved P3 acceleration. NOT evolution. */
#define GFS_MATERIAL_P3_HYBRID_LIBRARY
#include "test_material_p3_hybrid_v1.cpp"
#include "material_p3_finite_batch.h"
#include <omp.h>
static void cuda_check(int code){if(code){fprintf(stderr,"CUDA backend: %s\n",p3_finite_cuda_error(code));MPI_Abort(MPI_COMM_WORLD,code);}}
static P3Solution read_seed(const char *path,const P3Assembly &a){
    FILE *f=fopen(path,"rb");need(f,"preserved P3 seed state");uint64_t header[3];need(fread(header,sizeof(header),1,f)==1 && header[0]==0x4746503348594231ULL && header[1]==(uint64_t)a.dofs,"P3 seed identity/dimensions");P3Solution s;
    for(auto *v:{&s.x,&s.y,&s.vx,&s.vy,&s.ax,&s.ay,&s.lx,&s.ly}){size_t n=(v==&s.lx || v==&s.ly)?header[2]:header[1];std::vector<double> raw(n);need(fread(raw.data(),8,n,f)==n,"complete P3 seed arrays");v->assign(raw.begin(),raw.end());}need(fgetc(f)==EOF,"P3 seed exact length");fclose(f);
    // Only change basis units if the independent integration rule changes
    // residual normalization. Never re-evaluate an IC to replace the state.
    for(const auto &c:a.cell)for(int d=0;d<2;d++){
        int k=c.offset+c.n-2+d;Real old=d?s.vy[k]:s.vx[k];need(old>0,"preserved input residual normalization");Real ratio=a.scale[c.component][d]/old;
        for(auto *v:{&s.x,&s.y,&s.vx,&s.vy,&s.ax,&s.ay})(*v)[k]*=ratio;
    }return s;
}
static P3Solution control_state(const P3Solution &seed,const P3Assembly &a,double h,double bx,double by,bool rigid){
    P3Solution s=seed;for(int k=0;k<a.dofs;k++){s.x[k]=seed.x[k]+h*seed.vx[k]+.5L*h*h*seed.ax[k];s.y[k]=seed.y[k]+h*seed.vy[k]+.5L*h*h*seed.ay[k];s.vx[k]=seed.vx[k]+h*seed.ax[k];s.vy[k]=seed.vy[k]+h*seed.ay[k];}
    for(const auto &c:a.cell)for(int j=0;j<c.n;j++){
        int k=c.offset+j;Real one=j<c.n-2?1:0;if(rigid){s.vx[k]=.2L*one-1.7L*s.y[k];s.vy[k]=-.1L*one+1.7L*s.x[k];}else{s.vx[k]+=bx*one;s.vy[k]+=by*one;}
    }return s;
}
static std::vector<P3FinitePoint> finite_points(const P3Mesh &mesh,const P3Assembly &a){
    std::vector<P3FinitePoint> out;out.reserve(a.owned_point.size());for(const auto &p:a.owned_point){P3FinitePoint q={};q.triangle=p.triangle;q.rho0=p.rho;q.weight=p.w;int comp=mesh.linear.tri[p.triangle].component;for(int j=0;j<12;j++){double norm=j<10?1:a.scale[comp][j-10];q.N[j]=p.N[j]/norm;q.Nx[j]=p.Nx[j]/norm;q.Ny[j]=p.Ny[j]/norm;}out.push_back(q);}return out;
}
static std::vector<P3FiniteMap> finite_maps(const P3Mesh &mesh,const P3Assembly &a,const P3Solution &s,const std::vector<P3FinitePoint> &point,const std::vector<P3FiniteMap> *fixed_mass=nullptr){
    std::vector<P3FiniteMap> out=fixed_mass?*fixed_mass:std::vector<P3FiniteMap>(mesh.triangle.size());for(size_t it=0;it<out.size();it++){
        auto &map=out[it];map.map.n=12;for(int d=0;d<2;d++)map.origin[d]=mesh.triangle[it].geometry.x[0][d];int ids[12];p3_local_ids(mesh,a.cell,it,ids);const auto &c=a.cell[mesh.linear.tri[it].cell];
        for(int j=0;j<12;j++){int k=c.offset+ids[j];map.map.x[0][j]=(double)(s.x[k]-(j<10?map.origin[0]:0));map.map.x[1][j]=(double)(s.y[k]-(j<10?map.origin[1]:0));map.map.v[0][j]=(double)s.vx[k];map.map.v[1][j]=(double)s.vy[k];}
    }
    if(!fixed_mass)for(const auto &p:point)need(gfs_material12_mass_point(12,p.N,p.weight*p.rho0,out[p.triangle].mass),"finite reference triangle mass");return out;
}
struct FiniteLedger {
    double mass=0,area=0,Vdot=0,U=0,Udot=0,gpower=0,gtorque=0,gforce[2]={0},q[2]={0},p[2]={0},K=0,L=0,Kint=0,I=0,spin=0,Idot=0;
    double physical_work=0,physical_torque=0,physical_force[2]={0},external_work=0,external_torque=0,external_force[2]={0},boundary_Vdot=0;
    double matrix_mass=0,matrix_q[2]={0},matrix_p[2]={0},matrix_K=0,matrix_Kint=0,matrix_I=0,matrix_spin=0;
};
struct FiniteQuality {double minJ=DBL_MAX,minP=DBL_MAX,maxJ=0,energy_relative=0,angular_relative=0;int invalid=0;};
static void finite_volume_aggregate(const P3Mesh &mesh,P3Assembly &a,const std::vector<P3FinitePoint> &point,const std::vector<P3FiniteMap> &maps,const std::vector<P3FiniteOutput> &output,const std::vector<P3FiniteMoment> &moments,std::vector<FiniteLedger> &ledger,FiniteQuality &quality,int rank,int ranks){
    std::vector<std::vector<size_t>> bycell(a.cell.size()),triangles(a.cell.size());for(size_t i=0;i<point.size();i++)bycell[mesh.linear.tri[point[i].triangle].cell].push_back(i);for(size_t it=0;it<maps.size();it++)if(mesh.linear.tri[it].cell%ranks==rank)triangles[mesh.linear.tri[it].cell].push_back(it);
    std::vector<std::array<int,12>> ids(mesh.triangle.size());for(size_t it=0;it<ids.size();it++)p3_local_ids(mesh,a.cell,it,ids[it].data());
    for(auto &c:a.cell)for(auto &v:c.force[2])std::fill(v.begin(),v.end(),0);
    for(size_t i=0;i<point.size();i++){
        const auto &o=output[i];if(!o.valid){quality.invalid++;continue;}quality.minJ=std::min(quality.minJ,o.point.J);quality.maxJ=std::max(quality.maxJ,o.point.J);quality.minP=std::min(quality.minP,o.pressure);quality.energy_relative=std::max(quality.energy_relative,fabs(o.energy_identity)/(o.energy_scale+1.e-300));quality.angular_relative=std::max(quality.angular_relative,fabs(o.angular_identity)/(o.angular_scale+1.e-300));
    }
    if(quality.invalid){fprintf(stderr,"rank%d invalid finite volume points:%d\n",rank,quality.invalid);MPI_Abort(MPI_COMM_WORLD,61);}
    #pragma omp parallel for schedule(static)
    for(int ic=0;ic<(int)a.cell.size();ic++)if(ic%ranks==rank){
        auto &z=ledger[ic];auto &c=a.cell[ic];for(size_t i:bycell[ic]){
            const auto &q=point[i];const auto &o=output[i];if(!o.valid)continue;const auto &map=maps[q.triangle];double X=o.point.position[0]+map.origin[0],Y=o.point.position[1]+map.origin[1],vx=o.point.velocity[0],vy=o.point.velocity[1],dm=q.weight*q.rho0;
            for(int j=0;j<12;j++){c.force[2][0][ids[q.triangle][j]]+=o.force[0][j];c.force[2][1][ids[q.triangle][j]]+=o.force[1][j];}
            z.mass+=dm;z.area+=q.weight*o.point.J;z.Vdot+=q.weight*o.point.J*(o.point.grad_v[0]+o.point.grad_v[3]);z.U+=o.thermal;z.Udot+=o.Udot;z.gpower+=o.gravity_power;z.gtorque+=o.gravity_torque;for(int d=0;d<2;d++)z.gforce[d]+=o.gravity_force[d];z.q[0]+=dm*X;z.q[1]+=dm*Y;z.p[0]+=dm*vx;z.p[1]+=dm*vy;z.K+=.5*dm*(vx*vx+vy*vy);z.L+=dm*(X*vy-Y*vx);
        }
        double center[2]={z.q[0]/z.mass,z.q[1]/z.mass},vcenter[2]={z.p[0]/z.mass,z.p[1]/z.mass};
        for(size_t i:bycell[ic]){
            const auto &q=point[i];const auto &o=output[i];if(!o.valid)continue;const auto &map=maps[q.triangle];double dx=o.point.position[0]+map.origin[0]-center[0],dy=o.point.position[1]+map.origin[1]-center[1],dvx=o.point.velocity[0]-vcenter[0],dvy=o.point.velocity[1]-vcenter[1],dm=q.weight*q.rho0;
            z.Kint+=.5*dm*(dvx*dvx+dvy*dvy);z.I+=dm*(dx*dx+dy*dy);z.spin+=dm*(dx*dvy-dy*dvx);z.Idot+=2*dm*(dx*dvx+dy*dvy);
        }
        for(size_t it:triangles[ic]){
            const auto &o=moments[it];if(!o.valid)continue;const auto &m=o.moment;double x=m.center[0]+maps[it].origin[0],y=m.center[1]+maps[it].origin[1];z.matrix_mass+=m.mass;z.matrix_q[0]+=m.mass*x;z.matrix_q[1]+=m.mass*y;z.matrix_p[0]+=m.momentum[0];z.matrix_p[1]+=m.momentum[1];z.matrix_K+=m.kinetic;
        }
        double cm[2]={z.matrix_q[0]/z.matrix_mass,z.matrix_q[1]/z.matrix_mass},vm[2]={z.matrix_p[0]/z.matrix_mass,z.matrix_p[1]/z.matrix_mass};
        for(size_t it:triangles[ic]){
            const auto &o=moments[it];if(!o.valid)continue;const auto &m=o.moment;double dx=m.center[0]+maps[it].origin[0]-cm[0],dy=m.center[1]+maps[it].origin[1]-cm[1],dvx=m.momentum[0]/m.mass-vm[0],dvy=m.momentum[1]/m.mass-vm[1];z.matrix_Kint+=m.internal_kinetic+.5*m.mass*(dvx*dvx+dvy*dvy);z.matrix_I+=m.polar_inertia+m.mass*(dx*dx+dy*dy);z.matrix_spin+=m.intrinsic_angular+m.mass*(dx*dvy-dy*dvx);
        }
    }
    for(size_t it=0;it<maps.size();it++)if(mesh.linear.tri[it].cell%ranks==rank && !moments[it].valid)quality.invalid++;
}
static void finite_boundaries(const P3Mesh &mesh,const GfsEnrichConfig &cfg,const P3Assembly &a,const std::vector<P3FiniteMap> &maps,std::vector<FiniteLedger> &ledger,std::vector<std::array<Vec,2>> &force,int rank,int ranks){
    for(const auto &edge:mesh.linear.edge){
        if(edge.second.size()==2 && mesh.linear.tri[edge.second[0]].cell==mesh.linear.tri[edge.second[1]].cell)continue;
        for(int it:edge.second){
            const auto &t=mesh.linear.tri[it];if(t.cell%ranks!=rank)continue;auto &z=ledger[t.cell];int va=-1,vb=-1;for(int j=0;j<3;j++){int x=t.vertex[j],y=t.vertex[(j+1)%3];if(std::min(x,y)==edge.first.first && std::max(x,y)==edge.first.second){va=x;vb=y;}}
            need(va>=0,"finite physical boundary orientation");auto x=mesh.linear.node[va],y=mesh.linear.node[vb];double len=hypot(y[0]-x[0],y[1]-x[1]),nx=(y[1]-x[1])/len,ny=-(y[0]-x[0])/len;int ids[12];p3_local_ids(mesh,a.cell,it,ids);
            for(int q=0;q<cfg.na;q++){
                double r=cfg.an[q],w=len*cfg.aw[q],X=(1-r)*x[0]+r*y[0],Y=(1-r)*x[1]+r*y[1],N[12],Nx[12],Ny[12];p3_shape_stable(mesh.triangle[it],cfg,X,Y,a.scale[t.component].data(),N,Nx,Ny);GfsMaterialPoint12 p;need(gfs_material12_point(&maps[it].map,N,Nx,Ny,&p),"finite boundary positive map");double P=cfg.field.ic.pressure*pow(p.J,-cfg.field.ic.gamma),sigma[3]={-P,0,-P},traction[2],ratio;need(gfs_material12_piola_traction(&p,sigma,nx,ny,traction,&ratio),"finite physical Piola traction");double px=p.position[0]+maps[it].origin[0],py=p.position[1]+maps[it].origin[1],work=w*(p.velocity[0]*traction[0]+p.velocity[1]*traction[1]),torque=w*(px*traction[1]-py*traction[0]);
                z.physical_work+=work;z.physical_torque+=torque;for(int d=0;d<2;d++)z.physical_force[d]+=w*traction[d];double cnx=p.F[3]*nx-p.F[2]*ny,cny=-p.F[1]*nx+p.F[0]*ny;z.boundary_Vdot+=w*(p.velocity[0]*cnx+p.velocity[1]*cny);
                if(edge.second.size()==1){z.external_work+=work;z.external_torque+=torque;for(int d=0;d<2;d++)z.external_force[d]+=w*traction[d];for(int j=0;j<12;j++){force[t.cell][0][ids[j]]+=w*N[j]*traction[0];force[t.cell][1][ids[j]]+=w*N[j]*traction[1];}}
            }
        }
    }
}
static void export_control(const P3Mesh &mesh,const GfsEnrichConfig &cfg,const P3Assembly &a,const std::vector<P3FiniteMap> &maps,const P3Solution &s,const char *dest,const char *name,double h){
    std::string state_path=std::string(dest)+"/"+name+"_state.bin";FILE *state=fopen(state_path.c_str(),"wb");need(state,"new finite-control state");uint64_t hs[3]={0x4746503346495331ULL,(uint64_t)a.dofs,(uint64_t)s.lx.size()};fwrite(hs,sizeof(hs),1,state);
    for(const auto *v:{&s.x,&s.y,&s.vx,&s.vy,&s.ax,&s.ay,&s.lx,&s.ly})for(auto value:*v){double z=(double)value;need(fwrite(&z,8,1,state)==1,"finite-control state arrays");}fclose(state);
    if(std::string(name)!="curved_Kepler")return;
    std::string path=std::string(dest)+"/curved_map_control.bin";FILE *f=fopen(path.c_str(),"wb");need(f,"new curved-map control packet");uint64_t header[8]={0x4746503346494e31ULL,1,(uint64_t)maps.size(),(uint64_t)a.cell.size(),(uint64_t)a.dofs,12,10,864};need(fwrite(header,sizeof(header),1,f)==1,"curved control header");
    double parameters[12]={h,cfg.field.ic.gamma,cfg.field.ic.pressure,cfg.field.ic.eps,cfg.field.ic.rin,cfg.field.ic.rout,.001,1.e-6,.25,2.,0,0};fwrite(parameters,sizeof(parameters),1,f);
    for(size_t it=0;it<maps.size();it++){
        const auto &t=mesh.triangle[it];int cell=mesh.linear.tri[it].cell,component=mesh.linear.tri[it].component;uint64_t tag[3]={mesh.linear.cell[cell].id,(uint64_t)component,12};fwrite(tag,sizeof(tag),1,f);
        double geometry[13];for(int j=0;j<3;j++)for(int d=0;d<2;d++){geometry[2*j+d]=t.geometry.x[j][d];geometry[6+2*j+d]=t.geometry.grad[j][d];}geometry[12]=t.geometry.area;fwrite(geometry,sizeof(geometry),1,f);fwrite(t.v0,sizeof(t.v0),1,f);fwrite(t.remainder,sizeof(t.remainder),1,f);fwrite(maps[it].map.x,sizeof(maps[it].map.x),1,f);fwrite(maps[it].map.v,sizeof(maps[it].map.v),1,f);fwrite(maps[it].origin,sizeof(maps[it].origin),1,f);fwrite(a.scale[component].data(),sizeof(double),2,f);
    }fclose(f);
}
static void finite_report(const P3Mesh &mesh,const GfsEnrichConfig &cfg,const P3Assembly &a,const P3System &sys,const P3Solution &input,const P3Solution &seed,const std::vector<P3FiniteMap> &maps,const std::vector<FiniteLedger> &ledger,const FiniteQuality &quality,const char *name,const char *dest,double h,double eval_seconds,double kernel_seconds,double aggregate_seconds,bool rigid){
    auto begin=std::chrono::steady_clock::now();auto s=p3_solve(mesh,cfg,a,sys,2,0,0,&input);auto faces=p3_measure_faces(mesh,cfg,a,sys,s,dest,name,1);Real internal_jump=p3_internal_jump(mesh,cfg,a,s);
    Real maxE=0,maxL=0,maxspin=0,totalE=0,totalL=0,totalU=0,totalK=0,mass_error=0,moment_error=0,K_error=0,volume_error=0,max_heat=0,rigid_rest=0,seed_accel_difference=0;int nc=a.cell.size();
    printf("{\"case\":\"%s\",\"control_parameter_h\":%.17g,\"rigid_control\":%s,\"invalid_points\":%d,\"sampled_min_J\":%.17g,\"sampled_max_J\":%.17g,\"sampled_min_P\":%.17g,\"point_energy_identity_relative\":%.17g,\"point_angular_identity_relative\":%.17g,\"momentum_equation_relative\":%.17Lg,\"constraint_max\":%.17Lg,\"direct_reference_traction_load_relative\":%.17Lg,\"shared_position_jump\":%.17Lg,\"shared_velocity_jump\":%.17Lg,\"shared_acceleration_jump\":%.17Lg,\"internal_triangle_acceleration_jump\":%.17Lg,\"cells\":[",name,h,rigid?"true":"false",quality.invalid,quality.minJ,quality.maxJ,quality.minP,quality.energy_relative,quality.angular_relative,s.equation/s.force_scale,s.constraint,faces.load_error/s.force_scale,faces.position_jump,faces.velocity_jump,faces.acceleration_jump,internal_jump);
    for(int ic=0;ic<nc;ic++){
        const auto &c=a.cell[ic];const auto &z=ledger[ic];Real K=0,Kdot=0,Ldot=0,pdot[2]={0},numerical_force[2]={0};
        for(int i=0;i<c.n;i++)for(int j=0;j<c.n;j++){
            int I=c.offset+i,J=c.offset+j;Real M=c.M[i*c.n+j];K+=.5L*M*(s.vx[I]*s.vx[J]+s.vy[I]*s.vy[J]);Kdot+=M*(s.vx[I]*s.ax[J]+s.vy[I]*s.ay[J]);Ldot+=M*(s.x[I]*s.ay[J]-s.y[I]*s.ax[J]);if(i<c.n-2){pdot[0]+=M*s.ax[J];pdot[1]+=M*s.ay[J];}
        }
        for(int i=0;i<c.n-2;i++){numerical_force[0]+=faces.directx[c.offset+i];numerical_force[1]+=faces.directy[c.offset+i];}
        Real center[2]={z.q[0]/z.mass,z.q[1]/z.mass},velocity[2]={z.p[0]/z.mass,z.p[1]/z.mass};Real E=Kdot+z.Udot-faces.work[ic]-z.external_work-z.gpower,L=Ldot-faces.torque[ic]-z.external_torque-z.gtorque;
        Real intrinsic_rate=Ldot-center[0]*pdot[1]+center[1]*pdot[0],forcex=numerical_force[0]+z.external_force[0]+z.gforce[0],forcey=numerical_force[1]+z.external_force[1]+z.gforce[1];Real intrinsic_torque=faces.torque[ic]+z.external_torque+z.gtorque-center[0]*forcey+center[1]*forcex;
        Real omega=z.spin/z.I,Krot=.5L*z.spin*z.spin/z.I,Krest=z.Kint-Krot,Krot_dot=omega*intrinsic_rate-.5L*omega*omega*z.Idot,Kint_dot=Kdot-velocity[0]*pdot[0]-velocity[1]*pdot[1];
        Real mr=fabsl(z.matrix_mass-z.mass)/z.mass,me=std::max({fabsl(z.matrix_Kint-z.Kint)/z.Kint,fabsl(z.matrix_I-z.I)/z.I,fabsl(z.matrix_spin-z.spin)/sqrtl(2*z.I*z.Kint)});mass_error=std::max(mass_error,mr);moment_error=std::max(moment_error,me);K_error=std::max(K_error,fabsl(K-z.K)/z.U);volume_error=std::max(volume_error,(Real)fabs(z.Vdot-z.boundary_Vdot)/z.area);max_heat=std::max(max_heat,(Real)fabs(z.Udot)/z.U);if(rigid)rigid_rest=std::max(rigid_rest,fabsl(Krest)/z.Kint);
        maxE=std::max(maxE,fabsl(E)/z.U);maxL=std::max(maxL,fabsl(L)/(2*K));maxspin=std::max(maxspin,fabsl(intrinsic_rate-intrinsic_torque)/(2*K));totalE+=E;totalL+=L;totalU+=z.U;totalK+=K;
        if(h==0)for(int j=0;j<c.n;j++){seed_accel_difference=std::max({seed_accel_difference,fabsl(s.ax[c.offset+j]-seed.ax[c.offset+j]),fabsl(s.ay[c.offset+j]-seed.ay[c.offset+j])});}
        if(ic)printf(",");printf("{\"id\":%llu,\"mass\":%.17g,\"area\":%.17g,\"mean_density\":%.17g,\"U\":%.17g,\"K\":%.17Lg,\"energy_rate_residual_over_U\":%.17Lg,\"angular_rate_residual_over_2K\":%.17Lg,\"spin_rate_residual_over_2K\":%.17Lg,\"physical_pressure_energy_residual_over_U\":%.17Lg,\"intrinsic_angular\":%.17g,\"intrinsic_angular_rate\":%.17Lg,\"polar_inertia\":%.17g,\"internal_kinetic\":%.17g,\"rigid_rotation_kinetic\":%.17Lg,\"resolved_remainder\":%.17Lg,\"rotation_kinetic_rate\":%.17Lg,\"remainder_kinetic_rate\":%.17Lg,\"physical_thermal_rate\":%.17g}",(unsigned long long)mesh.linear.cell[ic].id,z.mass,z.area,z.mass/z.area,z.U,K,E/z.U,L/(2*K),(intrinsic_rate-intrinsic_torque)/(2*K),(Kdot+z.Udot-z.physical_work-z.gpower)/z.U,z.spin,intrinsic_rate,z.I,z.Kint,Krot,Krest,Krot_dot,Kint_dot-Krot_dot,z.Udot);
    }
    double solve_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    printf("],\"max_energy_rate_residual_over_U\":%.17Lg,\"max_angular_rate_residual_over_2K\":%.17Lg,\"max_spin_rate_residual_over_2K\":%.17Lg,\"global_energy_rate_residual_over_U\":%.17Lg,\"global_angular_rate_residual_over_2K\":%.17Lg,\"mass_moment_relative_error\":%.17Lg,\"centered_matrix_moment_relative_error\":%.17Lg,\"mass_matrix_kinetic_error_over_U\":%.17Lg,\"volume_rate_identity_over_area\":%.17Lg,\"thermal_rate_over_U_max\":%.17Lg,\"rigid_remainder_over_internalK_max\":%.17Lg,\"initial_acceleration_vs_preserved_seed_max\":%.17Lg,\"backend_eval_seconds_max\":%.17g,\"CUDA_kernel_seconds_max\":%.17g,\"aggregate_seconds_max\":%.17g,\"reference_solve_and_ledgers_seconds\":%.17g}",maxE,maxL,maxspin,totalE/totalU,totalL/(2*totalK),mass_error,moment_error,K_error,volume_error,max_heat,rigid_rest,seed_accel_difference,eval_seconds,kernel_seconds,aggregate_seconds,solve_time);fflush(stdout);export_control(mesh,cfg,a,maps,s,dest,name,h);
    std::string path=std::string(dest)+"/"+name+"_force.bin";FILE *f=fopen(path.c_str(),"wb");need(f,"finite force export");uint64_t header[2]={0x47465033464f5231ULL,(uint64_t)a.dofs};fwrite(header,sizeof(header),1,f);for(const auto &c:a.cell)for(int d=0;d<2;d++)need(fwrite(c.force[2][d].data(),8,c.n,f)==(size_t)c.n,"finite force vector");fclose(f);
}
static void finite_case(const P3Mesh &mesh,const GfsEnrichConfig &cfg,P3Assembly &a,const P3System *system,const P3Solution &seed,const std::vector<P3FinitePoint> &point,const std::vector<P3FiniteMap> &fixed_mass,void *gpu,bool use_gpu,int test,int rank,int ranks,const char *dest){
    const char *names[6]={"initial_Kepler","curved_Kepler","curved_boost","curved_gravity_off","curved_linear_gravity","curved_rigid_rotation"};double h=test==0?0:2.e-4,bx=test==2?.013:0,by=test==2?-.007:0;bool rigid=test==5;
    auto input=control_state(seed,a,h,bx,by,rigid);auto maps=finite_maps(mesh,a,input,point,&fixed_mass);P3FinitePhysics physics={cfg.field.ic.gamma,cfg.field.ic.pressure,cfg.field.ic.eps,test==3 || rigid?0:test==4?1:2};std::vector<P3FiniteOutput> output(point.size());std::vector<P3FiniteMoment> moments(maps.size());double kernel_time=0,eval_time=0;
    auto begin=std::chrono::steady_clock::now();
    if(use_gpu)cuda_check(p3_finite_cuda_eval(gpu,maps.data(),&physics,output.data(),moments.data(),&kernel_time,&eval_time));
    else{
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<point.size();i++)output[i]=p3_finite_evaluate(&maps[point[i].triangle],&point[i],&physics);
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<maps.size();i++)if(maps[i].mass[0][0]>0){double constant[12];for(int j=0;j<12;j++)constant[j]=j<10?1:0;moments[i].valid=gfs_material12_moments(&maps[i].map,maps[i].mass,constant,&moments[i].moment);}
        eval_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    }
    auto aggregate_start=std::chrono::steady_clock::now();std::vector<FiniteLedger> ledger(a.cell.size());FiniteQuality quality;finite_volume_aggregate(mesh,a,point,maps,output,moments,ledger,quality,rank,ranks);
    std::vector<std::array<Vec,2>> exterior(a.cell.size());for(size_t c=0;c<a.cell.size();c++)for(auto &f:exterior[c])f.resize(a.cell[c].n);finite_boundaries(mesh,cfg,a,maps,ledger,exterior,rank,ranks);
    for(size_t c=0;c<a.cell.size();c++)for(int d=0;d<2;d++){for(int j=0;j<a.cell[c].n;j++)a.cell[c].force[2][d][j]+=(double)exterior[c][d][j];all_sum(a.cell[c].force[2][d]);}
    static_assert(sizeof(FiniteLedger)%sizeof(double)==0,"homogeneous double ledger");MPI_Allreduce(MPI_IN_PLACE,ledger.data(),ledger.size()*sizeof(FiniteLedger)/sizeof(double),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    double minima[2]={quality.minJ,quality.minP},maxima[3]={quality.maxJ,quality.energy_relative,quality.angular_relative};MPI_Allreduce(MPI_IN_PLACE,minima,2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,maxima,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);quality.minJ=minima[0];quality.minP=minima[1];quality.maxJ=maxima[0];quality.energy_relative=maxima[1];quality.angular_relative=maxima[2];MPI_Allreduce(MPI_IN_PLACE,&quality.invalid,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);need(quality.invalid==0,"finite point/moment validity");
    double timings[3]={eval_time,kernel_time,std::chrono::duration<double>(std::chrono::steady_clock::now()-aggregate_start).count()};MPI_Allreduce(MPI_IN_PLACE,timings,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(rank==0)finite_report(mesh,cfg,a,*system,input,seed,maps,ledger,quality,names[test],dest,h,timings[0],timings[1],timings[2],rigid);
}
int main(int argc,char **argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);need(provided>=MPI_THREAD_FUNNELED,"MPI main-thread communication with OpenMP workers");int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==5 || argc==6,"polygon packet, seed state, new output, cpu|cuda, optional fine_quadrature");bool use_gpu=std::string(argv[4])=="cuda";need(use_gpu || std::string(argv[4])=="cpu","known native finite backend");bool fine=argc==6;if(fine)need(std::string(argv[5])=="fine_quadrature","known finite quadrature option");
    std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};auto linear=p3_refine_once(mesh_from(load(argv[1],ids)));auto cfg=config(fine?64:32,fine?28:20);auto mesh=p3_build(linear,cfg);auto assembly=p3_assemble(mesh,cfg,rank,ranks);auto seed=read_seed(argv[2],assembly);auto point=finite_points(mesh,assembly);auto fixed_mass=finite_maps(mesh,assembly,seed,point);P3System *system=nullptr;if(rank==0)system=new P3System(p3_system(mesh,cfg,assembly));
    void *gpu=nullptr;double setup=0;if(use_gpu)cuda_check(p3_finite_cuda_create(point.data(),point.size(),fixed_mass.size(),&gpu,&setup));MPI_Allreduce(MPI_IN_PLACE,&setup,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);uint64_t points=point.size();MPI_Allreduce(MPI_IN_PLACE,&points,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    if(rank==0)printf("{\"scope\":\"finite material-map RHS control from preserved P3 acceleration; adiabatic geometry-control thermal state; no timestep, restart, remap or long-evolution claim\",\"backend\":\"%s\",\"ranks\":%d,\"OMP_threads\":%d,\"quadrature\":[%d,%d],\"positive_volume_points\":%llu,\"triangles\":%zu,\"physical_cells\":%zu,\"scalar_dofs\":%d,\"conforming_scalar_dofs\":%d,\"CUDA_setup_and_static_upload_seconds_max\":%.17g,\"cases\":[",argv[4],ranks,omp_get_max_threads(),cfg.na,cfg.nr,(unsigned long long)points,mesh.triangle.size(),assembly.cell.size(),assembly.dofs,system->scalar_nullity,setup);
    for(int test=0;test<6;test++){if(rank==0 && test)printf(",\n");finite_case(mesh,cfg,assembly,system,seed,point,fixed_mass,gpu,use_gpu,test,rank,ranks,argv[3]);}
    if(use_gpu)cuda_check(p3_finite_cuda_destroy(gpu));delete system;if(rank==0)printf("]}\n");MPI_Finalize();return 0;
}
