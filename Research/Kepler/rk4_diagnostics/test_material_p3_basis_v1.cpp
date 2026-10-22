/* Necessary P3 representation/mass preflight. No force solve or time evolution. */
#define GFS_MATERIAL_HYBRID_LIBRARY
#include "test_material_hybrid_v1.cpp"
#include "material_p3_mesh.h"
#include <chrono>
#include <sys/resource.h>
struct P3Mass {int n;std::vector<double> matrix;double mass=0,K=0,U=0;};
static void basis_case(const SharedMesh &linear,int level,const GfsEnrichConfig &cfg,int rank,int ranks,const char *dest){
    auto start=std::chrono::steady_clock::now();auto mesh=p3_build(linear,cfg);int nc=linear.cell.size();std::vector<P3Mass> cell(nc);size_t storage=0;int broken=0;
    for(int c=0;c<nc;c++){cell[c].n=mesh.cell_nodes[c].size()+2;cell[c].matrix.resize(cell[c].n*cell[c].n);storage+=cell[c].matrix.size()*sizeof(double);broken+=cell[c].n;}
    double partition=0,gradient=0,nodal=0,polynomial=0,polynomial_gradient=0,verror=0,Berror=0;size_t points=0;
    for(size_t it=0;it<linear.tri.size();it++){
        const auto &tri=linear.tri[it];if(tri.cell%ranks!=rank)continue;const auto &p3=mesh.triangle[it];auto &out=cell[tri.cell];const auto &ids=mesh.cell_nodes[tri.cell];int local[12];
        for(int a=0;a<10;a++){auto pos=std::lower_bound(ids.begin(),ids.end(),mesh.connectivity[it][a]);need(pos!=ids.end() && *pos==mesh.connectivity[it][a],"P3 local ownership");local[a]=pos-ids.begin();}local[10]=out.n-2;local[11]=out.n-1;
        double length=sqrt(2*p3.geometry.area);
        for(int a=0;a<10;a++){
            const auto &x=mesh.node[mesh.connectivity[it][a]];double N[10],Nx[10],Ny[10];gfs_material_p3_polynomial(&p3,x[0],x[1],N,Nx,Ny);for(int b=0;b<10;b++)nodal=fmax(nodal,fabs(N[b]-(a==b)));
        }
        double px[3],py[3];for(int a=0;a<3;a++){px[a]=tri.t.x[a][0];py[a]=tri.t.x[a][1];}GfsMomentPolygon poly={px,py,3,0};GfsEnrichBasis qb={};qb.generator[0]=qb.generator[1]=12;GfsEnrichQuadrature qr;need(gfs_enrich_quadrature_init(poly,&cfg,&qb,&qr),"P3 mass positive quadrature");
        for(int q=0;q<gfs_enrich_nquad(poly,&cfg,&qr);q++){
            double X,Y,w,u[4],v[2],B[4],N[12],Nx[12],Ny[12],one[2]={1,1};need(gfs_enrich_quad(poly,&cfg,&qr,q,&X,&Y,&w),"P3 quadrature point");if(w==0)continue;points++;
            velocity(cfg,X+12,Y+12,v,B);need(gfs_enrich_reference(&cfg,X+12,Y+12,u,0,0),"original P3 input rho/energy");p3_shape_stable(p3,cfg,X,Y,one,N,Nx,Ny);
            double sum=0,sx=0,sy=0,vr[2]={N[10],N[11]},Br[4]={Nx[10],Ny[10],Nx[11],Ny[11]};
            for(int a=0;a<10;a++){sum+=N[a];sx+=Nx[a];sy+=Ny[a];for(int k=0;k<2;k++){vr[k]+=N[a]*p3.v0[a][k];Br[2*k]+=Nx[a]*p3.v0[a][k];Br[2*k+1]+=Ny[a]*p3.v0[a][k];}}
            partition=fmax(partition,fabs(sum-1));gradient=fmax(gradient,length*fmax(fabs(sx),fabs(sy)));for(int k=0;k<2;k++)verror=fmax(verror,fabs(vr[k]-v[k]));for(int k=0;k<4;k++)Berror=fmax(Berror,length*fabs(Br[k]-B[k]));
            // Independent monomial reproduction in barycentric coordinates.
            double rx=X-p3.geometry.x[0][0],ry=Y-p3.geometry.x[0][1],l1=p3.geometry.grad[1][0]*rx+p3.geometry.grad[1][1]*ry,l2=p3.geometry.grad[2][0]*rx+p3.geometry.grad[2][1]*ry;
            for(int i=0;i<=3;i++)for(int j=0;j<=3-i;j++){
                double value=0,dx=0,dy=0;for(int a=0;a<10;a++){int m[3];gfs_material_p3_index(a,m);double z=pow(m[1]/3.,i)*pow(m[2]/3.,j);value+=N[a]*z;dx+=Nx[a]*z;dy+=Ny[a]*z;}
                double exact=pow(l1,i)*pow(l2,j),d1=i?i*pow(l1,i-1)*pow(l2,j):0,d2=j?j*pow(l1,i)*pow(l2,j-1):0;
                polynomial=fmax(polynomial,fabs(value-exact));polynomial_gradient=fmax(polynomial_gradient,length*fmax(fabs(dx-d1*p3.geometry.grad[1][0]-d2*p3.geometry.grad[2][0]),fabs(dy-d1*p3.geometry.grad[1][1]-d2*p3.geometry.grad[2][1])));
            }
            double dm=w*u[0];out.mass+=dm;out.K+=.5*dm*(v[0]*v[0]+v[1]*v[1]);out.U+=w*cfg.field.ic.pressure/(cfg.field.ic.gamma-1);
            for(int a=0;a<12;a++)for(int b=a;b<12;b++){double z=dm*N[a]*N[b];out.matrix[local[a]*out.n+local[b]]+=z;if(a!=b)out.matrix[local[b]*out.n+local[a]]+=z;}
        }
    }
    std::vector<double> errors={partition,gradient,nodal,polynomial,polynomial_gradient,verror,Berror};MPI_Allreduce(MPI_IN_PLACE,errors.data(),errors.size(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    uint64_t npoints=points;MPI_Allreduce(MPI_IN_PLACE,&npoints,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    for(auto &c:cell){all_sum(c.matrix);std::vector<double> v={c.mass,c.K,c.U};all_sum(v);c.mass=v[0];c.K=v[1];c.U=v[2];}
    std::vector<Pair> scale(linear.components,{0,0});std::vector<double> compmass(linear.components);std::vector<int> component(nc,-1);
    for(const auto &t:linear.tri)component[t.cell]=t.component;
    for(int c=0;c<nc;c++){int n=cell[c].n,k=component[c];compmass[k]+=cell[c].mass;for(int d=0;d<2;d++)scale[k][d]+=cell[c].matrix[(n-2+d)*n+n-2+d];}
    for(int k=0;k<linear.components;k++)for(int d=0;d<2;d++){need(scale[k][d]>0,"P3 initial-velocity residual independent norm");scale[k][d]=sqrt(scale[k][d]/compmass[k]);}
    double assemble=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if(rank==0){
        auto solve_start=std::chrono::steady_clock::now();Real condition=0,minpivot=LDBL_MAX,masserror=0,Kerror=0,minscale=LDBL_MAX;std::string path=std::string(dest)+"/level"+std::to_string(level)+"_mass.bin";FILE *file=fopen(path.c_str(),"wb");need(file,"new P3 mass packet");uint64_t header[4]={0x474650334d415331ULL,(uint64_t)nc,(uint64_t)level,(uint64_t)broken};fwrite(header,sizeof(header),1,file);
        for(int c=0;c<nc;c++){
            int n=cell[c].n,k=component[c];Vec M(n*n),v[2]={Vec(n),Vec(n)};std::vector<double> scale_dof(n,1);
            for(int d=0;d<2;d++){scale_dof[n-2+d]=scale[k][d];minscale=std::min(minscale,(Real)scale[k][d]);}
            for(int i=0;i<n;i++)for(int j=0;j<n;j++)M[i*n+j]=(Real)cell[c].matrix[i*n+j]/scale_dof[i]/scale_dof[j];Factor factor(M,n);auto eig=extrema(M,n);need(eig[0]>0,"P3 cell mass eigenvalue positive");condition=std::max(condition,eig[1]/eig[0]);minpivot=std::min(minpivot,factor.minpivot);
            for(int i=0;i<n-2;i++){const auto &x=mesh.node[mesh.cell_nodes[c][i]];double velocity0[2],B[4];velocity(cfg,x[0]+12,x[1]+12,velocity0,B);for(int d=0;d<2;d++)v[d][i]=velocity0[d];}v[0][n-2]=scale[k][0];v[1][n-1]=scale[k][1];
            Real mass=0,K=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(i<n-2 && j<n-2)mass+=M[i*n+j];K+=.5L*M[i*n+j]*(v[0][i]*v[0][j]+v[1][i]*v[1][j]);}
            masserror=std::max(masserror,fabsl(mass-cell[c].mass)/cell[c].mass);Kerror=std::max(Kerror,fabsl(K-cell[c].K)/cell[c].U);
            uint64_t row[3]={linear.cell[c].id,(uint64_t)n,(uint64_t)k};need(fwrite(row,sizeof(row),1,file)==1,"P3 cell mass header");for(auto z:M){double d=(double)z;need(fwrite(&d,8,1,file)==1,"P3 cell mass matrix");}
        }
        fclose(file);double factor_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-solve_start).count();struct rusage usage;getrusage(RUSAGE_SELF,&usage);
        printf("{\"level\":%d,\"triangles\":%zu,\"physical_cells\":%d,\"global_polynomial_nodes\":%zu,\"broken_scalar_dofs\":%d,\"points\":%llu,\"partition_error\":%.17g,\"constant_gradient_error_scaled\":%.17g,\"nodal_error\":%.17g,\"polynomial_reproduction_error\":%.17g,\"polynomial_gradient_error_scaled\":%.17g,\"initial_velocity_error\":%.17g,\"initial_gradient_error_scaled\":%.17g,\"minimum_residual_scale\":%.17Lg,\"cell_mass_condition_max\":%.17Lg,\"cell_mass_min_pivot\":%.17Lg,\"cell_mass_relative_error\":%.17Lg,\"initial_kinetic_integral_error_over_U\":%.17Lg,\"mass_matrix_storage_bytes\":%zu,\"rank0_assembly_seconds\":%.17g,\"rank0_factor_and_spectrum_seconds\":%.17g,\"rank0_maxrss_kib\":%ld}",level,linear.tri.size(),nc,mesh.node.size(),broken,(unsigned long long)npoints,errors[0],errors[1],errors[2],errors[3],errors[4],errors[5],errors[6],minscale,condition,minpivot,masserror,Kerror,storage,assemble,factor_time,usage.ru_maxrss);fflush(stdout);
    }
}
int main(int argc,char **argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==3,"native polygon packet and new output");
    std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};auto mesh=mesh_from(load(argv[1],ids));auto cfg=config(32,20);
    if(rank==0)printf("{\"scope\":\"P3 basis/initial-field/positive-mass prerequisite only; no hybrid force solve, curvature convergence, remap or time evolution claim\",\"ranks\":%d,\"levels\":[",ranks);
    basis_case(mesh,0,cfg,rank,ranks,argv[2]);auto refined=p3_refine_once(mesh);if(rank==0)printf(",\n");basis_case(refined,1,cfg,rank,ranks,argv[2]);if(rank==0)printf("]}\n");MPI_Finalize();return 0;
}
