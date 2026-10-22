// Native binary80 acceptance of the one positive-support proposal. Reads only
// saved matrices/coefficients/M; no quadrature query and no field evolution.
#include "material_full_initial_driver.cpp"
static std::string json_string(const std::string&s,const char*key){std::string t="\""+std::string(key)+"\": \"";auto at=s.find(t);need(at!=std::string::npos,"support record path");at+=t.size();auto end=s.find('"',at);need(end!=std::string::npos,"complete support path");return s.substr(at,end-at);}
int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{{
        need(argc==7&&ranks==4,"saved_initial immutable_manifest accepted_prefix prepared_k1_prefix proposal_directory NEW_output");double start=MPI_Wtime();
        full_evolution::Material m;m.saved=argv[1];m.immutable=full_evolution::textfile(argv[2]);m.gamma=double(1.666667f);m.P0=1e-6;m.eps=.047;
        m.owned=material_full::load_basis(m.saved+"/full_conditioned_basis_r"+std::to_string(rank)+".bin");md::Ownership layout(m.owned.global,m.owned.owner,MPI_COMM_WORLD);full_evolution::thermal_import(m,rank);full_evolution::checkpoint(m,layout,std::string(argv[3])+std::to_string(rank)+".bin",false);full_evolution::sync(m,layout);
        FILE*saved=fopen((std::string(argv[4])+std::to_string(rank)+".bin").c_str(),"rb");need(saved,"actual next high-rule acceleration k1");uint64_t kh[8];material_full::read(saved,kh,8);need(kh[0]==0x4d464e4558544b31ULL&&kh[1]==1&&kh[2]==uint64_t(rank)&&kh[4]==m.owned.global.size()&&kh[6]==full_evolution::hash(m.owned.f.coefficient.data(),m.owned.f.coefficient.size()*sizeof(STCoefficient))&&kh[7]==full_evolution::hash(m.immutable.data(),m.immutable.size()),"actual k1 same accepted state and immutable binding");double time;R dt,rates[3];material_full::read(saved,&time);material_full::read(saved,&dt);material_full::read(saved,rates,3);need(time==m.ledger.time,"actual saved anchor time");std::vector<STCoefficient>k1(kh[4]);material_full::read(saved,k1.data(),k1.size());need(fgetc(saved)==EOF&&fclose(saved)==0,"complete saved actual k1");
        std::ifstream proposals(std::string(argv[5])+"/progress.jsonl");need(bool(proposals),"fixed support proposal records");std::string line;uint64_t sources=0,passed=0,old_points=0,new_points=0;
        FILE*out=fopen((std::string(argv[6])+"/audit_rank"+std::to_string(rank)+".jsonl").c_str(),"wx");need(out,"new native positive rule physical contraction audit");
        while(std::getline(proposals,line)){
            if(int(full_evolution::json_number(line,"rank"))!=rank)continue;
            need(line.find("\"positive_support_proposed\": true")!=std::string::npos,"single fixed proposal exists");
            auto input=json_string(line,"input"),support=json_string(line,"support_file");FILE*f=fopen(input.c_str(),"rb");need(f,"retained native constraints");uint64_t h[8];material_full::read(f,h,8);size_t rows=h[3],n=h[4],dofs=h[6];need(h[0]==0x4d46435542434e31ULL&&h[5]==sizeof(R)&&dofs<=ST_MAX,"original full-row constraint ABI");std::vector<uint32_t>slots(n);std::vector<double>w(n);std::vector<R>A(rows*n),target(rows),absolute(rows);material_full::read(f,slots.data(),n);material_full::read(f,w.data(),n);material_full::read(f,A.data(),A.size());material_full::read(f,target.data(),rows);material_full::read(f,absolute.data(),rows);need(fgetc(f)==EOF&&fclose(f)==0,"complete unchanged constraint matrix");
            f=fopen((std::string(argv[5])+"/"+support).c_str(),"rb");need(f,"single candidate positive weights");uint64_t ph[8];material_full::read(f,ph,8);need(ph[0]==0x4d46504f53535531ULL&&ph[2]==rows&&ph[3]==n&&ph[5]==h[2]&&ph[6]==uint64_t(rank),"positive candidate source binding");size_t nn=ph[4];std::vector<uint32_t>index(nn);std::vector<double>weight(nn);material_full::read(f,index.data(),nn);material_full::read(f,weight.data(),nn);need(fgetc(f)==EOF&&fclose(f)==0,"complete candidate support");
            for(size_t j=0;j<nn;j++)need(index[j]<n&&(j==0||index[j-1]<index[j])&&weight[j]>0&&std::isfinite(weight[j]),"strictly positive finite weights, unique original nodes; no clamp");
            std::vector<R>difference(rows),sum_bound(rows);R maxentry=0;
            for(size_t r=0;r<rows;r++){R sum=0,abs_sum=0;for(size_t q=0;q<nn;q++){R z=A[r*n+index[q]]*weight[q];sum+=z;abs_sum+=fabsl(z);}difference[r]=sum-target[r];sum_bound[r]=(n+nn+2)*LDBL_EPSILON*(absolute[r]+abs_sum);if(absolute[r]>0)maxentry=std::max(maxentry,fabsl(difference[r])/absolute[r]);else need(difference[r]==0,"exact zero constraint retained");}
            auto found=std::lower_bound(m.owned.source_triangle.begin(),m.owned.source_triangle.end(),h[2]);need(found!=m.owned.source_triangle.end()&&*found==h[2],"saved same canonical triangle");size_t s=found-m.owned.source_triangle.begin();const auto&t=m.owned.f.triangle[s];need(size_t(t.n)==dofs,"same actual conditioned basis");
            R Mdiff[ST_MAX][ST_MAX]={},Cdiff[ST_MAX][ST_MAX]={},Mquad[ST_MAX][ST_MAX]={},Cbound[ST_MAX][ST_MAX]={};size_t r=1;
            for(int a=0;a<t.n;a++)for(int b=a;b<t.n;b++,r++){Mdiff[a][b]=Mdiff[b][a]=difference[r];Mquad[a][b]=Mquad[b][a]=target[r]+difference[r];}
            for(int a=0;a<t.n;a++)for(int b=a+1;b<t.n;b++,r++){Cdiff[a][b]=difference[r];Cdiff[b][a]=-difference[r];Cbound[a][b]=Cbound[b][a]=fabsl(difference[r])+sum_bound[r];}
            R U=target[r],dU=difference[r],meanP_old=target[r+1]/target[0],meanP_new=(target[r+1]+difference[r+1])/(target[0]+difference[0]);need(U>0&&target[0]>0,"original positive cold heat/measure");size_t gradient=r+2,gravity=gradient+2*t.n,potential=gravity+2*t.n;need(potential+1==rows,"every original constraint retained");
            R K=0,dK=0,dL=0,dWork=0,dMass=0,dQJ=0,dQJdot=0,Jbound=0,Gerror=0,Gscale=0,force_error=0,force_scale=0;
            R G0[2][ST_MAX]={},Gd[2][ST_MAX]={},centered_x[ST_MAX][2]={},centered_v[ST_MAX][2]={};
            const auto&origin=m.owned.f.coefficient[t.id[0]];
            for(int a=0;a<t.n;a++)for(int d=0;d<2;d++){const auto&c=m.owned.f.coefficient[t.id[a]];centered_x[a][d]=double(c.x[d]-(a<10?origin.x[d]:0));centered_v[a][d]=double(c.v[d]-(a<10?origin.v[d]:0));}
            // Original determinant moments reconstruct the same derivative span.
            size_t cstart=1+t.n*(t.n+1)/2;
            for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++){
                const auto&ca=m.owned.f.coefficient[t.id[a]];const auto&cb=m.owned.f.coefficient[t.id[b]];R dm=Mdiff[a][b],dc=Cdiff[a][b];
                for(int d=0;d<2;d++){K+=.5L*ca.v[d]*cb.v[d]*(Mquad[a][b]-dm);dK+=.5L*ca.v[d]*cb.v[d]*dm;dWork+=R(ca.v[d])*k1[t.id[b]].v[d]*dm;}
                dL+=(R(ca.x[0])*cb.v[1]-R(ca.x[1])*cb.v[0])*dm;
                if(m.owned.f.dof[t.id[a]].kind<3&&m.owned.f.dof[t.id[b]].kind<3)dMass+=dm;
                dQJ+=centered_x[a][0]*centered_x[b][1]*dc;dQJdot+=(centered_v[a][0]*centered_x[b][1]+centered_x[a][0]*centered_v[b][1])*dc;Jbound+=fabsl(centered_x[a][0]*centered_x[b][1])*Cbound[a][b];
                Gd[0][a]+=centered_x[b][1]*dc;Gd[1][b]+=centered_x[a][0]*dc;
            }
            size_t row=cstart;for(int a=0;a<t.n;a++)for(int b=a+1;b<t.n;b++,row++){R c=target[row];G0[0][a]+=centered_x[b][1]*c;G0[0][b]-=centered_x[a][1]*c;G0[1][b]+=centered_x[a][0]*c;G0[1][a]-=centered_x[b][0]*c;}
            for(int d=0;d<2;d++)for(int a=0;a<t.n;a++){
                Gerror=std::max(Gerror,fabsl(Gd[d][a]));Gscale=std::max(Gscale,fabsl(G0[d][a]));
                // Volume force after the mean projector. The source boundary
                // term needs the later full-cell/full-stage oracle comparison.
                R old=-target[gradient+d*t.n+a]-meanP_old*G0[d][a]+target[gravity+d*t.n+a];
                R change=-difference[gradient+d*t.n+a]-(meanP_new-meanP_old)*G0[d][a]-meanP_new*Gd[d][a]+difference[gravity+d*t.n+a];
                force_scale=std::max(force_scale,fabsl(old));force_error=std::max(force_error,fabsl(change));
            }
            // Read, never replace, the authoritative saved binary80 mass block.
            f=fopen((m.saved+"/full_fixed_mass_r"+std::to_string(rank)+".bin").c_str(),"rb");need(f,"authoritative fixed mass for physical contraction");uint64_t mh[4];material_full::read(f,mh,4);need(mh[0]==0x4d464d4153533031ULL&&mh[3]==sizeof(R),"saved fixed metric ABI");need(fseeko(f,off_t(sizeof(mh)+s*ST_MAX*ST_MAX*sizeof(R)),SEEK_SET)==0,"direct saved local mass block");R metric[ST_MAX*ST_MAX];material_full::read(f,metric,ST_MAX*ST_MAX);fclose(f);
            R mass_metric_K=0,new_K=0;for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++)for(int d=0;d<2;d++){R vv=R(m.owned.f.coefficient[t.id[a]].v[d])*m.owned.f.coefficient[t.id[b]].v[d];mass_metric_K+=.5L*vv*metric[a*ST_MAX+b];new_K+=.5L*vv*Mquad[a][b];}
            bool pass=nn*4<=n&&fabsl(dK)/U<1e-8L&&fabsl(dWork)/U<1e-8L&&fabsl(dL)/(2*K)<1e-11L&&fabsl(new_K-mass_metric_K)/U<1e-8L&&fabsl(dU)/U<1e-8L&&Jbound/target[0]<1e-10L&&fabsl(dQJdot)/target[0]<1e-8L&&Gscale>0&&Gerror/Gscale<1e-11L&&fabsl(meanP_new-meanP_old)/m.P0<1e-11L&&force_scale>0&&force_error/force_scale<1e-11L;
            fprintf(out,"{\"triangle\":%llu,\"original_points\":%zu,\"support_points\":%zu,\"strict_positive_original_subset\":true,\"all_constraint_rows_checked\":%zu,\"max_entry_residual_over_abs_sum\":%.21Lg,\"delta_K_over_U\":%.21Lg,\"new_point_vs_fixed_metric_K_over_U\":%.21Lg,\"delta_mass\":%.21Lg,\"delta_L_over_2K\":%.21Lg,\"delta_Kdot_over_U\":%.21Lg,\"delta_U_over_U\":%.21Lg,\"delta_QJ_over_A\":%.21Lg,\"delta_QJdot_over_A\":%.21Lg,\"conservative_delta_QJ_bound_over_A\":%.21Lg,\"delta_gradient_relative\":%.21Lg,\"delta_mean_pressure_over_P0\":%.21Lg,\"delta_volume_force_relative\":%.21Lg,\"candidate_physical_contractions_pass\":%s,\"full_cell_oracle_still_required\":true,\"rule_committed\":false}\n",(unsigned long long)h[2],n,nn,rows,maxentry,dK/U,(new_K-mass_metric_K)/U,dMass,dL/(2*K),dWork/U,dU/U,dQJ/target[0],dQJdot/target[0],Jbound/target[0],Gerror/Gscale,(meanP_new-meanP_old)/m.P0,force_error/force_scale,pass?"true":"false");fflush(out);sources++;passed+=pass;old_points+=n;new_points+=pass?nn:n;
        }
        material_full::durable(out);uint64_t totals[4]={sources,passed,old_points,new_points};MPI_Allreduce(MPI_IN_PLACE,totals,4,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);if(!rank)printf("{\"scope\":\"saved positive subset: native binary80 all-row and actual-state/M contractions\",\"representatives\":%llu,\"passed\":%llu,\"original_points\":%llu,\"proposed_points_with_failures_uncompressed\":%llu,\"seconds\":%.9g,\"full_rule_still_authoritative\":true,\"production_rule_committed\":false}\n",(unsigned long long)totals[0],(unsigned long long)totals[1],(unsigned long long)totals[2],(unsigned long long)totals[3],MPI_Wtime()-start);
    }MPI_Finalize();return 0;}catch(const std::exception&e){fprintf(stderr,"rank%d POSITIVE_RULE_AUDIT_FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}return 2;
}
