#define main archived_main
#include "source/Hydro1DExam/laguerre_sod.c"
#undef main
int main(void) {
    const char *names[]={"gfs","hllc","ns"};
    printf("[\n");
    for(int k=0;k<3;k++) {
        setup_problem(PROB_CONTACT); mode=VORONOI;
        if(k==0)set_method_gfs();else if(k==1)set_method_hllc_cd10();else set_method_ns_av();
        init(200);compute_rhs();
        double amax=0,emax=0,vmax=0;
        for(int i=0;i<N;i++) {amax=fmax(amax,fabs(P[i].av));emax=fmax(emax,fabs(P[i].ae));}
        for(int i=0;i<=N;i++)vmax=fmax(vmax,fabs(vf[i]));
        printf("%s{\"method\":\"%s\",\"max_acceleration\":%.17g,\"max_de_dt\":%.17g,\"max_work_face_speed\":%.17g,\"pair_calls\":%ld}\n",k?",":"",names[k],amax,emax,vmax,n_pair_1d);
    }
    printf("]\n");return 0;
}
