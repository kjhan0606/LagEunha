/* Independent read of the original problem registry, not the new harness. */
#define main archived_original_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main
int main(void){
    const int ids[]={PROB_SOD,PROB_BLAST,PROB_COLLISION,PROB_SHUOSHER,PROB_NOH,PROB_LAX,PROB_DBLFAN,PROB_CONTACT,PROB_CONTACT};
    const char *names[]={"sod","blast","collision","shu","noh","lax","doublefan","contact_open","moving_contact_open"};
    printf("[");
    for(int k=0;k<9;k++){
        setup_problem(ids[k]);if(k==8)vL_ic=vR_ic=1;
        printf("%s{\"problem\":\"%s\",\"domain\":[%.17g,%.17g],\"gamma\":%.17g,\"endpoint\":%.17g,\"regions\":%d,\"cut\":%.17g,\"cut2\":%.17g,\"lower_reflect\":%s,\"upper_reflect\":%s,\"states\":[[%.17g,%.17g,%.17g]",k?",":"",names[k],x0_dom,x1_dom,GAM,tend,n_regions,xmid1,xmid2,bc_left==BC_REFLECT?"true":"false",bc_right==BC_REFLECT?"true":"false",rhoL_ic,vL_ic,pL_ic);
        if(n_regions==3)printf(",[%.17g,%.17g,%.17g]",rhoM_ic,vM_ic,pM_ic);
        if(n_regions>=2)printf(",[%.17g,%.17g,%.17g]",rhoR_ic,vR_ic,pR_ic);
        printf("]}");
    }printf("]\n");return 0;
}
