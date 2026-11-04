/* Direct fixtures for a guard rarely exercised by valid meshes. */
#include "../../Voro/voro.c"
int main(void){
    if(sizeof(postype)!=sizeof(double)){
        int rejected=Voro3D_FindVCStrict(NULL,NULL,NULL,0,NULL,0,0,0,NULL)==-1;
        printf("{\"scope\":\"reject unsupported float geometry\",\"pass\":%s}\n",rejected?"true":"false");return rejected?0:1;
    }
    int failures=0;
    for(int test=0;test<5;test++){
        Voro3D_Vertex a={0},b={0};Voro3D_point normal={0};normal.x=1;
        if(test==0){a.x=0;b.x=1;}
        if(test==1){a.x=.5+DBL_EPSILON;b.x=1;}
        if(test==2){a.x=0;b.x=.5-DBL_EPSILON;}
        if(test==3){a.x=.6;b.x=1.6;}
        if(test==4){a.x=b.x=.5;a.y=0;b.y=1;}
        voro3d_strict_geometry=1;voro3d_strict_failure=0;memset(&voro3d_strict_stats,0,sizeof(voro3d_strict_stats));
        Voro3D_Vertex out=Voro3D_FindIntersection(&a,&b,&normal,.5);
        int expected_failure=test>=3;
        int pass=voro3d_strict_failure==expected_failure;
        if(!expected_failure)pass&=fabs(out.x-.5)<=DBL_EPSILON;
        if(test==1||test==2)pass&=voro3d_strict_stats.endpoint_clamps==1&&voro3d_strict_stats.max_endpoint_residual_over_bound<=1;
        if(expected_failure)pass&=voro3d_strict_stats.rejected_intersections==1;
        printf("{\"fixture\":%d,\"failed\":%d,\"clamps\":%lld,\"endpoint_residual_over_bound\":%.17g,\"pass\":%s}\n",test,voro3d_strict_failure,voro3d_strict_stats.endpoint_clamps,voro3d_strict_stats.max_endpoint_residual_over_bound,pass?"true":"false");failures+=!pass;
    }
    return failures?1:0;
}
