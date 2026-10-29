#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Exam/gfs_moment_rk4_mpi.h"
#include "Exam/gfs_weight_mass_exchange.h"
static void ok(int c,const char *msg){if(!c){fprintf(stderr,"%s\n",msg);MPI_Abort(MPI_COMM_WORLD,1);}}
static void near(double x,double y){ok(isfinite(x)&&fabs(x-y)<1e-12*(1+fabs(y)),"numeric mismatch");}
int main(int argc,char **argv){
    MPI_Init(&argc,&argv);int rank,nrank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    ok(nrank==2,"expected two ranks");
    GfsMomentState accepted;gfs_moment_zero(&accepted);accepted.a[GFM_M]=rank?7:5;
    accepted.a[GFM_QXX]=accepted.a[GFM_QYY]=accepted.a[GFM_M]/12.;accepted.a[GFM_E]=10*accepted.a[GFM_M];
    GfsMomentState original=accepted,da,db,rhs;double flux[4];
    gfs_moment_zero(&da);gfs_moment_zero(&db);
    ok(gfs_weight_exchange_flux(.01,-.02,1,1,.1,-.2,1.5,flux),"flux");
    ok(gfs_weight_exchange_segment(&da,&db,.5,-.2,.5,.2,0,0,1,0,flux),"face");
    rhs=rank?db:da;GfsMomentRK4 scratch={0};
    ok(gfs_moment_rk4_begin(&scratch,&accepted,.01,0),"begin");
    for(int k=0;k<4;k++)ok(gfs_moment_rk4_add_rhs(&scratch,&rhs),"staging");
    ok(gfs_moment_rk4_collective_commit(MPI_COMM_WORLD,&scratch,&accepted,1,rank==0)==0,"expected global rejection");
    ok(!memcmp(&accepted,&original,sizeof(accepted)),"rejected state changed");
    if(rank==1)scratch.dt=.02;
    ok(gfs_moment_rk4_collective_commit(MPI_COMM_WORLD,&scratch,&accepted,1,1)==0,"expected timestep mismatch rejection");
    ok(!memcmp(&accepted,&original,sizeof(accepted)),"timestep rejection modified state");
    scratch.dt=.01;
    ok(gfs_moment_rk4_collective_commit(MPI_COMM_WORLD,&scratch,&accepted,1,1)==1,"expected commit");
    GfsMomentState candidate;ok(gfs_moment_rk4_candidate(&scratch,&candidate),"candidate");
    ok(!memcmp(&accepted,&candidate,sizeof(accepted)),"candidate not committed");
    double local[5]={accepted.a[GFM_M],accepted.a[GFM_PX],accepted.a[GFM_PY],accepted.a[GFM_E],gfs_moment_Lz(&accepted,rank,0)};
    double total[5];MPI_Allreduce(local,total,5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    near(total[0],12);near(total[1],0);near(total[2],0);near(total[3],120);near(total[4],0);
    /* A stale base on one rank must prevent partial writes everywhere. */
    GfsMomentState frozen=accepted;
    ok(gfs_moment_rk4_collective_commit(MPI_COMM_WORLD,&scratch,&accepted,1,1)==0,"expected stale-base rejection");
    ok(!memcmp(&accepted,&frozen,sizeof(accepted)),"stale rejection modified state");
    /* Empty ranks still join the vote and do not change the timestep range. */
    GfsMomentRK4 empty_check={0};GfsMomentState zero_rhs;gfs_moment_zero(&zero_rhs);
    if(rank==0){
        ok(gfs_moment_rk4_begin(&empty_check,&accepted,.01,1),"empty-rank begin");
        for(int k=0;k<4;k++)ok(gfs_moment_rk4_add_rhs(&empty_check,&zero_rhs),"empty-rank stage");
    }
    ok(gfs_moment_rk4_collective_commit(MPI_COMM_WORLD,rank?NULL:&empty_check,rank?NULL:&accepted,rank?0:1,1)==1,"empty-rank collective");
    if(rank==0)puts("MPI two-rank reject/commit/stale-base and conserved pair checks pass");
    MPI_Finalize();return 0;
}
