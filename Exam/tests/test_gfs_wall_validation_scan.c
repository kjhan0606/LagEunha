#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    size_t id;
    double mass,ie,volume,w2,x,y,vx,vy,ax,ay,die;
} TestParticle;
typedef struct { double x,y,m; unsigned int present,seen; } GfsWallBase;

static int test_state_valid(const TestParticle *p,int rhs){
    return p->ie>0 && p->volume>0 && isfinite(p->mass) && p->mass>0 && p->w2==0 &&
        isfinite(p->ie) && isfinite(p->volume) && isfinite(p->x) && isfinite(p->y) &&
        isfinite(p->vx) && isfinite(p->vy) &&
        (!rhs || (isfinite(p->ax) && isfinite(p->ay) && isfinite(p->die)));
}

#define GFS_WALL_SCAN_PARTICLE_TYPE TestParticle
#define GFS_WALL_SCAN_ID(p) ((p)->id)
#define GFS_WALL_SCAN_MASS(p) ((p)->mass)
#define GFS_WALL_SCAN_STATE_VALID(p,rhs) test_state_valid((p),(rhs))
#include "../gfs_wall_validation_scan.h"

static void prepare(GfsWallBase *base,TestParticle *p,size_t n){
    for(size_t i=0;i<n;i++){
        base[i]=(GfsWallBase){0,0,1.0,1,base[i].seen};
        p[i]=(TestParticle){n-1-i,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    }
}

static int expect(const char *name,int got,int wanted){
    if(got==wanted)return 1;
    fprintf(stderr,"FAIL %s: got=%d expected=%d\n",name,got,wanted);
    return 0;
}

int main(void){
    const size_t n=100003;
    GfsWallBase *base=calloc(n,sizeof(*base));
    TestParticle *p=malloc(n*sizeof(*p));
    if(!base || !p){free(base);free(p);return 2;}

    int ok=1;
    prepare(base,p,n);
    ok &= expect("unique reversed IDs",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,1,base,p,0),1);

    prepare(base,p,n);p[1].id=p[0].id;
    ok &= expect("duplicate ID",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,2,base,p,0),0);

    prepare(base,p,n);p[0].id=n;
    ok &= expect("out-of-range ID",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,3,base,p,0),0);

    prepare(base,p,n);p[17].mass=2.0;
    ok &= expect("changed mass",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,4,base,p,0),0);

    prepare(base,p,n);p[19].ie=NAN;
    ok &= expect("invalid thermal state",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,5,base,p,0),0);

    prepare(base,p,n);p[21].ax=NAN;
    ok &= expect("invalid RHS state",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,6,base,p,1),0);

    prepare(base,p,n);p[23].w2=1.0;
    ok &= expect("nonzero W0 weight",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,7,base,p,0),0);

    prepare(base,p,n);
    ok &= expect("new generation after rejects",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,8,base,p,0),1);

    prepare(base,p,n);base[n-1].present=0;
    ok &= expect("missing accepted ID",gfs_wall_parallel_validation_scan((ptrdiff_t)n,n,9,base,p,0),0);

    free(base);free(p);
    if(!ok)return 1;
    puts("PASS: OpenMP wall-validation scan cases");
    return 0;
}
