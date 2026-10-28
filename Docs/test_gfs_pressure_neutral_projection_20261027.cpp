#include "gfs_pressure_neutral_projection_20261027.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

static int require(bool ok, const char *message) {
    if (!ok) std::fprintf(stderr,"FAIL: %s\n",message);
    return ok ? 0 : 1;
}

int main() {
    int failed=0;
    const unsigned char active[4]={1,1,0,1};
    const double original[4]={1.25,-2.5,9.0,3.75};
    const double zero[4]={0.,0.,0.,0.};
    double r[4],alpha=-1.;
    std::memcpy(r,original,sizeof r);
    int status=gfs_pressure_neutral_projection(r,zero,active,4,&alpha);
    failed+=require(status==0 && alpha==0. && std::memcmp(r,original,sizeof r)==0,
                    "zero gradient leaves finite active rates bitwise unchanged");

    const double h[4]={.5,-.25,123.,.125};
    double legacy[4];
    std::memcpy(legacy,original,sizeof legacy);
    long double dot=0.,norm=0.;
    for(int i=0;i<4;++i) if(active[i]) {
        dot+=(long double)h[i]*legacy[i];
        norm+=(long double)h[i]*h[i];
    }
    double legacy_alpha=(double)(dot/norm);
    for(int i=0;i<4;++i) if(active[i]) legacy[i]-=legacy_alpha*h[i];
    std::memcpy(r,original,sizeof r);
    status=gfs_pressure_neutral_projection(r,h,active,4,&alpha);
    failed+=require(status==0 && alpha==legacy_alpha && std::memcmp(r,legacy,sizeof r)==0,
                    "nonzero-gradient projection preserves legacy arithmetic");
    long double remaining=0.;
    for(int i=0;i<4;++i) if(active[i]) remaining+=(long double)h[i]*r[i];
    failed+=require(std::fabs((double)remaining)<1e-14,
                    "nonzero-gradient projected rate is orthogonal");

    const double tiny[4]={.5e-30,-.25e-30,123.e-30,.125e-30};
    std::memcpy(r,original,sizeof r);
    status=gfs_pressure_neutral_projection(r,tiny,active,4,&alpha);
    double max_difference=0.;
    for(int i=0;i<4;++i) max_difference=std::fmax(max_difference,std::fabs(r[i]-legacy[i]));
    failed+=require(status==0 && max_difference<1e-14,
                    "small nonzero gradient is projected without arbitrary cutoff");

    const unsigned char none[4]={0,0,0,0};
    std::memcpy(r,original,sizeof r);
    status=gfs_pressure_neutral_projection(r,zero,none,4,&alpha);
    failed+=require(status==-3 && std::memcmp(r,original,sizeof r)==0,
                    "empty active mask remains an error");
    double bad[4]={0.,std::numeric_limits<double>::quiet_NaN(),0.,0.};
    std::memcpy(r,original,sizeof r);
    status=gfs_pressure_neutral_projection(r,bad,active,4,&alpha);
    failed+=require(status==-2 && std::memcmp(r,original,sizeof r)==0,
                    "nonfinite active gradient remains an error");
    const unsigned char one[4]={1,0,0,0};
    const double subnormal[4]={1e-320,0.,0.,0.};
    const double extreme[4]={1e300,0.,0.,0.};
    std::memcpy(r,extreme,sizeof r);
    status=gfs_pressure_neutral_projection(r,subnormal,one,4,&alpha);
    failed+=require(status==-4 && std::memcmp(r,extreme,sizeof r)==0,
                    "overflowing coefficient rejects without changing rates");
    if(failed) return 1;
    std::printf("PASS: zero, nonzero, tiny, inactive, nonfinite and overflow projection cases\n");
    return 0;
}
