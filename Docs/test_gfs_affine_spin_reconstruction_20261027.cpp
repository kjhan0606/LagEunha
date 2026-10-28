#include "gfs_affine_spin_reconstruction_20261027.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>

static bool near(double x,double y) {return std::fabs(x-y)<1e-12*(1.+std::fabs(y));}

int main() {
    gfs_affine_spin_input x{2.,.3,1.,.2,-1.1,.4,.6,0.};
    gfs_affine_spin_output ref;
    if(!gfs_affine_spin_reconstruct(x,&ref))return 1;
    if(!near(ref.spin_energy,0.))return 2;
    const double nonrigid_spin=ref.ayx*x.qxx-ref.axy*x.qyy+(ref.ayy-ref.axx)*x.qxy;
    if(!near(nonrigid_spin,0.))return 3;
    x.spin=ref.reference_spin;
    gfs_affine_spin_output same;
    if(!gfs_affine_spin_reconstruct(x,&same) || !near(same.axx,x.axx) ||
       !near(same.axy,x.axy) || !near(same.ayx,x.ayx) || !near(same.ayy,x.ayy))return 4;
    const double full=gfs_affine_energy(same.axx,same.axy,same.ayx,same.ayy,x.qxx,x.qxy,x.qyy);
    if(!near(full,same.nonrigid_energy+same.spin_energy))return 5;
    x.spin=3.1;
    gfs_affine_spin_output changed;
    if(!gfs_affine_spin_reconstruct(x,&changed))return 6;
    const double actual_spin=changed.ayx*x.qxx-changed.axy*x.qyy+
                             (changed.ayy-changed.axx)*x.qxy;
    const double actual_energy=gfs_affine_energy(changed.axx,changed.axy,changed.ayx,changed.ayy,
                                                  x.qxx,x.qxy,x.qyy);
    if(!near(actual_spin,x.spin) || !near(actual_energy,changed.nonrigid_energy+changed.spin_energy))return 7;
    x.spin=0.;
    gfs_affine_spin_output zero_excess;
    if(!gfs_affine_spin_reconstruct_excess(x,&zero_excess) ||
       zero_excess.axx!=x.axx || zero_excess.axy!=x.axy ||
       zero_excess.ayx!=x.ayx || zero_excess.ayy!=x.ayy)return 11;
    x.spin=.7;
    gfs_affine_spin_output excess;
    if(!gfs_affine_spin_reconstruct_excess(x,&excess))return 12;
    const double total_spin=excess.ayx*x.qxx-excess.axy*x.qyy+(excess.ayy-excess.axx)*x.qxy;
    if(!near(total_spin,excess.reference_spin+x.spin))return 13;
    gfs_affine_spin_output original=changed;
    x.qxy=2.;
    if(gfs_affine_spin_reconstruct(x,&changed) || std::memcmp(&original,&changed,sizeof(changed)))return 8;
    std::mt19937_64 rng(20261027);
    std::uniform_real_distribution<double> scale(.05,2.), value(-2.,2.);
    for(int i=0;i<100000;i++) {
        const double l00=scale(rng),l10=value(rng),l11=scale(rng);
        gfs_affine_spin_input q{l00*l00,l00*l10,l10*l10+l11*l11,
                                value(rng),value(rng),value(rng),value(rng),value(rng)};
        gfs_affine_spin_output o;
        if(!gfs_affine_spin_reconstruct(q,&o))return 9;
        const double spin=o.ayx*q.qxx-o.axy*q.qyy+(o.ayy-o.axx)*q.qxy;
        const double kinetic=gfs_affine_energy(o.axx,o.axy,o.ayx,o.ayy,q.qxx,q.qxy,q.qyy);
        if(!near(spin,q.spin) || !near(kinetic,o.nonrigid_energy+o.spin_energy))return 10;
    }
    std::puts("PASS: anisotropic projection, 100000 moment/energy invariants and invalid rollback");
}
