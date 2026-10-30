#include "../gfs_boundary_kind.h"
#include <assert.h>

int main(void)
{
    assert(gfs_particle_kind_from_flags(0, 0) == GFS_PARTICLE_FLUID);
    assert(gfs_particle_kind_from_flags(0, 1) == GFS_PARTICLE_FLUID);
    assert(gfs_particle_kind_from_flags(1, 0) == GFS_PARTICLE_MIRROR_GHOST);
    assert(gfs_particle_kind_from_flags(1, 1) == GFS_PARTICLE_VACUUM_GENERATOR);

    assert(gfs_face_kind_from_particle_kind(GFS_PARTICLE_FLUID)
           == GFS_FACE_INTERIOR);
    assert(gfs_face_kind_from_particle_kind(GFS_PARTICLE_MIRROR_GHOST)
           == GFS_FACE_MIRROR_GHOST);
    assert(gfs_face_kind_from_particle_kind(GFS_PARTICLE_VACUUM_GENERATOR)
           == GFS_FACE_VACUUM);
    assert(GFS_FACE_FIXED_WALL != GFS_FACE_VACUUM);
    return 0;
}
