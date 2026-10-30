#ifndef GFS_BOUNDARY_KIND_H
#define GFS_BOUNDARY_KIND_H

#if defined(__CUDACC__)
#define GFS_BOUNDARY_HD __host__ __device__
#else
#define GFS_BOUNDARY_HD
#endif

/* Particle kinds carried through the host/device SoA. */
enum GfsParticleKind {
    GFS_PARTICLE_FLUID = 0,
    GFS_PARTICLE_MIRROR_GHOST = 1,
    GFS_PARTICLE_VACUUM_GENERATOR = 2
};

/* Face kinds carried through the host/device CSR. */
enum GfsFaceKind {
    GFS_FACE_INTERIOR = 0,
    GFS_FACE_MIRROR_GHOST = 1,
    GFS_FACE_VACUUM = 2,
    GFS_FACE_FIXED_WALL = 3
};

/* Keep this independent of eunha.h so the classification contract is easy to
 * test without the MPI/CUDA application headers. */
static GFS_BOUNDARY_HD inline int
gfs_particle_kind_from_flags(int is_ghost, int is_vacuum)
{
    if (!is_ghost) return GFS_PARTICLE_FLUID;
    return is_vacuum ? GFS_PARTICLE_VACUUM_GENERATOR
                     : GFS_PARTICLE_MIRROR_GHOST;
}

static GFS_BOUNDARY_HD inline int
gfs_face_kind_from_particle_kind(int particle_kind)
{
    if (particle_kind == GFS_PARTICLE_VACUUM_GENERATOR)
        return GFS_FACE_VACUUM;
    if (particle_kind == GFS_PARTICLE_MIRROR_GHOST)
        return GFS_FACE_MIRROR_GHOST;
    return GFS_FACE_INTERIOR;
}

#endif
