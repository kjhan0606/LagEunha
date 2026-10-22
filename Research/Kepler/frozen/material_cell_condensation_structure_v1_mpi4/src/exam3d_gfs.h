#ifndef EXAM3D_GFS_H
#define EXAM3D_GFS_H
/*
 * 3D geometric-face-scheme (GFS) hydro driver, CPU/OpenMP, single MPI rank.
 * Marker string embedded in any binary that links it: LAGEUNHA_3D_GFS_V1.
 *
 * Entry points
 *   exam3d_gfs_is_hydro3d(file)  1 if the params file has
 *                                "define Simulation Model = Hydro3D"
 *   exam3d_gfs_run(file, myid, nranks)
 *                                run the test described by the params file.
 *                                Only myid == 0 computes; other ranks return 0.
 * Used by eunha2.c (dispatch before the 2D parameter reader) and by the
 * standalone Exam/exam3d_main.c (lag3d.exe). See the header of
 * exam3d_gfs.c for the physics and the flags.
 */
int exam3d_gfs_is_hydro3d(const char *paramfile);
int exam3d_gfs_run(const char *paramfile, int myid, int nranks);
#endif
