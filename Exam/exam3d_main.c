/* Standalone 3D GFS driver (lag3d.exe). Same source as the eunha2.exe
 * Hydro3D dispatch; no MPI, no FFTW. Usage: lag3d.exe params3d.dat */
#include <stdio.h>
#include "exam3d_gfs.h"
int main(int argc, char **argv){
	if(argc < 2){
		fprintf(stderr, "usage: %s params3d.dat\n", argv[0]);
		return 100;
	}
	if(!exam3d_gfs_is_hydro3d(argv[1])){
		fprintf(stderr, "%s: %s is not a Hydro3D parameter file\n", argv[0], argv[1]);
		return 101;
	}
	return exam3d_gfs_run(argv[1], 0, 1);
}
