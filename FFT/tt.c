#include "mpi.h"
#include <stdio.h>
 
int main( int argc, char *argv[] )
{
    int errs = 0, len;
    int provided, flag, claimed,eclass;
	char estring[MPI_MAX_ERROR_STRING];
 
    int error = MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided );
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
 




	MPI_Error_class(error,&eclass);
	MPI_Error_string(error, estring, &len);
	printf("Error %d: %s :::: %d :: %d\n", eclass, estring, MPI_THREAD_FUNNELED,provided); fflush(stdout);

    error = MPI_Is_thread_main( &flag );
    if (!flag) {
        errs++;
        printf( "This thread called init_thread but Is_thread_main gave false\n" );fflush(stdout);
    }
    MPI_Query_thread( &claimed );
    if (claimed != provided) {
        errs++;
        printf( "Query thread gave thread level %d but Init_thread gave %d\n", claimed, provided );fflush(stdout);
    }
	printf("Errors %d\n",errs);
 
    MPI_Finalize();
    return errs;
}
