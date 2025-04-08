/*
#ifndef EUNHA_H
#include "eunha.h"
#endif
*/
void my_MPI_Send(void *, long , MPI_Datatype , int , int , MPI_Comm );
void my_MPI_Recv(void *, long , MPI_Datatype , int , int , MPI_Comm , MPI_Status *);
void my_MPI_Sendrecv(void *, long, MPI_Datatype , int , int , void *, long , MPI_Datatype , int , int , MPI_Comm , MPI_Status *);
void Mpi_Basic_Set(SimParameters *, MPI_Comm );
void StartParallelRW(int , int , MPI_Comm); 
void CloseParallelRW(int , int , int, MPI_Comm ); 
