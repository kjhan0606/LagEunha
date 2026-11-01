#include <stddef.h>
#include <stdio.h>
#include <mpi.h>
#include "eunha.h"
#ifdef GFS_WEIGHT_DYNAMICS
_Static_assert(sizeof(treevorostressrk4particletype)==696,"Native record size");
_Static_assert(offsetof(treevorostressrk4particletype,weight)==576,"Legacy prefix");
_Static_assert(offsetof(treevorostressrk4particletype,weight_rk4)==608,"Weight RK4 offset");
_Static_assert(offsetof(treevorostressrk4particletype,weight_thermal)==672,"Explicit thermal offset");
#else
_Static_assert(sizeof(treevorostressrk4particletype)==576,"Default legacy ABI unchanged");
#endif
int main(void){printf("native_record_bytes=%zu\n",sizeof(treevorostressrk4particletype));return 0;}
