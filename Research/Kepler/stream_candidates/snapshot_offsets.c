#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <mpi.h>
#include "eunha.h"
#define T treevorostressrk4particletype
#define FIELD(name) printf("\"" #name "\":{\"offset\":%zu,\"size\":%zu},", offsetof(T,name),sizeof(((T*)0)->name))
int main(void) {
 printf("{\"particle_size\":%zu,\"fields\":{",sizeof(T));
 FIELD(u4if); FIELD(x); FIELD(y); FIELD(mass); FIELD(vx); FIELD(vy);
 FIELD(den); FIELD(pressure); FIELD(volume); FIELD(ie); FIELD(w2);
 FIELD(stress.gUxx); FIELD(stress.gUxy); FIELD(stress.gUyx); FIELD(stress.gUyy);
 printf("\"stress.dRhodx\":{\"offset\":%zu,\"size\":%zu}}}\n",offsetof(T,stress.dRhodx),sizeof(((T*)0)->stress.dRhodx));
 return 0;
}
