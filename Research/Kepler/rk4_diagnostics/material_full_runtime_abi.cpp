// Native resource introspection/build unit, not a hydro test or time solver.
#include <cstdio>
#include "material_distributed_mass.h"
#include "material_geometric_pressure_batch.h"
int main(){
    printf("{\"scope\":\"compiled full-state resource ABI; no physical queries\",\"sizeof_binary80_storage\":%zu,\"sizeof_STTriangle\":%zu,\"sizeof_STCoefficient\":%zu,\"sizeof_STPoint\":%zu,\"sizeof_GPInput\":%zu,\"sizeof_GPOutput\":%zu,\"sizeof_local_mass_factor\":%zu,\"local_mass_block_entries\":%d,\"sparse_halo_code_builds\":true,\"full_mass_solve_validated\":false}\n",sizeof(material_distributed::R),sizeof(STTriangle),sizeof(STCoefficient),sizeof(STPoint),sizeof(GPInput),sizeof(GPOutput),sizeof(material_distributed::Factor),ST_MAX*ST_MAX);
}
