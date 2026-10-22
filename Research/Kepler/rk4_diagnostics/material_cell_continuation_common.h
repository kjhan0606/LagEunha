// First actual selected-rule RHS at retained step5; never a fresh IC.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
#include <cuda_runtime.h>
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
using R=long double;using U=uint64_t;
#include "material_full_owned.h"
#include "material_full_initial_batch.h"
#include "material_distributed_mass.h"
namespace st=shared_trace;namespace md=material_distributed;
namespace shared_trace_rhs {using Vec=std::vector<R>;}
#include "material_geometric_pressure_geometry.h"
namespace gp=geometric_pressure;
static R value(MFDoubleDouble a){return R(a.hi)+a.lo;}
struct Device {void*ptr=nullptr;~Device(){mf_initial_destroy(ptr);}};
static void config_rule(GfsEnrichConfig&cfg,double gamma,double pressure,double eps){
    cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.kind=2;cfg.field.ic.gamma=gamma;cfg.field.ic.pressure=pressure;cfg.field.ic.cx=cfg.field.ic.cy=12;cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=eps;
    std::vector<double>x,w;st::gauss(64,x,w);for(int i=0;i<64;i++){cfg.an[i]=x[i];cfg.aw[i]=w[i];}st::gauss(28,x,w);for(int i=0;i<28;i++){cfg.rn[i]=x[i];cfg.rw[i]=w[i];}
}
#include "material_cell_condensation_solver.h"
// This isolated frozen header has one explicit solve hook; all original
// physical force, energy/torque, per-source spin and reaction gates are retained.
#include "material_reference_cache.h"
#include "material_full_evolution_host.h"
namespace fe=full_evolution;
#include "material_reference_cache_host.h"
