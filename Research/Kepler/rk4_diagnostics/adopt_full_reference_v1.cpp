#include <mpi.h>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <unistd.h>
using R=long double;
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
#include "gfs_enriched_euler.h"
#include "material_shared_trace_factory.h"
#include "material_full_reference_packet.h"
#include "material_geometric_pressure_batch.h"
// Instantiate the same canonical factory interface for the compact full
// geometry view. Do not build a globally replicated full mass/DOF operator.
template shared_trace::Factory shared_trace::topology<full_reference::Packet>(const full_reference::Packet&,int);
int main(int argc,char**argv){try{
    need(argc==2,"retained full reference artifact");auto m=full_reference::load(argv[1]);auto p=full_reference::packet(m,double(1.666667f),1.e-6,.047,2,8);
    need(p.record.size()==m.triangle.size()&&m.external_edges==2048,"native full reference packet/context dimensions");
    printf("{\"scope\":\"actual compact full-reference import and canonical factory C++ interface instantiation; no mass/field queries or evolution\",\"cells\":%zu,\"triangles\":%zu,\"native_compact_record_bytes\":%zu,\"compact_reference_packet_bytes\":%llu,\"STPoint_bytes\":%zu,\"GPOutput_bytes\":%zu,\"point_and_two_output_bytes\":%zu,\"stream_65536_point_and_two_output_bytes\":%zu,\"full_reference_import_pass\":true,\"full_mass_or_native_evolution_implemented\":false}\n",m.cell.size(),m.triangle.size(),sizeof(full_reference::Record),(unsigned long long)(p.record.size()*sizeof(full_reference::Record)),sizeof(STPoint),sizeof(GPOutput),sizeof(STPoint)+2*sizeof(GPOutput),size_t(65536)*(sizeof(STPoint)+2*sizeof(GPOutput)));return 0;
}catch(const std::exception&e){fprintf(stderr,"FULL_ADOPTION_FAIL %s\n",e.what());return 2;}}
