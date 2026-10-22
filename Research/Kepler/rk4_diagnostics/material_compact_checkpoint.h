#ifndef MATERIAL_COMPACT_CHECKPOINT_H
#define MATERIAL_COMPACT_CHECKPOINT_H
// Read the actual accepted v3 payload without initial-force/k1 replay.
// This imports coefficients and both ledgers, not a completed rezoning.
#include "material_compact_p3.h"
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <memory>
#include <string>

namespace compact_p3 {
struct PhysicalLedger {
    double time=0,h=0;uint64_t steps=0,epoch=0,stage=0;
    R external_work=0,external_torque=0,thermal_work=0;
};
struct Endpoint {R mass=0,K=0,U=0,potential=0,L=0;};
struct NumericalLedger {
    uint64_t generation=0,paired_steps=0,reduced=0;
    R U_bridge=0,Phi_bridge=0,mass_bridge=0,pressure_error=0,density_error=0,velocity_error=0;
};
struct AnchorExpectation {
    uint64_t rank=0,ranks=0,local_coefficients=0,steps=0,generation=0;
    double time=0;std::string exact_binding;
};
struct ImportedCheckpoint {
    PhysicalLedger physical;Endpoint initial,current;NumericalLedger numerical;
    std::vector<STCoefficient>coefficient;uint64_t payload_hash=0,coefficient_hash=0;
    std::string binding;
};
inline uint64_t bytes_hash(const void*ptr,size_t count,uint64_t h=14695981039346656037ULL){
    const auto*p=static_cast<const unsigned char*>(ptr);
    for(size_t i=0;i<count;i++){h^=p[i];h*=1099511628211ULL;}return h;
}
inline ImportedCheckpoint read_checkpoint(const std::string&path,const AnchorExpectation&e){
    static_assert(sizeof(PhysicalLedger)==96&&sizeof(Endpoint)==80&&sizeof(NumericalLedger)==128,"native accepted ledger ABI");
    static_assert(sizeof(STCoefficient)==32&&sizeof(R)==16&&LDBL_MANT_DIG==64,"original binary80/native coefficient ABI");
    require(e.ranks==4&&e.rank<e.ranks&&e.local_coefficients>0&&e.steps>=5&&e.time>0&&std::isfinite(e.time),"explicit actual accepted owner/time expectation");
    require(!e.exact_binding.empty()&&e.exact_binding.size()<1048576,"exact immutable physical/numerical reference binding supplied");
    std::unique_ptr<FILE,int(*)(FILE*)>f(fopen(path.c_str(),"rb"),fclose);
    require(bool(f),"actual accepted checkpoint readable; no initial-state substitute");
    uint64_t checksum=14695981039346656037ULL;
    auto read=[&](void*p,size_t n){require(fread(p,1,n,f.get())==n,"complete accepted checkpoint payload");checksum=bytes_hash(p,n,checksum);};
    uint64_t header[8];read(header,sizeof(header));
    require(header[0]==0x4d4645564f4c3031ULL&&header[1]==3&&header[2]==e.rank&&header[3]==e.ranks&&header[4]==e.local_coefficients&&header[5]==sizeof(STCoefficient)&&header[6]==sizeof(PhysicalLedger)&&header[7]==sizeof(R),"same accepted v3 owner/count/physical ABI");
    uint64_t size=0;read(&size,sizeof(size));require(size==e.exact_binding.size(),"same immutable binding length");
    ImportedCheckpoint out;out.binding.resize(size);read(out.binding.data(),size);
    require(out.binding==e.exact_binding,"same actual generating and positive-rule epoch; no refit/rebase");
    read(&out.physical,sizeof(out.physical));read(&out.initial,sizeof(out.initial));read(&out.current,sizeof(out.current));
    require(out.physical.stage==0&&out.physical.steps==e.steps&&out.physical.time==e.time&&out.physical.h>0&&std::isfinite(out.physical.h),"actual accepted state, not an intermediate/failed stage");
    uint64_t numerical_header[4];read(numerical_header,sizeof(numerical_header));
    require(numerical_header[0]==0x4d464e554d4c4431ULL&&numerical_header[1]==1&&numerical_header[2]==sizeof(NumericalLedger)&&numerical_header[3]==1,"actual active positive functional and numerical ledger ABI");
    read(&out.numerical,sizeof(out.numerical));
    require(out.numerical.generation==e.generation&&out.numerical.reduced==1&&out.numerical.paired_steps>0,"same independently verified positive-rule generation");
    for(R value:{out.physical.external_work,out.physical.external_torque,out.physical.thermal_work,out.numerical.U_bridge,out.numerical.Phi_bridge,out.numerical.mass_bridge,out.numerical.pressure_error,out.numerical.density_error,out.numerical.velocity_error})require(std::isfinite(value),"finite imported physical and numerical ledgers");
    require(out.numerical.pressure_error>=0&&out.numerical.density_error>=0&&out.numerical.velocity_error>=0,"retain accumulated nonnegative numerical budgets");
    for(const auto*p:{&out.initial,&out.current}){
        for(R value:{p->mass,p->K,p->U,p->potential,p->L})require(std::isfinite(value),"finite accepted physical endpoint");
        require(p->mass>0&&p->K>=0&&p->U>0,"positive actual accepted mass/thermal state");
    }
    out.coefficient.resize(e.local_coefficients);
    read(out.coefficient.data(),out.coefficient.size()*sizeof(STCoefficient));
    uint64_t stored=0;require(fread(&stored,sizeof(stored),1,f.get())==1&&stored==checksum&&fgetc(f.get())==EOF,"accepted checksum covers immutable binding, both ledgers and all coefficients");
    require(fclose(f.release())==0,"accepted input closed without modification");
    for(const auto&c:out.coefficient)for(int k=0;k<2;k++)require(std::isfinite(c.x[k])&&std::isfinite(c.v[k]),"finite full accepted physical field coefficients");
    out.payload_hash=checksum;out.coefficient_hash=bytes_hash(out.coefficient.data(),out.coefficient.size()*sizeof(STCoefficient));
    return out;
}
} // namespace compact_p3
#endif
