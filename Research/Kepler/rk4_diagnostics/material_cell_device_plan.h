#ifndef MATERIAL_CELL_DEVICE_PLAN_H
#define MATERIAL_CELL_DEVICE_PLAN_H
// Metadata-only planning of exactly the allocations made by Solver. No force,
// factorization, state conversion or GPU solver allocation is performed here.
namespace cell_native {
struct DevicePlan {
    U local=0,aggregate=0;
    explicit DevicePlan(const md::Mass&M,const std::string&directory){
        cc::Partition p(directory+"/cell_partition_r"+std::to_string(M.layout.rank)+".bin",M.layout.rank,M.layout.ranks);
        md::Ownership l(p.global_G,p.owner_G,M.layout.comm);
        need(p.nf==M.layout.global.size()&&l.global_owned==p.total,"device plan same original full/skeleton identity");
        U ni=0,ngsum=0,nL=0,nIG=0;
        for(const auto&c:p.cell){ni+=c.I.size();ngsum+=c.G.size();nL+=c.I.size()*(c.I.size()+1)/2;nIG+=c.I.size()*c.G.size();}
        need(ni==p.ni&&ni+p.ng==p.nf,"device plan complete original supports");
        U ns=0,na=0;for(int rank=0;rank<l.ranks;rank++){ns+=l.serve[rank].size();na+=l.ask[rank].size();}
        U n=l.global.size(),no=l.owned.size();
        U csr=(n+1)*8+p.column.size()*20+n*16+no*8+(ns+na)*8+5*n*16+
            std::max<U>(1,std::max(ns,na))*16+((no+255)/256)*16+16;
        local=p.cell.size()*sizeof(MCCell)+ni*4+ngsum*4+p.ng*4+p.ng+(p.ng+1)*8+ngsum*8+
            (nL+nIG)*16+4*p.nf*16+2*ngsum*16+2*p.ng*16+csr;
        MPI_Allreduce(&local,&aggregate,1,MPI_UINT64_T,MPI_SUM,M.layout.comm);
        need(aggregate==19244169964ULL,"device plan equals verified original-factor/skeleton allocation, without a changed trial space");
    }
};
}
#endif
