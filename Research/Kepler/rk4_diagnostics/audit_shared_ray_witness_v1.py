"""Exact saved witness audit and same-ray polynomial branch certificate."""
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import sha
ROOT=Path(__file__).resolve().parent
def ball(b,at):return quad(b,at),quad(b,at+16)
def main():
    run=ROOT/'shared_ray_witness_v1_mpi1';dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_ray_dependency_v1');assert dest.parent==ROOT and not dest.exists()
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('original_source_integrity.txt','original_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert (run/'native_exit_code.txt').read_text().strip()=='0';c=json.loads((run/'contract.json').read_text());assert c['failed_panels_reproduced']==2
    packet=ROOT/'material_p3_finite_v2/mpi1_cpu/curved_map_control.bin';raw=packet.read_bytes();width=(len(raw)-160)//384;xy=struct.unpack_from('<6d',raw,160+133*width+24)
    vertices=[tuple(F(v) for v in xy[2*j:2*j+2]) for j in range(3)]
    assert vertices==[(F(-3,1024),F(-3,1024)),(F(0),F(0)),(F(-3,512),F(0))]
    # Exact standard chart2: A=vertex1, B=(vertex1+vertex2)/2,
    # C=(vertex0+vertex1+vertex2)/3. Therefore R(t)=(-3,-t)/1024.
    A=vertices[1];B=tuple((vertices[1][k]+vertices[2][k])/2 for k in range(2));C=tuple(sum(v[k] for v in vertices)/3 for k in range(2))
    assert A==(0,0) and B==(F(-3,1024),0) and C==(F(-3,1024),F(-1,1024))
    h,args=wire(run/'sqrt_arguments_rank0.bin');assert h[3]==496 and len(args)==2
    ht,traces=wire(run/'regular_q_correlations_rank0.bin');assert ht[3]==1072 and len(traces)==2
    hp,panels=wire(run/'original_failed_panels_rank0.bin');assert hp[3]==704 and len(panels)==2
    oldpath=ROOT/'shared_origin_atlas_v2_mpi1/graded_source133_atlas_rank0.bin';old=oldpath.read_bytes();oh=struct.unpack_from('<8Q',old);oldrows={struct.unpack_from('<4Q',old,64+j*704)[::2]:old[64+j*704:64+(j+1)*704] for j in range(oh[4])}
    proof=[]
    for a,t,p in zip(args,traces,panels):
        slot,source,target,index,chart,axis=struct.unpack_from('<6Q',a);assert (slot,index) in ((3200,94),(3203,2)) and source==133 and chart==2 and axis==0
        assert p==oldrows[slot,index] and struct.unpack_from('<3Q',t)==(slot,target,chart)
        assert ball(t,32)==ball(a,112) and ball(t,64)==ball(a,144) and ball(t,96)==ball(a,208)
        radial=ball(a,176);angular=ball(a,208);Dx=ball(a,400);Dy=ball(a,432);square=ball(a,464)
        assert square[0]-square[1]<=0 and Dx[0]==F(-3,1024)
        n,e=angular;center=(F(9)+n*n)/1048576;error=(2*abs(n)*e+e*e)/1048576;lower=center-error
        assert lower>0 and abs(center-square[0])<=square[1]
        assert Dx[0]-Dx[1]<=F(-3,1024)<=Dx[0]+Dx[1]
        # The real-coefficient complex disk centered at n has image
        # (9+t²)/1048576 within [center ± error] in complex modulus.
        proof.append({'slot_id':slot,'target':target,'original_panel':index,'axis':0,'angular_center':str(n),'angular_disk_radius':str(e),'old_Dx_disk_radius':float(Dx[1]),'exact_Dx':str(F(-3,1024)),'ray_square_formula':'(9+t*t)/1048576','collected_square_center_exact':str(center),'collected_square_disk_radius_exact':str(error),'positive_half_plane_lower_exact':str(lower),'positive_half_plane_lower_display':float(lower),'old_square_lower_display':float(square[0]-square[1]),'same_actual_ray_branch_certified':True})
    assert c['new_inverse_solves']==c['new_canonical_root_jobs']==c['new_physical_quadrature_queries']==c['normal_panel_replays']==0 and c['state_unchanged']
    dest.mkdir();paths=[Path(__file__).resolve(),packet,run/'contract.json',run/'sqrt_arguments_rank0.bin',run/'regular_q_correlations_rank0.bin',run/'original_failed_panels_rank0.bin',oldpath]
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths))
    result={'scope':'exact saved two-failure ABI/identity and algebraic same-ray complex disk certificate; retains native canonical-field and implicit graph premises','original_witness_accounting_pass':True,'failure_cause':'independent barycentric direction intervals lose cancellation of exactly constant ray x component','both_original_failure_disks_in_positive_ray_square_half_plane_after_collection':True,'same_function_no_subdivision_or_precision_change':True,'proofs':proof,'new_physical_quadrature_queries':0,'remaining_incomplete_atlas_slots':31,'full_quadrature_allowed':False}
    out=dest/'contract.json';out.write_text(json.dumps(result,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(out)}  {out}\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
