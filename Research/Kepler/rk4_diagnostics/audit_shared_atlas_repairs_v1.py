"""Exact arithmetic audit of limited saved proofs; not a field recomputation."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import Q,rows,sha

ROOT=Path(__file__).resolve().parent
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def ball(raw,at):return quad(raw,at),quad(raw,at+16)
def contains(b,x):return b[0]-b[1]<=x<=b[0]+b[1]

def main():
    run=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_repairs_v1_mpi1')
    assert run.parent==ROOT
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('original_witness_source_integrity.txt','original_witness_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert (run/'native_exit_code.txt').read_text().strip()=='0'
    c=json.loads((run/'contract.json').read_text())
    assert json.loads((ROOT/'SHARED_ATLAS_WITNESS_20261021.json').read_text())['independent_saved_witness_accounting_pass']
    packet=(ROOT/'material_p3_finite_v2/mpi1_cpu/curved_map_control.bin').read_bytes()
    width=(len(packet)-160)//384;assert len(packet)==160+384*width
    def triangle(s):
        a=struct.unpack_from('<6d',packet,160+s*width+24)
        return [tuple(F(x) for x in a[j*2:j*2+2]) for j in range(3)]
    costs=rows(run/'cost_repairs.csv');assert {int(r['slot_id']) for r in costs}=={3248,4499}
    byslot={int(r['slot_id']):r for r in costs}
    _,panels=wire(run/'cost_repaired_panels_rank0.bin');seen=Counter();grids=defaultdict(set);grid_sizes={}
    for r in costs:
        tri=triangle(int(r['source']));minimum=None
        for j,a in enumerate(tri):
            b=tri[(j+1)%3];d=sub(b,a);parameter=max(F(),min(F(1),-dot(a,d)/dot(d,d)));closest=tuple(a[k]+parameter*d[k] for k in range(2));r2=dot(closest,closest)
            minimum=r2 if minimum is None else min(minimum,r2)
        assert minimum==F(9,524288)==F(r['exact_source_distance_squared'])
        lower=Q(r['new_source_lower']);assert 0<lower and lower*lower<=minimum
        assert Q(r['new_reference_lower'])==max(Q(r['old_reference_lower']),lower)
        assert int(r['panels'])*256==int(r['queries'])<=817152 and r['slot_cost_pass']=='1'
    for raw in panels:
        slot,source,index,target,policy,dependent,lower,upper=struct.unpack_from('<4Q4i',raw)
        assert slot in byslot and source==slot//24 and target==int(byslot[slot]['target'])
        assert index==seen[slot];seen[slot]+=1
        assert policy==0 and dependent==lower==upper==-1 and raw[192]==1
        chart,ir,it,nr,nt=struct.unpack_from('<5i',raw,208)
        assert 0<=chart<6 and 0<=ir<nr and 0<=it<nt and nt>=4
        key=slot,chart
        if key in grid_sizes:assert grid_sizes[key]==(nr,nt)
        grid_sizes[key]=nr,nt
        assert (ir,it) not in grids[key];grids[key].add((ir,it))
        for at,value in zip((240,272,304,336),(F(ir,nr),F(ir+1,nr),F(it,nt),F(it+1,nt))):assert contains(ball(raw,at),value)
        apex,part=struct.unpack_from('<2i',raw,368);assert 2*apex+part==chart
        ideal=[[F(0)]*3 for _ in range(3)];ideal[0][apex]=1;other=(apex+(2 if part else 1))%3
        ideal[2 if part else 1]=[F(int(k in (apex,other)),2) for k in range(3)]
        ideal[1 if part else 2]=[F(1,3)]*3
        for a in range(3):
            for j in range(3):assert contains(ball(raw,384+32*(3*a+j)),ideal[a][j])
        assert quad(raw,672)>quad(raw,688)>=0
        # Independent exact geometric widths may be slightly less than the
        # original outward metric. Require sufficient same-rule resolution.
        tri=triangle(source)
        pts=[tuple(sum(ideal[a][j]*tri[j][k] for j in range(3)) for k in range(2)) for a in range(3)]
        D=max(abs(pts[a][k]-pts[0][k]) for a in (1,2) for k in range(2))
        E=max(abs(pts[2][k]-pts[1][k]) for k in range(2))
        scale=min(Q(byslot[slot]['new_reference_lower']),F(float(.047)))/8
        assert nr*scale>=D and nt*scale>=E
    for (slot,ch),(nr,nt) in grid_sizes.items():assert grids[slot,ch]=={(i,j) for i in range(nr) for j in range(nt)}
    assert len(grid_sizes)==12 and dict(seen)=={slot:int(r['panels']) for slot,r in byslot.items()}
    directions=rows(run/'oblique_direction_proofs.csv');success=defaultdict(list)
    for r in directions:
        d=[Q(r['direction_x']),Q(r['direction_y'])];L=Q(r['L'])
        normals=[(Q(r[f'normal_{k}']),Q(r[f'normal_{k}_radius'])) for k in ('x','y')]
        gradient=[(Q(r[f'cone_g{k}']),Q(r[f'cone_g{k}_radius'])) for k in ('x','y')]
        nominal=sum(normals[k][0]*d[k] for k in range(2));uncertainty=sum(normals[k][1]*abs(d[k]) for k in range(2))
        threshold=sum(abs(a)+b for a,b in normals)*L*max(map(abs,d))
        stored=Q(r['sufficient_L_threshold']);assert stored>=threshold
        sufficient=abs(nominal)-uncertainty>stored
        assert r['sufficient_L_pass']==str(int(sufficient))
        center=sum(gradient[k][0]*d[k] for k in range(2));radius=sum(gradient[k][1]*abs(d[k]) for k in range(2))
        n,e=Q(r['derivative']),Q(r['derivative_radius']);assert n-e<=center-radius<=center+radius<=n+e
        passed=abs(n)>e;assert r['actual_gradient_pass']==str(int(passed))
        success[int(r['slot_id']),int(r['chart'])].append(passed)
    assert len(success)==14 and sum(all(v) for v in success.values())==c['oblique_charts_with_proved_direction']
    h,inv=wire(run/'vertex_seed_inverse_rank0.bin');assert len(inv)==1 and h[3]==608
    b=inv[0];inside=struct.unpack_from('<Q',b,32)[0]
    ir=rows(run/'vertex_seed_inverse.csv');assert len(ir)==1;ir=ir[0]
    assert bool(inside)==c['shared_vertex_inverse_pass'] and ir['inside']==str(inside)
    if inside:
        X,Y=quad(b,528),quad(b,544);L,R,res=quad(b,560),quad(b,576),quad(b,592)
        assert 0<=L<1 and R>0 and R>=2*res/(1-L)
        for k in range(2):
            inp=ball(b,48+32*k);pos=ball(b,176+32*k)
            assert res>=abs(pos[0]-inp[0])+pos[1]+inp[1]
            assert ball(b,112+32*k)==((X,Y)[k],R)
        tri=triangle(133);a,bv,cv=tri;D=sub(bv,a);E=sub(cv,a);det=D[0]*E[1]-D[1]*E[0]
        gx=[None,E[1]/det,-D[1]/det];gy=[None,-E[0]/det,D[0]/det];gx[0]=-gx[1]-gx[2];gy[0]=-gy[1]-gy[2]
        values=[None,((X-a[0])*E[1]-(Y-a[1])*E[0])/det,(D[0]*(Y-a[1])-D[1]*(X-a[0]))/det];values[0]=1-values[1]-values[2]
        for j in range(3):
            n,e=ball(inv[0],400+32*j);rad=(abs(gx[j])+abs(gy[j]))*R
            assert n-e>0 and n-e<=values[j]-rad<=values[j]+rad<=n+e
        assert quad(inv[0],496)==quad(inv[0],512)==0 # exact original reference vertex seed, no shift
        assert ir['seed_global_vertex']=='118'
    assert c['repaired_cost_slots']==2 and c['remaining_incomplete_atlas_slots']==33
    assert c['new_quadrature_queries']==c['new_canonical_root_jobs']==c['normal_chart_reconstructions']==0
    assert c['state_unchanged'] and not c['atlas_complete'] and not c['full_quadrature_allowed']
    print(json.dumps({'scope':'saved proof exact inequality/ABI/whole-chart coverage accounting; retained map/C0/L/gradient and engineering field/libm premises not independently reevaluated','run':str(run),'contract':c,'limited_saved_repair_audit_pass':True,'two_cost_slot_plans_complete':True,'cost_slots':costs,'direction_charts_with_all_active_faces_proved':sum(all(v) for v in success.values()),'shared_node_inverse_interior_and_contraction_inequalities_pass':bool(inside),'remaining_incomplete_atlas_slots':33,'full_quadrature_allowed':False},indent=2))

if __name__=='__main__':main()
