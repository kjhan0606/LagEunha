"""Baseline input identity and selected-gradient recovery audit; no PDE replay."""
from collections import defaultdict
from pathlib import Path
from fractions import Fraction as F
import contextlib,io,json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,encloses
from audit_shared_regular_q_transfer_v1 import rows,sha
import audit_shared_pending_clearance_v1 as regular

ROOT=Path(__file__).resolve().parent
RUN='shared_baseline_fields_v2_mpi1'

def verify_baseline(run):
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    assert (run/'native_exit_code.txt').read_text().strip()=='0'
    for name in ('previous_source_integrity.txt','registry_input_integrity.txt',
                 'registry_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    return json.loads((run/'contract.json').read_text())

def main():
    run=ROOT/RUN
    c=verify_baseline(run)
    # The same full field/boundary/branch accounting used for the new29
    # applies to these unchanged baseline panels, including failed prefixes.
    saved=io.StringIO()
    with contextlib.redirect_stdout(saved):
        regular.main(RUN,RUN,'unchanged_baseline_panels_rank0.bin',34955,500,verify_baseline)
    account=json.loads(saved.getvalue())
    jobs=rows(ROOT/'shared_baseline_field_inputs_v1/slot_jobs.csv')
    _,panels=wire(run/'unchanged_baseline_panels_rank0.bin')
    pm={struct.unpack_from('<Q',b)[0:1]+struct.unpack_from('<Q',b,16):b for b in panels}
    raw={};expected=set()
    for r in jobs:
        path=r['typed_packet']
        if path not in raw:raw[path]=Path(path).read_bytes()
        first,n,slot=map(int,(r['record_first'],r['record_count'],r['slot_id']))
        for j in range(n):
            original=raw[path][64+704*(first+j):64+704*(first+j+1)]
            assert pm[slot,j]==original;expected.add((slot,j))
    assert set(pm)==expected and len(expected)==34955
    _,derivatives=wire(run/'selected_derivative_proofs_rank0.bin')
    bychart={};first_by_slot={};selected=defaultdict(set)
    for (slot,j),p in pm.items():
        first_by_slot.setdefault(slot,p)
        ci=struct.unpack_from('<i',p,48)[0]
        ref=p[368:704]
        key=slot,ci
        if key in bychart:assert bychart[key]==ref
        else:bychart[key]=ref
        if not p[192]:
            for face in struct.unpack_from('<2i',p,40):
                if face>=0:selected[key].add(face)
    state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes()
    seen=set()
    for d in derivatives:
        slot,source,target,ch,face,dependent=struct.unpack_from('<6Q',d)
        assert len(d)==544 and d[48:384]==bychart[slot,ch]
        assert face in selected[slot,ch] and dependent in (0,1)
        key=slot,ch,face;assert key not in seen;seen.add(key)
        p=first_by_slot[slot]
        assert (source,target)==(struct.unpack_from('<Q',p,8)[0],struct.unpack_from('<Q',p,24)[0])
        v=[tuple(F(x) for x in struct.unpack_from('<2d',state,64+1064*source+16*j)) for j in range(3)]
        apex=struct.unpack_from('<i',d,48)[0]
        for k in range(2):
            B=(F(0),F(0));C=B
            for j in range(3):
                B=add(B,mul(ball(d,48+16+32*(3+j)),(v[j][k],F(0))))
                C=add(C,mul(ball(d,48+16+32*(6+j)),(v[j][k],F(0))))
            direction=add(C,(-B[0],B[1]))
            if not dependent:
                # D(t)=B-A+t*(C-B), over exactly t in [0,1].
                direction=add(add(B,(-v[apex][k],F(0))),mul((F(1,2),F(1,2)),direction))
            encloses(ball(d,448+32*k),direction)
        dot=add(mul(ball(d,384),ball(d,448)),mul(ball(d,416),ball(d,480)))
        encloses(ball(d,512),dot)
        n,e=ball(d,512);assert n-e>0 or n+e<0
    assert seen=={(slot,ch,face) for (slot,ch),faces in selected.items() for face in faces}
    assert len(derivatives)==c['selected_derivative_records']
    account.update(scope='exact frozen baseline panels, selected derivative direction/cone dot accounting, full field/measure/cache and positive origin branch accounting; native analytic/cone/libm premises retained; no all12 or trajectory claim',
                   independent_baseline_field_audit_pass=True,original_typed_panel_bytes_preserved=34955,
                   selected_derivative_records=len(derivatives),
                   previous110332_full_fields_replayed=0,
                   historical_source133_pass_data_recovery_scope_axes=7952,
                   first_v1_input_failure='CSV CRLF final column; terminal native exit2 before any field query; preserved',
                   original800_oblique_fields_already_recovered=True)
    print(json.dumps(account,indent=2))

if __name__=='__main__':main()
