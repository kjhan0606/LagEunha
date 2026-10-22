"""Saved-prefix accounting and exact geometric correction-domain checks.

This is not a trajectory or physical-heating test. Native analytic field and
elementary-function enclosure premises remain those of the frozen evaluator.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire, quad
from audit_shared_pending29_projection_v1 import ball, add, mul, neg, encloses, verify
from audit_shared_regular_q_transfer_v1 import binary_round, sha

ROOT = Path(__file__).resolve().parent

def mag(x):
    return abs(x[0]) + x[1]

def gradients(state, source):
    at = 64 + 1064 * source
    v = [tuple(F(x) for x in struct.unpack_from('<2d', state, at+16*j)) for j in range(3)]
    ax, ay = v[0]; bx, by = v[1]; cx, cy = v[2]
    det = (bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
    g1 = ((cy-ay)/det, -(cx-ax)/det)
    g2 = (-(by-ay)/det, (bx-ax)/det)
    return v, [tuple(-g1[k]-g2[k] for k in range(2)), g1, g2]

def binomial_cap(anchor, X, Y, direction, eps):
    ax, ay = anchor
    dx, dy = add(X,(-ax,F(0))),add(Y,(-ay,F(0)))
    z = ax*ax+ay*ay+eps*eps
    numerator = add(mul((F(2),F(0)),add(mul((ax,F(0)),dx),mul((ay,F(0)),dy))),
                    add(mul(dx,dx),mul(dy,dy)))
    argument = numerator[0]/z,numerator[1]/z
    gap = F(3,4)-mag(argument)
    assert gap > 0
    Dx,Dy = direction
    first = mul((F(2)/z,F(0)),add(mul(X,(Dx,F(0))),mul(Y,(Dy,F(0)))))
    derivative = mag(first)+(Dx*Dx+Dy*Dy)/z
    return min(F(1),gap/(2*derivative)) if derivative else F(1)

def audit_domains(domains, fm, state):
    by_key = defaultdict(list)
    for z in domains:
        slot,source,target,branch,panel,axis,stage,complete = struct.unpack_from('<8Q',z)
        assert stage == 6 and complete == 1
        frame = fm[slot]
        assert (source,target) == struct.unpack_from('<2Q',frame,8)
        assert 0 <= branch < 23
        by_key[slot,panel,axis].append(z)
        u,alpha,beta = ball(z,64),ball(z,96),quad(z,736)
        radius = quad(z,752)
        assert 0 < radius <= 1
        r,dr = ball(z,128),ball(z,160)
        assert r[0]-r[1] > 0 and 2*mag(dr)*radius <= r[0]-r[1]
        assert radius <= quad(z,768) <= 1
        v,g = gradients(state,source)
        eps = F(struct.unpack_from('<d',state,64+1064*source+256)[0])
        e = [ball(frame,512+32*k)[0] for k in range(2)]
        d = [ball(frame,448+32*k)[0] for k in range(2)]
        gu = [sum(g[a][k]*e[k] for k in range(2)) for a in range(3)]
        gv = [sum(g[a][k]*d[k] for k in range(2)) for a in range(3)]
        apex = struct.unpack_from('<i',frame,48)[0]
        def combination(weights):
            return add((weights[apex],F(0)),add(
                mul((sum(weights[a]*(gu[a]+beta*gv[a]) for a in range(3)),F(0)),u),
                mul((sum(weights[a]*gv[a] for a in range(3)),F(0)),alpha)))
        for a in range(3):
            w = [F(0)]*3; w[a] = F(1)
            encloses(ball(z,192+32*a),combination(w))
            encloses(ball(z,288+32*a),(gv[a],F(0)))
        for k in range(2):
            encloses(ball(z,384+32*k),add(ball(frame,384+32*k),add(
                mul((e[k]+beta*d[k],F(0)),u),mul((d[k],F(0)),alpha))))
        assert quad(z,768) <= binomial_cap(v[apex],ball(z,384),ball(z,416),d,eps)
        for edge in range(3):
            a,b = struct.unpack_from('<2i',state,64+1064*source+944+8*edge)
            w = [F(0)]*3; w[a] = w[b] = F(1)
            D = ball(z,448+32*edge)
            encloses(D,combination(w))
            floor = D[0]-D[1]
            derivative = add(ball(z,288+32*a),ball(z,288+32*b))
            assert floor > 0 and 2*mag(derivative)*radius <= floor
            N = ball(z,192+32*b)
            center = binary_round(N[0]/D[0],113,-16382)[0]
            w[a],w[b] = -center,1-center
            # Verify the exact correlated rational parameter at every u disk;
            # saved ball includes the same uncertainty and native rounding.
            from audit_shared_pending29_projection_v1 import div
            encloses(ball(z,544+32*edge),add((center,F(0)),div(combination(w),D)))
            w[a],w[b] = gv[b],-gv[a]
            determinant = ball(z,640+32*edge)
            encloses(determinant,combination(w))
            lipschitz = quad(z,832+16*edge)
            assert lipschitz >= 2*mag(determinant)/(floor*floor)
            assert 0 < quad(z,784+16*edge) <= 1
            parameter = ball(z,544+32*edge)
            edge_direction = tuple(v[b][k]-v[a][k] for k in range(2))
            edge_X = [add((v[a][k],F(0)),mul(parameter,(edge_direction[k],F(0)))) for k in range(2)]
            anchor = v[b] if parameter[0] > F(1,2) else v[a]
            assert quad(z,784+16*edge) <= binomial_cap(anchor,*edge_X,edge_direction,eps)
            assert radius*lipschitz <= quad(z,784+16*edge)
    return by_key


def main():
    old = ROOT/'shared_oblique_panels_v1_mpi1'
    run = ROOT/'shared_oblique_v_suffix_v1_mpi1'
    c = verify(run)
    def read(name, width):
        h, records = wire(run/name)
        assert h[3] == width
        return records
    frames = read('unchanged_frames_rank0.bin', 2064)
    panels = read('unchanged_candidate_panels_rank0.bin', 192)
    slabs = read('unchanged_actual_slabs_rank0.bin', 304)
    assert frames == wire(old/'unchanged_oblique_frames_rank0.bin')[1]
    assert panels == wire(old/'candidate_oblique_panels_rank0.bin')[1]
    assert slabs == wire(old/'unchanged_actual_slabs_rank0.bin')[1]
    retained = read('retained800_axis_proofs_rank0.bin', 656)
    assert retained == [b for b in wire(old/'candidate_axis_proofs_rank0.bin')[1]
                        if struct.unpack_from('<Q', b, 24)[0]]
    axes = read('new_geometric_v_axis_proofs_rank0.bin', 656)
    points = read('new_passed_full_field_points_rank0.bin', 1488)
    domains = read('geometric_v_prior_domains_rank0.bin', 880)
    fm = {struct.unpack_from('<Q', f)[0]: f for f in frames}
    pb = defaultdict(list)
    for b in panels:
        pb[struct.unpack_from('<Q', b)[0]].append(b)
    kept = {struct.unpack_from('<3Q', b) for b in retained}
    passed = [b for b in axes if struct.unpack_from('<Q', b, 24)[0]]
    assert len(passed) == len(points) == c['new_complex_axis_passes'] == 2005
    ab = defaultdict(list)
    for b in axes:
        slot, index, axis, ok = struct.unpack_from('<4Q', b)
        assert (slot, index, axis) not in kept and axis in (0, 1)
        ab[slot].append(b)
        panel = pb[slot][index]
        U = (F(1,2), F(1,2) if axis else F(5,4))
        E = (F(1,2), F(5,4) if axis else F(1,2))
        encloses(ball(b, 32), add(ball(panel,64), mul(add(ball(panel,96),neg(ball(panel,64))),U)))
        encloses(ball(b, 64), add(ball(panel,128),mul(add(ball(panel,160),neg(ball(panel,128))),E)))
        if ok:
            assert ball(b,224)[0] > ball(b,224)[1] >= 0
            assert 0 <= quad(b,352) < 1 and 0 <= quad(b,368) < 1
            assert 0 <= quad(b,400) <= quad(b,416)
            measure = mul(mul(mul(ball(fm[slot],576),ball(b,192)),
                              add(ball(panel,96),neg(ball(panel,64)))),
                          add(ball(panel,160),neg(ball(panel,128))))
            encloses(ball(b,288),measure)
    for slot, records in ab.items():
        expected = [(slot,j,a) for j in range(len(pb[slot])) for a in (0,1)
                    if (slot,j,a) not in kept]
        assert [struct.unpack_from('<3Q',b) for b in records] == expected[:len(records)]
        assert all(struct.unpack_from('<Q',b,24)[0] for b in records[:-1])
        assert not struct.unpack_from('<Q',records[-1],24)[0]
    for b, p in zip(passed,points):
        assert b[:24] == p[:24] and b[256:320] == p[1424:1488]
    assert len(axes) == 2019 and len(retained) == 800 and len(panels) == 4384
    state_path = ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin'
    assert sha(state_path) == '74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    state = state_path.read_bytes()
    by_key = audit_domains(domains, fm, state)
    margins = []
    for b in axes:
        if struct.unpack_from('<Q',b,24)[0]:
            continue
        key = struct.unpack_from('<3Q',b)
        phase,edge = struct.unpack_from('<2i',b,448)
        if phase != 66 or edge < 0:
            continue
        checks = []
        for z in by_key[key]:
            source = struct.unpack_from('<Q',z,8)[0]
            a,b0 = struct.unpack_from('<2i',state,64+1064*source+944+8*edge)
            D = ball(z,448+32*edge)
            derivative = add(ball(z,288+32*a),ball(z,288+32*b0))
            lower = D[0]-D[1]-mag(derivative)*quad(z,752)
            assert lower > 0
            checks.append(str(lower))
        assert checks
        margins.append(dict(slot=key[0],panel=key[1],axis=key[2],
                            correlated_denominator_positive_lower_bounds=checks))
    assert len(margins) == 13
    assert c['passing_axis_replays'] == c['new_inverse_solves'] == c['new_canonical_root_jobs'] == c['new_physical_quadrature_queries'] == 0
    assert c['state_unchanged'] and not c['full_quadrature_allowed']
    print(json.dumps(dict(scope='exact frozen plan/prefix,139 geometric domain algebra and inequality checks; retained native analytic enclosure premises; no physical-heating validation',
        run=str(run),independent_geometric_v_suffix_audit_pass=True,
        retained_passing_axes=800,new_passing_axes=2005,total_passing_axes=2805,
        complete_geometric_prior_domains=len(domains),first_failures=14,
        failed_Hessian_cases_with_positive_correlated_edge_denominator=margins,
        remaining_incomplete_atlas_slots=29,full_quadrature_allowed=False),indent=2))

if __name__ == '__main__':
    main()
