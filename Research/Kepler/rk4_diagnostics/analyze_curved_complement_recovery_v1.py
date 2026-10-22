#!/usr/bin/env python3
"""Audit completed native CSVs/positive packets after wrapper corruption.

No native JSON/exit status is invented. No hydrodynamics or time evolution.
"""
import csv
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import struct
import sys

FIELDS = ('mass', 'px', 'py', 'angular', 'K', 'U', 'volume', 'qx', 'qy', 'polar', 'thermal_rate', 'potential')


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    recovery = Path(sys.argv[1])
    original = recovery.parent / 'curved_complement_v2'
    table = rows(recovery / 'complement_rows.csv')
    assert sorted(int(r['triangle']) for r in table) == list(range(384))
    flux = {(int(r['source']), int(r['target'])): r for r in rows(recovery / 'source_outgoing_moments.csv')}
    target = {int(r['id']): r for r in rows(recovery / 'target_moments.csv')}
    guard = Path((recovery / 'event_guard.txt').read_text().strip())
    ids = [int(r['id']) for r in rows(guard / 'site_weights.csv')]
    assert len(target) == len(ids) == 24 and len(flux) == 576
    source = {j: flux[j, j] for j in ids}
    with localcontext() as context:
        context.prec = 50
        assembly = {}
        global_error = {}
        for k in FIELDS:
            assembly[k] = float(max(abs(Decimal(target[j][k]) - Decimal(source[j][k]) - sum((Decimal(flux[i, j][k]) - Decimal(flux[j, i][k]) for i in ids if i != j), Decimal(0))) for j in ids))
            global_error[k] = float(abs(sum((Decimal(target[j][k]) - Decimal(source[j][k]) for j in ids), Decimal(0))))
    ref = json.loads((recovery.parent / 'material_p3_finite_v2/mpi1_cpu/contract.json').read_text())['cases'][1]['cells']
    ref = {r['id']: r for r in ref}
    source_error = {
        'mass_relative_max': max(abs(float(source[i]['mass']) / ref[i]['mass'] - 1) for i in ids),
        'K_over_cold_U_max': max(abs(float(source[i]['K']) - ref[i]['K']) / ref[i]['U'] for i in ids),
        'U_relative_max': max(abs(float(source[i]['U']) / ref[i]['U'] - 1) for i in ids),
    }
    packet = (recovery.parent / 'material_p3_finite_v2/mpi1_cpu/curved_map_control.bin').read_bytes()
    h = struct.unpack_from('=8Q', packet)
    assert h[:2] == (0x4746503346494e31, 1) and h[2] == 384 and h[7] == 864
    source_ids = [struct.unpack_from('=Q', packet, 160 + 864 * i)[0] for i in range(384)]
    sums = {(i, j): [0.0, 0.0] for i in ids for j in ids if i != j}
    packets = []
    count = 0
    for rank in range(4):
        path = original / f'positive_outgoing_rank{rank}.bin'
        digest = hashlib.sha256()
        with path.open('rb') as f:
            head = f.read(48)
            digest.update(head)
            qh = struct.unpack('=6Q', head)
            assert qh[:4] == (0x474650334f564c31, 1, rank, 4) and qh[5] == 56
            nr = 0
            while True:
                data = f.read(56 * 16384)
                if not data:
                    break
                digest.update(data)
                assert len(data) % 56 == 0
                for it, dest, flags, x, y, weight, rho in struct.iter_unpack('=3Q4d', data):
                    assert it < 384 and dest < 24 and flags in (0, 1, 2)
                    assert all(math.isfinite(v) for v in (x, y, weight, rho)) and weight > 0 and rho > 0
                    i, j = source_ids[it], ids[dest]
                    assert i != j
                    pair = sums[i, j]
                    increment = weight * rho - pair[1]
                    next_sum = pair[0] + increment
                    pair[1] = (next_sum - pair[0]) - increment
                    pair[0] = next_sum
                    nr += 1
            assert nr == qh[4]
            count += nr
        packets.append({'path': str(path), 'points': nr, 'bytes': 48 + 56 * nr, 'sha256': digest.hexdigest()})
    assert count == sum(int(r['positive_points']) for r in table)
    mass_roundtrip = max(abs(sums[i, j][0] - float(flux[i, j]['mass'])) / ref[i]['mass'] for i, j in sums)
    totalU = math.fsum(float(source[i]['U']) for i in ids)
    failed = [r for r in table if r['quadrature_complete'] != '1']
    failed.sort(key=lambda r: float(r['quadrature_ratio']), reverse=True)
    report = {
        'scope': 'recovered complete first-execution native moment artifacts; wrapper later duplicated execution and truncated native JSON; original exit code unknown; no cold-transfer/evolution acceptance',
        'triangles': len(table),
        'quadrature_complete_triangles': 384 - len(failed),
        'worst_quadrature_rows': failed[:12],
        'max_quadrature_estimate_ratio': max(float(r['quadrature_ratio']) for r in table),
        'source_vs_independent_native_finite_control': source_error,
        'assembly_absolute_error': assembly,
        'global_conservation_absolute_error': global_error,
        'angular_mass_bound_at_radial_nodes': math.fsum(float(r['angular_mass_bound']) for r in table),
        'angular_K_bound_at_radial_nodes_over_total_U': math.fsum(float(r['angular_K_bound']) for r in table) / totalU,
        'positive_outgoing_points': count,
        'positive_packets': packets,
        'packet_outgoing_mass_roundtrip_relative_max': mass_roundtrip,
        'max_target_mass_drift': max(abs(float(target[i]['mass']) / float(source[i]['mass']) - 1) for i in ids),
        'positive_target_m_U_I': all(float(r['mass']) > 0 and float(r['U']) > 0 and float(r['I']) > 0 for r in target.values()),
        'minimum_resolved_Krest': min(float(r['Krest']) for r in target.values()),
        'first_native_exit_code_preserved': False,
        'complement_mass_matrix_positivity_proved': False,
        'full_cold_transfer_proved': False,
    }
    with (recovery / 'recovered_analysis.json').open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('positive_packets', 'assembly_absolute_error', 'global_conservation_absolute_error')}, indent=2))


if __name__ == '__main__':
    main()
