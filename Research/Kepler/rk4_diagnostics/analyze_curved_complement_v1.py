#!/usr/bin/env python3
"""Read and independently audit native artifacts; no fluid/time evolution."""
import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import struct

FIELDS = ('mass', 'px', 'py', 'angular', 'K', 'U', 'volume', 'qx', 'qy', 'polar', 'thermal_rate', 'potential')


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.directory
    native = json.loads((root / 'contract.json').read_text())
    assert not native['full_cold_transfer_proved']
    assert not native['complement_mass_matrix_positivity_proved']
    guard = Path((root / 'event_guard.txt').read_text().strip())
    assert (guard / 'status.txt').read_text().strip() == 'AUDIT_COMPLETE'
    target = {int(r['id']): r for r in table(root / 'target_moments.csv')}
    flux = {(int(r['source']), int(r['target'])): r for r in table(root / 'source_outgoing_moments.csv')}
    ids = [int(r['id']) for r in table(guard / 'site_weights.csv')]
    assert len(ids) == len(target) == 24
    source = {i: flux[i, i] for i in ids}
    with localcontext() as context:
        context.prec = 50
        assembly = {}
        global_residual = {}
        for field in FIELDS:
            largest = Decimal(0)
            for j in ids:
                expected = Decimal(source[j][field]) + sum((Decimal(flux[i, j][field]) - Decimal(flux[j, i][field]) for i in ids if i != j), Decimal(0))
                largest = max(largest, abs(expected - Decimal(target[j][field])))
            assembly[field] = float(largest)
            global_residual[field] = float(abs(sum((Decimal(source[j][field]) - Decimal(target[j][field]) for j in ids), Decimal(0))))
    reference = json.loads((root.parent / 'material_p3_finite_v2/mpi1_cpu/contract.json').read_text())['cases'][1]['cells']
    reference = {r['id']: r for r in reference}
    comparison = {
        'source_mass_relative_difference_max': max(abs(float(source[i]['mass']) / reference[i]['mass'] - 1) for i in ids),
        'source_K_difference_over_U_max': max(abs(float(source[i]['K']) - reference[i]['K']) / reference[i]['U'] for i in ids),
        'source_U_relative_difference_max': max(abs(float(source[i]['U']) / reference[i]['U'] - 1) for i in ids),
    }
    packet = root.parent / 'material_p3_finite_v2/mpi1_cpu/curved_map_control.bin'
    raw = packet.read_bytes()
    header = struct.unpack_from('=8Q', raw)
    assert header[:2] == (0x4746503346494E31, 1) and header[7] == 864
    source_id = [struct.unpack_from('=Q', raw, 160 + it * 864)[0] for it in range(header[2])]
    assert len(raw) == 160 + header[2] * 864
    saved_mass = {(i, j): [0.0, 0.0] for i in ids for j in ids if i != j}
    packets = []
    total_points = 0
    for rank in range(native['ranks']):
        path = root / f'positive_outgoing_rank{rank}.bin'
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            hraw = stream.read(48)
            digest.update(hraw)
            h = struct.unpack('=6Q', hraw)
            assert h[:4] == (0x474650334F564C31, 1, rank, native['ranks']) and h[5] == 56
            points = 0
            while True:
                block = stream.read(56 * 16384)
                if not block:
                    break
                digest.update(block)
                assert len(block) % 56 == 0
                for tri, target_index, flags, x, y, weight, rho in struct.iter_unpack('=3Q4d', block):
                    assert tri < len(source_id) and target_index < len(ids) and flags in (0, 1, 2)
                    assert all(math.isfinite(v) for v in (x, y, weight, rho)) and weight > 0 and rho > 0
                    i, j = source_id[tri], ids[target_index]
                    assert i != j, 'own fragment must be an implicit positive complement'
                    pair = saved_mass[i, j]
                    increment = weight * rho - pair[1]
                    updated = pair[0] + increment
                    pair[1] = (updated - pair[0]) - increment
                    pair[0] = updated
                    points += 1
            assert points == h[4]
            total_points += points
            packets.append({'rank': rank, 'points': points, 'bytes': 48 + 56 * points, 'sha256': digest.hexdigest()})
    assert total_points == native['positive_outgoing_points']
    mass_roundtrip = max(abs(saved_mass[i, j][0] - float(flux[i, j]['mass'])) / reference[i]['mass'] for i, j in saved_mass)
    report = {
        'scope': 'independent native-source/positive-outgoing artifact audit; not cold-remap or time-evolution acceptance',
        'native': native,
        'assembly_absolute_residual_by_moment': assembly,
        'global_conservation_absolute_residual_by_moment': global_residual,
        'source_vs_independent_finite_control': comparison,
        'positive_outgoing_packets': packets,
        'positive_outgoing_points': total_points,
        'saved_outgoing_mass_roundtrip_relative_max': mass_roundtrip,
        'full_cold_transfer_proved': False,
    }
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'source_vs_finite_control': comparison, 'outgoing_points': total_points, 'outgoing_mass_roundtrip_relative_max': mass_roundtrip, 'full_cold_transfer_proved': False}, indent=2))


if __name__ == '__main__':
    main()
