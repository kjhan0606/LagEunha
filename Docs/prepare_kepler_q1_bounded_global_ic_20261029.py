"""Create an independently audited bounded Q1 IC without changing disk data.

This starts a new boundary experiment at t=0. It is never a repair/restart of
the contaminated periodic t=3 state. Run geometry work on a compute server.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import time
import numpy as np
from scipy.spatial import cKDTree
from audit_kepler_native_radial_20261026 import ordered, GAMMA_MINUS_ONE
from make_kepler_bounded_polygon_ic_20261027 import bounded_areas

SOURCE_SHA = 'edb1f9d30700a95093a4f442afd3461238dd4e93ad9b075302557c15dcd693c9'


def prepare(source, output):
    start = time.monotonic()
    if output.exists():
        raise FileExistsError('Refusing to overwrite a Q1 boundary experiment')
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError('Not the canonical saved Q1 t=0 source')
    p, meta = ordered(source)
    if meta['t'] != 0 or len(p) != 384**2 or np.any(p['w2'] != 0):
        raise ValueError('Unexpected original Q1 state')
    r = np.hypot(p['x'] - 3., p['y'] - 3.)
    kept = np.flatnonzero(r < 2.65)
    disk = (r[kept] >= .5) & (r[kept] < 2.)
    if len(kept) != 90352 or np.count_nonzero(disk) != 48680:
        raise ValueError('Canonical Q1 retained/disk population changed')
    xy = np.column_stack((p['x'][kept], p['y'][kept]))
    estimates = [bounded_areas(xy, radius) for radius in (4., 5.)]
    area = estimates[-1][0]
    guard_difference = float(np.max(abs(estimates[0][0] / area - 1)))
    target_area = .5 * 1024 * 2.75**2 * np.sin(2 * np.pi / 1024)
    partition_error = abs(float(np.sum(area, dtype=np.longdouble)) - target_area)
    disk_area_error = float(np.max(abs(p['volume'][kept][disk] / area[disk] - 1)))
    clearance = float(2.75 * np.cos(np.pi / 1024) - np.max(r[kept]))
    separation = float(cKDTree(xy).query(xy, k=2, workers=1)[0][:, 1].min())
    if not (guard_difference < 1e-8 and partition_error < 1e-7
            and disk_area_error < 1e-6 and clearance > .099 and separation > 0):
        raise ValueError('Independent geometry gate failed')
    # Read raw records in original file order, then select by their preserved
    # global IDs. Never reconstruct a full restart from a partial dtype.
    records = np.frombuffer(raw, dtype=np.uint8, count=len(p)*576, offset=16).reshape(-1, 576)
    raw_ids = np.ndarray((len(p),), dtype='<u8', buffer=raw, offset=24, strides=(576,)) & np.uint64((1<<56)-1)
    records = records[np.argsort(raw_ids)][kept].copy()
    state = np.ndarray((len(kept),), dtype=p.dtype, buffer=records)
    exterior = r[kept] >= 2.
    # Retain every core/disk record byte exactly. Reinitialize the exterior
    # area/mass/thermal energy once, using native float volume storage and the
    # original Hopkins density profile, then keep these masses fixed.
    re = np.sqrt((state['x'][exterior]-3.)**2 + (state['y'][exterior]-3.)**2)
    rho = .01 + (1. + (re-2.)/.1)**-3
    state['volume'][exterior] = area[exterior]
    state['mass'][exterior] = rho * state['volume'][exterior].astype(float)
    state['den'][exterior] = rho
    state['ie'][exterior] = state['pressure'][exterior].astype(float)*state['volume'][exterior].astype(float)/GAMMA_MINUS_ONE
    preserved_fields = ['x', 'y', 'vx', 'vy', 'mass', 'ie', 'u4if']
    field_checks = {name: bool(np.array_equal(state[name][disk], p[name][kept][disk])) for name in preserved_fields}
    original_sorted = np.frombuffer(raw, dtype=np.uint8, count=len(p)*576, offset=16).reshape(-1, 576)[np.argsort(raw_ids)]
    interior_exact = bool(np.array_equal(records[~exterior], original_sorted[kept][~exterior]))
    thermal_pressure = GAMMA_MINUS_ONE*state['ie'].astype(float)/state['volume'].astype(float)
    pressure_error = float(np.max(abs(thermal_pressure/state['pressure'].astype(float)-1)))
    if not (all(field_checks.values()) and interior_exact and pressure_error < 2e-7):
        raise ValueError('Disk preservation or uniform-pressure gate failed')
    output.mkdir(parents=True)
    # Positive restart argument selects native reading. Step 1 labels the new
    # t=0 initial record; it represents no prior physical evolution.
    snapshot = output/'khout.000001.dat'
    result_raw = struct.pack('<4i', -1, len(kept), 576, 1)+records.tobytes()+raw[-52:]
    snapshot.write_bytes(result_raw)
    np.save(output/'independent_areas.npy', area)
    manifest = {'passed': True, 'source': str(source), 'source_sha256': SOURCE_SHA,
                'initial': str(snapshot), 'initial_sha256': hashlib.sha256(result_raw).hexdigest(),
                'initial_time': 0., 'native_restart_argument': 1,
                'sites': len(kept), 'disk_sites': int(disk.sum()), 'original_ID_extent': int(kept[-1]+1),
                'disk_field_exact': field_checks, 'all_core_and_disk_record_bytes_exact': interior_exact,
                'new_exterior_mass_and_energy_once': int(exterior.sum()),
                'boundary_radius': 2.75, 'sides': 1024, 'retain_radius': 2.65,
                'conservative_minimum_plane_clearance': clearance, 'minimum_pair_distance': separation,
                'domain_area': float(target_area), 'area_partition_absolute_error': partition_error,
                'independent_guard_area_relative_difference': guard_difference,
                'disk_original_float_volume_relative_error': disk_area_error,
                'maximum_thermal_pressure_relative_error': pressure_error,
                'analysis_seconds': time.monotonic()-start,
                'scope': 'Independent new bounded Q1 t=0 IC, global disk IDs/data preserved. Not a Kepler evolution or long-time accuracy result.'}
    (output/'initial_audit.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path);parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.output), indent=2))
