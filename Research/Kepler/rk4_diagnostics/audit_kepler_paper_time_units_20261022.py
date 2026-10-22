#!/usr/bin/env python3
"""Saved-input arithmetic only; never launches or advances a simulation."""
import datetime
import hashlib
import json
import math
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parent
EXPECTED = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics')
if ROOT != EXPECTED:
    raise RuntimeError('Unexpected research directory')
EXAM = ROOT / 'material_full_positive_transition_v2_ada_mpi4/Exam'
THERMAL = ROOT / 'material_full_initial_v2_h200_mpi4/run/full_generating_thermal_r0.bin'
DRIVER = EXAM / 'material_full_initial_driver.cpp'
FIELD = EXAM / 'gfs_enriched_euler.h'
driver, field = DRIVER.read_text(), FIELD.read_text()
if 'cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=eps;' not in driver:
    raise RuntimeError('Reference radii require reinspection')
if 'rho=.001+ic->rin/den*inner*outer;' not in field:
    raise RuntimeError('Reference density requires reinspection')
if 'double z=r*r+ic->eps*ic->eps,omega=pow(z,-.75);' not in field:
    raise RuntimeError('Reference circular law requires reinspection')
with THERMAL.open('rb') as f:
    header = f.read(80)
magic, version, cells, real_bytes, rank, context = struct.unpack('<6Q', header[:48])
if (magic, version, real_bytes, rank, context) != (0x4d46544845524D31, 1, 16, 0, 1):
    raise RuntimeError('Generating-law header mismatch')
gamma, pressure, eps, reserved = struct.unpack('<4d', header[48:])
if not (gamma > 1 and pressure > 0 and eps == .047 and reserved == 0):
    raise RuntimeError('Unexpected generating-law values')
GM = 1.0  # unit prefactor in the inspected circular velocity/potential source
def period(radius):
    return 2 * math.pi * (radius * radius + eps * eps) ** .75 / math.sqrt(GM)
def taper(z):
    e = math.exp(-2 * abs(z))
    return (1 if z >= 0 else e) / (1 + e)
def rho(radius):
    return .001 + 2 / (radius + 1e-6) * taper((radius-2)/.25) * taper((8-radius)/.25)

result = {
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'scope': 'saved generating header, inspected frozen sources, and analytical arithmetic only; no evolution',
    'root': str(ROOT),
    'sources_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [DRIVER, FIELD]},
    'thermal_header_path': str(THERMAL),
    'thermal_header_first_80_bytes_sha256': hashlib.sha256(header).hexdigest(),
    'native_law': {'GM': GM, 'gamma': gamma, 'pressure': pressure, 'eps': eps,
                   'reference_rin': 2, 'reference_rout': 8, 'rank0_cells': cells,
                   'density': '0.001 + 2/(r+1e-6) * logistic(8*(r-2)) * logistic(8*(8-r))'},
    'period_formula': '2*pi*(r*r+eps*eps)**(3/4)/sqrt(GM)',
    'native_periods': {str(r): period(r) for r in [0., .5, 1., 2., 4., 8.]},
    'native_historical_times': {str(t): {'orbits_r2': t/period(2), 'orbits_r4': t/period(4)} for t in [6, 8, 10]},
    'native_1000_orbit_durations': {'r2': 1000*period(2), 'r4': 1000*period(4)},
    'native_rho_samples_in_scaled_paper_annulus': {str(r): rho(r) for r in [3.2, 4., 4.8]},
    'paper_source': 'https://arxiv.org/pdf/1409.7395',
    'paper_locations': 'section 4.2.4, equation 34, figures 6-8, printed pages 18-22',
    'paper_reference_eps_used_here': None,
    'paper_eps_zero_caption_arithmetic_only': {
        str(t): {'orbits_r1': t/(2*math.pi), 'orbits_rhalf': t/(2*math.pi*.5**1.5)}
        for t in [120, 600]},
    'paper_reported_1000_orbit_reference_radius': 'not uniquely established by inspected statement',
    'direct_matched_paper_run': False,
    'endpoint_selected_by_this_audit': False,
    'long_Kepler_validated': False,
}
output = ROOT / 'KEPLER_PAPER_TIME_UNITS_AUDIT_20261022.json'
with output.open('x') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps(result, indent=2))
