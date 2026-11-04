"""Evaluate the existing exact20pi endpoint; do not launch or resume a solver."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029')
BASE=ROOT/'native_u64_kdk_individual_20261104_v1/kdk_block2_wallkdk_v1'
RUN=BASE/'run_wall_lts_kdk_t4030_to20pi_v2'
PKG=Path(__file__).resolve().parent
DOCS=Path('/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs')
INITIAL=ROOT/'native_u64_w0_initial_20261101_v1/initial696.bin'
ADAPTER=BASE/'analysis_tools/full_reference_audit_20261104_v1/audit.py'
TARGET=20*math.pi


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path,data):
    with path.open('x') as stream:
        json.dump(data,stream,indent=2,allow_nan=False)
        stream.write('\n')


def main():
    if not (RUN/'KDK_LTS_ADIABATIC_T20PI_COMPLETE').is_file():
        raise RuntimeError('Existing run has not yet completed its native endpoint verification')
    if (RUN/'native_exit_code.txt').read_text().strip()!='0':
        raise RuntimeError('Native solver did not exit0')
    receipt=dict(line.split('=',1) for line in (RUN/'completion_receipt.txt').read_text().splitlines())
    endpoint=RUN/f"khout.{int(receipt['final_step'])}.dat"
    a=module('audit_adapter',ADAPTER)
    actual=a.gate.snapshot_time(endpoint)
    if abs(actual-TARGET)>1e-12 or abs(float(receipt['final_time'])-TARGET)>1e-12:
        raise ValueError('Requested20pi and actual footer/receipt differ')
    binary=BASE/'build_kdk_wall_v8/weight_native_U64_KDK.exe'
    if sha(binary)!='b26c3a1140178c973ef29c1af3b3437a87a60af4ac5ff44e9e6283fe6830e667':
        raise ValueError('Authoritative baseline binary changed')
    out=PKG/'final20pi'
    out.mkdir(exist_ok=False)
    full=a.compare(INITIAL,endpoint)
    restart=a.compare(RUN/'khout.165209.dat',endpoint)
    if not all(full['state_checks'].values()):
        raise ValueError('Endpoint native state-integrity check failed')
    save(out/'original_t0_metrics.json',full)
    save(out/'continuation_metrics.json',restart)
    subprocess.run([sys.executable,str(PKG/'audit_lineage_ledgers.py'),str(endpoint),
                    str(out/'lineage_ledgers.json')],check=True)
    ledger=json.loads((out/'lineage_ledgers.json').read_text())
    if not ledger['all_accepted_steps_guarded']:
        raise ValueError('Accepted-lineage nonzero correction or invalid-state guard')
    renderer=module('native_density_renderer',DOCS/'render_kepler_u64_gizmo_density_20261102.py')
    renderer.render([INITIAL,endpoint],out/'density_t0_t20pi_gizmo_scale.png')
    # The pre-existing orbit helper has an undefined circular diagnostic at
    # the exact r=0 reference site. Keep its return code and qualification.
    orbit=BASE/'analysis_tools/orbit_audit_20261104_v1/kepler_orbit_snapshot'
    proc=subprocess.run([str(orbit),str(INITIAL),str(endpoint)],capture_output=True,text=True)
    (out/'orbit_from_t0.txt').write_text(proc.stdout)
    (out/'orbit_stderr.txt').write_text(proc.stderr)
    if proc.returncode!=1 or 'mass_mismatch=0 invalid=1' not in proc.stdout or proc.stderr:
        raise ValueError('Unexpected orbital-helper error')
    p0,_=a.gate.ordered(INITIAL)
    p,_=a.gate.ordered(endpoint)
    r0=np.hypot(p0['x']-3,p0['y']-3)
    r=np.hypot(p['x']-3,p['y']-3)
    eos=a.gate.GM1*p['ie']/p['volume']
    eos_error=float(np.max(abs(p['pressure']/eos-1)))
    area_error=float(np.sum(p['volume'],dtype=np.longdouble)-a.gate.DOMAIN_AREA)
    dx,dy=p['x']-3,p['y']-3
    delta=2*np.pi/1024
    normal=(np.floor(np.arctan2(dy,dx)/delta)+.5)*delta
    excess=dx*np.cos(normal)+dy*np.sin(normal)-2.75*np.cos(np.pi/1024)
    tolerance=64*np.finfo(float).eps*np.maximum(2.75,abs(dx)+abs(dy))
    geometry_checks=dict(pressure_EOS=bool(eos_error<1e-6),
                         area_tiles_polygon=bool(abs(area_error)<1e-5),
                         sites_inside_polygon=bool(np.all(excess<=tolerance)))
    save(out/'EOS_geometry.json',dict(checks=geometry_checks,pressure_EOS_relative_max=eos_error,
                                     area_error=area_error,site_polygon_excess_max=float(excess.max()),
                                     pressure_min=float(p['pressure'].min())))
    if not all(geometry_checks.values()):
        raise ValueError('Endpoint EOS or native geometric closure failed')
    if int(np.count_nonzero(r0==0))!=1:
        raise ValueError('Unexpected center-site population')
    disk=(r0>=.5)&(r0<=2)
    mass=p0['mass']
    j0=(p0['x']-3)*p0['vy']-(p0['y']-3)*p0['vx']
    j=(p['x']-3)*p['vy']-(p['y']-3)*p['vx']
    L0=float(np.sum(mass*j0,dtype=np.longdouble))
    L=float(np.sum(mass*j,dtype=np.longdouble))
    h=4/256
    density_ratio=p['den'][disk]/p0['den'][disk]
    density_change=abs(density_ratio-1)
    disk_mass=mass[disk]
    morphology=dict(original_disk_cells=int(disk.sum()),
                    current_radius_min_max=[float(r[disk].min()),float(r[disk].max())],
                    absolute_radial_departure_quantiles={str(q):float(np.quantile(abs(r[disk]-r0[disk]),q))
                                                        for q in (.5,.95,.99,1.)},
                    disk_mass_fraction_with_radial_departure_above_initial_spacing=float(
                        np.sum(mass[disk & (abs(r-r0)>h)],dtype=np.longdouble)
                        /np.sum(mass[disk],dtype=np.longdouble)),
                    absolute_density_change_quantiles={str(q):float(np.quantile(density_change,q))
                                                       for q in (.5,.95,.99,.999,1.)},
                    density_ratio_min_max=[float(density_ratio.min()),float(density_ratio.max())],
                    original_disk_density_min_max=[float(p['den'][disk].min()),float(p['den'][disk].max())],
                    density_change_tail_mass_fraction={str(level):float(
                        np.sum(disk_mass[density_change>level],dtype=np.longdouble)
                        /np.sum(disk_mass,dtype=np.longdouble)) for level in (.1,.5,1.)},
                    density_change_tail_cell_count={str(level):int(np.count_nonzero(density_change>level))
                                                    for level in (.1,.5,1.)},
                    initial_spacing_reference=h,
                    note='Descriptive geometry/orbit measurements; no additional rejection threshold invented.')
    save(out/'disk_morphology_metrics.json',morphology)
    source_contract=json.loads((RUN/'run_contract.json').read_text())
    summary=dict(job='449508',time=actual,endpoint=str(endpoint),endpoint_sha256=sha(endpoint),
                 binary=str(binary),binary_sha256=sha(binary),
                 original_t0_energy_drift=full['total_energy_relative_from_reference'],
                 continuation_energy_drift=restart['total_energy_relative_from_reference'],
                 original_disk_J_drift=full['disk_angular_momentum_relative'],
                 full_domain_J_initial=L0,full_domain_J_final=L,full_domain_J_relative=L/L0-1,
                 global_J_change_plus_known_wall_torques=(L-L0+ledger['known_collision_torque']
                                                        +ledger['known_RK4_pressure_torque']),
                 angular_accounting_qualification='KDK pressure-wall torque is not separately logged; '
                                                 'this is not a complete discrete angular-balance residual.',
                 annular_velocity_L1=full['gizmo_L1_velocity'],
                 annular_density_L1_literal=full['gizmo_L1_density_literal'],
                 annular_density_L1_initial_plateau=full['gizmo_L1_density_plateau'],
                 disk_density_change=full['disk_density_abs_mass_weighted'],
                 disk_entropy_log_abs=full['disk_entropy_log_abs_mass_weighted'],
                 published_GIZMO_velocity_L1_approx=6e-3,
                 velocity_ratio_to_approx_paper=full['gizmo_L1_velocity']/6e-3,
                 velocity_magnitude_not_worse_than_paper_graph=full['gizmo_L1_velocity']<=6e-3,
                 paper_reference='https://arxiv.org/pdf/1409.7395#page=20',
                 published_exact_time_density_L1=None,
                 density_first_crossing_threshold_screen=full['gizmo_L1_density_literal']<.25,
                 comparison_qualification='Same metric, physical time and effective resolution; '
                                          'different IC particle placement and boundary realization. '
                                          'Graph ordinate is approximate. .25 is a first-crossing threshold, '
                                          'not actual GIZMO density error at20pi. No density-error ratio claimed.',
                 original_cohort_1pct_flags_diagnostic_only=full['original_t0_physics_checks'],
                 all_accepted_steps_guarded=ledger['all_accepted_steps_guarded'],
                 accepted_steps=ledger['accepted_steps'],morphology=morphology,
                 manual_density_image_review='pending',source_contract=source_contract,
                 model=dict(W=0,fixed_mass=True,spin=False,viscosity='off; Re=infinity',steering=False,
                            pressure_operator='GFS_VARPRESS_W0 variational geometric pressure force; exact adiabatic-volume thermal map',
                            pressure_operator_source='exam_gpu.cu:1246-1270,1290-1315; build_kdk_wall_v8.sh defines GFS_VARPRESS_W0 for CUDA',
                            MPI_ranks=1,OMP_threads=8,GPU='H200:1',
                            lineage='RK4to40.000091; globalKDKto40.2; block2KDKto20pi',
                            LTS_scope='two-level prototype; global geometry at fine events; no multi-rank scaling claim'))
    save(out/'summary.json',summary)
    manifest={str(f.relative_to(out)):dict(bytes=f.stat().st_size,sha256=sha(f))
              for f in sorted(out.iterdir()) if f.is_file()}
    save(out/'analysis_manifest.json',manifest)
    print(json.dumps({k:summary[k] for k in ['time','annular_velocity_L1','annular_density_L1_literal',
                     'velocity_ratio_to_approx_paper','original_t0_energy_drift','original_disk_J_drift',
                     'disk_density_change','disk_entropy_log_abs','all_accepted_steps_guarded']},indent=2))


if __name__=='__main__':
    main()
