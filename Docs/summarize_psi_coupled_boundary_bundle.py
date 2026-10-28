"""Collect the completed coupled physical/energy/cost and arc-wall evidence."""
import argparse,hashlib,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--docs',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
docs=args.docs;base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_psi_coupled_midpoint_once_20261019';modules=base/'Codex_psi_coupled_modules_20261019'
files=dict(density_api='canonical_psi_coupled_density_api_2026-10-19.json',coupled='canonical_psi_coupled_midpoint_2026-10-19.json',connection='canonical_psi_coupled_connection_audit_2026-10-19.json',saved_density='canonical_psi_coupled_saved_density_2026-10-19.json',physical='canonical_psi_coupled_physical_comparison_2026-10-19.json',circular='canonical_circular_power_reference_2026-10-19.json',shared_flux='canonical_circular_shared_flux_2026-10-19.json')
evidence={name:json.loads((docs/file).read_text()) for name,file in files.items()}
assert evidence['coupled']['status']=='completed_connection_gate' and evidence['coupled']['provider_calls']==2
assert evidence['connection']['status']=='passed_saved_connection_audit'
assert all(evidence[k]['status']=='passed' for k in ('density_api','saved_density','circular','shared_flux'))
assert evidence['physical']['status']=='completed_comparison'
manifest=json.loads((modules/'manifest.json').read_text())
for name,digest in manifest['python_modules'].items():assert hashlib.sha256((modules/(name+'.py')).read_bytes()).hexdigest()==digest
for name,digest in manifest['native_files'].items():assert hashlib.sha256((modules/name).read_bytes()).hexdigest()==digest
energy=evidence['connection']['energy'];quadrature=evidence['saved_density']['states'];q48_uncertainty=energy['measured_initial_and_final_q48_constraint_uncertainty'];qdelta_uncertainty=sum(quadrature[k]['active_q24_minus_q48_absolute_sum'] for k in ('initial','final'))
old=json.loads((docs/'canonical_numerical_Psi_endpoint_bundle_2026-10-18.json').read_text())
report=dict(status='completed_bundle_for_Astra_evaluation',provider_calls=2,evolution_steps=1,accepted_time=1e-4,evidence=evidence,
 energy_assessment=dict(reconstructed_residual=energy['reconstructed_q48_residual'],transported_residual=energy['transported_residual'],measured_q48_constraint_uncertainty=q48_uncertainty,empirical_endpoint_q24_q48_uncertainty=qdelta_uncertainty,sum_constraint_and_empirical_quadrature_uncertainty=q48_uncertainty+qdelta_uncertainty,stored_sum_roundoff_bound=energy['stored_sum_roundoff_bound'],old_saved_step_reconstruction_transport_defect=old['regional_defect']['all_active']['signed_reconstruction_minus_transport'],limits='Different actual coupled trajectory. Fit/quadrature diagnostics are empirical, not certified bounds. Old/new comparison changes density closure and volume quadrature together; do not attribute all changes to Psi alone.'),
 cost_interpretation=dict(stage0_seconds=json.loads((root/'stage0/report.json').read_text())['seconds'],actual_half_seconds=json.loads((root/'midpoint/report.json').read_text())['seconds'],timing_scope='Stage times include constrained reconstruction and compiled shadows/reference/capture, exclude independent q48. Aggregate timers span three reconstructions and six RHS outer iterations. velocity_prepare nests inside rhs_outer for six calls; do not sum all timers.',concurrency='Old endpoint reference and small circular static fixtures overlapped some stage0 work; not exclusive wall-time benchmarking.',production_total_measured=False),
 current_scope='Existing .60 anchored patch only, fixed h=1e-4. True circular wall static reference/flux contract only. No new full-domain IC, timestep increase, long-time run, steering, viscosity, or fallback.',
 source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*[docs/file for file in files.values()],docs/'canonical_numerical_Psi_endpoint_bundle_2026-10-18.json',modules/'manifest.json',root/'manifest.json',base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so']})
args.output.write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],energy=report['energy_assessment'],physical_active={k:v['active'] for k,v in evidence['physical']['comparison'].items()}),indent=2))
