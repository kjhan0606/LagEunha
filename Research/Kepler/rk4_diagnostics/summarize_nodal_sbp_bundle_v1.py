"""Summarize terminal native spatial checks without rerunning physics."""
import hashlib,json
from pathlib import Path

root=Path(__file__).resolve().parent
names=['v1','v2','v3','v4','angular_v1','angular_v2','upwind_A_v1']
audits={}
for name in names:
    directory=root/('nodal_sbp_'+name)
    assert (directory/'status.txt').read_text().strip()=='AUDIT_COMPLETE'
    analysis=directory/'analysis.json'
    audits[name]=dict(path=str(directory),native_calculation='terminal',
        physics_approval=False,
        analysis_sha256=hashlib.sha256(analysis.read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256((directory/'contract.exe').read_bytes()).hexdigest())
split=json.loads((root/'nodal_sbp_v3/analysis.json').read_text())
central=json.loads((root/'nodal_sbp_v4/analysis.json').read_text())
upwind=json.loads((root/'nodal_sbp_upwind_A_v1/analysis.json').read_text())
bundle=dict(scope='completed native spatial rejection bundle; no evolved nodal state',
    previous_goal_turn='no progress: read-only opinion response; current turn revalidated source and completed new native checks',
    current_goal_turn='progress: native operators, exact Jacobians, independent controls and angular failures change next choice',
    final_goal='native accelerated Kepler stability through historical t6-10 and longest GIZMO duration; not achieved',
    evolution_authorized_for_these_candidates=False,
    original_runtime_sources='enriched_v5_build/source.sha256 independently rechecked; all25 entries OK',
    audits=audits,
    rejected_candidates=[
      dict(name='positive P2 nodal KEP split',reason='advected density-wave local stability, cold divergence accuracy and angular weak balance',
           exact_advected_Jacobian_max_real=split['exact_conservative_spectrum']['max_real'],
           max_pressure_rate_over_P0=split['maximum_pressure_rate_over_P0']),
      dict(name='full Euler arithmetic central diagnostic control',reason='cold thermal kinetic cancellation fails',
           exact_advected_Jacobian_max_real=central['exact_conservative_spectrum']['max_real'],
           max_pressure_rate_over_P0=central['maximum_pressure_rate_over_P0']),
      dict(name='midpoint-reconstructed Euler HLL volume correction',reason='first angular gate fails; same RHS has huge cold thermal damage',
           max_angular_residual_over_kinetic_scale=upwind['maximum_abs_angular_residual_over_kinetic_scale'],
           max_pressure_rate_over_P0=upwind['maximum_pressure_rate_over_P0'],
           downstream_P4_or_Jacobian_or_runtime='not run')],
    next_planning_input=str(root/'VOLUME_UPWIND_A_RESULTS.md'),
    planning_delegate='/root/astra_planning_review',
    current_review='one complete review delivered ASTRA_VOLUME_UPWIND_PLAN.md; its subsequent first gate is now measured and rejected',
    live_native_handles=[])
destination=root/'NODAL_SBP_BUNDLE_20261021.json'
with destination.open('x') as stream:json.dump(bundle,stream,indent=2);stream.write('\n')
print(json.dumps(dict(summary=str(destination),audits=len(audits),physics_approved=False,live_native_handles=[])))
