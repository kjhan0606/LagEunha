"""Small saved-array aggregation of the completed numerical Psi audit."""
import hashlib
import json
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
root=base/'Codex_numerical_Psi_endpoints_v2_20261018'
out=Path('Docs/canonical_numerical_Psi_endpoint_bundle_2026-10-18.json')
assert not out.exists(),'Preserve existing result'
report=json.loads((root/'report.json').read_text())
ledger_root=base/'Codex_saved_midpoint_Phi_flux_20261018'
ledger=json.loads((ledger_root/'report.json').read_text())
initial=np.load(root/'initial_diagnostics.npz');end=np.load(root/'accepted_diagnostics.npz')
first=np.load(root/'initial_reconstruction.npz');last=np.load(root/'accepted_reconstruction.npz')
half=np.load(ledger_root/'midpoint.npz');h=ledger['dt'];active=initial['active']
assert np.array_equal(initial['ids'],end['ids']) and np.array_equal(initial['ids'],half['ids'])
dPsi=end['q48'][:,1]-initial['q48'][:,1]
dHydro=last['state'][:-4].reshape(-1,20)[:,19]-first['state'][:-4].reshape(-1,20)[:,19]
transport=h*half['Psi_evolution_masked']
uncertainty=abs(end['q48'][:,1]-end['q24'][:,1])+abs(initial['q48'][:,1]-initial['q24'][:,1])
expected=h*half['combined_evolution_rate']
def region(mask):
    v=dPsi-transport;combined=dHydro+dPsi-expected
    return dict(cells=int(mask.sum()),reconstructed_Psi_increment=float(dPsi[mask].sum()),
        transported_Psi_increment=float(transport[mask].sum()),
        signed_reconstruction_minus_transport=float(v[mask].sum()),
        absolute_reconstruction_minus_transport_sum=float(abs(v[mask]).sum()),
        maximum_cell_reconstruction_minus_transport=float(abs(v[mask]).max(initial=0)),
        hydro_plus_reconstruction_minus_local_face_budget=float(combined[mask].sum()),
        quadrature_difference_absolute_sum=float(uncertainty[mask].sum()))
regions={}
for key in initial.files:
    if key.startswith('region_'):
        name=key[len('region_'):]
        regions[name+'_active']=region(initial[key]&active)
regions['frozen_diagnostic_only']=region(~active)
result=dict(status='completed_diagnostic_not_conservation_claim',endpoint_report=report,
    flux_ledger_report=ledger,regional_defect=regions,
    face_statistics_note='Face RMS uses Gauss8 length weights. Minimum/median/p90/p99/maximum are sampled-point statistics, with unweighted quantiles; these are not length-weighted quantiles. Interior faces only; analytical exterior mirror states are not synthesized.',
    regional_note='Regions fixed by initial diagnostic masks; intersected with active. Frozen diagnostic has masked transport zero and is not an evolved Psi law.',
    aggregation_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (root/'report.json',ledger_root/'report.json',root/'initial_diagnostics.npz',
         root/'accepted_diagnostics.npz',ledger_root/'midpoint.npz')})
out.write_text(json.dumps(result,indent=2))
print(json.dumps(regions,indent=2))
