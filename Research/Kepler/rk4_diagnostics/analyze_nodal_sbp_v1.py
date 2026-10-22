"""Read native spatial audit and diagonalize its recorded Jacobian; no evolution."""
import json, sys
from pathlib import Path
import numpy as np

dest=Path(sys.argv[1]);packet=json.loads((dest/'contract.json').read_text())
def read_matrix(path):
    with path.open('rb') as f:
        header=np.fromfile(f,dtype='<u8',count=2)
        assert header[0]==0x4746535342504a31
        n=int(header[1]);return np.fromfile(f,dtype='<f8').reshape(n,n)
eigen=[];previous=None
for row in packet['density_wave']['jacobians']:
    matrix=read_matrix(dest/row['file']);n=len(matrix)
    spectrum=np.linalg.eigvals(matrix)
    item=dict(relative_FD_step=row['relative_FD_step'],variables=n,
              max_real=float(spectrum.real.max()),min_real=float(spectrum.real.min()),
              max_abs_imaginary=float(abs(spectrum.imag).max()),
              positive_real_gt_1e_7=int((spectrum.real>1.e-7).sum()),
              matrix_max_abs=float(abs(matrix).max()))
    if previous is not None:item['matrix_max_change']=float(abs(matrix-previous).max())
    previous=matrix;eigen.append(item)
    np.save(dest/(Path(row['file']).stem+'_eigenvalues.npy'),spectrum)
analytic=None
exact_spectrum=None
if (dest/'density_jacobian_exact_conservative.bin').exists():
    exact=read_matrix(dest/'density_jacobian_exact_conservative.bin')
    exact_eig=np.linalg.eigvals(exact)
    np.save(dest/'density_exact_conservative_eigenvalues.npy',exact_eig)
    exact_spectrum=dict(max_real=float(exact_eig.real.max()),min_real=float(exact_eig.real.min()),
                        max_abs_imag=float(abs(exact_eig.imag).max()))
if (dest/'density_jacobian_exact_primitive.bin').exists():
    matrix=read_matrix(dest/'density_jacobian_exact_primitive.bin')
    acoustic=np.array([i for i in range(len(matrix)) if i%4])
    entropy=np.array([i for i in range(len(matrix)) if not i%4])
    block=matrix[np.ix_(acoustic,acoustic)]
    eig=np.linalg.eigvals(block)
    metric=np.fromfile(dest/'density_acoustic_metric.bin',dtype='<f8')
    assert len(metric)==len(block) and np.all(metric>0)
    normalized=np.sqrt(metric[:,None])*block/np.sqrt(metric[None,:])
    energy_eig=np.linalg.eigvalsh((normalized+normalized.T)*.5)
    analytic=dict(scope='exact first variation of recorded SBP volume and matching-trace HLLC; stationary translating density wave',
        acoustic_variables=len(block),entropy_column_max=float(abs(matrix[:,entropy]).max()),
        max_acoustic_real=float(eig.real.max()),min_acoustic_real=float(eig.real.min()),
        max_acoustic_abs_imag=float(abs(eig.imag).max()),
        max_acoustic_energy_growth=float(energy_eig.max()),min_acoustic_energy_growth=float(energy_eig.min()))
    np.save(dest/'density_exact_acoustic_eigenvalues.npy',eig)
cells=packet['actual_cells']
states=[s for c in cells for s in c['states']]
report=dict(scope='analysis of recorded native C++ matrices and actual polygon spatial tests; no nodal timestep',
    central_Euler_diagnostic_control=packet.get('central_Euler_diagnostic_control',False),
    reconstructed_HLL_volume_candidate=packet.get('reconstructed_HLL_volume_candidate',False),
    upwind_affine_contract=packet.get('upwind_affine_contract'),
    cells=len(cells),nodes=sorted(set(c['nodes'] for c in cells)),
    algebra=packet['algebra'],
    maximum_SBP_error=max(c['SBP'] for c in cells),
    maximum_row_sum_error=max(c['row_sum'] for c in cells),
    maximum_P2_exactness_error=max(c['P2_exactness'] for c in cells),
    maximum_cubature_compatibility=max(c['cubature_compatibility'] for c in cells),
    minimum_nodal_mass_fraction=min(c['minimum_mass_fraction'] for c in cells),
    maximum_initial_pressure_error=max(c['initial_P_error'] for c in cells),
    maximum_abs_initial_mass_quadrature_relative=max(abs(c['initial_mass_quadrature_relative']) for c in cells),
    maximum_abs_initial_energy_quadrature_over_thermal=max(abs(c['initial_energy_quadrature_over_thermal']) for c in cells),
    maximum_pressure_chain_error_over_P0=max(s['chain_rule_error_over_P0'] for s in states),
    maximum_pressure_rate_over_P0=max(s['max_Pdot_over_P0'] for s in states),
    maximum_divergence_error=max(s['max_divergence_error'] for s in states),
    maximum_gravity_thermal_residual_over_P0=max(s['gravity_thermal_residual_over_P0'] for s in states),
    maximum_affine_GCL_error=max(s['affine_metric_GCL'] for s in states),
    maximum_abs_angular_residual=max(abs(s.get('angular_volume_plus_boundary_residual',0)) for s in states),
    maximum_abs_angular_residual_over_kinetic_scale=max(abs(s.get('angular_residual_over_kinetic_scale',0)) for s in states),
    maximum_abs_upwind_angular_correction_over_kinetic_scale=max(abs(s.get('upwind_angular_correction_over_kinetic_scale',0)) for s in states),
    maximum_upwind_conservative_sum_relative=max(s.get('upwind_conservative_sum_relative',0) for s in states),
    maximum_upwind_abs_pressure_rate_over_P0=max(s.get('upwind_max_abs_pressure_rate_over_P0',0) for s in states),
    density_wave=packet['density_wave'],eigenvalues=eigen,analytic_density_wave=analytic,
    exact_conservative_spectrum=exact_spectrum,
    limitations=['Whole-polygon operator supports exact affine mesh GCL only; actual Laguerre vertex motion not certified.',
      'Actual initial matching traces are tested locally; runtime MPI/weight solve is not yet connected.',
      'Positive P2 cubature may not resolve the original nonpolynomial IC; errors remain explicit.',
      'The density-wave Jacobian is a spatial subsystem on periodic squares with frozen translating metrics.',
      'No Sod/Sedov, nodal restart or long Kepler evolution has been run.'])
(dest/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
