"""One deterministic positive support proposal from saved native constraints.

All original constraint rows remain present. No singular direction is deleted;
this creates candidates only. Native binary80 refinement/physical audits govern
adoption. It does not evolve fields or replace the original quadrature.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '8'
os.environ['MKL_NUM_THREADS'] = '8'
os.environ['OMP_NUM_THREADS'] = '8'
import hashlib
import json
import pathlib
import time
import numpy as np
import scipy
from scipy.optimize import nnls

root = pathlib.Path(__file__).resolve().parent
source = root/'material_cubature_feasibility_v1_mpi4'
dest = root/'material_positive_support_v1'
dest.mkdir()
reports = []
for rank in range(4):
    for line in (source/'run'/f'constraints_rank{rank}.jsonl').read_text().splitlines():
        meta = json.loads(line); p = source/'run'/meta['file']; start = time.monotonic()
        with p.open('rb') as f:
            h = np.fromfile(f,'<u8',8); rows,n = map(int,h[3:5])
            slots = np.fromfile(f,'<u4',n); weight = np.fromfile(f,'<f8',n)
            A = np.fromfile(f,np.longdouble,rows*n).reshape(rows,n)
            target = np.fromfile(f,np.longdouble,rows); absolute = np.fromfile(f,np.longdouble,rows)
            assert not f.read()
        # Even identically zero rows stay in the proposed least-squares system.
        scale = np.where(absolute>0,absolute,np.longdouble(1))
        B = np.asarray(A*weight.astype(np.longdouble)/scale[:,None],dtype=np.float64)
        b = np.asarray(target/scale,dtype=np.float64)
        report = {'rank':rank,'triangle':meta['triangle'],'constraints':rows,'original_points':n,
                  'all_constraint_rows_retained':True,'small_singular_directions_deleted':False,
                  'rank_certified':False,'physical_rule_committed':False}
        try:
            proposal,resnorm = nnls(B,b,maxiter=3*n)
            assert np.all(np.isfinite(proposal)) and np.all(proposal>=0)
            active = np.flatnonzero(proposal>0); w = proposal[active]*weight[active]
            assert np.all(w>0) and len(active)>0
            residual = np.sum(A[:,active]*w.astype(np.longdouble),axis=1,dtype=np.longdouble)-target
            relative = np.abs(residual)/scale
            name = f'support_r{rank}_triangle{meta["triangle"]}.bin'
            with (dest/name).open('xb') as f:
                np.asarray([0x4d46504f53535531,1,rows,n,len(active),meta['triangle'],rank,16],dtype='<u8').tofile(f)
                np.asarray(active,dtype='<u4').tofile(f)
                np.asarray(w,dtype='<f8').tofile(f)
            report.update(support_points=len(active),point_reduction_factor=n/len(active),
                          nnls_scaled_residual_norm=float(resnorm),
                          binary80_max_all_constraint_residual_over_abs_sum=float(relative.max()),
                          input=str(p),input_sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),
                          support_file=name,positive_support_proposed=True)
        except Exception as exc:
            report.update(positive_support_proposed=False,error=str(exc))
        report['seconds']=time.monotonic()-start
        reports.append(report)
        with (dest/'progress.jsonl').open('a') as f: f.write(json.dumps(report)+'\n')
output={'scope':'one fixed positive-support algorithm on all 24 retained matrices; no physical evolution',
        'scipy_version':scipy.__version__,'reference':'https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html',
        'sources':len(reports),'proposals':sum(r['positive_support_proposed'] for r in reports),
        'original_points':sum(r['original_points'] for r in reports),
        'proposed_points_including_uncompressed_failures':sum(r.get('support_points',r['original_points']) for r in reports),
        'maximum_binary80_raw_proposal_constraint_residual':max((r.get('binary80_max_all_constraint_residual_over_abs_sum',float('inf')) for r in reports)),
        'all_rows_retained':True,'physical_rule_committed':False,'native_refinement_and_physical_contractions_required':True}
(dest/'contract.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
