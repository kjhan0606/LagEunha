"""Nested independent audits of saved lifted paths; never reintegrate geometry."""
import json,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from adaptive_coupled_mass import atomic
from smooth_mass_path import initialize_worker,sample_record_file
from predictive_evolution_batch import record_digest
from smooth_path_quadrature import node_indices,assess_values


def run():
    cfg=json.loads(Path('lift_refinement_config.json').read_text());rows=[]
    with ProcessPoolExecutor(max_workers=8,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
        for n in [8,16,32]:
            root=Path(cfg['source']);path=root/('n'+str(n)+'_lifted_record.json');record=json.loads(path.read_text())['record'];digest=record_digest(record)
            source=json.loads((root/'preserved'/('n'+str(n)+'_initial_state.json')).read_text());mass=np.asarray(source['mass'])
            cache_path=Path('n'+str(n)+'_samples.json');cache=json.loads(cache_path.read_text()) if cache_path.exists() else {};history=[]
            for count in [9,17,33,65,129]:
                missing=[j for j in node_indices(count) if str(j) not in cache]
                futures={executor.submit(sample_record_file,str(path),digest,mass,j/128):j for j in missing}
                for f in as_completed(futures):cache[str(futures[f])]=f.result();atomic(cache_path,cache)
                checked=assess_values([cache[str(j)] for j in node_indices(count)],mass,record['dt'],.05,cfg['budgets'],True)
                row=dict(nodes=count,ratio=checked['ratio'],physical_ratio=checked['physical_ratio'],mass_ratio=checked['mass_ratio'],quadrature_ratio=checked['quadrature_ratio'],generator_defect=checked['generator'])
                history.append(row);print(n,row,flush=True)
                if checked['ratio']<=1:break
            rows.append(dict(n=n,history=history,passed=bool(checked['ratio']<=1)))
    atomic('lift_refinement_report.json',dict(completed=True,rows=rows,note='Read-only reuse of saved lifted polynomials; no new geometry/root integration. Independent actual polynomial derivatives and original budgets.'))


if __name__=='__main__':run()
