"""Compare corrected source/loss integration; no overwrite of prior runs."""
import json,os
from pathlib import Path
import numpy as np
import factor_moment_evolve as short
import factor_2d_suite as spatial
from factor_moment_evolve import advance_factors_sylvester
from kinetic_factor_rhs import rhs
from factor_followup_batch import gate,PREVIOUS

short.rhs=rhs;spatial.rhs=rhs
short.advance_factors=advance_factors_sylvester
spatial.advance_factors=advance_factors_sylvester

def main():
    passed,detail=gate(PREVIOUS['kinetic'])
    if not passed:
        Path('gate.json').write_text(json.dumps({'status':'blocked_by_predecessor','detail':detail},indent=2));return
    rows=[];orders={}
    for kind,a in [('wave',.1),('shear',0.)]:
        states=[]
        for dt in (.001,.0005,.00025,.000125):
            r=short.run(kind,a,dt);r['method']='sylvester_loss';rows.append(r)
            with open('results.jsonl','a') as f:f.write(json.dumps(r)+'\n')
            if r['status']!='completed':break
            z=json.loads(Path(f'{kind}_a{a:g}_dt{dt:g}.checkpoint.json').read_text());states.append(np.array(z['state'])[:-4].reshape(-1,20))
        if len(states)==4:
            orders[kind]={}
            for name,sl in [('p',slice(9,11)),('Pi',slice(15,19)),('E',slice(19,20))]:
                errors=[float(np.linalg.norm(x[:,sl]-y[:,sl])) for x,y in zip(states[:-1],states[1:])]
                orders[kind][name]={'successive_L2_differences':errors,'orders':[float(np.log2(x/y)) if min(x,y)>1e-13 else None for x,y in zip(errors[:-1],errors[1:])]}
        Path('convergence.json').write_text(json.dumps(orders,indent=2))
    if len(rows)==8 and all(r['status']=='completed' and r['min_pressure']>0 and r['min_factor_R']>=-1e-20 and r['max_mass_flux']<1e-10 and r['final']['ledger_error']<1e-10 for r in rows):
        states=[]
        for dt in (.002,.001,.0005):
            r=spatial.run('gresho',4,.02,dt);r['method']='sylvester_loss'
            with open('results.jsonl','a') as f:f.write(json.dumps(r)+'\n')
            if r['status']!='completed':break
            z=json.loads(Path(f'gresho_n4_dt{dt:g}.checkpoint.json').read_text());states.append(np.array(z['state'])[:-4].reshape(-1,20))
        if len(states)==3:
            errors=[float(np.linalg.norm(x[:,15:19]-y[:,15:19])) for x,y in zip(states[:-1],states[1:])]
            orders['gresho_Pi']={'differences':errors,'order':float(np.log2(errors[0]/errors[1])) if min(errors)>1e-13 else None}
        Path('convergence.json').write_text(json.dumps(orders,indent=2))
if __name__=='__main__':main()
