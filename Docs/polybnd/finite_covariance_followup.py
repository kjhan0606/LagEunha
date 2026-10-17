"""Finite residual covariance time-order discriminator and corrected spatial check."""
import json,os
from pathlib import Path
import numpy as np
import factor_2d_suite as suite
from factor_moment_evolve import advance_factors_exponential,advance_factors_sylvester
from kinetic_factor_rhs import rhs
ROOT=Path.cwd();previous=Path('/gpfs/kjhan/LagForce/Codex_B_sylvester_2d_20261010_004842')
old=Path('/gpfs/kjhan/LagForce/Codex_B_followup_2d_20261010_004623/kinetic')

def good(r):
    return r.get('status')=='completed' and r.get('min_factor_R',-1)>=-1e-20 and r.get('min_pressure',0)>0 and r.get('max_mass_flux',1)<1e-10 and r.get('final',{}).get('ledger_error',1)<1e-10

def main():
    p=previous/'results.jsonl';rows=[json.loads(s) for s in p.read_text().splitlines()] if p.exists() else []
    if len(rows)!=11 or not all(good(r) for r in rows):
        (ROOT/'gate.json').write_text(json.dumps(dict(status='blocked_by_predecessor',rows=rows),indent=2));return
    suite.rhs=rhs;report={}
    for mode,advance in [('old_loss',advance_factors_exponential),('corrected_loss',advance_factors_sylvester)]:
        folder=ROOT/mode;folder.mkdir();os.chdir(folder);suite.advance_factors=advance
        cases=[(4,.002),(4,.0005)] if mode=='old_loss' else [(6,.001),(8,.001)]
        for n,dt in cases:
            r=suite.run('gresho',n,.02,dt);r['method']=mode
            with (ROOT/'results.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
            if not good(r):
                (ROOT/'gate.json').write_text(json.dumps(dict(status='stopped_after_failure',record=r),indent=2));return
        if mode=='old_loss':
            values=[]
            for dt in (.002,.001,.0005):
                where=old if dt==.001 else folder
                z=json.loads((where/f'gresho_n4_dt{dt:g}.checkpoint.json').read_text());values.append(np.array(z['state'])[:-4].reshape(-1,20))
            for name,sl in [('p',slice(9,11)),('Pi',slice(15,19)),('E',slice(19,20))]:
                errors=[float(np.linalg.norm(a[:,sl]-b[:,sl])) for a,b in zip(values[:-1],values[1:])]
                report[name]=dict(differences=errors,order=float(np.log2(errors[0]/errors[1])) if min(errors)>1e-13 else None)
            (ROOT/'old_loss_time_order.json').write_text(json.dumps(report,indent=2))
        os.chdir(ROOT)
if __name__=='__main__':main()
