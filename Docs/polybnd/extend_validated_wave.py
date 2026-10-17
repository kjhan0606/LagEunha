"""Extend a completed inspection interval within its original validation horizon.

Only the execution endpoint and configuration identities change. Original
initial arrays, accepted state, counters and cumulative budgets are retained.
The caller verifies project ownership and terminal process state before use.
"""
import argparse,hashlib,json,math
from pathlib import Path


def prepare(source,destination,end,receipt='execution.txt'):
    source=Path(source).resolve();destination=Path(destination)
    names=('simple_wave_config.json','initial_state.json','checkpoint.json','simple_wave_report.json')
    payload={name:(source/name).read_bytes() for name in names}
    cfg,initial,cp,report=(json.loads(payload[name]) for name in names)
    old_identity=hashlib.sha256(payload['simple_wave_config.json']).hexdigest()
    if (source/receipt).read_text().splitlines()[-1]!='exit=0':raise ValueError('source receipt is not successful')
    if any((source/name).exists() for name in ('pending.json','pending_samples.json')):
        raise ValueError('source has a pending transaction')
    if initial['config_sha256']!=old_identity or cp['config_sha256']!=old_identity:
        raise ValueError('source configuration identity mismatch')
    if report['config']!=cfg or not all(report.get(k,False) for k in ('completed','numerical_passed','physical_accuracy_passed')):
        raise ValueError('source inspection gates not passed')
    if cp['mass']!=initial['mass'] or cp['t']!=cfg['end'] or report['t']!=cp['t']:
        raise ValueError('source endpoint or mass mismatch')
    horizon=cfg.get('validation_end')
    if horizon is None or not all(math.isfinite(x) for x in (end,horizon)) or not cp['t']<end<=horizon:
        raise ValueError('extension must stay within the original validation horizon')
    if report.get('path_budget_duration')!=horizon-initial['t']:
        raise ValueError('source did not use the full-horizon budget')
    new_cfg=dict(cfg,end=end);config_bytes=(json.dumps(new_cfg,indent=2)+'\n').encode()
    identity=hashlib.sha256(config_bytes).hexdigest()
    new_initial=dict(initial,config_sha256=identity);new_cp=dict(cp,config_sha256=identity)
    provenance=dict(source=str(source),source_hashes={k:hashlib.sha256(v).hexdigest() for k,v in payload.items()},
        source_receipt=receipt,source_t=cp['t'],source_steps=cp['steps'],new_end=end,
        validation_end=horizon,config_sha256=identity,physical_and_cumulative_state_unchanged=True,
        history_note='Earlier accepted paths remain in source; only subsequent paths belong to this continuation.',full_validation=False)
    destination.mkdir(parents=False,exist_ok=False)
    (destination/'simple_wave_config.json').write_bytes(config_bytes)
    for name,value in [('initial_state.json',new_initial),('checkpoint.json',new_cp),('extension_provenance.json',provenance)]:
        (destination/name).write_text(json.dumps(value)+'\n')
    return provenance


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('destination');p.add_argument('--end',required=True,type=float);p.add_argument('--receipt',default='execution.txt')
    a=p.parse_args();print(json.dumps(prepare(a.source,a.destination,a.end,a.receipt)))
