"""One-shot bounded log/process check on the execution host. Never launches jobs.

Call no more than once per minute. Missing/ambiguous process identity is not
authorization to restart. This checks progress, not physical correctness.
"""
import argparse,datetime,json,os,re
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pid',required=True,type=int)
parser.add_argument('--root',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_reversible_memory_physical_to_t10_20261018'))
args=parser.parse_args();root=args.root;case='memory_from_t0';folder=root/case
manifest=json.loads((root/'manifest.json').read_text());meta=manifest['cases'][case]
process={'pid':args.pid,'identity_verified':False}
try:
    proc=Path('/proc')/str(args.pid)
    cmd=[s.decode() for s in (proc/'cmdline').read_bytes().split(b'\0') if s]
    cwd=os.readlink(proc/'cwd');stat=(proc/'stat').read_text();state=stat[stat.rfind(')')+2:].split()[0]
    process.update(command=cmd,cwd=cwd,state=state,
        identity_verified=bool(cmd and cmd[0]==meta['binary'] and Path(cwd)==folder and state!='Z'))
except OSError as error:
    process['observation_error']=str(error)
with (folder/'log').open('rb') as f:
    f.seek(0,2);end=f.tell();start=max(0,end-1048576);f.seek(start);raw=f.read(end-start)
lines=raw.split(b'\n')[1:-1] if start else raw.split(b'\n')[:-1]
pending=None;latest=None;paired=None
for rawline in lines:
    line=rawline.decode(errors='replace')
    if line.startswith('[ACCEPTED_STATE]'):
        pending=dict(re.findall(r'(\w+)=([^\s]+)',line));latest=pending
    elif line.startswith('[RMEM_ACCEPTED]') and pending:
        values=dict(re.findall(r'(\w+)=([^\s]+)',line))
        if values.get('dt')==pending.get('dt') and values.get('retries')==pending.get('retries'):
            paired={'accepted':pending,'memory':values}
        pending=None
exit_path=folder/'exit.code'
out={'observed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'process':process,'runner_exit_code':exit_path.read_text().strip() if exit_path.exists() else None,
    'log_bytes_read':len(raw),'latest_accepted':latest,'latest_paired':paired,
    'scope':__doc__,'restart_performed':False}
print(json.dumps(out,indent=2))
