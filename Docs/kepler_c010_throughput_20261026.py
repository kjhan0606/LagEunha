"""Extract useful native timestep and stage-time distributions from C010/to10."""

from pathlib import Path
import json
import re

import numpy as np

SOURCE = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/C010/to10/log')
TIME = re.compile(r'^Time=\s*([0-9.eE+-]+)\s+step=\s*(\d+)\s+dt=\s*([0-9.eE+-]+)')
CPU = re.compile(r'^\[CPU step\s+(\d+)\].*total=([0-9.eE+-]+)s')
TESS = re.compile(r'^\[GPU tess\].*time=([0-9.eE+-]+) ms')


def stats(a):
    a = np.asarray(a, float)
    return dict(n=len(a), min=float(np.min(a)), p05=float(np.quantile(a, .05)),
                median=float(np.median(a)), p95=float(np.quantile(a, .95)),
                max=float(np.max(a)))


def main():
    times, steps, dt, cpu, tess = [], [], [], [], []
    with SOURCE.open(errors='replace') as stream:
        for line in stream:
            if m := TIME.match(line):
                times.append(float(m.group(1)))
                steps.append(int(m.group(2)))
                dt.append(float(m.group(3)))
            elif m := CPU.match(line):
                cpu.append(float(m.group(2)))
            elif m := TESS.match(line):
                tess.append(float(m.group(1)))
    times=np.asarray(times);steps=np.asarray(steps);dt=np.asarray(dt)
    if len(times)<100 or not np.all(np.diff(times)>=0) or not np.all(np.diff(steps)==1):
        raise ValueError('Unexpected incomplete or nonmonotonic C010 log')
    bins=[]
    for lo,hi in ((6.,6.25),(6.25,6.5),(6.5,6.6)):
        take=(times>=lo)&(times<hi)
        bins.append(dict(tlo=lo,thi=hi,steps=int(np.sum(take)),dt=stats(dt[take]),
                         count_dt_below_1e_minus_5=int(np.sum(dt[take]<1e-5)),
                         count_dt_below_1e_minus_6=int(np.sum(dt[take]<1e-6))))
    result=dict(source=str(SOURCE), start_t=float(times[0]), final_t=float(times[-1]),
                start_step=int(steps[0]), final_step=int(steps[-1]),
                time_advance=float(times[-1]-times[0]), records=len(times),
                global_dt=stats(dt), cpu_step_seconds=stats(cpu),
                gpu_tess_milliseconds=stats(tess), bins=bins,
                limitation='Log line timings only; CPU and GPU stage measurements can overlap or omit communication. Throughput from wall clock must be measured separately.')
    output=Path(__file__).with_name('C010_TO10_THROUGHPUT_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
