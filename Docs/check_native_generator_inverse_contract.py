"""Check saved inverse bounds against a conservative FTZ/DAZ extension.

No geometry, fluid RHS or trajectory. The original RN interval remains an
exact RN enclosure; this checks extra contracts without assuming process MXCSR.
"""
from fractions import Fraction as F
import argparse,json,sys
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--geometry-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_full_disk_native_volume_map_modules_20261019')
    from audit_saved_native_volume_map import inverse,rn_interval,floor_float,ceil_float
    tiny=F(float(np.finfo(float).tiny));h=F(.005)
    def expanded_inverse(a,b):
        lo,hi=rn_interval(b)
        if b==0:lo=min(lo,-tiny);hi=max(hi,tiny)
        choices=[F(a)]
        if 0<abs(F(a))<tiny:choices.append(F(0))
        intervals=[]
        for original in choices:
            pl=ceil_float(lo-original);pu=floor_float(hi-original)
            if pl>pu:continue
            low=rn_interval(pl)[0]/h;high=rn_interval(pu)[1]/h
            low=min(low,(lo-original)/h);high=max(high,(hi-original)/h)
            if pl<=0<=pu:low=min(low,-tiny/h);high=max(high,tiny/h)
            if low<=0<=high:low=min(low,-tiny);high=max(high,tiny)
            intervals.append((low,high))
        assert intervals
        return min(x[0] for x in intervals),max(x[1] for x in intervals)
    # Exact rational forward rounding, and explicit input/product/output flush.
    vals=(0.,-0.,.5,-.5,1.,-1.,1e-310,-1e-310,float(np.finfo(float).tiny),-float(np.finfo(float).tiny))
    sub=float(np.nextafter(0.,1.));rates=(0.,-0.,sub,-sub,1e-307,-1e-307,.03,-.03,.99,-.99)
    cases=0
    for a in vals:
        for g in rates:
            for daz in (False,True):
                aa=0. if daz and 0<abs(a)<float(tiny) else a
                gg=0. if daz and 0<abs(g)<float(tiny) else g
                product=float(h*F(gg));flushed=0. if abs(product)<float(tiny) else product
                values=[float(F(aa)+F(product)),float(F(aa)+h*F(gg)),float(F(aa)+F(flushed))]
                values += [0. if abs(v)<float(tiny) else v for v in values]
                for b in values:
                    lo,hi=expanded_inverse(a,b);assert lo<=F(g)<=hi,(a,g,b,daz);cases+=1
    def coords(name):
        with np.load(args.geometry_root/name/'state.npz') as f:return f['state'][:-4].reshape(-1,20)[:,:2]
    a=coords('base');b=coords('next');different=[];widened=0
    for i in range(len(a)):
        for j in (0,1):
            old= inverse(float(a[i,j]),float(b[i,j]),h)[:2]
            new=expanded_inverse(float(a[i,j]),float(b[i,j]));assert new[0]<=old[0] and new[1]>=old[1]
            if new!=old:different.append([i,j]);widened+=new!=old
    report=dict(status='passed_extended_inverse_contract_and_saved_input_comparison',forward_cases=cases,
        saved_coordinate_cases=a.size,saved_bounds_identical_to_original=len(different)==0,
        changed_coordinate_indices=different,expanded_contracts=['separate_RN','FMA_RN','subnormal_input_DAZ','product_FTZ','addition_output_FTZ'],
        physical_providers=0,accepted_time_steps=0,geometry_calls=0,
        limits='This validates the additional generic flush contracts. If saved bounds differ, original map certificates require re-evaluation; it does not silently overwrite or certify them.')
    args.output.mkdir();(args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
