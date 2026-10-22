"""Read-only evaluations of actual native GPU face inputs and both branches."""
from pathlib import Path
import csv, ctypes, hashlib, json, subprocess
import numpy as np

root=Path(__file__).resolve().parents[1]
base=root/'rk4_diagnostics/branch_trace'; run=base/'quarter_native'
replay=json.loads((run/'manifest.json').read_text())
if not replay['instrumentation_equivalent']:
    raise RuntimeError('instrumented replay differs; no physics interpretation permitted')
out=base/'analysis';out.mkdir(exist_ok=False)
native=root/'raw_vertex_dedup/v6/exam_gpu.after.cu'
source=native.read_text()

def extract(name):
    start=source.index('void '+name+'('); opening=source.index('{',start)
    depth=1;end=opening+1
    while depth:
        if source[end]=='{':depth+=1
        elif source[end]=='}':depth-=1
        end+=1
    return source[start:end]

functions='\n\n'.join(extract(n) for n in ('dev_hll_star_state','dev_hllc_face_2d','dev_hllc_face_2d_rest_frame'))
(out/'native_riemann_functions.c').write_text(functions+'\n')
columns=['cell_extreme','reconstructed_extreme','floor_mask','P_L_raw','P_R_raw','rho_L_raw','rho_R_raw',
         'P_L','P_R','rho_L','rho_R','vn_L','vn_R','geom_W','HC_pressure','HC_vstar',
         'HLL_pressure','HLL_vstar','pair_base','pair_HC','pair_HLL',
         'selected_pressure','selected_W','selected_Udot','HC_Udot','HLL_Udot',
         'pressure_jump_contribution','velocity_jump_contribution','pressure_ratio',
         'selected_work_is_vstar','geom_vx','geom_vy']
wrapper=r'''
#include <math.h>
#include <string.h>
#include "native_riemann_functions.c"
#include "/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src/Exam/gfs_pair.h"
/* Cell rows are exact GPU inputs: index,id,x,y,vx,vy,m,V,rho,P,cs,U,W,Wold,
 * Bxx,Bxy,Byx,Byy,dPdx,dPdy,dRhodx,dRhody. No fluid evolution occurs here. */
int eval_face(const double *a,const double *b,const double *g,const double *par,double *out){
    const double Gamma=par[0],dt=par[1]; const int geom_fv=(int)par[2];
    const double dx=b[2]-a[2],dy=b[3]-a[3],dist=hypot(dx,dy),d2=dx*dx+dy*dy;
    const double ex=dx/(dist+1.e-30),ey=dy/(dist+1.e-30),sx=g[4],sy=g[5];
    const double area=hypot(sx,sy),nx=sx/(area+1.e-30),ny=sy/(area+1.e-30);
    const double cx=.5*(g[0]+g[2]),cy=.5*(g[1]+g[3]);
    const double dvx=b[4]-a[4],dvy=b[5]-a[5];
    const double delta=-(dvx*(cx-.5*dx)+dvy*(cy-.5*dy))/dist;
    const double wx=a[4]+.5*dvx+delta*ex,wy=a[5]+.5*dvy+delta*ey;
    const double geomW=wx*sx+wy*sy,wn=wx*nx+wy*ny;
    const double lx=cx,ly=cy,rx=cx-dx,ry=cy-dy;
    const double pL0=a[9]+a[18]*lx+a[19]*ly,pR0=b[9]+b[18]*rx+b[19]*ry;
    const double rL0=a[8]+a[20]*lx+a[21]*ly,rR0=b[8]+b[20]*rx+b[21]*ry;
    const double pL=fmax(pL0,1.e-10),pR=fmax(pR0,1.e-10);
    const double rL=fmax(rL0,1.e-10),rR=fmax(rR0,1.e-10);
    const int floors=(pL0<1.e-10?1:0)|(pR0<1.e-10?2:0)|(rL0<1.e-10?4:0)|(rR0<1.e-10?8:0);
    const double vnL=a[4]*nx+a[5]*ny+nx*(a[14]*lx+a[15]*ly)+ny*(a[16]*lx+a[17]*ly);
    const double vnR=b[4]*nx+b[5]*ny+nx*(b[14]*rx+b[15]*ry)+ny*(b[16]*rx+b[17]*ry);
    double ph,uh,pc,uc;
    dev_hllc_face_2d_rest_frame(rL,pL,vnL,a[10],rR,pR,vnR,b[10],wn,Gamma,1,&pc,&uc);
    dev_hll_star_state(a[8],a[9],a[4]*ex+a[5]*ey,a[10],
                       b[8],b[9],b[4]*ex+b[5]*ey,b[10],Gamma,&ph,&uh);
    const int extreme=fmax(a[9],b[9])>100.*fmax(fmin(a[9],b[9]),1.e-30);
    const int inner=fmax(pL,pR)>100.*fmax(fmin(pL,pR),1.e-30);
    const double viS=a[4]*sx+a[5]*sy,qA=dvx*sx+dvy*sy,q=dvx*ex+dvy*ey;
    const double hW=geom_fv?geomW:uh*(ex*sx+ey*sy),sc=geomW-viS,sh=hW-viS;
    double pair0=gfs_pair_pressure(dist,a[7],b[7],a[8],b[8],a[10],b[10],fmax(-q,0));
    double pw=gfs_pair_work_limit(pair0,area,dt,a[6],b[6],q,a[11],b[11],6.);
    const double pairc=gfs_pair_charge_limit(pw,sc,qA-sc,dt,a[11],b[11],6.);
    const double pairh=gfs_pair_charge_limit(pw,sh,qA-sh,dt,a[11],b[11],6.);
    pc+=pairc;ph+=pairh;
    const double actualP=extreme?ph:pc,actualW=extreme?hW:geomW;
    const double uHC=-pc*sc,uH=-ph*sh;
    int k=0;
#define O(v) out[k++]=(v)
    O(extreme);O(inner);O(floors);O(pL0);O(pR0);O(rL0);O(rR0);
    O(pL);O(pR);O(rL);O(rR);O(vnL);O(vnR);O(geomW);O(pc);O(uc);
    O(ph);O(uh);O(pair0);O(pairc);O(pairh);O(actualP);O(actualW);
    O(-actualP*(actualW-viS));O(uHC);O(uH);O(-(ph-pc)*sc);O(-ph*(hW-geomW));
    O(fmax(a[9],b[9])/fmax(fmin(a[9],b[9]),1.e-300));O(extreme&&!geom_fv);O(wx);O(wy);
#undef O
    return k;
}
'''
(out/'face_wrapper.c').write_text('#define HLLC_VACUUM_PRATIO 100.0\n#define HLLC_VACUUM_PMIN_GUARD 0.01\n'+wrapper)
subprocess.run(['gcc','-shared','-fPIC','-O2','-o',str(out/'face_wrapper.so'),str(out/'face_wrapper.c'),'-lm'],check=True)
lib=ctypes.CDLL(str(out/'face_wrapper.so'))
array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
lib.eval_face.argtypes=[array,array,array,array,array];lib.eval_face.restype=ctypes.c_int

def evaluate(a,b,g,par):
    result=np.zeros(64)
    count=lib.eval_face(np.ascontiguousarray(a),np.ascontiguousarray(b),np.ascontiguousarray(g),np.ascontiguousarray(par),result)
    assert count==len(columns)
    return dict(zip(columns,result[:count]))

with (run/'rk4_stage_r0.tsv').open() as fh:
    stages={(int(r['step']),r['stage'],int(r['id'])):r for r in csv.DictReader(fh,delimiter='\t')}
focus={67935,68001,68526,68592};allrows=[];percall=[];boundary=[]
for call in range(1,17):
    meta=json.loads((run/f'face_input_meta_r0_call{call}.json').read_text())
    assert meta['av_mode']==5 and meta['use_muscl']==1 and meta['phase1']==1
    assert all(meta[k]==0 for k in ('alphavis','betavis','nu_phys','cd_amax','hyperv_alpha','face_charge','w2lin'))
    par=np.array([meta['Gamma'],meta['dtold'],meta['geom_fv']])
    inputs=np.loadtxt(run/f'face_input_r0_call{call}.tsv',skiprows=1,ndmin=2)
    faces=np.loadtxt(run/f'face_dump_r0_call{call}.tsv',skiprows=1,ndmin=2)
    parts={int(a[0]):a for a in inputs}
    assert all(np.all(a[12:14]==0) for a in inputs)
    digest=hashlib.sha256()
    for prefix,suffix in (('face_input','tsv'),('face_dump','tsv'),('face_input_meta','json')):
        digest.update((run/f'{prefix}_r0_call{call}.{suffix}').read_bytes())
    input_hash=digest.hexdigest()
    stage=('k1_rhs','k2_rhs','k3_rhs','k4_rhs')[(call-1)%4];step=(call-1)//4+1
    cellrows={};errors=[]
    for face in faces:
        ai,bi=int(face[2]),int(face[3]);a,b=parts[ai],parts[bi]
        assert face[6]==0
        # Dumped world vertices round at world-coordinate precision. Native
        # dS and unwrapped particle positions are captured without that loss.
        g=np.r_[face[12:14]-a[2:4],face[14:16]-a[2:4],face[10:12]]
        z=evaluate(a,b,g,par)
        viS=np.dot(a[4:6],face[10:12]);native_udot=face[9]+face[7]*viS
        errp=abs(z['selected_pressure']-face[7]);errw=abs(z['selected_W']-face[8])
        errors.append((errp,errw,abs(z['selected_Udot']-native_udot)))
        if not np.isclose(z['selected_pressure'],face[7],rtol=2e-8,atol=1e-16) or errw>2e-12:
            raise RuntimeError(f'face reproduction failed call{call} IDs{a[1]},{b[1]}: {errors[-1]}')
        row=dict(call=call,step=step,stage=stage,id_i=int(a[1]),id_j=int(b[1]),
            P_i=float(a[9]),P_j=float(b[9]),U_i=float(a[11]),input_sha256=input_hash,
            native_pressure=float(face[7]),native_W=float(face[8]),native_Udot=float(native_udot),**z)
        allrows.append(row);cellrows.setdefault(int(a[1]),[]).append(row)
        if call==1 and int(a[1]) in focus:
            # Partial positive-U probes retain captured gradients. At U->0 all
            # positive-pressure neighbors select cell-centered HLL, so these
            # gradients/floors are bypassed. The exact limit below therefore
            # does not rely on a frozen gradient approximation for target flux.
            for factor in (1.,1e-2,1e-4,0.):
                aa=a.copy();aa[11]*=factor
                aa[9]=float(np.float32((meta['Gamma']-1)*aa[11]/aa[7]))
                aa[10]=float(np.float32(np.sqrt(meta['Gamma']*aa[9]/aa[8])))
                zz=evaluate(aa,b,g,par)
                boundary.append(dict(id_i=int(a[1]),id_j=int(b[1]),factor=factor,U=float(aa[11]),P=float(aa[9]),
                    cell_extreme=bool(zz['cell_extreme']),pair_base=zz['pair_base'],Udot=zz['selected_Udot'],
                    work_vstar=bool(zz['selected_work_is_vstar']),
                    scope='one-sided cell-HLL limit; gradients bypassed' if factor==0 else 'partial positive-input probe, captured gradients retained'))
    summary=[]
    for ident in sorted(focus):
        rr=cellrows[ident];r=stages[step,stage,ident]
        summary.append(dict(id=ident,U=float(r['U_raw']),P=float(r['P_stored']),native_Udot=float(r['Udot']),
            reconstructed_from_faces_Udot=sum(a['native_Udot'] for a in rr),
            extreme_faces=int(sum(a['cell_extreme'] for a in rr)),faces=len(rr),
            pressure_ratio_max=max(a['pressure_ratio'] for a in rr),
            HC_Udot=sum(a['HC_Udot'] for a in rr),HLL_Udot=sum(a['HLL_Udot'] for a in rr),
            branch_pressure_jump=sum(a['pressure_jump_contribution'] for a in rr),
            branch_velocity_jump=sum(a['velocity_jump_contribution'] for a in rr)))
    percall.append(dict(call=call,step=step,stage=stage,input_sha256=input_hash,
        maximum_face_reproduction_errors=np.max(errors,axis=0).tolist(),focus=summary))

with (out/'faces.tsv').open('w') as fh:
    w=csv.DictWriter(fh,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
limit=[]
for ident in sorted(focus):
    data=[dict(factor=factor,Udot=sum(r['Udot'] for r in boundary if r['id_i']==ident and r['factor']==factor),
               all_faces_cell_HLL=all(r['cell_extreme'] for r in boundary if r['id_i']==ident and r['factor']==factor),
               pair_base_max=max(r['pair_base'] for r in boundary if r['id_i']==ident and r['factor']==factor))
          for factor in (1.,1e-2,1e-4,0.)]
    limit.append(dict(id=ident,probes=data))
result=dict(scope='bounded branch and one-sided pressure boundary diagnosis; no evolved operator modification',
    replay_equivalent=replay['instrumentation_equivalent'],
    native_source_sha256=hashlib.sha256(native.read_bytes()).hexdigest(),
    extracted_functions_sha256=hashlib.sha256(functions.encode()).hexdigest(),
    actual_calls=percall,one_sided_U_boundary=limit,
    boundary_faces=boundary,
    boundary_contract='At Pi->0 with positive neighbor P, all target faces take cell-centered HLL. Target pressure/density/velocity gradients are then unused. Geometry/rho/mass/velocity fixed, EOS Pi and cs_i recomputed, pair contribution checked. Positive probes outside that limit retain captured gradients and are partial. Zero is an analytic helper limit, never an evolved invalid state or pressure-floor repair.',
    conclusions_pending='Interpret percall branch/force/work identities before choosing a physics or time integrator change.')
(out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print('limits',json.dumps(limit,indent=2),flush=True)
print('focus68526',json.dumps([dict(call=c['call'],stage=c['stage'],**next(r for r in c['focus'] if r['id']==68526)) for c in percall],indent=2),flush=True)
