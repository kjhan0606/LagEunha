"""Standalone figures of completed physical accuracy experiments only.

No state is advanced or corrected. Probe arrays and source hashes are saved
alongside the figures so that local errors remain inspectable.
"""
import hashlib,json,os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from adaptive_coupled_mass import atomic
from simple_wave_evolution import install_domain
from analyze_simple_wave_time import sampled_fields
from physical_field_accuracy_audit import environment


def load(entry):
    root=Path(entry['root']);cfg=json.loads((root/'simple_wave_config.json').read_text());r=json.loads((root/'simple_wave_report.json').read_text());cp=json.loads((root/'checkpoint.json').read_text());initial=json.loads((root/'initial_state.json').read_text())
    identity=hashlib.sha256((root/'simple_wave_config.json').read_bytes()).hexdigest()
    if not r['completed'] or not r['numerical_passed'] or abs(cp['t']-cfg['end'])>1e-14 or cp['config_sha256']!=identity or initial['config_sha256']!=identity or r['config']!=cfg:raise ValueError('figure source incomplete or inconsistent')
    np.testing.assert_array_equal(cp['mass'],initial['mass'])
    return cfg,r,cp


def save(fig,name):
    fig.savefig(name+'.pdf',bbox_inches='tight');fig.savefig(name+'.svg',bbox_inches='tight');fig.savefig(name+'.png',dpi=180,bbox_inches='tight');plt.close(fig)


def convergence(entries):
    fig,axes=plt.subplots(2,3,figsize=(11,6.5),layout='constrained');colors={1:'#777777',2:'#006c9d',4:'#b53d27'};markers={1:'x',2:'o',4:'s'};rows=[]
    for problem,axrow in zip(['nonlinear_simple_wave','steady_vortex'],axes):
        selected=[]
        for entry in entries:
            cfg,r,cp=load(entry)
            if cfg.get('problem','nonlinear_simple_wave')!=problem:continue
            row=dict(label=entry['label'],root=entry['root'],problem=problem,n=cfg['n'],t=cp['t'],degree=entry['degree'],errors=r['final_score']['errors_over_amplitude'],numerical_passed=r['numerical_passed'],physical_passed=r['physical_accuracy_passed'],checkpoint_sha256=hashlib.sha256((Path(entry['root'])/'checkpoint.json').read_bytes()).hexdigest())
            rows.append(row);selected.append(row)
        for degree in [1,2,4]:
            data=sorted([v for v in selected if v['degree']==degree],key=lambda v:v['n'])
            if not data:continue
            for ax,key in zip(axrow,['density_field','velocity_field','pressure_field']):
                ax.loglog([v['n'] for v in data],[v['errors'][key] for v in data],color=colors[degree],marker=markers[degree],label='pressure degree '+str(degree))
        for ax,key in zip(axrow,['density_field','velocity_field','pressure_field']):
            goal=.001 if problem=='steady_vortex' and key=='density_field' else .01
            ax.axhline(goal,color='black',linestyle='--',linewidth=.9,label='L1 target');ax.grid(True,which='both',alpha=.25);ax.set_xlabel('cells per direction N');ax.set_xticks([16,32,64],labels=['16','32','64']);ax.xaxis.set_minor_formatter(NullFormatter());ax.set_title(('Wave' if problem=='nonlinear_simple_wave' else 'Vortex')+' / '+key.replace('_field',''));ax.set_ylabel('L1 / benchmark scale');ax.legend(fontsize=7)
    fig.suptitle('Smooth Euler benchmarks at T=0.05; positive midpoint RK2 time integration',fontsize=12)
    save(fig,'physical_spatial_convergence');return rows


def probes(entry,n):
    cfg,r,cp=load(entry);os.environ.update(environment(cfg));x,y=np.meshgrid((np.arange(n)+.3819660112501051)/n,(np.arange(n)+.6180339887498949)/n);points=np.column_stack([x.ravel(),y.ravel()]);fields=sampled_fields(cp,points)
    if cfg.get('problem','nonlinear_simple_wave')=='steady_vortex':from steady_vortex_exact import fields as reference
    else:from simple_wave_exact import fields as reference
    rho,v,P=reference(points,cp['t'],cfg['amplitude']);density_scale=1. if cfg.get('problem')=='steady_vortex' else cfg['amplitude']
    from moment_ale_evolve import GAMMA
    result=dict(density=(fields['density']-rho)/density_scale,velocity=np.linalg.norm(fields['velocity']-v,axis=1)/cfg['amplitude'],entropy=(fields['pressure']/fields['density']**GAMMA)/(P/rho**GAMMA)-1)
    np.savez_compressed(entry['label']+'_probe_fields.npz',points=points,density_numerical=fields['density'],density_reference=rho,velocity_numerical=fields['velocity'],velocity_reference=v,pressure_numerical=fields['pressure'],pressure_reference=P,**result)
    return {k:a.reshape(n,n) for k,a in result.items()}


def maps(entries,n):
    rows=[]
    for problem in ['nonlinear_simple_wave','steady_vortex']:
        selected=[v for v in entries if v.get('map',False) and load(v)[0].get('problem','nonlinear_simple_wave')==problem]
        if not selected:continue
        data=[probes(v,n) for v in selected];fig,axes=plt.subplots(3,len(selected),figsize=(3.2*len(selected)+.8,8.5),layout='constrained',squeeze=False)
        for i,(key,title) in enumerate([('density','density error / benchmark scale'),('velocity','velocity error / peak speed'),('entropy','relative entropy error')]):
            positive=key=='velocity';limit=max(float(np.max(abs(d[key]))) for d in data)
            for j,(entry,d) in enumerate(zip(selected,data)):
                ax=axes[i,j];im=ax.imshow(d[key],origin='lower',extent=[0,1,0,1],interpolation='nearest',cmap='viridis' if positive else 'RdBu_r',vmin=0 if positive else -limit,vmax=limit);ax.set_aspect('equal');ax.set_xlabel('x');ax.set_ylabel('y');ax.set_title(entry['label'] if i==0 else '')
            fig.colorbar(im,ax=list(axes[i]),shrink=.8,label=title)
        fig.suptitle(('Wave' if problem=='nonlinear_simple_wave' else 'Vortex')+' at T=.05; '+str(n)+'x'+str(n)+' independent probe points\nSampled errors; maxima are not certified continuum bounds',fontsize=11)
        save(fig,('wave' if problem=='nonlinear_simple_wave' else 'vortex')+'_physical_error_maps')
        rows.extend(dict(label=e['label'],maximum_sampled={k:float(np.max(abs(a))) for k,a in d.items()}) for e,d in zip(selected,data))
    return rows


def sections(entries,n):
    for problem in ['nonlinear_simple_wave','steady_vortex']:
        selected=[v for v in entries if v.get('map',False) and load(v)[0].get('problem','nonlinear_simple_wave')==problem]
        if not selected:continue
        cols=min(2,len(selected));rows=(len(selected)+cols-1)//cols
        fig,axes=plt.subplots(rows,cols,figsize=(5*cols,3.3*rows),layout='constrained',squeeze=False)
        for ax,entry in zip(axes.ravel(),selected):
            with np.load(entry['label']+'_probe_fields.npz') as data:
                points=data['points'].reshape(n,n,2);rho=data['density_numerical'].reshape(n,n);reference=data['density_reference'].reshape(n,n);k=n//2;x=points[k,:,0]
                ax.plot(x,reference[k],color='black',linewidth=1.6,label='exact Euler density');ax.plot(x,rho[k],color='#006c9d',linewidth=1.,label='numerical density')
                ax.set_title(entry['label']+' / y='+format(points[k,0,1],'.4f'));ax.set_xlabel('x');ax.set_ylabel('density');ax.grid(alpha=.25);ax.legend(fontsize=8)
        for ax in axes.ravel()[len(selected):]:ax.set_visible(False)
        fig.suptitle('Unfiltered reconstructed density along one independent probe row at T=.05',fontsize=11)
        save(fig,('wave' if problem=='nonlinear_simple_wave' else 'vortex')+'_density_sections')


def run():
    install_domain();cfg=json.loads(Path('physical_accuracy_figures_config.json').read_text());plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
    n=cfg.get('probe_grid',512);rows=convergence(cfg['cases']);local=maps(cfg['cases'],n);sections(cfg['cases'],n);atomic('physical_accuracy_figures_report.json',dict(completed=True,rows=rows,local=local,note='Read-only completed states only. L1 targets and numerical models unchanged. All maps and sections use actual reconstructed numerical fields. No filtering, density refitting or state correction. Figures are standalone PDF/SVG/PNG; independent probe arrays saved.'))


if __name__=='__main__':run()
