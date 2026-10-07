"""Read completed reference histories and generate the accelerator-note figure."""
from pathlib import Path
import os
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    cases={'BB':'261001-192802-radial-bb-diagnostics-0plus-T030-Fs54-radial-reference',
           'Polyak':'261001-233335-radial-polyak-scratch-0plus-T030-Fs54-radial-reference',
           'Anderson (pmax=5)':'261002-094542-radial-anderson-aa20-0plus-pmax5-tight-continuation-radial-reference',
           'Anderson (pmax=10)':'261002-111709-radial-anderson-aa20-0plus-pmax10-scratch-T030-Fs54-radial-reference',
           'Anderson (pmax=3)':'261002-130520-radial-anderson-aa20-0plus-pmax3-scratch-T030-Fs54-radial-reference'}
    result={};initial=None
    fig,axes=plt.subplots(2,1,figsize=(9,7),layout='constrained')
    for label,name in cases.items():
        folder=ROOT/'runs'/name
        manifest=json.loads((folder/'manifest.json').read_text())
        parent=Path(manifest['history_parent_run']) if manifest.get('history_parent_run') else None
        h=np.genfromtxt(folder/('combined_iteration_history.dat' if parent else 'iteration_history.dat'),names=True)
        metrics=dict(line.split('=',1) for line in (folder/'metrics.txt').read_text().splitlines() if '=' in line)
        state=np.loadtxt((parent or folder)/'input_fields_2d.dat')
        if initial is None: initial=state
        else: np.testing.assert_array_equal(initial,state)
        color,style={'BB':('#2166ac','-'),'Polyak':('#b2182b','--'),
                     'Anderson (pmax=5)':('#333333','-.'),
                     'Anderson (pmax=10)':('#762a83',':'),
                     'Anderson (pmax=3)':('#555555',(0,(5,1,1,1,1,1)))}[label]
        for ax,key in zip(axes,('map_rms','map_max')):
            ax.semilogy(h['iteration'],h[key],color=color,ls=style,label=label)
            ax.set_ylabel(key.replace('_',' '));ax.legend();ax.grid(alpha=.15)
        record={'run':name,'maps':len(h),'terminal_status':metrics['terminal_status'],
                'final_rms':float(h['map_rms'][-1]),'final_max':float(h['map_max'][-1]),
                'median_map_seconds':float(np.median(h['map_seconds'])),
                'summed_map_seconds':float(sum(h['map_seconds'])),
                'thresholds':{},'rms_increases':int(sum(np.diff(h['map_rms'])>0)),
                'late_rms_iterations_per_decade':float(-1/np.polyfit(h['iteration'][-200:],np.log10(h['map_rms'][-200:]),1)[0]),
                'regions':{key:float(metrics[key]) for key in ('fine_residual_rms','medium_residual_rms','outer_residual_rms')},
                'launcher_error':manifest['launcher_error']}
        if label.startswith('Anderson'):
            arrays=np.concatenate([np.stack([np.loadtxt(f) for f in sorted(p.glob('metrics.txt.signed.map*.dat'))])
                                   for p in ((parent,folder) if parent else (folder,))])
            x,r,u=(arrays[:,:,i] for i in (1,2,3))
            np.testing.assert_allclose(np.sqrt(np.mean(r*r,axis=1)),h['map_rms'],rtol=1e-12)
            np.testing.assert_allclose(np.max(abs(r),axis=1),h['map_max'],rtol=1e-12)
            np.testing.assert_allclose(x[1:],x[:-1]+u[:-1],atol=1e-14,rtol=1e-12)
            nr=np.linalg.norm(r,axis=1)
            cosine=np.sum(r[1:]*r[:-1],axis=1)/(nr[1:]*nr[:-1])
            layout=np.loadtxt(folder/'metrics.txt.signed.layout.dat');n=len(layout)
            rr=r[-1].reshape(21,n);zones=layout[:,3].astype(int)
            worst=np.argmax(abs(r[-1]))
            names=[f'{part} A{a}{b}' for a in 'xyz' for b in 'xyz' for part in ('Re','Im')]+['mf_x','mf_y','mf_z']
            applied=np.linalg.norm(u,axis=1)>0
            cap=float(metrics['anderson_maximum_mixing'])
            record.update(parent_run=str(parent) if parent else None,signed_integrity='norms and state continuity verified',
                negative_residual_alignments=int(sum(cosine<0)),
                capped_applied_steps=int(sum((h['p']>=cap-1e-12)&applied)),
                applied_steps=int(sum(applied)),
                rms_doublings=int(sum(nr[1:]>2*nr[:-1])),
                stochastic_fallbacks=int(sum(h['stochastic_fallback'])),
                final_zone_squared_norm_fraction=[float(np.sum(rr[:,zones==z]**2)/nr[-1]**2) for z in (1,2,3)],
                final_worst_radius=float(layout[worst%n,2]),final_worst_component=names[worst//n],
                parent_launcher_error=json.loads((parent/'manifest.json').read_text())['launcher_error'] if parent else None)
            if parent:
                for ax in axes:
                    ax.axvline(80.5,color='0.6',lw=.8,ls=':')
        for tol in (2e-4,2e-5,1e-5,2e-6):
            ix=np.flatnonzero(h['map_max']<=tol)
            record['thresholds'][str(tol)]={'iteration':int(h['iteration'][ix[0]]),
                 'map_seconds_to_crossing':float(sum(h['map_seconds'][:ix[0]+1]))} if len(ix) else None
        result[label]=record
    axes[1].axhline(2e-5,color='0.4',ls=':',label='Daytime screening tolerance')
    axes[1].axhline(2e-6,color='0.6',ls='-.',label='Original tolerance')
    axes[1].legend();axes[1].set_xlabel('Map evaluation')
    target=ROOT/'docs/technical_note/figures/accelerator_comparison.png'
    fig.savefig(target,dpi=180);plt.close(fig)
    (ROOT/'work/signed-bb-analysis/accelerator_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
