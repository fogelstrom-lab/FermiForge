"""Analyse saved signed BB diagnostics without changing a run or its solver."""
import argparse
import json
import os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT/'work/signed-bb-analysis')
    args = parser.parse_args()
    run = args.run.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    h = np.genfromtxt(run/'iteration_history.dat', names=True)
    d = np.genfromtxt(run/'metrics.txt.bb_diagnostics.dat', names=True)
    layout = np.loadtxt(run/'metrics.txt.signed.layout.dat')
    files = sorted(run.glob('metrics.txt.signed.map*.dat'))
    assert len(files) == len(h) == len(d)
    arrays = np.stack([np.loadtxt(f) for f in files])
    x,r,u = (arrays[:,:,i] for i in (1,2,3))
    n = len(layout)
    assert x.shape[1] == 21*n
    np.testing.assert_array_equal(arrays[:,:,0],np.broadcast_to(np.arange(1,21*n+1),x.shape))
    nr = np.linalg.norm(r,axis=1)
    np.testing.assert_allclose(np.sqrt(np.mean(r*r,axis=1)),h['map_rms'],rtol=1e-12)
    np.testing.assert_allclose(np.max(abs(r),axis=1),h['map_max'],rtol=1e-12)
    np.testing.assert_allclose(u,h['p'][:,None]*r,atol=1e-15,rtol=1e-8)
    np.testing.assert_allclose(x[1:],x[:-1]+u[:-1],atol=1e-14,rtol=1e-12)
    s = x[1:]-x[:-1]
    y = r[:-1]-r[1:]
    sy = np.sum(s*y,axis=1)
    bb1 = np.sum(s*s,axis=1)/sy
    bb2 = sy/np.sum(y*y,axis=1)
    valid=d['secant_valid'][1:]>0
    np.testing.assert_allclose(bb1[valid],d['raw_BB1'][1:][valid],rtol=1e-10)
    np.testing.assert_allclose(bb2[valid],d['raw_BB2'][1:][valid],rtol=1e-10)
    cosine=np.sum(r[1:]*r[:-1],axis=1)/(nr[1:]*nr[:-1])
    np.testing.assert_allclose(cosine,d['residual_alignment'][1:],atol=1e-12)
    radius=layout[:,2]
    zones=layout[:,3].astype(int)
    rr=r.reshape(-1,21,n)
    xx=x.reshape(-1,21,n)
    energy=np.stack([np.sum(rr[:,:,zones==z]**2,axis=(1,2)) for z in (1,2,3)],axis=1)
    fraction=energy/np.sum(energy,axis=1)[:,None]
    worst=np.argmax(abs(r),axis=1)
    worst_radius=radius[worst%n]
    names=[f'{part} A{a}{b}' for a in 'xyz' for b in 'xyz' for part in ('Re','Im')]+['mf_x','mf_y','mf_z']
    growth=nr[1:]/nr[:-1]
    mxgrowth=h['map_max'][1:]/h['map_max'][:-1]
    alpha=h['p'][:-1]
    report={'source':str(run),'maps':len(h),'integrity':'signed vectors agree with history, secants and applied steps',
            'final_rms':float(h['map_rms'][-1]),'final_max':float(h['map_max'][-1]),
            'negative_secants':int(np.sum(sy<0)),
            'negative_residual_alignment':int(np.sum(cosine<0)),
            'alpha_cap_count':int(np.sum(alpha>=99.999)),
            'formula_counts':{str(int(f)):int(np.sum(d['formula']==f)) for f in np.unique(d['formula'])},
            'status_counts':{str(int(f)):int(np.sum(d['status']==f)) for f in np.unique(d['status'])},
            'final_zone_energy_fraction':fraction[-1].tolist(),
            'final_worst_radius':float(worst_radius[-1]),'final_worst_component':names[worst[-1]//n],
            'threshold_first_crossings':{str(t):int(h['iteration'][np.flatnonzero(h['map_max']<=t)[0]])
                  if np.any(h['map_max']<=t) else None for t in (2e-4,2e-5,1e-5,2e-6)},
            'windows':[], 'steps':{}}
    for label,mask in [('all',np.ones(len(alpha),dtype=bool)),('alpha_ge50',alpha>=50),
                       ('alpha_lt50',alpha<50),('BB1',d['formula'][:-1]==1),('BB2',d['formula'][:-1]==2)]:
        report['steps'][label]={'count':int(sum(mask)),'median_rms_ratio':float(np.median(growth[mask])),
             'rms_doubles':int(sum(growth[mask]>2)),'max_doubles':int(sum(mxgrowth[mask]>2)),
             'mean_log10_rms_ratio':float(np.mean(np.log10(growth[mask])))}
    for lo,hi in [(1,100),(101,250),(251,len(h))]:
        mask=(h['iteration']>=lo)&(h['iteration']<=hi)
        slope=np.polyfit(h['iteration'][mask],np.log10(h['map_rms'][mask]),1)[0]
        counts=np.bincount(worst[mask]//n,minlength=21)
        report['windows'].append({'range':[lo,hi],'fitted_iterations_per_decade':float(-1/slope),
              'median_zone_energy_fraction':np.median(fraction[mask],axis=0).tolist(),
              'worst_radius_median':float(np.median(worst_radius[mask])),
              'max_component_counts':{names[j]:int(v) for j,v in enumerate(counts) if v},
              'median_raw_BB1':float(np.median(d['raw_BB1'][mask])),
              'median_raw_BB2':float(np.median(d['raw_BB2'][mask]))})
    # Small singular subspace of late unit residuals, not a Jacobian spectrum.
    late=h['iteration']>=251
    unit=r[late]/nr[late,None]
    eig,vec=np.linalg.eigh(unit@unit.T)
    eig=eig[::-1]; vec=vec[:,::-1]
    modes=vec[:,:3].T@unit/np.sqrt(eig[:3,None])
    report['late_unit_residual_svd_energy']=(eig[:8]/sum(eig)).tolist()
    report['late_modes_zone_energy']=[np.sum(mode.reshape(21,n)**2,axis=0)[zones==z].sum().item()
                                      for mode in modes for z in (1,2,3)]
    report['late_modes_components']=[{names[j]:float(value) for j,value in enumerate(
        np.sum(mode.reshape(21,n)**2,axis=1)) if value>.02} for mode in modes]
    # Compare the leading pattern with an infinitesimal spin rotation about y.
    # This is a geometric overlap, not proof of a symmetry or zero eigenvalue.
    a=(xx[-1,:18:2]+1j*xx[-1,1:18:2]).reshape(3,3,n)
    da=np.stack([a[2],np.zeros_like(a[1]),-a[0]])
    tangent=np.zeros((21,n))
    tangent[:18:2]=da.reshape(9,n).real
    tangent[1:18:2]=da.reshape(9,n).imag
    lead=modes[0].reshape(21,n)
    outside=zones==3
    report['leading_pattern_spin_y_rotation_overlap']={
        'absolute_cosine_all':float(abs(np.sum(tangent*lead)/np.linalg.norm(tangent))),
        'absolute_cosine_outer':float(abs(np.sum(tangent[:,outside]*lead[:,outside])/
                  np.linalg.norm(tangent[:,outside])/np.linalg.norm(lead[:,outside]))),
        'qualification':'Geometric overlap with infinitesimal spin rotation of final A; not a proven zero mode.'}
    region_alignment={}
    for z,label in enumerate(('core','transition','outer'),1):
        rz=rr[:,:,zones==z].reshape(len(h),-1)
        cz=np.sum(rz[1:]*rz[:-1],axis=1)/(np.linalg.norm(rz[1:],axis=1)*np.linalg.norm(rz[:-1],axis=1))
        region_alignment[label]={'negative_count':int(sum(cz<0)),
            'late_median':float(np.median(cz[250:])),
            'late_negative_count':int(sum(cz[250:]<0))}
    report['regional_residual_alignment']=region_alignment
    logres=np.log10(h['map_rms'][late]); t=np.arange(len(logres))
    signal=logres-np.polyval(np.polyfit(t,logres,1),t)
    frequency=np.fft.rfftfreq(len(signal)); power=abs(np.fft.rfft(signal*np.hanning(len(signal))))**2
    keep=(frequency>=1/30)&(frequency<=.5)
    peaks=np.argsort(power[keep])[-3:][::-1]
    report['late_detrended_rms_spectral_periods']=(1/frequency[keep][peaks]).tolist()
    report['largest_rms_growth_events']=[{'map_before':int(h['iteration'][i]),
         'map_after':int(h['iteration'][i+1]),'alpha':float(alpha[i]),
         'formula':int(d['formula'][i]),'rms_ratio':float(growth[i]),
         'residual_cosine':float(cosine[i]),'worst_radius_after':float(worst_radius[i+1])}
         for i in np.argsort(growth)[-5:][::-1]]
    # Effective response along each observed secant: Rayleigh quotient of I-J_F.
    # Local, direction-dependent and not an eigenvalue estimate.
    report['secant_rayleigh_range']=np.quantile(sy/np.sum(s*s,axis=1),[0,.25,.5,.75,1]).tolist()
    previous=ROOT/'runs/261001-130107-radial-bb-scratch-0plus-T030-Fs54-radial-reference'
    if (previous/'iteration_history.dat').exists():
        old=np.genfromtxt(previous/'iteration_history.dat',names=True)
        report['previous_scratch_comparison']={'maps':len(old),'same_length':len(old)==len(h)}
        if len(old)==len(h):
            report['previous_scratch_comparison']['max_absolute_residual_history_difference']=float(np.max(abs(old['map_max']-h['map_max'])))
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,ax=plt.subplots(3,1,figsize=(12,10),layout='constrained')
    it=h['iteration']
    for key,label,color,style in [('map_rms','RMS','#2166ac','-'),('map_max','Maximum','#b2182b','--')]:
        ax[0].semilogy(it,h[key],label=label,color=color,ls=style)
    ax[0].set_ylabel('Fixed-point residual');ax[0].legend()
    ax[1].plot(it,d['residual_alignment'],label='Residual alignment',color='#2166ac')
    ax[1].plot(it,d['signed_cosine'],label='Signed secant cosine',color='#b2182b',alpha=.6)
    ax[1].axhline(0,color='0.5',lw=.6);ax[1].set_ylim(-1.05,1.05);ax[1].legend()
    ax[2].semilogy(it[1:],abs(d['raw_BB1'][1:]),label='|raw BB1|',alpha=.7)
    ax[2].semilogy(it[1:],abs(d['raw_BB2'][1:]),label='|raw BB2|',color='#b2182b',ls='--',alpha=.7)
    ax[2].semilogy(it[:-1],alpha,label='Applied alpha',color='black',lw=.8)
    ax[2].axhline(100,color='0.5',ls='--');ax[2].legend();ax[2].set_ylabel('Step estimate');ax[2].set_xlabel('Map iteration')
    fig.savefig(out/'signed_bb_dynamics.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(3,1,figsize=(12,10),layout='constrained')
    for j,z in enumerate(('core','transition','outer')):
        color=('#2166ac','#b2182b','black')[j]; style=('-','--',':')[j]
        ax[0].semilogy(it,d[z+'_rms'],label=z,color=color,ls=style)
        ax[1].plot(it,fraction[:,j],label=z,color=color,ls=style)
    ax[0].set_ylabel('Regional RMS');ax[0].legend()
    ax[1].set_ylabel('Fraction of residual norm squared');ax[1].legend()
    ax[2].plot(it,worst_radius,'.',color='black',ms=3);ax[2].set_ylabel('Radius of largest component')
    ax[2].axhline(3,color='0.5',ls=':');ax[2].axhline(15,color='0.5',ls='--');ax[2].set_xlabel('Map iteration')
    fig.savefig(out/'signed_bb_regions.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(2,1,figsize=(11,7),layout='constrained')
    for j,mode in enumerate(modes[:2]):
        profile=np.sqrt(np.sum(mode.reshape(21,n)**2,axis=0))
        ax[0].plot(radius,profile,label=f'Pattern {j+1}: {100*eig[j]/sum(eig):.1f}% of unit-residual content',
                   color=('#2166ac','#b2182b')[j],ls=('-','--')[j])
        ax[1].plot(it[late],unit@mode,label=f'Pattern {j+1}',color=('#2166ac','#b2182b')[j],ls=('-','--')[j])
    ax[0].set_xlabel('Radius / xi0');ax[0].set_ylabel('Pattern magnitude per node');ax[0].legend()
    ax[1].set_xlabel('Map iteration');ax[1].set_ylabel('Signed projection of unit residual');ax[1].legend()
    fig.suptitle('Late residual patterns (maps 251–433); not Jacobian eigenmodes')
    fig.savefig(out/'signed_bb_patterns.png',dpi=160)
    fig.savefig(ROOT/'docs/technical_note/figures/signed_bb_patterns.png',dpi=160)
    plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
