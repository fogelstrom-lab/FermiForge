"""Read-only analysis of completed BB histories; write separate diagnostic artifacts."""
from pathlib import Path
import json
import os
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES=[
 ('2D 1–100','260928-232830-single-harmonic-0plus-bb-bb-radial-style'),
 ('2D 101–198','260929-161751-single-harmonic-0plus-bb-continued-bb-radial-style-101-200'),
 ('2D 199–298','260930-092558-single-harmonic-0plus-bb-tight-tolerance-2e-8'),
 ('Radial 1–266','261001-082221-a-phase-radial-symmetry-radial-reference')]
OUT=ROOT/'work/bb-spike-analysis'
OUT.mkdir(parents=True,exist_ok=True)

def spectral(y):
    x=np.arange(len(y))
    trend=np.polyval(np.polyfit(x,y,1),x)
    f=np.fft.rfftfreq(len(y))
    power=abs(np.fft.rfft((y-trend)*np.hanning(len(y))))**2
    keep=(f>=1/30)&(f<=.5)
    return f[keep],power[keep]/sum(power[keep])

def cosine(a,b):
    return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))

reports=[]
fig,axes=plt.subplots(4,3,figsize=(15,15),layout='constrained')
for row,(label,name) in enumerate(CASES):
    folder=ROOT/'runs'/name
    d=np.genfromtxt(folder/'iteration_history.dat',names=True)
    n=len(d); y=np.log10(d['map_max']); rms=np.log10(d['map_rms'])
    med=np.array([np.median(y[max(0,i-10):min(n,i+11)]) for i in range(n)])
    z=y-med
    peaks=np.where((z[1:-1]>z[:-2])&(z[1:-1]>=z[2:])&(z[1:-1]>.3))[0]+1
    f,p=spectral(y)
    ac=np.correlate(z-z.mean(),z-z.mean(),'full')[n-1:]/sum((z-z.mean())**2)
    rise=np.diff(y); loga=np.log10(np.maximum(d['p'],1e-30))
    before=d['p'][:-1]
    spike=rise>.3
    report=dict(label=label,source=str(folder),n=n,
                spike_definition='local maximum >0.3 dex above centered 21-point running median of log10(max residual)',
                peaks_local=(peaks+1).tolist(),peak_intervals=np.diff(peaks).tolist(),
                spectral_peaks=[dict(period=float(1/f[i]),power_fraction=float(p[i]))
                                for i in np.argsort(p)[-5:][::-1]],
                lag6_acf=float(ac[6]),lag12_acf=float(ac[12]),
                alpha_next_log_growth_correlation=float(np.corrcoef(loga[:-1],rise)[0,1]),
                doubled_error_steps=int(sum(spike)),doubled_after_alpha_ge50=int(sum(spike&(before>=50))),
                steps_alpha_ge50=int(sum(before>=50)),alpha_cap_hits=int(sum(before>=99.999)),
                median_six_step_max_ratio=float(np.median(10**(y[6:]-y[:-6]))))
    # Every snapshot contains unsigned component-norm diagnostics, not signed residual vectors.
    frames=[]
    for k in range(1,n+1):
        file=folder/f'sampled_residuals.dat.map{k:06d}.dat'
        a=np.loadtxt(file,ndmin=2)
        a=a[a[:,10]==1]
        if frames and not np.array_equal(a[:,:2],frames[0][:,:2]):
            raise ValueError('grid changed within segment')
        frames.append(a)
    radius=np.linalg.norm(frames[0][:,:2],axis=1)
    pattern=np.array([a[:,8] for a in frames])
    pointmax=np.array([a[:,9] for a in frames])
    rmax=radius[np.argmax(pointmax,axis=1)]
    fraction=(pattern[:,radius<3]**2).sum(axis=1)/(pattern**2).sum(axis=1)
    events=[]
    for k in peaks:
        if k<2 or k+3>=n: continue
        events.append(dict(iteration=int(k+1),pre_alpha=float(d['p'][k-1]),
                           spike_factor=float(10**(y[k]-y[k-1])),
                           after3_vs_before=float(10**(y[k+3]-y[k-1])),
                           pre_spike_shape_cosine=cosine(pattern[k-1],pattern[k]),
                           rmax_before=float(rmax[k-1]),rmax_spike=float(rmax[k]),
                           rmax_after3=float(rmax[k+3]),
                           core_fraction_before=float(fraction[k-1]),
                           core_fraction_spike=float(fraction[k])))
    report['events']=events
    report['median_adjacent_peak_shape_cosine']=float(np.median([
        cosine(pattern[a],pattern[b]) for a,b in zip(peaks[:-1],peaks[1:])]))
    report['median_peak_vs_pre_shape_cosine']=float(np.median([e['pre_spike_shape_cosine'] for e in events]))
    report['fraction_recovered_below_pre_spike_after3']=float(np.mean([e['after3_vs_before']<1 for e in events]))
    report['median_rmax_peak']=float(np.median(rmax[peaks]))
    report['median_rmax_nonpeak']=float(np.median(np.delete(rmax,peaks)))
    report['median_core_fraction_peak']=float(np.median(fraction[peaks]))
    report['median_core_fraction_nonpeak']=float(np.median(np.delete(fraction,peaks)))
    report['six_step_ratio_halves']=[float(np.median(10**(s[6:]-s[:-6]))) for s in np.array_split(y,2)]
    report['median_post3_vs_pre_spike']=float(np.median([e['after3_vs_before'] for e in events]))
    if label.startswith('Radial'):
        late=np.arange(n)>=n-100
        ispeak=np.isin(np.arange(n),peaks)
        report['late_radial']={}
        for tag,mask in [('spike',late&ispeak),('between',late&~ispeak)]:
            outer=(pattern[:,radius>=20]**2).sum(axis=1)/(pattern**2).sum(axis=1)
            report['late_radial'][tag]=dict(count=int(sum(mask)),
                median_maximum_radius=float(np.median(rmax[mask])),
                median_outer_nodal_squared_residual_fraction=float(np.median(outer[mask])))
    ax=axes[row,0]
    ax.plot(np.arange(1,n+1),d['map_max'],color='black',label='max')
    ax.plot(np.arange(1,n+1),d['map_rms'],color='#4b0082',label='RMS')
    ax.scatter(peaks+1,d['map_max'][peaks],marker='x',color='tab:orange',s=22)
    ax.set(yscale='log',xlabel='Local iteration',ylabel='Residual',title=label)
    ax.legend()
    ax=axes[row,1]
    ax.plot(1/f,p,'k-',label='max')
    fr,pr=spectral(rms);ax.plot(1/fr,pr,color='#4b0082',ls='--',label='RMS')
    ax.axvline(6,color='gray',ls=':')
    ax.set(xlim=(2,20),xlabel='Period [iterations]',ylabel='Fraction of band power',
           title='Detrended log-error spectrum (Hann window)')
    ax=axes[row,2]
    ax.plot(np.arange(1,n+1),rmax,'k-',lw=.8)
    ax.scatter(peaks+1,rmax[peaks],marker='x',color='tab:orange',s=22)
    ax.set(xlabel='Local iteration',ylabel='Radius of maximum residual [xi0]',
           title='Spatial motion of residual maximum')
    reports.append(report)
fig.savefig(OUT/'bb_spike_comparison.png',dpi=150)
plt.close(fig)
(OUT/'analysis.json').write_text(json.dumps(reports,indent=2))
for r in reports:
    print(r['label'],json.dumps({k:v for k,v in r.items() if k not in
        ('source','events','peaks_local','peak_intervals','spike_definition')},indent=2))
print('Artifacts:',OUT)
