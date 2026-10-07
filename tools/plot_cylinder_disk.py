#!/usr/bin/env python3
"""Reconstruct a radial cylinder solution on a disk; no new physics solve."""
import argparse
import json
from pathlib import Path
import numpy as np
from plot_2d_fields import (FieldMap, cartesian_to_harmonics, plot_amplitudes,
    plot_phase_cosines, plot_basis_amplitudes, plot_basis_phase_cosines,
    plot_density_and_current, HARMONIC_COMPONENTS)
import matplotlib.pyplot as plt


def interpolate(radius, values, query):
    """Continuous average of bracketing quadratics, matching sample_radial_state."""
    if len(radius)<3 or np.any(np.diff(radius)<=0):
        raise ValueError('need at least three strictly increasing radial nodes')
    cell=np.clip(np.searchsorted(radius,query,side='right')-1,0,len(radius)-2)
    left=np.maximum(0,cell-1); right=np.minimum(cell,len(radius)-3)
    result=np.zeros(query.shape+values.shape[1:],dtype=values.dtype)
    extra=(None,)*(values.ndim-1)
    for start in (left,right):
        for offset in range(3):
            i=start+offset; weight=np.ones(query.shape)
            for other in range(3):
                if other!=offset:
                    j=start+other
                    weight *= (query-radius[j])/(radius[i]-radius[j])
            result += .5*weight[(...,)+extra]*values[i]
    return result


def reconstruct(radius, gap, current, x, y, winding):
    xx,yy=np.meshgrid(x,y); r=np.hypot(xx,yy); theta=np.arctan2(yy,xx)
    inside=r<=radius[-1]+1e-12
    axis=interpolate(radius,gap,np.minimum(r,radius[-1]))
    flow=interpolate(radius,current,np.minimum(r,radius[-1]))
    c=np.cos(theta); s=np.sin(theta)
    rotation=np.zeros(r.shape+(3,3)); rotation[...,0,0]=c; rotation[...,1,1]=c
    rotation[...,0,1]=-s; rotation[...,1,0]=s; rotation[...,2,2]=1
    # Equivalent to C_sk(phi)=exp[i(m-s-k)phi] C_sk(0).
    a=np.einsum('...ij,...jk,...lk->...il',rotation,axis,rotation)
    a *= np.exp(1j*winding*theta)[...,None,None]
    flow=np.einsum('...ij,...j->...i',rotation,flow)
    a[~inside]=np.nan; flow[~inside]=np.nan
    return FieldMap(x,y,np.moveaxis(a,(-2,-1),(0,1)),
                    np.sum(abs(a)**2,axis=(-2,-1))/3,np.moveaxis(flow,-1,0),
                    {'current_kind':'current_related_mean_field'})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    parser.add_argument('--points',type=int,default=401)
    args=parser.parse_args(); run=args.run.resolve()
    if args.points<3 or args.points%2==0:
        parser.error('--points must be odd and at least 3')
    meta=json.loads((run/'run_metadata.json').read_text())
    winding=1 if meta['initialization']==-2 else meta['winding']
    data=np.loadtxt(run/'cylinder_final_op_xyz'); cur=np.loadtxt(run/'cylinder_final_curr')
    if not np.allclose(data[:,0],cur[:,0]): raise ValueError('profile radii differ')
    gap=(data[:,1::2]+1j*data[:,2::2]).reshape(-1,3,3)
    x=np.linspace(-data[-1,0],data[-1,0],args.points)
    field=reconstruct(data[:,0],gap,cur[:,2:5],x,x,winding)
    h=cartesian_to_harmonics(field.order_parameter)
    dest=run/'plots'/'disk'; dest.mkdir(parents=True,exist_ok=True)
    plots=[('order_parameter_amplitude',plot_amplitudes(field)),
           ('order_parameter_phase_cosine',plot_phase_cosines(field,1e-3)),
           ('order_parameter_harmonic_amplitude',plot_basis_amplitudes(x,x,h,HARMONIC_COMPONENTS,'Harmonic amplitudes: symmetry-reconstructed disk')),
           ('order_parameter_harmonic_phase_cosine',plot_basis_phase_cosines(x,x,h,HARMONIC_COMPONENTS,1e-3,'Harmonic phase cosines: symmetry-reconstructed disk')),
           ('density_and_current',plot_density_and_current(field,23,1e-10))]
    for name,figure in plots:
        for ax in figure.axes:
            if ax.get_xlabel()==r'$x/\xi_0$' or ax.get_aspect()==1:
                ax.add_patch(plt.Circle((0,0),data[-1,0],fill=False,color='black',lw=.5))
        figure.savefig(dest/(name+'.png'),dpi=180); plt.close(figure)
    mid=len(x)//2
    for name,a,labels in [('axis_profiles',field.order_parameter,('x','y','z')),
                         ('harmonic_axis_profiles',h,('+','0','-'))]:
        fig,axes=plt.subplots(3,3,figsize=(13,10),constrained_layout=True)
        for i in range(3):
            for j in range(3):
                ax=axes[i,j]
                for v,label,color in [(a[i,j,mid,:],'x axis','#332288'),(a[i,j,:,mid],'y axis','#aa6600')]:
                    ax.plot(x,v.real,color=color,label=label+' Re')
                    ax.plot(x,v.imag,color=color,ls='--',label=label+' Im')
                ax.set(title=f'{labels[i]}{labels[j]}',xlabel=r'coordinate / $\xi_0$')
                ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=8)
        fig.suptitle('Symmetry-reconstructed cylinder: '+name.replace('_',' '))
        fig.savefig(dest/(name+'.png'),dpi=180); plt.close(fig)
    (dest/'reconstruction.json').write_text(json.dumps(dict(source=str(run),
        symmetry_winding=winding,plotting_points_per_axis=args.points,
        interpolation='blended quadratic, same as radial sampler',
        outside_disk='masked; no extrapolation',independent_2d_solve=False),indent=2)+'\n')
    print(dest)


if __name__=='__main__': main()
