#!/usr/bin/env python3
"""Reproducible implementation schematics; not converged solution data."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np

OUT = ROOT / 'docs/technical_note/figures'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13})


def smooth(width, cells, stretch=2.3):
    return width*np.sinh(stretch*np.linspace(-1, 1, cells+1))/np.sinh(stretch)


def extend(core, width, growth=1.15):
    n = 1
    while (core[-1]-core[-2])*(growth**n-1)/(growth-1) < width-core[-1]:
        n += 1
    step = (width-core[-1]) / ((growth**n-1)/(growth-1))
    outer = core[-1] + np.cumsum(step*growth**np.arange(n))
    outer[-1] = width
    return np.r_[-outer[::-1], core, outer]


def grids():
    fig, axes = plt.subplots(2, 1, figsize=(7, 10), layout='constrained')
    for ax, coords, radii, title in [
        (axes[0], smooth(22,24), (14,17,20), 'Smooth grid: compact core calculation'),
        (axes[1], extend(smooth(20,24),60), (40,52,58), 'Extended grid: preserve inner coordinates')]:
        x,y=np.meshgrid(coords,coords); active=x*x+y*y<=radii[2]**2
        ax.scatter(x[~active],y[~active],s=5,facecolors='none',edgecolors='.6',lw=.4)
        ax.scatter(x[active],y[active],s=5,c='#3f007d')
        phi=np.linspace(0,2*np.pi,500)
        for r,style in zip(radii,(':','--','-')):
            ax.plot(r*np.cos(phi),r*np.sin(phi),style,c='.2',lw=1.4)
        if ax is axes[1]:
            ax.plot([-20,20,20,-20,-20],[-20,-20,20,20,-20],c='#0072b2',lw=1.7)
        ax.set(aspect='equal',title=title,xlabel=r'$x/\xi_0$',ylabel=r'$y/\xi_0$')
        ax.set_title(title+'\nFit radii: '+str(radii[:2])+r' $\xi_0$; active: '+str(radii[2])+r' $\xi_0$',fontsize=11)
    handles=[Line2D([],[],marker='o',ls='',color='#3f007d',label='Updated nodes'),
             Line2D([],[],marker='o',ls='',markerfacecolor='none',color='.6',label='Dependent halo'),
             Line2D([],[],color='.2',ls=':',label='Inner fit'),
             Line2D([],[],color='.2',ls='--',label='Matching'),
             Line2D([],[],color='.2',ls='-',label='Active boundary')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Every dot is a stored node; circles are masks, not mesh lines\n'
                 'Reduced 24-cell base grids (production uses 64)',fontsize=12)
    fig.savefig(OUT/'grid_construction.pdf'); fig.savefig(OUT/'grid_construction.png',dpi=180)
    plt.close(fig)


def seeds():
    x=np.linspace(-8,8,161); X,Y=np.meshgrid(x,x)
    r=np.hypot(X,Y); phi=np.arctan2(Y,X); s=r/(3/np.sqrt(1-.3**2))
    f=np.ones_like(s); np.divide(np.tanh(s),s,out=f,where=s!=0)
    zero=np.zeros_like(f)
    fields=[(zero,zero),(f*np.cos(phi)**2,f*np.sin(phi)**2),
            (f*np.cos(phi)**2,zero),(f,-f)]
    names=['Normal: nop','A-phase: aop','Double-core: dop','Quadrupole: qop']
    cmap=LinearSegmentedColormap.from_list('bwr_seed',['#2166ac','white','#b2182b'])
    fig,axes=plt.subplots(2,4,figsize=(12,7),layout='constrained')
    for col,(name,pair) in enumerate(zip(names,fields)):
        for row,values in enumerate(pair):
            ax=axes[row,col]
            im=ax.pcolormesh(x,x,values,cmap=cmap,vmin=-1,vmax=1,shading='nearest',rasterized=True)
            ax.set(aspect='equal',xlabel=r'$x/\xi_0$',ylabel=r'$y/\xi_0$')
            ax.set_title(name if row==0 else '')
            ax.set_xticks([-8,0,8]);ax.set_yticks([-8,0,8])
    fig.suptitle(r'Initial core filling: top $\mathrm{Re}\,A_{zx}/\Delta_B$; bottom $\mathrm{Im}\,A_{zy}/\Delta_B$' '\n'
                 r'Exact seed formulas, before halo regeneration; $T/T_c=0.30$, $F_1^s=0$',fontsize=14)
    fig.colorbar(im,ax=axes,location='bottom',shrink=.6,label='Signed component: blue < 0, white = 0, red > 0')
    fig.savefig(OUT/'core_initializations.pdf');fig.savefig(OUT/'core_initializations.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    grids();seeds()
