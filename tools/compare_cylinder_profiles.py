#!/usr/bin/env python3
"""Compare modern and new_src cylinder profiles on the same uniform radial grid."""
import argparse
import json
from pathlib import Path
import numpy as np

def compare(modern, reference):
    a=np.loadtxt(modern/'cylinder_final_op_xyz')
    b=np.loadtxt(reference/'op_xyz')
    if a.shape!=b.shape or a.shape[1]!=19:
        raise ValueError('incompatible radial profile shapes')
    # new_src prints coordinates with three decimal places. Compare fields at
    # their common grid indices, not at falsely distinct rounded positions.
    if not np.allclose(a[:,0],b[:,0],rtol=0,atol=5.1e-4):
        raise ValueError('radial coordinates differ')
    if not np.allclose(a[:,0],np.linspace(a[0,0],a[-1,0],len(a)),rtol=0,atol=1e-10):
        raise ValueError('comparison expects a uniform cylinder grid')
    x=a[:,1::2]+1j*a[:,2::2]; y=b[:,1::2]+1j*b[:,2::2]
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError('nonfinite profile')
    difference=np.abs(x-y); idx=np.unravel_index(np.argmax(difference),difference.shape)
    names=['xx','xy','xz','yx','yy','yz','zx','zy','zz']
    report=dict(modern=str(modern),reference=str(reference),
                maximum_absolute_gap_difference=float(difference[idx]),
                relative_l2_gap_difference=float(np.linalg.norm(x-y)/np.linalg.norm(y)),
                maximum_component=names[idx[1]],maximum_radius=float(a[idx[0],0]))
    return a[:,0],x,y,report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('modern',type=Path); p.add_argument('reference',type=Path)
    p.add_argument('--previous',type=Path)
    args=p.parse_args()
    modern=args.modern.resolve(); reference=args.reference.resolve()
    r,x,y,report=compare(modern,reference)
    if args.previous:
        _,_,_,report['previous']=compare(args.previous.resolve(),reference)
    dest=modern/'comparison'; dest.mkdir(exist_ok=True)
    (dest/'radial_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for ax,bx,k,name in zip(axes[0],axes[1],(0,4,8),('xx','yy','zz')):
        ax.plot(r,y[:,k].real,color='black',label='new_src')
        ax.plot(r,x[:,k].real,color='#332288',ls='--',label='modern continued')
        ax.set(title=f'Re A{name}',xlabel=r'$r/\xi_0$'); ax.legend(); ax.grid(alpha=.2)
        bx.semilogy(r,np.maximum(abs(x[:,k]-y[:,k]),1e-18),color='#332288',label='continued')
        if args.previous:
            _,prior,_,_=compare(args.previous.resolve(),reference)
            bx.semilogy(r,np.maximum(abs(prior[:,k]-y[:,k]),1e-18),color='black',ls='--',label='original')
        bx.set(xlabel=r'$r/\xi_0$',ylabel='Absolute complex-component difference')
        bx.legend(); bx.grid(alpha=.2)
    fig.savefig(dest/'radial_comparison.png',dpi=160); plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
