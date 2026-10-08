#!/usr/bin/env python3
"""Archive a radially constrained cylinder run using the modern Fortran stack."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import math
from run_and_plot import find_program, RunnerError

ROOT = Path(__file__).resolve().parents[1]

def restart_files(directory):
    directory = Path(directory).expanduser().resolve()
    for names in (('cylinder_final_op_xyz','cylinder_final_curr'),('op_xyz','curr')):
        paths=tuple(directory/name for name in names)
        if all(path.is_file() for path in paths):
            return paths
    raise RunnerError(f'no complete modern or new_src restart profile pair in {directory}')

def resolve_toolchain():
    """Respect PATH first, then find common Mac installations without shell edits."""
    prefixes = (Path('/opt/homebrew/bin'), Path('/opt/local/bin'), Path('/usr/local/bin'))
    cmake = find_program('cmake', tuple(prefix/'cmake' for prefix in prefixes))
    mpirun = find_program('mpirun', tuple(prefix/'mpirun' for prefix in prefixes))
    env = os.environ.copy()
    # CMake and MPI wrappers must also be able to locate compiler executables.
    # Preserve the user's PATH precedence and wrapper symlink names.
    entries = env.get('PATH', '').split(os.pathsep)
    for prefix in (cmake.parent, mpirun.parent, *prefixes):
        if prefix.is_dir() and str(prefix) not in entries:
            entries.append(str(prefix))
    env['PATH'] = os.pathsep.join(entries)
    return str(cmake), str(mpirun), env

def plot(run):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    meta=json.loads((run/'run_metadata.json').read_text())
    if meta.get('spatial_mode')=='full_2d':
        plot_disk(run,meta)
        if meta.get('texture','none')!='none':
            from plot_cylinder_texture import plot_texture
            plot_texture(run,meta)
        return
    from run_and_plot import (write_complex_grid_matplotlib, write_current_matplotlib,
                              CARTESIAN_COMPONENTS, HARMONIC_COMPONENTS)
    from plot_2d_fields import cartesian_to_harmonics
    dest=run/'plots'; dest.mkdir(exist_ok=True)
    final=np.loadtxt(run/'cylinder_final_op_xyz')
    a=(final[:,1::2]+1j*final[:,2::2]).reshape(-1,3,3)
    h=np.stack([cartesian_to_harmonics(value) for value in a]).reshape(-1,9)
    table=np.empty((len(final),19)); table[:,0]=final[:,0]
    table[:,1::2]=h.real; table[:,2::2]=h.imag
    np.savetxt(run/'cylinder_final_op_harm',table)
    write_complex_grid_matplotlib(run/'cylinder_final_op_xyz',dest/'cartesian_profiles.png',
                                  CARTESIAN_COMPONENTS,'Cylinder: Cartesian components')
    write_complex_grid_matplotlib(run/'cylinder_final_op_harm',dest/'harmonic_profiles.png',
                                  HARMONIC_COMPONENTS,'Cylinder: harmonic components')
    write_current_matplotlib(run/'cylinder_final_curr',dest/'amplitude_and_current_field.png')
    history=np.atleast_2d(np.loadtxt(run/'cylinder_history.dat'))
    fig,ax=plt.subplots()
    ax.semilogy(history[:,0],np.maximum(history[:,1],1e-30),label='RMS')
    ax.semilogy(history[:,0],np.maximum(history[:,2],1e-30),label='maximum')
    ax.set(xlabel='Evaluated state (0 = initial)',ylabel='Absolute fixed-point residual',
           title='Modern cylinder: Anderson iteration')
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(dest/'convergence.png',dpi=160); plt.close(fig)
    fig,ax=plt.subplots()
    for file,label,style in [('cylinder_initial_op_xyz','initial','--'),
                            ('cylinder_final_op_xyz','final','-'),('cylinder_mapped_op_xyz','map(final)',':')]:
        data=np.loadtxt(run/file)
        for c,name in enumerate(('xx','yy','zz')):
            ax.plot(data[:,0],data[:,1+8*c],style,color=('#332288','#000000','#AA6600')[c],
                    label=f'{label} Re A{name}')
    ax.set(xlabel=r'$r/\xi_0$',ylabel='Order parameter',title='Cylinder diagonal profiles')
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(dest/'initial_final.png',dpi=160); plt.close(fig)

def plot_disk(run,meta):
    """Actual independently iterated 2D fields, not a radial reconstruction."""
    if meta.get('disk_grid')=='annular':
        from plot_annular_disk import plot_annular_disk
        plot_annular_disk(run,meta)
        return
    import numpy as np
    from dataclasses import replace
    from plot_2d_fields import (load_field_map,cartesian_to_harmonics,plot_amplitudes,
        plot_phase_cosines,plot_basis_amplitudes,plot_basis_phase_cosines,
        plot_density_and_current,HARMONIC_COMPONENTS)
    import matplotlib.pyplot as plt
    field=load_field_map(run/'cylinder_final_fields_2d.dat')
    xx,yy=np.meshgrid(field.x,field.y); outside=np.hypot(xx,yy)>meta['radius']+1e-12
    field.order_parameter[:,:,outside]=np.nan
    field.current[:,outside]=np.nan; field.pair_density[outside]=np.nan
    field=replace(field,metadata={**field.metadata,'current_kind':'current_related_mean_field'})
    h=cartesian_to_harmonics(field.order_parameter)
    dest=run/'plots'; dest.mkdir(exist_ok=True)
    initial=load_field_map(run/'cylinder_initial_fields_2d.dat')
    mapped=load_field_map(run/'cylinder_mapped_fields_2d.dat')
    delta=field.order_parameter-initial.order_parameter
    reference_norm=np.linalg.norm(initial.order_parameter[:,:,~outside])
    gap_res=mapped.order_parameter-field.order_parameter
    cur_res=mapped.current-field.current
    rms=np.sqrt((np.sum(abs(gap_res)**2,axis=(0,1))+np.sum(cur_res**2,axis=0))/21)
    maximum=np.maximum(np.max(np.maximum(abs(gap_res.real),abs(gap_res.imag)),axis=(0,1)),
                       np.max(abs(cur_res),axis=0))
    peak=np.unravel_index(np.nanargmax(maximum),maximum.shape)
    comparison={'relative_l2_gap_change_from_initial':float(np.linalg.norm(delta[:,:,~outside])/reference_norm),
        'maximum_gap_change_from_initial':float(np.nanmax(abs(delta))),
        'maximum_residual_location':[float(field.x[peak[1]]),float(field.y[peak[0]])],
        'active_points':int(np.sum(~outside)),
        'initial_is_imported_radial_reference':bool(meta.get('initial_run'))}
    (run/'disk_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    figs=[('order_parameter_amplitude',plot_amplitudes(field)),
          ('order_parameter_phase_cosine',plot_phase_cosines(field,1e-3)),
          ('harmonic_amplitude',plot_basis_amplitudes(field.x,field.y,h,HARMONIC_COMPONENTS,'Independent 2D disk: harmonics')),
          ('harmonic_phase',plot_basis_phase_cosines(field.x,field.y,h,HARMONIC_COMPONENTS,1e-3,'Independent 2D disk: harmonic phases')),
          ('density_and_current',plot_density_and_current(field,20,1e-10))]
    history=np.atleast_2d(np.loadtxt(run/'cylinder_history.dat'))
    fig,ax=plt.subplots()
    ax.semilogy(history[:,0],np.maximum(history[:,1],1e-30),label='RMS')
    ax.semilogy(history[:,0],np.maximum(history[:,2],1e-30),label='maximum')
    ax.set(xlabel='Evaluated state',ylabel='Absolute residual',title='Unconstrained 2D cylinder')
    ax.legend(); ax.grid(alpha=.2); fig.tight_layout(); figs.append(('convergence',fig))
    fig,axes=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
    for ax,values,title in zip(axes,(rms,maximum),('Point RMS residual','Maximum component residual')):
        image=ax.pcolormesh(field.x,field.y,values,shading='auto',cmap='Purples')
        ax.set(title=title,xlabel='x / xi0',ylabel='y / xi0',aspect='equal')
        fig.colorbar(image,ax=ax)
    figs.append(('residual_maps',fig))
    mid=len(field.x)//2
    for basis,final_values,start_values in [('cartesian',field.order_parameter,initial.order_parameter),
            ('harmonic',h,cartesian_to_harmonics(initial.order_parameter))]:
        fig,axes=plt.subplots(3,3,figsize=(13,10),constrained_layout=True)
        for i in range(3):
            for j in range(3):
                ax=axes[i,j]
                for values,tag,alpha in [(start_values,'initial',.35),(final_values,'final',1.)]:
                    for a,label,color in [(values[i,j,mid,:],'x','#332288'),(values[i,j,:,mid],'y','#aa6600')]:
                        ax.plot(field.x,a.real,color=color,alpha=alpha,label=f'{tag} {label} Re')
                        ax.plot(field.x,a.imag,color=color,alpha=alpha,ls='--',label=f'{tag} {label} Im')
                ax.set(title=f'{basis} ({i+1},{j+1})',xlabel='coordinate / xi0'); ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=6,ncol=2)
        figs.append((basis+'_axis_profiles',fig))
    for name,fig in figs:
        fig.savefig(dest/(name+'.png'),dpi=170); plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ranks',type=int,default=10)
    p.add_argument('--iterations',type=int,default=0,help='AA updates; 0 is one map')
    p.add_argument('--tolerance',type=float,default=2e-6,help='absolute maximum residual threshold')
    p.add_argument('--initialization',type=int,choices=(-1,-2,0,1),default=-1,
                   help='-1 bulk B, -2 bulk A (no vortex), 0 normal-core, 1 A-core (unit winding)')
    p.add_argument('--fs1',type=float,default=0.,help='Landau F_1^s; converted internally to A_1^s')
    p.add_argument('--texture',choices=('none','a_mermin_ho','a_planar','a_panam'),default='none',
                   help='full-2D A seed; a_planar is the original mixed A/polar trial, not the standard planar phase')
    p.add_argument('--spatial-mode',choices=('radial_symmetry','full_2d'),default='radial_symmetry')
    p.add_argument('--disk-cells',type=int,default=32,help='uniform cells across the full disk diameter')
    p.add_argument('--disk-grid',choices=('cartesian','annular'),default='cartesian')
    p.add_argument('--disk-radial-layout',choices=('wall','core_wall'),default='wall',
                   help='annular clustering at wall only, or at both core and wall')
    p.add_argument('--disk-rings',type=int,default=20,help='annular radial intervals, excluding centre node')
    p.add_argument('--disk-stretch',type=float,default=2.,help='0 uniform radii; positive increases selected clustering')
    p.add_argument('--disk-tangent-spacing',type=float,default=1.,help='target arc spacing in xi0')
    p.add_argument('--perturbation',type=float,default=0.,help='2D Axx perturbation amplitude relative to bulk gap')
    p.add_argument('--restart-2d',type=Path,help='same-grid full 2D cylinder checkpoint')
    p.add_argument('--azimuths',type=int,default=32)
    p.add_argument('--radial-points',type=int,default=100)
    p.add_argument('--radius',type=float,default=10.)
    p.add_argument('--step',type=float,default=10/99)
    p.add_argument('--half-length',type=float,default=8000/99)
    p.add_argument('--collision-aligned',action='store_true')
    p.add_argument('--initial-run',type=Path,help='modern cylinder or new_src archive; fresh Anderson history')
    p.add_argument('--build-dir',type=Path,default=ROOT/'work/modern-cylinder-build')
    p.add_argument('--plot-only',type=Path)
    a=p.parse_args()
    if a.plot_only:
        plot(a.plot_only.resolve()); return
    if a.texture!='none' and (a.spatial_mode!='full_2d' or a.initialization!=-2 or a.initial_run or a.restart_2d):
        p.error('texture requires full_2d, initialization=-2, and no restart/reference')
    if min(a.ranks,a.azimuths,a.radial_points-2)<=0 or min(a.radius,a.step,a.half_length)<=0 or a.iterations<0:
        p.error('positive counts/lengths required; iterations must be nonnegative')
    if not math.isfinite(a.tolerance) or a.tolerance<=0:
        p.error('--tolerance must be finite and positive')
    if not math.isfinite(a.fs1) or a.fs1<0:
        p.error('--fs1 must be finite and nonnegative')
    if a.disk_cells<8 or a.disk_cells%2 or not math.isfinite(a.perturbation):
        p.error('--disk-cells must be even and >=8; perturbation must be finite')
    if (a.disk_rings<4 or not math.isfinite(a.disk_stretch) or not 0<=a.disk_stretch<=4
            or not math.isfinite(a.disk_tangent_spacing) or a.disk_tangent_spacing<=0):
        p.error('annular grid requires >=4 rings, stretch in [0,4], positive finite tangent spacing')
    if a.disk_grid!='cartesian' and a.spatial_mode!='full_2d':
        p.error('--disk-grid annular requires --spatial-mode full_2d')
    if a.spatial_mode!='full_2d' and (a.perturbation or a.restart_2d):
        p.error('perturbation and 2D restart require --spatial-mode full_2d')
    if a.restart_2d and (not a.restart_2d.is_file() or a.initial_run):
        p.error('2D restart must exist and cannot be combined with --initial-run')
    build=a.build_dir.resolve()
    if a.initial_run:
        try:
            sources=restart_files(a.initial_run)
        except RunnerError as exc:
            p.error(str(exc))
    try:
        cmake, mpirun, env = resolve_toolchain()
    except RunnerError as exc:
        p.error(str(exc))
    subprocess.run([cmake,'-S',str(ROOT),'-B',str(build),'-DFERMIFORGE_REQUIRE_MPI=ON',
                    '-DCMAKE_BUILD_TYPE=Release'],check=True,env=env)
    subprocess.run([cmake,'--build',str(build),'--target','benchmark_radial_cylinder_2d',
                    '--parallel',str(a.ranks)],check=True,env=env)
    stamp=dt.datetime.now().strftime('%y%m%d-%H%M%S-%f')
    phase={-2:'bulk-a',-1:'bulk-b',0:'normal-core',1:'a-core'}[a.initialization]
    if a.texture!='none':
        phase=a.texture
    winding=1 if a.initialization>=0 else 0
    mode_tag='-full2d' if a.spatial_mode=='full_2d' else ''
    run=ROOT/'runs'/f'{stamp}-modern-cylinder{mode_tag}-R{a.radius:g}-n{winding}-{phase}-Fs{a.fs1:g}'
    run.mkdir(parents=True)
    shutil.copy2(build/'benchmark_radial_cylinder_2d',run/'solver')
    shutil.copy2(ROOT/'new_src/gauss11.dat',run/'gauss11.dat')
    shutil.copy2(ROOT/'new_src/ozaki_T=0.3.dat',run/'ozaki.dat')
    initial=''
    if a.restart_2d:
        shutil.copy2(a.restart_2d,run/'restart_fields_2d.dat')
        initial="restart_2d_file='restart_fields_2d.dat',"
    if a.initial_run:
        for src,dst in zip(sources,('initial_op_xyz','initial_curr')):
            shutil.copy2(src,run/dst)
        initial="initial_gap_file='initial_op_xyz', initial_current_file='initial_curr',"
    config=f'''&radial_cylinder
 disk_grid='{a.disk_grid}', disk_radial_layout='{a.disk_radial_layout}', disk_rings={a.disk_rings},
 disk_stretch={a.disk_stretch:.17g}, disk_tangent_spacing={a.disk_tangent_spacing:.17g},
 spatial_mode='{a.spatial_mode}', disk_cells={a.disk_cells}, perturbation={a.perturbation:.17g},
 radius={a.radius:.17g}, radial_points={a.radial_points}, azimuths={a.azimuths},
 temperature=0.30, fs1={a.fs1:.17g}, winding={winding},
 initialization={a.initialization}, texture='{a.texture}',
 maximum_step={a.step:.17g}, half_length={a.half_length:.17g},
 boundary_relaxation_distance={10*a.step:.17g},
 collision_aligned={'.true.' if a.collision_aligned else '.false.'},
 max_iterations={a.iterations}, tolerance={a.tolerance:.17g}, p_max=3.0,
 gauss_file='gauss11.dat', ozaki_file='ozaki.dat', output_prefix='cylinder',
 {initial}
/
'''
    (run/'input.nml').write_text(config)
    meta={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    meta.update(temperature=0.30,fs1=a.fs1,winding=winding,units='xi_0',status='running')
    if a.texture!='none':
        meta.update(winding_meaning='additional imposed winding, not total texture circulation',
                    seed_wall_circulation_quanta=1 if a.texture in ('a_mermin_ho','a_panam') else None,
                    texture_amplitude='sqrt(2) * legacy B bulk gap: initial scale only')
    if a.initial_run:
        meta.update(restart_sources=[str(path) for path in sources],anderson_history='fresh')
    (run/'run_metadata.json').write_text(json.dumps(meta,indent=2))
    env.setdefault('OMP_NUM_THREADS','1'); env.setdefault('OPENBLAS_NUM_THREADS','1')
    if sys.platform=='darwin':
        env.setdefault('OMPI_MCA_pml','ob1'); env.setdefault('OMPI_MCA_btl','self,sm')
    command=[mpirun,'--map-by','slot:OVERSUBSCRIBE','--bind-to','none','-n',str(a.ranks),
             './solver','input.nml']
    print(f'Run directory: {run}',flush=True)
    with (run/'run.log').open('w') as log,(run/'mpi.log').open('w') as err:
        proc=subprocess.Popen(command,cwd=run,env=env,text=True,stdout=subprocess.PIPE,stderr=err)
        for line in proc.stdout:
            print(line,end='',flush=True); log.write(line); log.flush()
        code=proc.wait()
    meta.update(status='completed' if code==0 else 'failed',return_code=code)
    (run/'run_metadata.json').write_text(json.dumps(meta,indent=2))
    if code:
        raise SystemExit(f'Run failed ({code}); see {run / "mpi.log"}')
    plot(run)
    print(f'Plots: {run / "plots"}')

if __name__=='__main__':
    main()
