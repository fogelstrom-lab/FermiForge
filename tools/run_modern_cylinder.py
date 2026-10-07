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

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ranks',type=int,default=10)
    p.add_argument('--iterations',type=int,default=0,help='AA updates; 0 is one map')
    p.add_argument('--tolerance',type=float,default=2e-6,help='absolute maximum residual threshold')
    p.add_argument('--initialization',type=int,choices=(-1,-2,0,1),default=-1,
                   help='-1 bulk B, -2 bulk A (no vortex), 0 normal-core, 1 A-core (unit winding)')
    p.add_argument('--fs1',type=float,default=0.,help='Landau F_1^s; converted internally to A_1^s')
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
    if min(a.ranks,a.azimuths,a.radial_points-2)<=0 or min(a.radius,a.step,a.half_length)<=0 or a.iterations<0:
        p.error('positive counts/lengths required; iterations must be nonnegative')
    if not math.isfinite(a.tolerance) or a.tolerance<=0:
        p.error('--tolerance must be finite and positive')
    if not math.isfinite(a.fs1) or a.fs1<0:
        p.error('--fs1 must be finite and nonnegative')
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
    winding=1 if a.initialization>=0 else 0
    run=ROOT/'runs'/f'{stamp}-modern-cylinder-R{a.radius:g}-n{winding}-{phase}-Fs{a.fs1:g}'
    run.mkdir(parents=True)
    shutil.copy2(build/'benchmark_radial_cylinder_2d',run/'solver')
    shutil.copy2(ROOT/'new_src/gauss11.dat',run/'gauss11.dat')
    shutil.copy2(ROOT/'new_src/ozaki_T=0.3.dat',run/'ozaki.dat')
    initial=''
    if a.initial_run:
        for src,dst in zip(sources,('initial_op_xyz','initial_curr')):
            shutil.copy2(src,run/dst)
        initial="initial_gap_file='initial_op_xyz', initial_current_file='initial_curr',"
    config=f'''&radial_cylinder
 radius={a.radius:.17g}, radial_points={a.radial_points}, azimuths={a.azimuths},
 temperature=0.30, fs1={a.fs1:.17g}, winding={winding},
 initialization={a.initialization},
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
