#!/usr/bin/env python3
"""Independent reduced-pole radial cylinder map comparison; no production runs."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import shutil
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build-dir', type=Path, default=ROOT/'work/mac-linux-integration-build')
    args = parser.parse_args()
    modern = args.build_dir.resolve()/'benchmark_radial_cylinder_2d'
    with tempfile.TemporaryDirectory(prefix='fermiforge-cylinder-compare-') as tmp:
        work = Path(tmp)
        src = ROOT/'new_src'
        files = [src/n for n in ('global_dec.f90','bulkgap.f90','init_calc.f90','mpicalls.f90',
                  '3DFS_MPICodes/trajectories.f90','interpol.f90','riccati.f90','getnewses.f90')]
        files.append(ROOT/'tests/test_radial_cylinder_map.f90')
        subprocess.run(['mpifort','-std=f2018','-O2','-fdefault-real-8','-fdefault-double-8',
                        '-fcheck=all',*map(str,files),'-o','legacy_probe'],cwd=work,check=True)
        shutil.copy2(src/'gauss11.dat',work/'gauss11.dat')
        shutil.copy2(src/'ozaki_T=0.3.dat',work/'ozaki.dat')
        for radius,manufactured in ((4.,False),(10.,False),(10.,True)):
            fs1=5.4 if manufactured else 0.
            if manufactured:
                r=np.linspace(0,radius,100)
                seed=np.zeros((100,19)); seed[:,0]=r
                seed[:,1]=.28*(1-.8*(r/radius)**2)
                seed[:,9]=.28*(1+.1*(r/radius)**2)
                seed[:,17]=.28*(1+.05*(r/radius)**2)
                current=np.zeros((100,5)); current[:,0]=r; current[:,3]=.001*r
                np.savetxt(work/'op_xyz',seed); np.savetxt(work/'curr',current)
            start=3 if manufactured else -1
            inp=f'0.30\n0\n1\n{fs1}\n2\n1e-6\n{start}\n3\n0\n1\n{radius}\n'
            subprocess.run([str(work/'legacy_probe')],input=inp,text=True,cwd=work,check=True,
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=120)
            nml=f'''&radial_cylinder
 radius={radius}, fs1={fs1}, radial_points=100, azimuths=2, pole_limit=1,
 maximum_step={radius/99:.17g}, half_length={800*radius/99:.17g},
 boundary_relaxation_distance={10*radius/99:.17g}, collision_aligned=.false.,
 initial_gap_file='initial_op_xyz', initial_current_file='initial_curr',
 gauss_file='gauss11.dat', ozaki_file='ozaki.dat', output_prefix='modern'
/
'''
            (work/'input.nml').write_text(nml)
            subprocess.run([str(modern),'input.nml'],cwd=work,check=True,timeout=120)
            old=np.loadtxt(work/'cylinder_map.dat')
            new=np.loadtxt(work/'modern_mapped_op_xyz')
            current=np.loadtxt(work/'modern_mapped_curr')
            gap_error=np.max(np.abs(old[:,1:19]-new[:,1:19]))
            current_error=np.max(np.abs(old[:,[19,21,23]]-current[:,2:5]))
            print(f'R={radius:g}: gap difference={gap_error:.6e}; current difference={current_error:.6e}',flush=True)
            loc=np.unravel_index(np.argmax(np.abs(old[:,1:19]-new[:,1:19])),(100,18))
            print('maximum at row/component',loc,'radius',old[loc[0],0],flush=True)
            away_error=np.max(np.abs(old[1:,1:19]-new[1:,1:19]))
            print(f'away from exact centre: {away_error:.6e}',flush=True)
            # Centre transverse rays hit grid-coincident reflections: legacy
            # repeated additions choose sides by rounding. This is a bounded
            # discretization comparison, not a bitwise geometry oracle.
            assert np.isfinite(new).all() and gap_error<2e-6 and away_error<2e-7 and current_error<2e-7
            if manufactured:
                # Independent check of the collision-aligned branch and MPI
                # ownership/reduction, using the same nonuniform fields.
                (work/'input.nml').write_text(nml.replace('collision_aligned=.false.','collision_aligned=.true.'))
                subprocess.run([str(modern),'input.nml'],cwd=work,check=True,timeout=120)
                one=np.loadtxt(work/'modern_mapped_op_xyz')
                subprocess.run(['mpirun','--map-by','slot:OVERSUBSCRIBE','--bind-to','none','-n','10',
                                str(modern),'input.nml'],cwd=work,check=True,timeout=120)
                many=np.loadtxt(work/'modern_mapped_op_xyz')
                assert np.isfinite(many).all() and np.max(np.abs(one-many))<1e-13
                print('PASS collision-aligned map: one and ten MPI ranks agree',flush=True)
        print('PASS modern/legacy uniform-step bulk and nonuniform finite-Fs1 cylinder maps')

if __name__=='__main__':
    main()
