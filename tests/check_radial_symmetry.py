"""Bounded driver/MPI regression; never starts a production run."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from run_axisymmetric_core_benchmark import replace_namelist_value

def main():
    executable = Path(sys.argv[1]).resolve() / 'benchmark_axisymmetric_core_2d'
    template = (ROOT / 'examples/2d_normal_core_convergence_smoke.nml').read_text()
    template = template.replace('\n/\n', '\n  require_convergence = .false.\n/\n')
    changes = dict(bulk_gap_mode="'manual'", spatial_mode="'radial_symmetry'", half_width=3, number_of_cells=12,
                   active_radius=3, fine_region_half_width=0.75, medium_region_half_width=1.75,
                   endpoint_policy="'state_asymptotic'", asymptotic_matching_radius=2.5,
                   asymptotic_outer_radius=8, asymptotic_maximum_step=1,
                   trajectory_maximum_step=1, maximum_iterations=2, convergence_tolerance=1e-30)
    for key, value in changes.items():
        template = replace_namelist_value(template, key, str(value))
    env = dict(os.environ, OMPI_MCA_pml='ob1', OMPI_MCA_btl='self,sm')
    with tempfile.TemporaryDirectory(prefix='ff-radial-test-') as directory:
        root = Path(directory)
        results = {}
        for method, ranks, mode in [('bb', 1, 'file'), ('bb', 1, 'generate'),
                                    ('bb', 2, 'generate'), ('anderson', 2, 'generate'),
                                    ('polyak', 1, 'generate'), ('polyak', 2, 'generate')]:
            folder = root / f'{method}-{ranks}-{mode}'
            folder.mkdir()
            inp = folder / 'input.nml'
            case = replace_namelist_value(template, 'iteration_method', repr(method))
            case = replace_namelist_value(case, 'ozaki_mode', repr(mode))
            diagnostics = mode == 'generate'
            case = replace_namelist_value(case, 'save_iteration_diagnostics',
                                          '.true.' if diagnostics else '.false.')
            if mode == 'generate':
                case = replace_namelist_value(case, 'temperature', '0.3')
                case = replace_namelist_value(case, 'ozaki_file', "'intentionally-absent.dat'")
            inp.write_text(case)
            outputs = [folder / s for s in ['input.dat','mapped.dat','metrics.txt','final.dat','history.dat','checkpoint.dat']]
            command = ['mpiexec', '--host', f'localhost:{ranks}', '--map-by', f'ppr:{ranks}:node',
                       '--bind-to', 'none', '-n', str(ranks), str(executable),str(inp),*map(str,outputs)]
            run = subprocess.run(command,cwd=ROOT,env=env,text=True,capture_output=True,timeout=120)
            if run.returncode:
                raise RuntimeError(run.stdout+run.stderr)
            metrics = dict(line.split('=',1) for line in outputs[2].read_text().splitlines() if '=' in line)
            assert metrics['spatial_mode'].strip()=='radial_symmetry'
            assert int(metrics['independent_spatial_points'])==7
            history=np.genfromtxt(outputs[4],names=True)
            assert len(history)==2 and np.all(np.isfinite(history['map_rms']))
            if method == 'polyak':
                assert np.all(history['engine']==3) and np.all(history['history']==0)
                assert np.all(history['p']==2.0)
                assert abs(float(metrics['polyak_drag'])-0.5)<1e-14
                plot=subprocess.run([sys.executable,str(ROOT/'tools/plot_2d_iteration_history.py'),
                                     str(outputs[4]),'--output',str(folder/'history.png')],
                                    cwd=ROOT,env=env,text=True,capture_output=True,timeout=60)
                if plot.returncode:
                    raise RuntimeError(plot.stdout+plot.stderr)
            data=np.loadtxt(outputs[3])
            assert np.all(np.isfinite(data))
            results[method,ranks,mode]=data
            if diagnostics and method != 'bb':
                vectors=[np.loadtxt(str(outputs[2])+f'.signed.map{k:06d}.dat') for k in (1,2)]
                np.testing.assert_allclose(vectors[1][:,1],vectors[0][:,1]+vectors[0][:,3],atol=1e-14)
                for j,v in enumerate(vectors):
                    np.testing.assert_allclose(np.sqrt(np.mean(v[:,2]**2)),history['map_rms'][j],rtol=1e-12)
                assert Path(str(outputs[2])+'.regional_diagnostics.dat').is_file()
                assert not Path(str(outputs[2])+'.bb_diagnostics.dat').exists()
                if method=='polyak':
                    np.testing.assert_allclose(vectors[1][:,3],0.5*vectors[0][:,3]+2*vectors[1][:,2],atol=1e-14)
            if diagnostics and method == 'bb':
                diag=np.genfromtxt(str(outputs[2])+'.bb_diagnostics.dat',names=True)
                vectors=[np.loadtxt(str(outputs[2])+f'.signed.map{k:06d}.dat') for k in (1,2)]
                for j,v in enumerate(vectors):
                    np.testing.assert_allclose(v[:,3],history['p'][j]*v[:,2],atol=2e-16,rtol=1e-12)
                    assert abs(np.sqrt(np.mean(v[:,2]**2))-history['map_rms'][j])<1e-12
                np.testing.assert_allclose(vectors[1][:,1],vectors[0][:,1]+vectors[0][:,3],atol=1e-14)
                assert diag['secant_valid'][0]==0 and diag['secant_valid'][1]==1
                s=vectors[1][:,1]-vectors[0][:,1]
                y=vectors[0][:,2]-vectors[1][:,2]
                np.testing.assert_allclose(diag['raw_BB1'][1],np.dot(s,s)/np.dot(s,y),rtol=1e-10)
                np.testing.assert_allclose(diag['raw_BB2'][1],np.dot(s,y)/np.dot(y,y),rtol=1e-10)
            assert Path(str(outputs[2])+'.ozaki.dat').is_file()
            if (method,ranks)==('bb',1):
                for script, extra in [('plot_residual_dynamics.py', []),
                                      ('plot_residual_locations.py', ['--output',str(folder/'plots')])]:
                    plot = subprocess.run([sys.executable,str(ROOT/'tools'/script),str(folder),*extra],
                                          cwd=ROOT,env=env,text=True,capture_output=True,timeout=60)
                    if plot.returncode:
                        raise RuntimeError(plot.stdout+plot.stderr)
        np.testing.assert_allclose(results['bb',1,'generate'],results['bb',2,'generate'],atol=2e-12,rtol=2e-12)
        np.testing.assert_allclose(results['polyak',1,'generate'],results['polyak',2,'generate'],atol=2e-12,rtol=2e-12)
        np.testing.assert_allclose(results['bb',1,'file'],results['bb',1,'generate'],atol=2e-11,rtol=2e-11)
        # Validate the user-facing 60-point example without evaluating any maps.
        prepared = root/'prepared'
        prepared.mkdir()
        example=(ROOT/'examples/radial_symmetry_a_phase.nml').read_text()
        example=example.replace('\n/\n','\n  prepare_only = .true.\n/\n')
        inp=prepared/'input.nml'
        inp.write_text(example)
        outputs=[prepared/s for s in ['input.dat','mapped.dat','metrics.txt','final.dat','history.dat','checkpoint.dat']]
        command=['mpiexec','--host','localhost:1','--map-by','ppr:1:node','--bind-to','none','-n','1',
                 str(executable),str(inp),*map(str,outputs)]
        run=subprocess.run(command,cwd=ROOT,env=env,text=True,capture_output=True,timeout=60)
        if run.returncode:
            raise RuntimeError(run.stdout+run.stderr)
        import re
        assert re.search(r'PREPARED ONLY:.*14161\s+60\s+0',run.stdout),run.stdout
        match = re.search(r'Bulk gap \(ozaki\) =\s*(\S+)',run.stdout)
        assert match and abs(float(match[1])-0.279894824067945)<1e-12,run.stdout
    print('radial symmetry BB/AA/Polyak and 1/2-rank equivalence passed')

if __name__=='__main__':
    main()
