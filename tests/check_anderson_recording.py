"""Bounded AA20 reset test: signed recording must leave the trajectory unchanged."""
import os
import sys
import tempfile
import subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_axisymmetric_core_benchmark import replace_namelist_value

def main():
    exe=Path(sys.argv[1]).resolve()/'benchmark_axisymmetric_core_2d'
    text=(ROOT/'examples/2d_normal_core_convergence_smoke.nml').read_text()
    text=text.replace('\n/\n','\n  require_convergence = .false.\n/\n')
    for key,val in dict(spatial_mode="'radial_symmetry'",iteration_method="'anderson'",
        half_width=3,number_of_cells=12,active_radius=3,fine_region_half_width=.75,
        medium_region_half_width=1.75,endpoint_policy="'state_asymptotic'",
        asymptotic_matching_radius=2.5,asymptotic_outer_radius=8,asymptotic_maximum_step=1,
        trajectory_maximum_step=1,maximum_iterations=22,convergence_tolerance=1e-30,
        anderson_cycle_iterations=20,simple_cycle_iterations=0).items():
        text=replace_namelist_value(text,key,str(val))
    env=dict(os.environ,OMPI_MCA_pml='ob1',OMPI_MCA_btl='self,sm')
    with tempfile.TemporaryDirectory(prefix='ff-aa-recording-') as tmp:
        fields=[]; histories=[]
        for enabled in (False,True):
            folder=Path(tmp)/str(enabled);folder.mkdir()
            inp=folder/'input.nml'
            inp.write_text(replace_namelist_value(text,'save_iteration_diagnostics','.true.' if enabled else '.false.'))
            paths=[folder/s for s in ('initial.dat','mapped.dat','metrics.txt','final.dat','history.dat','checkpoint.dat')]
            run=subprocess.run(['mpiexec','--host','localhost:2','--map-by','ppr:2:node','--bind-to','none',
                '-n','2',str(exe),str(inp),*map(str,paths)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=180)
            if run.returncode:raise RuntimeError(run.stdout+run.stderr)
            fields.append(np.loadtxt(paths[3]));histories.append(np.genfromtxt(paths[4],names=True))
            if enabled:
                vectors=[np.loadtxt(str(paths[2])+f'.signed.map{k:06d}.dat') for k in range(1,23)]
                for old,new in zip(vectors[:-1],vectors[1:]):
                    np.testing.assert_allclose(old[:,1]+old[:,3],new[:,1],atol=1e-14)
                assert histories[-1]['history'][20]==1 and histories[-1]['cycle_position'][20]==1
                assert histories[-1]['p'][20]==.01
                assert not Path(str(paths[2])+'.bb_diagnostics.dat').exists()
                region=np.genfromtxt(str(paths[2])+'.regional_diagnostics.dat',names=True)
                assert len(region)==22
        np.testing.assert_array_equal(fields[0],fields[1])
        for name in histories[0].dtype.names:
            if name!='map_seconds':np.testing.assert_array_equal(histories[0][name],histories[1][name])
    print('PASS: AA20 reset, signed transitions and bitwise recording neutrality (two ranks)')
if __name__=='__main__':main()
