"""Prepare the real scratch-run inputs without doing production iterations."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
exe = Path(sys.argv[1]).resolve() / 'benchmark_axisymmetric_core_2d'
env = dict(os.environ, OMPI_MCA_pml='ob1', OMPI_MCA_btl='self,sm')
states = []
with tempfile.TemporaryDirectory(prefix='ff-bb-seeds-') as directory:
    for seed in ('0plus','plus0'):
        folder = Path(directory)/seed
        folder.mkdir()
        text = (ROOT/f'examples/radial_bb_localized_{seed}.nml').read_text()
        assert f"initialization_mode = 'localized_{seed}'" in text
        assert "restart_field_file = ''" in text
        text = text.replace('\n/\n','\n  prepare_only = .true.\n/\n')
        inp = folder/'input.nml'
        inp.write_text(text)
        outputs = [folder/s for s in ('input.dat','mapped.dat','metrics.txt','final.dat','history.dat','checkpoint.dat')]
        run = subprocess.run(['mpiexec','--host','localhost:1','--bind-to','none','-n','1',
                              str(exe),str(inp),*map(str,outputs)],
                             cwd=ROOT,env=env,text=True,capture_output=True,timeout=60)
        if run.returncode:
            raise RuntimeError(run.stdout+run.stderr)
        import re
        assert re.search(r'PREPARED ONLY:.*14161\s+60\s+0',run.stdout),run.stdout
        assert 'Bulk gap (ozaki)' in run.stdout
        assert 'Fs1, feedback_scale:' in run.stdout
        state = np.loadtxt(outputs[0])
        assert np.all(np.isfinite(state))
        states.append(state)
    assert states[0].shape == states[1].shape
    assert np.max(np.abs(states[0]-states[1])) > 0.01
print('Both 60-point scratch seeds prepared successfully, with distinct fields; no map evaluated.')
