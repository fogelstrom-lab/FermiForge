"""Check both historical positional layouts against Fortran and runner metadata."""
from pathlib import Path
import subprocess
import tempfile
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from run_and_plot import read_new_input_summary, RunnerError


def main():
    with tempfile.TemporaryDirectory(prefix='fermiforge-input-formats-') as directory:
        work = Path(directory)
        exe = work / 'input_probe'
        subprocess.run(['mpifort', '-fdefault-real-8', '-fdefault-double-8',
                        '-fcheck=all', '-std=f2018',
                        *[str(ROOT / 'new_src' / name) for name in
                          ('global_dec.f90', 'bulkgap.f90', 'init_calc.f90')],
                        str(ROOT / 'tests/test_radial_input_formats.f90'),
                        '-o', str(exe)], cwd=work, check=True)
        base = ['0.3','1','0','5.4','2','1e-6','1','3','4']
        cases = [(base,70,10,1), (base+['5'],70,10,5),
                 (base[:3]+['40','4']+base[3:],40,4,1),
                 (base[:3]+['40','4']+base[3:]+['5'],40,4,5)]
        cylinder = base.copy()
        cylinder[2]='1'
        cases += [(cylinder+['5','6'],70,6,5),
                  (cylinder[:3]+['40','4']+cylinder[3:]+['5','6'],40,6,5)]
        for records, extent, scale, pmax in cases:
            text = '! comment\n\n'+'\n'.join(v+' : record' for v in records)+'\n'
            path = work / 'input.inp'
            path.write_text(text)
            result = subprocess.run([str(exe)], input=text, text=True,
                                    capture_output=True, check=True, cwd=work)
            values = list(map(float,result.stdout.split()))
            expected = [.3,1,float(records[2]),extent,scale,5.4,2,1e-6,1,3,4,pmax]
            assert all(abs(a-b)<1e-12 for a,b in zip(values,expected)), values
            metadata = read_new_input_summary(path)
            assert metadata['fermi_liquid_f1s']==5.4
            assert metadata['anderson_p_max']==pmax
            assert metadata['maximum_iterations']==4
            assert metadata['free_grid_extent']==extent
            if records[2]=='1':
                assert metadata['radius_legacy_units']==scale
            else:
                assert metadata['free_radial_scale']==scale
        for records in (cylinder, cylinder+['5'], cylinder+['5','NaN'],
                        cylinder+['5','-1'], base+['0'], base+['5']*4):
            text='\n'.join(records)+'\n'
            path.write_text(text)
            result = subprocess.run([str(exe)], input=text, text=True,
                                    capture_output=True, cwd=work)
            assert result.returncode != 0
            try:
                read_new_input_summary(path)
            except RunnerError:
                pass
            else:
                raise AssertionError('Python parser accepted invalid input')
        metadata=read_new_input_summary(ROOT / 'new_src/qcv.inp')
        assert metadata['input_layout']=='explicit_grid'
        print('PASS: six radial layouts, invalid inputs, and current Mac input')


if __name__ == '__main__':
    main()
