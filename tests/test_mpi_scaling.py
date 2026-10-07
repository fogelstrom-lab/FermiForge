import importlib.util
import tempfile
import unittest
import os
import sys
from pathlib import Path

SPEC=importlib.util.spec_from_file_location('scaling',Path(__file__).resolve().parents[1]/'tools/benchmark_mpi_scaling.py')
SCALING=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCALING)

class ScalingTests(unittest.TestCase):
    def test_failed_launch_is_not_success(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(RuntimeError,'exited 7'):
                SCALING.launch([sys.executable,'-c',"import sys; print('diagnostic',file=sys.stderr); sys.exit(7)"],
                               Path(folder),os.environ.copy(),5,Path(folder))
            self.assertIn('diagnostic',(Path(folder)/'mpi.log').read_text())
    def test_namelist_override(self):
        text="&test\n maximum_iterations = 100 ! old\n/\n"
        self.assertIn('maximum_iterations = 0',SCALING.set_value(text,'maximum_iterations','0'))
        self.assertIn('prepare_only = .false.',SCALING.set_value(text,'prepare_only','.false.'))
        self.assertTrue(SCALING.set_value(text,'prepare_only','.false.').endswith('/\n'))

    def test_compare_fields(self):
        with tempfile.TemporaryDirectory() as folder:
            a=Path(folder)/'a';b=Path(folder)/'b'
            a.write_text('# header\n1 2 3\n');b.write_text('1 2 3.000000000001\n')
            self.assertLess(SCALING.compare_fields(a,b,1e-11,0),1e-11)
            b.write_text('1 2 nan\n')
            with self.assertRaises(ValueError): SCALING.compare_fields(a,b,1e-11,0)
            b.write_text('1 2 3\n4 5 6\n')
            with self.assertRaises(ValueError): SCALING.compare_fields(a,b,1e-11,0)

    def test_summary_excludes_warmup(self):
        rows=[dict(ranks=2,warmup=True,map_seconds=100,elapsed_seconds=100),
              dict(ranks=2,warmup=False,map_seconds=10,elapsed_seconds=12),
              dict(ranks=4,warmup=False,map_seconds=5,elapsed_seconds=7)]
        result=SCALING.summarize(rows)
        self.assertEqual(result[1]['relative_map_efficiency'],1)
        self.assertAlmostEqual(result[1]['relative_elapsed_efficiency'],6/7)
        self.assertEqual(result[0]['repeats'],1)

    def test_stage_dependencies(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); dest=root/'archive';dest.mkdir()
            text='&test\n'
            for key in ('order_parameter_file','current_field_file','gauss_file','ozaki_file'):
                (root/key).write_text('fixture')
                text+=f" {key} = '{key}'\n"
            text+=" restart_field_file = ''\n/\n"
            staged,records=SCALING.stage_inputs(text,dest,root)
            self.assertEqual(len(records),4)
            self.assertIn(str(dest/'ozaki_file'),staged)

if __name__=='__main__': unittest.main()
