import os
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_modern_cylinder import resolve_toolchain, restart_files, RunnerError

class ToolchainTests(unittest.TestCase):
    def test_modern_restart_pair(self):
        with patch.object(Path,'is_file',lambda path: path.name in ('cylinder_final_op_xyz','cylinder_final_curr')):
            paths=restart_files('/example')
            self.assertEqual(tuple(p.name for p in paths),('cylinder_final_op_xyz','cylinder_final_curr'))

    def test_legacy_restart_pair(self):
        with patch.object(Path,'is_file',lambda path: path.name in ('op_xyz','curr')):
            self.assertEqual(tuple(p.name for p in restart_files('/example')),('op_xyz','curr'))

    def test_incomplete_restart_rejected(self):
        with patch.object(Path,'is_file',lambda path: path.name=='cylinder_final_op_xyz'):
            with self.assertRaisesRegex(RunnerError,'no complete'):
                restart_files('/example')

    def test_fallbacks_and_child_path(self):
        with patch.dict(os.environ,{'PATH':'/usr/bin:/bin'}), \
             patch('run_modern_cylinder.find_program',side_effect=[Path('/opt/local/bin/cmake'),
                                                                 Path('/opt/homebrew/bin/mpirun')]) as find, \
             patch.object(Path,'is_dir',return_value=True):
            cmake,mpi,env=resolve_toolchain()
            self.assertEqual(cmake,'/opt/local/bin/cmake')
            self.assertEqual(mpi,'/opt/homebrew/bin/mpirun')
            self.assertTrue(env['PATH'].startswith('/usr/bin:/bin:'))
            self.assertIn('/opt/homebrew/bin',env['PATH'].split(os.pathsep))
            self.assertEqual(os.environ['PATH'],'/usr/bin:/bin')
            self.assertIn(Path('/opt/local/bin/cmake'),find.call_args_list[0].args[1])

    def test_selected_tools_not_replaced(self):
        with patch.dict(os.environ,{'PATH':'/custom/bin:/usr/bin'}), \
             patch('run_modern_cylinder.find_program',side_effect=[Path('/custom/bin/cmake'),
                                                                 Path('/custom/bin/mpirun')]):
            cmake,mpi,env=resolve_toolchain()
            self.assertEqual(cmake,'/custom/bin/cmake')
            self.assertEqual(mpi,'/custom/bin/mpirun')
            self.assertTrue(env['PATH'].startswith('/custom/bin:/usr/bin'))

    def test_missing_tool_reported(self):
        with patch('run_modern_cylinder.find_program',side_effect=RunnerError('missing cmake')):
            with self.assertRaisesRegex(RunnerError,'missing cmake'):
                resolve_toolchain()

if __name__=='__main__':
    unittest.main()
