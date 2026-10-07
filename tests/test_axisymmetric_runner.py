"""Regression for replacing an already-populated restart filename."""
import sys
import unittest
import os
import io
import tempfile
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from run_axisymmetric_core_benchmark import replace_namelist_value, parse_arguments, configure_mpi_shutdown
from run_axisymmetric_core_benchmark import run_logged, RunnerError
import run_axisymmetric_core_benchmark as runner


class NamelistReplacementTests(unittest.TestCase):
    def test_anderson_diagnostics_prepares_both_scratch_inputs(self):
        # Exercise the real runner through input preparation, stopping before
        # compilation/MPI. BB templates omit all three optional AA controls.
        class PreparedInput(Exception):
            pass

        for seed in ('0plus', 'plus0'):
            with self.subTest(seed=seed), tempfile.TemporaryDirectory() as directory:
                args = [
                    '--input', str(runner.PROJECT_ROOT / 'examples' /
                                   f'radial_bb_localized_{seed}.nml'),
                    '--run-root', directory, '--iteration-method', 'anderson',
                    '--anderson-cycle', '20', '--simple-cycle', '0',
                    '--anderson-history', '10', '--anderson-pmax', '5',
                    '--anderson-progress', '0.1', '--max-iterations', '600',
                    '--tolerance', '2e-5', '--checkpoint-interval', '1',
                    '--save-iteration-diagnostics', '--ranks', '10',
                ]
                with patch.object(runner, 'build_benchmark', side_effect=PreparedInput):
                    with self.assertRaises(PreparedInput):
                        runner.main(args)
                text = next(Path(directory).glob('*/input.nml')).read_text()
                for key, value in (
                    ('iteration_method', "'anderson'"),
                    ('anderson_history_limit', '10'),
                    ('anderson_maximum_mixing', '5.0'),
                    ('anderson_progress_threshold', '0.1'),
                    ('anderson_cycle_iterations', '20'),
                    ('simple_cycle_iterations', '0'),
                    ('save_iteration_diagnostics', '.true.'),
                ):
                    self.assertIn(f'{key} = {value}', text)
                    self.assertEqual(text.count(key + ' ='), 1)
                self.assertIn(f"initialization_mode = 'localized_{seed}'", text)

    def test_unknown_missing_control_is_rejected(self):
        with self.assertRaises(RunnerError):
            replace_namelist_value('&test\n/\n', 'anderson_histroy_limit', '10')

    def test_polyak_options(self):
        args = parse_arguments(['--iteration-method', 'polyak', '--polyak-step-size', '2', '--polyak-drag', '0.5'])
        self.assertEqual((args.iteration_method,args.polyak_step_size,args.polyak_drag),('polyak',2,0.5))
        text = '&test\n/\n'
        for key,value in [('polyak_step_size','2'),('polyak_drag','0.5')]:
            text=replace_namelist_value(text,key,value)
            self.assertIn(f'{key} = {value}',text)
    def test_signed_diagnostics_option(self):
        self.assertIsNone(parse_arguments([]).save_iteration_diagnostics)
        self.assertTrue(parse_arguments(['--save-iteration-diagnostics']).save_iteration_diagnostics)
        self.assertFalse(parse_arguments(['--no-save-iteration-diagnostics']).save_iteration_diagnostics)
        text=replace_namelist_value('&test\n/\n','save_iteration_diagnostics','.true.')
        self.assertIn('save_iteration_diagnostics = .true.',text)
    def test_ozaki_controls(self):
        self.assertEqual(parse_arguments(['--fs1','5.4']).fs1,5.4)
        self.assertIn('fs1 = 5.4', replace_namelist_value('&test\n/\n','fs1','5.4'))
        args = parse_arguments(['--temperature', '0.3', '--ozaki-cutoff', '50', '--bulk-gap-mode', 'legacy'])
        self.assertEqual(args.bulk_gap_mode, 'legacy')
        self.assertEqual(args.temperature, 0.3)
        self.assertEqual(args.ozaki_cutoff, 50)
        text = '&test\n/\n'
        for key, value in [('temperature', '0.3'), ('ozaki_cutoff', '50'), ('ozaki_mode', "'generate'")]:
            text = replace_namelist_value(text, key, value)
            self.assertIn(f'{key} = {value}', text)
    def test_radial_mode(self):
        args = parse_arguments(['--spatial-mode', 'radial_symmetry'])
        self.assertEqual(args.spatial_mode, 'radial_symmetry')
        self.assertIn("spatial_mode = 'radial_symmetry'",
                      replace_namelist_value('&test\n/\n', 'spatial_mode', "'radial_symmetry'"))
    def test_separate_runtime_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / 'run.log'
            mpi_log = Path(directory) / 'mpi.log'
            output = io.StringIO()
            with redirect_stdout(output):
                run_logged([sys.executable, '-c',
                            "import sys; print('iteration progress'); sys.stderr.write('PMIx diagnostic\\n' * 10000)"],
                           log, os.environ.copy(), echo=True, stderr_log_path=mpi_log)
            self.assertIn('iteration progress', output.getvalue())
            self.assertNotIn('PMIx diagnostic', output.getvalue())
            logged_lines = [line for line in log.read_text().splitlines()
                            if not line.startswith('command:')]
            self.assertNotIn('PMIx diagnostic', '\n'.join(logged_lines))
            self.assertEqual(mpi_log.read_text().count('PMIx diagnostic'), 10000)
            with self.assertRaises(RunnerError) as error:
                run_logged([sys.executable, '-c', "import sys; sys.stderr.write('failure'); sys.exit(7)"],
                           log, os.environ.copy(), stderr_log_path=mpi_log)
            self.assertIn('status 7', str(error.exception))
            self.assertIn(str(mpi_log), str(error.exception))
            self.assertEqual(mpi_log.read_text(), 'failure')

    def test_shutdown_environment(self):
        env = {}
        configure_mpi_shutdown(env, "darwin")
        self.assertEqual(env["PMIX_MCA_pmix_finalize_timeout"], "60")
        self.assertEqual(env["PMIX_MCA_pmix_client_base_verbose"], "2")
        env = {"PMIX_MCA_pmix_finalize_timeout": "120", "PMIX_MCA_pmix_client_base_verbose": "0"}
        configure_mpi_shutdown(env, "darwin")
        self.assertEqual(env["PMIX_MCA_pmix_finalize_timeout"], "120")
        self.assertEqual(env["PMIX_MCA_pmix_client_base_verbose"], "0")
        env = {}
        configure_mpi_shutdown(env, "linux")
        self.assertEqual(env, {})
    def test_bb_options(self):
        args = parse_arguments(["--iteration-method", "bb", "--bb-curvature", "absolute",
                                "--bb-maximum-mixing", "10"])
        self.assertEqual(args.iteration_method, "bb")
        self.assertEqual(args.bb_maximum_mixing, 10)
        text = replace_namelist_value("&test\n/\n", "iteration_method", "'bb'")
        self.assertIn("iteration_method = 'bb'", text)
    def test_iteration_cycle_options(self):
        restart_only = parse_arguments(["--anderson-cycle", "20", "--simple-cycle", "0"])
        self.assertEqual(restart_only.simple_cycle, 0)
        args = parse_arguments(["--anderson-cycle", "20", "--simple-cycle", "3",
                                "--simple-mixing", "0.1"])
        self.assertEqual((args.anderson_cycle, args.simple_cycle, args.simple_mixing),
                         (20, 3, 0.1))
        text = replace_namelist_value("&test\n/\n", "anderson_cycle_iterations", "20")
        self.assertIn("anderson_cycle_iterations = 20", text)
        text = replace_namelist_value(text, "anderson_cycle_iterations", "10")
        self.assertEqual(text.count("anderson_cycle_iterations"), 1)
        self.assertIn("anderson_cycle_iterations = 10", text)

    def test_localized_harmonic_modes(self):
        for mode in ('0plus', 'plus0', '0minus', 'minus0'):
            name = 'localized_' + mode
            args = parse_arguments(['--initialization-mode', name])
            self.assertEqual(args.initialization_mode, name)

    def test_repeated_restart(self):
        text = "&test\n restart_field_file = ''\n maximum_iterations = 0\n/\n"
        for path in ("'/old/run/checkpoint.dat'", "'/new/run/final.dat'"):
            text = replace_namelist_value(text, 'restart_field_file', path)
        text = replace_namelist_value(text, 'maximum_iterations', '20')
        self.assertEqual(text, "&test\n restart_field_file = '/new/run/final.dat'\n maximum_iterations = 20\n/\n")

    def test_quoted_delimiters_and_comments(self):
        for old in ("'/path/a,b!c.dat'", '"/path/a,b!c.dat"', "'/it''s/a.dat'"):
            text = ' restart_field_file = ' + old + ', ! keep comment\n/\n'
            self.assertEqual(replace_namelist_value(text, 'restart_field_file', "'/new/file'"),
                             " restart_field_file = '/new/file', ! keep comment\n/\n")


if __name__ == '__main__':
    unittest.main()
