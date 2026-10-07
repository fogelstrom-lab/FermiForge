#!/usr/bin/env python3
"""Test the modern map/iteration starting from archived converged vortex fields."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from run_double_core_from_scratch import PROJECT_ROOT, set_namelist_value
from run_double_core_radius_ladder import fortran_string


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core', choices=('double', 'normal', 'a-phase'), default='double')
    parser.add_argument('--iterations', type=int, default=0,
                        help='0 measures the initial map defect; positive values iterate')
    parser.add_argument('--halo', choices=('frozen', 'state_asymptotic'), default='frozen')
    parser.add_argument('--native-grid', action='store_true',
                        help='retain the full legacy 601x601 mesh (double core only)')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--probe-only', action='store_true',
                        help='double-core map only at center, half cores and outer probes')
    parser.add_argument('--ranks', type=int, default=10)
    parser.add_argument('--fine-spacing', type=float,
                        help='override central mesh spacing for a resolution study')
    parser.add_argument('--no-build', action='store_true')
    parser.add_argument('--build-dir', type=Path,
                        default=PROJECT_ROOT / 'work/double-core-scratch-build')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()
    if args.iterations < 0 or args.ranks < 1:
        parser.error('iterations must be nonnegative and ranks positive')
    if args.probe_only and (args.core != 'double' or args.iterations != 0):
        parser.error('--probe-only requires --core double --iterations 0')
    if args.core != 'double' and (args.native_grid or args.halo != 'frozen'):
        parser.error('--native-grid and --halo apply only to the double core')
    stamp = dt.datetime.now().strftime('%y%m%d-%H%M%S-%f')
    prepared = PROJECT_ROOT / 'work' / 'converged-reference-inputs' / stamp
    prepared.mkdir(parents=True)
    if args.core == 'double':
        source = PROJECT_ROOT / '2D_benchmarks/Double_core_vortex_T=0.30_Fs1=5.4'
        text = (PROJECT_ROOT / 'examples/2d_double_core_overnight.nml').read_text()
        changes = {
            'benchmark_directory': fortran_string(source),
            'initialization_mode': "'legacy_split'" if args.native_grid else "'legacy_resampled'",
            'gauss_file': fortran_string(source / 'gauss11.dat'),
            'ozaki_file': fortran_string(source / 'ozaki.dat'),
            'azimuth_count': '64',
            'asymptotic_bulk_gap': '0.2799462',
            'half_core_offset_y': '11.8',
            'scratch_half_width': '60.0',
            'scratch_fine_region_half_width': '16.0',
            'scratch_medium_region_half_width': '36.0',
            'scratch_fine_spacing': '0.4',
            'scratch_medium_spacing': '0.8',
            'scratch_coarse_spacing': '2.0',
            'update_region': "'half_core_stencils'" if args.probe_only else "'centered_ellipse'",
            'update_radius': '32.0',
            'update_radius_x': '28.0',
            'update_radius_y': '32.0',
            'global_probe_radius': '40.0',
            'asymptotic_fit_inner_radius': '24.0' if args.halo != 'frozen' else '50.0',
            'asymptotic_matching_radius': '27.0' if args.halo != 'frozen' else '58.0',
            'asymptotic_outer_radius': '100.0',
            'bulk_reference_fit_radius': '58.0',
            'inactive_halo_policy': repr(args.halo),
            'asymptotic_ray_probe_count': '0',
            'asymptotic_probe_iterations': '0',
            'zero_mode_projection': "'none'",
            'maximum_iterations': str(args.iterations),
        }
        source_files = ('op_x', 'op_y', 'op_z', 'curr', 'gauss11.dat', 'ozaki.dat', 'qcv.inp')
        command = [sys.executable, str(PROJECT_ROOT / 'tools/run_double_core_from_scratch.py'),
                   '--run-name', 'double-core-converged-reference', '--iterations', str(args.iterations)]
    else:
        stem = 'normal' if args.core == 'normal' else 'a_phase'
        source = PROJECT_ROOT / f'benchmarks/{stem}_core_2d/reference/current_new_src_T0.30_F1s5.4'
        text = (PROJECT_ROOT / f'examples/2d_{stem}_core_multiscale.nml').read_text()
        changes = {
            'maximum_iterations': str(args.iterations),
            'endpoint_policy': "'radial_reference'",
            'asymptotic_matching_radius': '0.0',
            'initialization_mode': "'radial_reference'",
            'initial_perturbation_amplitude': '0.0',
            'half_width': '32.0', 'active_radius': '9.0',
            'fine_region_half_width': '4.0', 'medium_region_half_width': '7.0',
            'fine_spacing': '0.25', 'medium_spacing': '0.5', 'coarse_spacing': '2.0',
            'checkpoint_interval': '1',
        }
        source_files = ('op_xyz', 'curr', 'gauss11.dat', 'ozaki.dat')
        command = [sys.executable, str(PROJECT_ROOT / 'tools/run_axisymmetric_core_benchmark.py'),
                   '--case-name', f'{args.core}-converged-reference',
                   '--initialization', 'converged-radial-reference']
    for key, value in changes.items():
        text = set_namelist_value(text, key, value)
    if args.fine_spacing is not None:
        if args.fine_spacing <= 0:
            parser.error('--fine-spacing must be positive')
        key = 'scratch_fine_spacing' if args.core == 'double' else 'fine_spacing'
        text = set_namelist_value(text, key, str(args.fine_spacing))
    hashes = {}
    for leaf in source_files:
        path = source / leaf
        with path.open('rb') as stream:
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
            hashes[str(path)] = digest.hexdigest()
    input_path = prepared / 'input.nml'
    input_path.write_text(text)
    command.extend(['--input', str(input_path), '--ranks', str(args.ranks)])
    command.extend(['--build-dir', str(args.build_dir.resolve())])
    if args.no_build:
        command.append('--no-build')
    if args.no_plot:
        command.append('--no-plot')
    (prepared / 'provenance.json').write_text(json.dumps({
        'source_sha256': hashes, 'command': command, 'controls': vars(args),
        'interpretation': 'Execution success is not fixed-point agreement. Inspect unprojected map residuals.'
    }, indent=2, default=str) + '\n')
    print(f'Prepared input: {input_path}', flush=True)
    if args.prepare_only:
        return 0
    return subprocess.run(command, cwd=PROJECT_ROOT, check=False).returncode


if __name__ == '__main__':
    raise SystemExit(main())
