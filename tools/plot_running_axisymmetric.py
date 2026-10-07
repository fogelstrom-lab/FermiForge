#!/usr/bin/env python3
"""Snapshot a running axisymmetric benchmark and plot it without running MPI.

Each invocation creates a new inspection-* directory. Residuals are pre-update
maps, whereas the checkpoint contains the field after the numbered update.
"""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def stable_read(path):
    before = path.stat()
    data = path.read_bytes()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (
            after.st_size, after.st_mtime_ns, after.st_ino) or len(data) != after.st_size:
        raise ValueError(f'{path} changed during reading; retry after checkpoint writing finishes')
    if not data or not data.endswith(b'\n'):
        raise ValueError(f'{path} is empty or incomplete; retry shortly')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path, nargs='?')
    parser.add_argument('--latest', metavar='NAME',
                        help='select newest run directory containing NAME')
    parser.add_argument('--no-radial-comparison', action='store_true',
                        help='plot the state and residuals without radial-reference overlays')
    args = parser.parse_args()
    if bool(args.run) == bool(args.latest):
        parser.error('provide either a run directory or --latest NAME')
    if args.latest:
        candidates = [p for p in (ROOT / 'runs').iterdir()
                      if p.is_dir() and args.latest in p.name and (p / 'input.nml').is_file()]
        if not candidates:
            parser.error(f'no run containing {args.latest!r} exists yet')
        run = max(candidates, key=lambda p: p.name)
    else:
        run = args.run.resolve()
    source = run / 'checkpoint_fields_2d.dat'
    if not source.is_file():
        parser.error(f'{run.name}: no saved checkpoint yet; retry after the first update')
    data = stable_read(source)
    match = re.search(rb'checkpoint_iteration_(\d+)', data)
    if not match:
        raise ValueError('Checkpoint has no recognized iteration label')
    iteration = int(match[1])
    # Retain immutable copies; never edit any solver-owned output.
    out = run / ('inspection-' + datetime.now().strftime('%y%m%d-%H%M%S-%f')
                 + f'-update{iteration:04d}')
    out.mkdir()
    field = out / 'final_fields_2d.dat'  # Compatibility name for the existing plotter.
    field.write_bytes(data)
    from plot_2d_fields import load_field_map
    loaded = load_field_map(field)
    for key, actual in [('nx_points', len(loaded.x)), ('ny_points', len(loaded.y))]:
        declared = re.search(rf'# {key}=(\d+)', data.decode())
        if not declared or int(declared[1]) != actual:
            raise ValueError('Checkpoint grid is incomplete; retry after writing finishes')
    reference = out / 'radial_reference.dat'
    if not args.no_radial_comparison:
        reference.write_bytes(stable_read(run / 'metrics.txt.radial_reference.dat'))
    inp = stable_read(run / 'input.nml')
    (out / 'input.nml').write_bytes(inp)
    clean = re.sub(r'!.*', '', inp.decode())
    radius = re.search(r'\bactive_radius\s*=\s*([\d.eEdD+\-]+)', clean, re.I)
    if not radius:
        raise ValueError('input.nml must explicitly specify active_radius')
    value = float(radius[1].lower().replace('d', 'e'))
    (out / 'metrics.txt').write_text(f'active_radius={value}\n')
    history = stable_read(run / 'iteration_history.dat').decode().splitlines()
    history = [line for line in history if line.startswith('#') or
               (line.strip() and int(line.split()[0]) <= iteration)]
    (out / 'iteration_history.dat').write_text('\n'.join(history) + '\n')
    snapshots = []
    for path in sorted(run.glob('sampled_residuals.dat.map*.dat')):
        number = int(re.search(r'map(\d+)', path.name)[1])
        if number <= iteration:
            (out / path.name).write_bytes(stable_read(path))
            snapshots.append(path.name)
    manifest = dict(source=str(source), checkpoint_iteration=iteration,
                    sha256=hashlib.sha256(data).hexdigest(), residual_snapshots=snapshots,
                    note='Checkpoint is post-update; residual maps are pre-update. '
                         'final_fields_2d.dat here is a checkpoint copy, not a converged result.')
    (out / 'snapshot_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    subprocess.run([
        sys.executable, str(ROOT / 'tools/plot_double_core_axis_profiles.py'),
        str(field), '--basis', 'harmonic', '--label', f'Checkpoint {iteration}',
        '--title', f'Harmonic axis profiles after update {iteration}',
        '--output', str(out / 'plots/harmonic_axis_profiles.png')], check=True)
    if args.no_radial_comparison:
        plots = out / 'plots'
        plots.mkdir(exist_ok=True)
        def call(script, *arguments):
            subprocess.run([sys.executable, str(ROOT / 'tools' / script),
                            *map(str, arguments)], check=True)
        call('plot_2d_fields.py', field, '--output-dir', plots,
             '--prefix', f'checkpoint_{iteration:04d}', '--formats', 'png')
        call('plot_double_core_axis_profiles.py', field,
             '--label', f'Checkpoint {iteration}',
             '--title', f'2D vortex: checkpoint after update {iteration}',
             '--output', plots / 'axis_profiles.png')
        if any(line.strip() and not line.startswith('#') for line in history):
            call('plot_2d_iteration_history.py', out / 'iteration_history.dat',
                 '--output', plots / 'convergence.png')
        if snapshots:
            call('plot_residual_locations.py', out, '--output', plots)
            call('plot_residual_dynamics.py', out)
        print(f'Plots: {plots}')
        return
    command = [sys.executable, str(ROOT / 'tools/recover_axisymmetric_plots.py'),
               str(out), '--radial-reference', str(reference),
               '--state-label', f'2D checkpoint after update {iteration}']
    if not snapshots:
        command.append('--comparison-only')
    print(f'Plotting saved update {iteration}; snapshot: {out}', flush=True)
    subprocess.run(command, check=True)
    print(f'Plots: {out / "plots"}')


if __name__ == '__main__':
    main()
