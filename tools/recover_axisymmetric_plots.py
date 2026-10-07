#!/usr/bin/env python3
"""Plot saved axisymmetric results without launching the solver or clearing errors.

The reference must be an unperturbed radial embedding on the same grid.
No registration, phase alignment, or renormalization is applied to differences.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'work/matplotlib'))
from plot_2d_fields import load_field_map, COMPONENTS
from plot_double_core_axis_profiles import zero_index
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--radial-reference', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--comparison-only', action='store_true')
    parser.add_argument('--state-label', default='2D final state')
    args = parser.parse_args()
    run = args.run.resolve()
    reference_path = args.radial_reference.resolve()
    final = load_field_map(run / 'final_fields_2d.dat')
    reference = load_field_map(reference_path)
    metrics = dict(line.split('=', 1) for line in (run / 'metrics.txt').read_text().splitlines()
                   if '=' in line)
    active_radius = float(metrics['active_radius'])
    if not (np.array_equal(final.x, reference.x) and np.array_equal(final.y, reference.y)):
        raise ValueError('Reference and final field must have identical coordinates')
    out = args.output_dir.resolve() if args.output_dir else run / 'plots'
    out.mkdir(exist_ok=True, parents=True)
    x, y = np.meshgrid(final.x, final.y)
    weights = np.outer(np.gradient(final.y), np.gradient(final.x)) * (x*x+y*y <= active_radius**2)
    def relative_l2(a, b):
        denominator = np.sum(np.abs(b)**2 * weights)
        return float(np.sqrt(np.sum(np.abs(a-b)**2 * weights) / denominator)) if denominator else None
    comparison = {
        'mask': f'active disk r <= {active_radius}; nodal-area weighted; no alignment',
        'gap_relative_l2': relative_l2(final.order_parameter, reference.order_parameter),
        'mean_field_relative_l2': relative_l2(final.current, reference.current),
        'component_relative_l2': {name: relative_l2(final.order_parameter[i,j], reference.order_parameter[i,j])
                                  for name, i, j in COMPONENTS},
    }
    (out / 'radial_comparison_metrics.json').write_text(json.dumps(comparison, indent=2)+'\n')
    print(json.dumps(comparison, indent=2), flush=True)

    def call(script, *arguments):
        subprocess.run([sys.executable, str(ROOT / 'tools' / script),
                        *map(str, arguments)], check=True, cwd=ROOT)

    for label, field in [('final', run / 'final_fields_2d.dat'),
                         ('radial_reference', reference_path)]:
        call('plot_2d_fields.py', field, '--output-dir', out,
             '--prefix', label, '--formats', 'png')
    call('plot_double_core_axis_profiles.py', run / 'final_fields_2d.dat',
         '--reference-field', reference_path, '--reference-label', 'preserved radial solution',
         '--label', args.state_label, '--title', 'Vortex core: radial reference versus 2D',
         '--output', out / 'axis_profiles_radial_comparison.png')
    call('plot_double_core_axis_profiles.py', run / 'final_fields_2d.dat',
         '--reference-field', reference_path, '--reference-label', 'preserved radial solution',
         '--label', args.state_label, '--title', 'Vortex harmonics: radial reference versus 2D',
         '--basis', 'harmonic', '--output', out / 'harmonic_axis_profiles_radial_comparison.png')
    if not args.comparison_only:
        call('plot_2d_iteration_history.py', run / 'iteration_history.dat',
             '--output', out / 'convergence.png')
        call('plot_residual_locations.py', run, '--output', out)
        call('plot_residual_dynamics.py', run)

    # Overlay each complex component separately, retaining small components.
    for direction in ('x', 'y'):
        coord = final.x if direction == 'x' else final.y
        index = zero_index(final.y if direction == 'x' else final.x, direction)
        def cut(array):
            return array[..., index, :] if direction == 'x' else array[..., :, index]
        fig, axes = plt.subplots(3, 3, figsize=(14, 10), layout='constrained')
        for name, i, j in COMPONENTS:
            ax = axes[i, j]
            for values, style, label in [(reference.order_parameter, '--', 'radial'),
                                         (final.order_parameter, '-', args.state_label)]:
                v = cut(values)[i, j]
                ax.plot(coord, v.real, style, color='#4b0082', label=label + ' Re')
                ax.plot(coord, v.imag, style, color='#0072B2', label=label + ' Im')
            ax.set(title=f'A_{name}', xlabel=direction + r'/$\xi_0$')
            ax.axvline(-active_radius, color='.7', lw=.6)
            ax.axvline(active_radius, color='.7', lw=.6)
            ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f'Radial reference versus {args.state_label}: {direction}-axis\n'
                     'Raw amplitudes; solid = 2D, dashed = radial; no phase alignment')
        fig.savefig(out / f'radial_component_comparison_{direction}.png', dpi=180)
        plt.close(fig)

    delta = np.abs(final.order_parameter - reference.order_parameter)
    cmap = LinearSegmentedColormap.from_list('white_indigo', ['white', '#4b0082'])
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), layout='constrained')
    for name, i, j in COMPONENTS:
        ax = axes[i, j]
        image = ax.pcolormesh(final.x, final.y, delta[i, j], shading='nearest',
                              cmap=cmap, vmin=0, vmax=max(delta.max(), 1e-30))
        ax.set(title=f'|A_{name}(2D) − A_{name}(radial)|', aspect='equal')
    fig.colorbar(image, ax=axes.ravel().tolist(), label='Absolute complex-component difference')
    fig.suptitle('Difference from the preserved radial solution (not the iteration residual)')
    fig.savefig(out / 'radial_component_difference_maps.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), layout='constrained')
    for row, direction in enumerate(('x', 'y')):
        coord = final.x if direction == 'x' else final.y
        index = zero_index(final.y if direction == 'x' else final.x, direction)
        for field, style, label in [(reference, '--', 'radial'), (final, '-', '2D final')]:
            for col, values in enumerate([field.pair_density, *field.current]):
                v = values[index, :] if direction == 'x' else values[:, index]
                axes[row, col].plot(coord, v, style, label=label)
                axes[row, col].set(xlabel=direction + r'/$\xi_0$', title=
                    ['Pair-density proxy', 'Current-related mean field x',
                     'Current-related mean field y', 'Current-related mean field z'][col])
                axes[row, col].grid(alpha=.2)
        axes[row, 0].legend()
    fig.suptitle('Radial versus 2D axis profiles: saved proxies, not superfluid-density response')
    fig.savefig(out / 'radial_density_current_comparison.png', dpi=180)
    plt.close(fig)
    (out / 'recovered_plot_provenance.json').write_text(json.dumps({
        'final': str(run / 'final_fields_2d.dat'), 'radial_reference': str(reference_path),
        'note': 'Postprocessing only. Does not clear the MPI launcher failure. No alignment applied.'
    }, indent=2) + '\n')
    print(f'Recovered plots: {out}')


if __name__ == '__main__':
    main()
