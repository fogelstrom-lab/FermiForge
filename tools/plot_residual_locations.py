#!/usr/bin/env python3
"""Plot saved raw map defects and track their maxima (no solver evaluation)."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs', type=Path, nargs='+', help='snapshot files or a run directory')
    p.add_argument('--labels', nargs='+')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    files = []
    for path in a.inputs:
        if path.is_dir():
            files.extend(sorted(path.glob('sampled_residuals.dat.map*.dat')) or
                         [path / 'sampled_residuals.dat'])
        else:
            files.append(path)
    labels = a.labels or [f'Map {i+1} (pre-update)' for i in range(len(files))]
    if len(labels) != len(files):
        p.error('provide one label per snapshot')
    arrays, peaks = [], []
    for path, label in zip(files, labels):
        data = np.loadtxt(path, ndmin=2)
        if data.shape[1] < 10 or not np.isfinite(data).all():
            raise ValueError(f'Invalid residual data: {path}')
        if data.shape[1] >= 11:
            active = data[data[:, 10] == 1]
        else:
            active = data[data[:, 2] == 1]
        if not len(active):
            raise ValueError(f'No interior points: {path}')
        arrays.append(active)
        r, m = active[np.argmax(active[:, 8])], active[np.argmax(active[:, 9])]
        peaks.append(dict(label=label, source=str(path.resolve()),
                          rms_xy=r[:2].tolist(), rms=float(r[8]),
                          max_xy=m[:2].tolist(), maximum=float(m[9])))
    a.output.mkdir(parents=True, exist_ok=True)
    cmap = LinearSegmentedColormap.from_list('white_indigo', ['white', '#4b0082'])
    selected = list(range(len(files))) if len(files) <= 3 else [0, len(files)-1]
    fig, axes = plt.subplots(2, len(selected), figsize=(5*len(selected)+1, 10),
                             squeeze=False, layout='constrained')
    for row, col in enumerate((8, 9)):
        norm = PowerNorm(.45, vmin=0, vmax=max(d[:, col].max() for d in arrays))
        for ax, index in zip(axes[row], selected):
            ax.set_facecolor('#e8e8e8')
            d = arrays[index]
            x, y = np.unique(d[:, 0]), np.unique(d[:, 1])
            if len(y) == 1:
                ax.plot(d[:, 0], d[:, col], color='#4b0082', linewidth=1)
                im = ax.scatter(d[:, 0], d[:, col], c=d[:, col], cmap=cmap, norm=norm)
                ax.set(title=labels[index]+' — independent radial ray',
                       xlabel=r'$r/\xi_0$', ylabel='Residual', ylim=(0,1.05*norm.vmax))
                ax.grid(alpha=.2)
                continue
            grid = np.full((len(y), len(x)), np.nan)
            grid[np.searchsorted(y, d[:, 1]), np.searchsorted(x, d[:, 0])] = d[:, col]
            im = ax.pcolormesh(x, y, np.ma.masked_invalid(grid), cmap=cmap,
                               norm=norm, shading='nearest', rasterized=True)
            peak = peaks[index]['rms_xy' if col == 8 else 'max_xy']
            ax.plot(*peak, 'kx', ms=10, mew=2)
            ax.set(title=labels[index] + f'\npeak at ({peak[0]:g}, {peak[1]:g})',
                   xlabel=r'$x/\xi_0$', ylabel=r'$y/\xi_0$', aspect='equal')
        fig.colorbar(im, ax=axes[row].tolist(), label=
                     'Point RMS (21 real components)' if col == 8 else 'Largest absolute real-component defect')
    fig.suptitle('Raw self-consistency residual; common scale within each row\n'
                 'White = zero; grey = not sampled; × = one maximum; intensity uses power 0.45')
    fig.savefig(a.output / 'residual_maps.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    for key, value, style, label in [('rms_xy', 'rms', 'o-', 'Point RMS peak'),
                                    ('max_xy', 'maximum', 'x--', 'Component maximum')]:
        xy = np.array([r[key] for r in peaks])
        axes[0].plot(xy[:, 0], xy[:, 1], style, label=label)
        for i, (x, y) in enumerate(xy):
            axes[0].annotate(str(i+1), (x, y), xytext=(4, 4+i*3), textcoords='offset points')
        axes[1].semilogy(np.arange(1, len(peaks)+1), [r[value] for r in peaks], style, label=label)
    axes[0].set(xlabel=r'$x/\xi_0$', ylabel=r'$y/\xi_0$', title='Peak locations (one tie representative)')
    axes[1].set(xlabel='Saved map index', ylabel='Residual', title='Peak residual history')
    for ax in axes:
        ax.legend(); ax.grid(alpha=.25)
    fig.savefig(a.output / 'residual_peak_track.png', dpi=160)
    plt.close(fig)
    (a.output / 'residual_peak_locations.json').write_text(json.dumps(peaks, indent=2)+'\n')
    print(json.dumps(peaks, indent=2))


if __name__ == '__main__':
    main()
