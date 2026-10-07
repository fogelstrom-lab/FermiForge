#!/usr/bin/env python3
"""Static iteration sheets and animated PNGs of spatial map defects."""
import argparse
import io
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
import numpy as np
from PIL import Image


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run', type=Path)
    p.add_argument('--duration-ms', type=int, default=900)
    a = p.parse_args()
    if a.duration_ms < 1:
        p.error('duration must be positive')
    files = sorted(a.run.glob('sampled_residuals.dat.map*.dat'))
    if not files:
        p.error('no per-map residual snapshots found')
    frames = []
    for f in files:
        d = np.loadtxt(f, ndmin=2)
        if d.shape[1] < 10 or not np.isfinite(d).all():
            raise ValueError(f'Invalid or incomplete snapshot: {f}')
        d = d[d[:, 10 if d.shape[1] >= 11 else 2] == 1]
        if not len(d):
            raise ValueError(f'No active points: {f}')
        if frames and not np.array_equal(d[:, :2], frames[0][1][:, :2]):
            raise ValueError('Mesh changed or a snapshot is incomplete; plot after map writing finishes')
        frames.append((int(re.search(r'map(\d+)', f.name)[1]), d))
    out = a.run / 'plots'
    out.mkdir(exist_ok=True)
    cmap = LinearSegmentedColormap.from_list('white_indigo', ['white', '#4b0082'])
    metadata = {'snapshots': [str(f.resolve()) for f in files],
                'meaning': 'Map k is evaluated before update k, after k-1 local updates.',
                'scale': 'fixed per metric across all frames; power 0.45', 'metrics': {}}
    for name, col, title in [('rms', 8, 'Point RMS residual'),
                              ('max', 9, 'Maximum real-component residual')]:
        vmax = max(d[:, col].max() for _, d in frames)
        norm = PowerNorm(.45, vmin=0, vmax=max(vmax, np.finfo(float).tiny))
        metadata['metrics'][name] = {'maximum': float(vmax)}

        def panel(ax, frame):
            number, d = frame
            x, y = np.unique(d[:, 0]), np.unique(d[:, 1])
            if len(y) == 1:
                ax.plot(d[:, 0], d[:, col], color='#4b0082', linewidth=1)
                artist = ax.scatter(d[:, 0], d[:, col], c=d[:, col], cmap=cmap, norm=norm)
                ax.set(xlabel=r'$r/\xi_0$', ylabel='Residual', ylim=(0, 1.05*norm.vmax),
                       title=f'Before update {number} — independent radial ray')
                ax.grid(alpha=.2)
                return artist
            z = np.full((len(y), len(x)), np.nan)
            z[np.searchsorted(y, d[:, 1]), np.searchsorted(x, d[:, 0])] = d[:, col]
            artist = ax.pcolormesh(x, y, np.ma.masked_invalid(z), shading='nearest', cmap=cmap, norm=norm)
            peak = d[np.argmax(d[:, col])]
            ax.plot(peak[0], peak[1], 'kx', ms=6)
            ax.set_facecolor('#e8e8e8')
            ax.set(aspect='equal', xlabel=r'$x/\xi_0$', ylabel=r'$y/\xi_0$',
                   title=f'Before update {number} | peak {peak[col]:.2e}\nat ({peak[0]:.2f}, {peak[1]:.2f})')
            return artist

        columns = min(4, len(frames))
        rows = (len(frames)+columns-1)//columns
        fig, axes = plt.subplots(rows, columns, figsize=(4*columns+1, 3.8*rows+0.8),
                                 squeeze=False, layout='constrained')
        for ax, frame in zip(axes.flat, frames):
            artist = panel(ax, frame)
        for ax in list(axes.flat)[len(frames):]:
            ax.set_visible(False)
        fig.suptitle(title+' — iteration dynamics\nFixed scale; white = zero, grey = unmeasured; × = one peak', fontsize=14)
        fig.colorbar(artist, ax=list(axes.flat), shrink=.75, label='Residual (power-0.45 intensity)')
        fig.savefig(out / f'residual_{name}_dynamics.png', dpi=150)
        plt.close(fig)

        images = []
        for frame in frames:
            fig, ax = plt.subplots(figsize=(7, 6), layout='constrained')
            artist = panel(ax, frame)
            fig.suptitle(title+' — fixed colour scale')
            fig.colorbar(artist, ax=ax, label='Residual (power-0.45 intensity)')
            buffer = io.BytesIO()
            fig.savefig(buffer, format='png', dpi=120)
            plt.close(fig)
            buffer.seek(0)
            images.append(Image.open(buffer).convert('RGBA'))
        path = out / f'residual_{name}_animation.png'
        images[0].save(path, save_all=True, append_images=images[1:],
                       duration=a.duration_ms, loop=0, disposal=0, blend=0)
        with Image.open(path) as check:
            assert check.n_frames == len(frames)
        print(out / f'residual_{name}_dynamics.png')
        print(path)
    (out / 'residual_dynamics_manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
