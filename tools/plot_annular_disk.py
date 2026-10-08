"""Plot native annular node data; triangulation is for display only."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap
from plot_2d_fields import cartesian_to_harmonics


def plot_annular_disk(run, meta):
    dest = run/'plots'
    dest.mkdir(exist_ok=True)
    tables = [np.loadtxt(run/f'cylinder_{tag}_fields_2d.dat', ndmin=2)
              for tag in ('initial', 'final', 'mapped')]
    initial, final, mapped = tables
    if any(t.shape != final.shape or not np.allclose(t[:, :2], final[:, :2], atol=1e-12, rtol=0)
           for t in tables):
        raise ValueError('annular checkpoint coordinates differ')
    x, y = final[:, 0], final[:, 1]
    tri = mtri.Triangulation(x, y)
    fields = [(t[:, 2:20:2]+1j*t[:, 3:20:2]).reshape(-1, 3, 3).transpose(1, 2, 0) for t in tables]
    start, gap, map_gap = fields
    harmonic = cartesian_to_harmonics(gap)
    indigo = LinearSegmentedColormap.from_list('white_indigo', ['white', '#332288'])
    def save(fig, name):
        fig.savefig(dest/f'{name}.png', dpi=170)
        plt.close(fig)
    def decorate(ax):
        ax.set(xlabel=r'$x/\xi_0$', ylabel=r'$y/\xi_0$', aspect='equal')
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    ax.scatter(x, y, s=5, color='#332288')
    ax.add_patch(plt.Circle((0, 0), meta['radius'], fill=False, color='black'))
    decorate(ax)
    ax.set_title(f'Independent annular grid: {len(x)} points')
    save(fig, 'grid_points')
    for name, data in [('order_parameter', gap), ('harmonic', harmonic)]:
        for phase in (False, True):
            fig, axes = plt.subplots(3, 3, figsize=(12, 11), constrained_layout=True)
            vmax = max(float(np.max(abs(data))), 1e-30)
            for i in range(3):
                for j in range(3):
                    ax = axes[i, j]
                    values = np.cos(np.angle(data[i, j])) if phase else abs(data[i, j])
                    local = mtri.Triangulation(x, y, triangles=tri.triangles)
                    if phase:
                        local.set_mask(np.any(abs(data[i, j][tri.triangles]) < vmax*1e-3, axis=1))
                    # Entirely vanishing components have no defined phase.
                    if local.mask is None or not np.all(local.mask):
                        im = ax.tripcolor(local, values, shading='gouraud',
                                          cmap='bwr' if phase else indigo,
                                          vmin=-1 if phase else 0, vmax=1 if phase else vmax)
                        fig.colorbar(im, ax=ax, shrink=.7)
                    decorate(ax)
                    indices = ('+', '0', '-') if name == 'harmonic' else ('x', 'y', 'z')
                    ax.set_title(f'{name} {indices[i]}{indices[j]}')
            save(fig, name+('_phase_cosine' if phase else '_amplitude'))
    residual = mapped[:, 2:20]-final[:, 2:20]
    residual = np.column_stack((residual, mapped[:, 21:24]-final[:, 21:24]))
    rms = np.sqrt(np.mean(residual**2, axis=1))
    maximum = np.max(abs(residual), axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for ax, val, title in zip(axes, (rms, maximum), ('Point RMS residual', 'Maximum component residual')):
        im = ax.tripcolor(tri, val, shading='gouraud', cmap=indigo)
        fig.colorbar(im, ax=ax); decorate(ax); ax.set_title(title)
    save(fig, 'residual_maps')
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for ax, val, title in zip(axes, (final[:, 20], np.hypot(final[:, 21], final[:, 22]), final[:, 23]),
                              ('Pair density', 'In-plane current-related mean field', 'Axial current-related mean field')):
        im = ax.tripcolor(tri, val, shading='gouraud', cmap=indigo)
        fig.colorbar(im, ax=ax); decorate(ax); ax.set_title(title)
    if np.max(np.hypot(final[:, 21], final[:, 22])) > 1e-30:
        stride = max(1, len(x)//200)
        axes[1].quiver(x[::stride], y[::stride], final[::stride, 21], final[::stride, 22])
    save(fig, 'density_and_current')
    for basis, data, old in [('cartesian', gap, start),
                             ('harmonic', harmonic, cartesian_to_harmonics(start))]:
        fig, axes = plt.subplots(3, 3, figsize=(13, 10), constrained_layout=True)
        for i in range(3):
            for j in range(3):
                ax = axes[i, j]
                for coord, cross, label, color in [(x, y, 'x', '#332288'), (y, x, 'y', '#aa6600')]:
                    ids = np.flatnonzero(abs(cross) < 1e-10*max(1., meta['radius']))
                    ids = ids[np.argsort(coord[ids])]
                    for values, tag, alpha in [(old, 'initial', .35), (data, 'final', 1.)]:
                        ax.plot(coord[ids], values[i, j, ids].real, color=color, alpha=alpha, label=f'{tag} {label} Re')
                        ax.plot(coord[ids], values[i, j, ids].imag, '--', color=color, alpha=alpha, label=f'{tag} {label} Im')
                indices = ('+', '0', '-') if basis == 'harmonic' else ('x', 'y', 'z')
                ax.set(title=f'{basis} {indices[i]}{indices[j]}', xlabel='coordinate / xi0'); ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=6, ncol=2)
        save(fig, basis+'_axis_profiles')
    history = np.loadtxt(run/'cylinder_history.dat', ndmin=2)
    fig, ax = plt.subplots(constrained_layout=True)
    for col, label in [(1, 'RMS'), (2, 'maximum')]:
        ax.semilogy(history[:, 0], np.maximum(history[:, col], 1e-30), label=label)
    ax.set(xlabel='Evaluated state', ylabel='Absolute fixed-point residual')
    ax.legend(); ax.grid(alpha=.2); save(fig, 'convergence')
    peak = int(np.argmax(maximum))
    report = {'active_points': len(x), 'maximum_residual_location': [float(x[peak]), float(y[peak])],
              'maximum_gap_change_from_initial': float(np.max(abs(gap-start))),
              'relative_l2_gap_change_from_initial': float(np.linalg.norm(gap-start)/max(np.linalg.norm(start), 1e-30)),
              'initial_is_imported_radial_reference': bool(meta.get('initial_run')),
              'norm_weighting': 'unweighted nodes, not area integrals'}
    (run/'disk_comparison.json').write_text(json.dumps(report, indent=2)+'\n')
