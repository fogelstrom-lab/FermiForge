"""Orbital chirality diagnostic: equals l only for a pure A-phase field."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def plot_texture(run, meta):
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    spinfig, spinaxes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    for row, tag in enumerate(('initial', 'final')):
        data = np.loadtxt(run/f'cylinder_{tag}_fields_2d.dat', ndmin=2)
        data = data[np.hypot(data[:, 0], data[:, 1]) <= meta['radius']+1e-12]
        x, y = data[:, 0], data[:, 1]
        a = (data[:, 2:20:2]+1j*data[:, 3:20:2]).reshape(-1, 3, 3)
        norm = np.sum(abs(a)**2, axis=(1, 2))
        orbital = (1j*np.cross(a, a.conj()).sum(axis=1)).real
        orbital = np.divide(orbital, norm[:, None], out=np.zeros_like(orbital), where=norm[:, None]>1e-30)
        tri = mtri.Triangulation(x, y)
        gram = np.einsum('nsi,nti->nst', a, a.conj()).real
        eigenvalues, eigenvectors = np.linalg.eigh(gram)
        director = eigenvectors[:, :, -1]
        # A director has no sign: choose a reproducible display convention.
        director *= np.where(director[:, 0]<0, -1., 1.)[:, None]
        fraction = np.divide(eigenvalues[:, -1], norm, out=np.zeros_like(norm), where=norm>1e-30)
        for col, val in enumerate((director[:, 2], fraction)):
            ax = spinaxes[row, col]
            im = ax.tripcolor(tri, val, shading='gouraud', cmap='bwr' if col==0 else 'Purples',
                              vmin=-1 if col==0 else 0, vmax=1)
            if col==0:
                ids = np.flatnonzero(fraction>.99)[::max(1, len(x)//200)]
                ax.quiver(x[ids], y[ids], director[ids, 0], director[ids, 1], scale=15)
            ax.set(title=f'{tag}: '+('principal spin axis z / in-plane arrows' if col==0 else 'largest spin eigenvalue / trace'),
                   xlabel='x / xi0', ylabel='y / xi0', aspect='equal')
            spinfig.colorbar(im, ax=ax)
        values = (orbital[:, 2], np.linalg.norm(orbital, axis=1), norm/3)
        for col, (val, title) in enumerate(zip(values, ('Orbital chirality z + in-plane arrows',
                                                       'Orbital chirality magnitude', 'Pair density'))):
            ax = axes[row, col]
            limits = dict(vmin=-1, vmax=1) if col == 0 else dict(vmin=0, vmax=1) if col == 1 else dict(vmin=0)
            im = ax.tripcolor(tri, val, shading='gouraud', cmap='bwr' if col==0 else 'Purples', **limits)
            if col == 0:
                stride = max(1, len(x)//200)
                ax.quiver(x[::stride], y[::stride], orbital[::stride, 0], orbital[::stride, 1], scale=15)
            ax.set(title=f'{tag}: {title}', xlabel='x / xi0', ylabel='y / xi0', aspect='equal')
            fig.colorbar(im, ax=ax)
    fig.suptitle('Normalized orbital chirality equals the unit l-vector only in pure A-phase')
    fig.savefig(run/'plots'/'orbital_texture.png', dpi=170)
    plt.close(fig)
    spinfig.suptitle('Spin director diagnostic: arrows shown only where rank-one fraction > 0.99')
    spinfig.savefig(run/'plots'/'spin_texture.png', dpi=170)
    plt.close(spinfig)
