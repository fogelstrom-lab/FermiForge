"""Compare saved signed residuals; no solver calls or changes to run data.

Use the two previously identified BB patterns as fixed probes in all runs.
Also compute each run's independent SVD over a common residual-threshold window.
Neither construction is a Jacobian eigenanalysis.
"""
from pathlib import Path
import json
import os
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = {
    'BB': ['261001-192802-radial-bb-diagnostics-0plus-T030-Fs54-radial-reference'],
    'AA cap 3': ['261002-130520-radial-anderson-aa20-0plus-pmax3-scratch-T030-Fs54-radial-reference'],
    'AA cap 5': ['261002-091852-radial-anderson-aa20-diagnostics-0plus-T030-Fs54-radial-reference',
                 '261002-094542-radial-anderson-aa20-0plus-pmax5-tight-continuation-radial-reference'],
    'AA cap 10': ['261002-111709-radial-anderson-aa20-0plus-pmax10-scratch-T030-Fs54-radial-reference'],
}


def modes_of(residual):
    unit = residual / np.linalg.norm(residual, axis=1)[:, None]
    val, vec = np.linalg.eigh(unit @ unit.T)
    val, vec = val[::-1], vec[:, ::-1]
    modes = vec[:, :3].T @ unit / np.sqrt(val[:3, None])
    return modes, val[:3] / val.sum()


def main():
    data = {}
    for label, names in CASES.items():
        folders = [ROOT / 'runs' / name for name in names]
        arr = np.concatenate([np.stack([np.loadtxt(f) for f in sorted(p.glob('metrics.txt.signed.map*.dat'))])
                              for p in folders])
        x, residual, update = (arr[:, :, i] for i in (1, 2, 3))
        history = np.concatenate([np.genfromtxt(p / 'iteration_history.dat', names=True) for p in folders])
        np.testing.assert_allclose(np.sqrt(np.mean(residual**2, axis=1)), history['map_rms'], rtol=1e-12)
        np.testing.assert_allclose(np.max(abs(residual), axis=1), history['map_max'], rtol=1e-12)
        np.testing.assert_allclose(x[1:], x[:-1] + update[:-1], atol=1e-14, rtol=1e-12)
        layout = np.loadtxt(folders[0] / 'metrics.txt.signed.layout.dat')
        if data:
            np.testing.assert_array_equal(layout, data['BB']['layout'])
        start = int(np.flatnonzero(history['map_max'] <= 2e-5)[0])
        data[label] = dict(x=x, r=residual, layout=layout, start=start)

    layout = data['BB']['layout']; radius = layout[:, 2]; n = len(radius)
    # Retain the exact reference window used in the earlier BB claim.
    probes, fractions = modes_of(data['BB']['r'][250:])
    for p in probes:
        if p[np.argmax(abs(p))] < 0:
            p *= -1
    figures = [plt.subplots(2, 2, figsize=(11, 7), layout='constrained', sharey='col') for _ in range(2)]
    # Keep identical column limits across both pages.
    projected = np.concatenate([d['r'][d['start']:] @ probes[:2].T for d in data.values()])
    result = {'reference_window': [251, len(data['BB']['r'])],
              'reference_unit_svd_fractions': fractions.tolist(), 'runs': {}}
    for row, (label, d) in enumerate(data.items()):
        residual = d['r']; start = d['start']; r = residual[start:]
        nr = np.linalg.norm(r, axis=1)
        shape = r.reshape(-1, 21, n)
        independent, energy = modes_of(r)
        a = (d['x'][-1].reshape(21, n)[:18:2] +
             1j * d['x'][-1].reshape(21, n)[1:18:2]).reshape(3, 3, n)
        da = np.stack([a[2], np.zeros_like(a[1]), -a[0]])
        tangent = np.zeros((21, n))
        tangent[:18:2] = da.reshape(9, n).real
        tangent[1:18:2] = da.reshape(9, n).imag
        outer = radius > 15
        inner = radius <= 15
        t = tangent[:, outer].ravel(); t /= np.linalg.norm(t)
        outer_r = shape[:, : , outer].reshape(len(r), -1)
        fl = shape[:, 18:, inner].reshape(len(r), -1)
        cosine = np.sum(fl[1:] * fl[:-1], axis=1) / (np.linalg.norm(fl[1:], axis=1) * np.linalg.norm(fl[:-1], axis=1))
        q = r @ probes[:2].T
        unit_q = q / nr[:, None]
        independent_info = []
        for mode, fraction in zip(independent, energy):
            m = mode.reshape(21, n)
            independent_info.append(dict(fraction=float(fraction),
                inner_fl_fraction=float(np.sum(m[18:, inner]**2)),
                outer_fraction=float(np.sum(m[:, outer]**2)),
                outer_spin_y_cosine=float(abs(m[:, outer].ravel() @ t) / np.linalg.norm(m[:, outer]))))
        result['runs'][label] = dict(window=[start+1, len(residual)],
            independent_patterns=independent_info,
            fixed_probe_mean_unit_squared_projection=np.mean(unit_q**2, axis=0).tolist(),
            fixed_probe_sign_changes=np.sum(q[1:] * q[:-1] < 0, axis=0).tolist(),
            inner_fl_reversals=int(sum(cosine < 0)), inner_fl_pairs=len(cosine),
            median_outer_rotation_squared_fraction=float(np.median((outer_r @ t)**2 / np.sum(outer_r**2, axis=1))),
            fixed_probe_first_last_absolute=np.abs(q[[0, -1]]).tolist())
        for col, (color, style) in enumerate([('#2166ac', '-'), ('#b2182b', '--')]):
            ax = figures[row//2][1][row%2, col]
            ax.plot(np.arange(start+1, len(residual)+1), q[:, col], color=color, ls=style, lw=1)
            ax.set_yscale('symlog', linthresh=1e-7)
            ax.set_ylim(min(projected[:,col].min()*1.3,-1e-7),max(projected[:,col].max()*1.3,1e-7))
            ax.axhline(0, color='0.5', lw=.5)
            ax.set_title(label)
            ax.set_xlabel('Recorded map'); ax.set_ylabel('Signed residual projection')
            ax.grid(alpha=.15)
            if label == 'AA cap 5':
                ax.axvline(80.5, color='0.5', ls=':', lw=1)
    for page,(fig,axes) in enumerate(figures,1):
        fig.suptitle('Fixed BB spatial probes in each method\nLeft: broad order-parameter pattern; right: inner Fermi-liquid pattern')
        fig.savefig(ROOT/f'docs/technical_note/figures/signed_patterns_comparison_{page}.png', dpi=180)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), layout='constrained')
    for j, (color, style) in enumerate([('#2166ac','-'),('#b2182b','--')]):
        mode = probes[j].reshape(21,n)
        for c in (4,12) if j == 0 else (18,):
            name = {4:'Re A_xz',12:'Re A_zx',18:'Fermi-liquid x'}[c]
            axes[j].plot(radius, mode[c], color=color, ls=style if c != 12 else ':', label=name)
        axes[j].axhline(0,color='0.5',lw=.5)
        axes[j].axvline(15,color='0.5',ls='--',lw=.7)
        axes[j].set_xlabel('Radius / xi0'); axes[j].set_ylabel('Unit-pattern component')
        axes[j].set_title(['Broad BB probe','Inner BB probe'][j]);axes[j].legend()
        if j == 1:
            axes[j].set_xlim(0,15)
    fig.savefig(ROOT/'docs/technical_note/figures/signed_pattern_shapes.png',dpi=180)
    plt.close(fig)
    (ROOT/'work/signed-bb-analysis/signed_pattern_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
