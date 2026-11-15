"""Plot the achieved S-branch configurations; lines only connect sampled nodes for display."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]

def norm(v):
    return (v-v[0])/(v[-1]-v[0])

def render(v, n=500):
    v = norm(v)
    t0 = np.linspace(0., 1., len(v))
    t = np.linspace(0., 1., n)
    return t, np.interp(t, t0, v.real) + 1j*np.interp(t, t0, v.imag)

z = np.load(ROOT/'.scratch/final_m100.npz')
old = norm(z['W'][z['Ds'].argmin()])
root = ROOT/'.ultramath/smwhdqmht/m1'
configs = [
    (100, old, '#555555', '$m=100$ stored best'),
    (400, norm(np.load(root/'s_m400.npy')), '#2ca02c', '$m=400$ S'),
    (800, norm(np.load(root/'s_m800_plain.npy')), '#1f77b4', '$m=800$ S'),
    (1600, norm(np.load(root/'s_m1600_plain.npy')), '#d62728', '$m=1600$ S'),
]
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), dpi=150)
for m, v, colour, label in configs:
    t, w = render(v)
    for ax in axes[:1]:
        ax.plot(w.real, w.imag, lw=.95 if m > 100 else .75, color=colour, alpha=.90, label=label)
    axes[1].plot(t, w.imag, lw=.9, color=colour, alpha=.9, label=label)
arc = render(norm(np.load(root/'arc_m1600.npy')))[1]
axes[0].plot(arc.real, arc.imag, color='#444444', lw=1.0, ls='--', alpha=.75,
             label='arc branch, $m=1600$ (lower $D$)')
axes[0].axhline(0, color='black', lw=.55, alpha=.6)
axes[0].set_aspect('equal', adjustable='datalim')
axes[0].set_title('S-shaped and arc-shaped node sets')
axes[0].set_xlabel('$(v-v_0)/(v_{m-1}-v_0)$, real part')
axes[0].set_ylabel('imaginary part')
axes[0].legend(fontsize=7.5, loc='best', frameon=False)
axes[0].grid(alpha=.18)
axes[1].axhline(0, color='black', lw=.55, alpha=.6)
axes[1].set_title('Imaginary coordinate versus index parameter')
axes[1].set_xlabel('$t=i/(m-1)$')
axes[1].set_ylabel('$\\operatorname{Im}v(t)$')
axes[1].grid(alpha=.18)
fig.text(.5, .008, 'Two achieved, nonstationary families; node joins are display only, not continuum limits.',
         ha='center', fontsize=7.5, color='#555555')
fig.tight_layout(rect=(0, .025, 1, 1))
out = ROOT/'figures/s-branch-convergence.png'
fig.savefig(out, bbox_inches='tight', pil_kwargs={'optimize': True})
print(out, out.stat().st_size)
