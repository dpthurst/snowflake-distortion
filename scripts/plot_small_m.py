"""The optimal configurations at m = 3, 4, 5, 6, with their tight pairs.

m=3 has D = 1 exactly (three points can carry all three ratios equal); m=4 is the closed-form
optimum of the note; m=5 and m=6 are the certified ones.  Drawn in the endpoint gauge
(v_0 = 0, v_{m-1} = 1) with the chord dotted and the nodes numbered.

The two rows show the two families that a minimax optimum ties: the pairs at the minimum ratio
(top) and the pairs at the maximum ratio (bottom), each taken as tight when its ratio is within
1e-6 (relative) of the extreme.  Every other pair of each configuration sits strictly between
the two families, which is what makes the shape rigid.

    python3 scripts/plot_small_m.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core      # noqa: E402

ALPHA = core.ALPHA
TOL = 1e-6                     # relative tolerance for calling a ratio tight
CMIN, CMAX = "#1f5fa8", "#c8342a"


def ratios(v):
    m = len(v)
    i, j = np.triu_indices(m, 1)
    return np.abs(v[i] - v[j]) / np.abs(i - j) ** ALPHA, i, j


def tight(v):
    """(pairs at the minimum ratio, pairs at the maximum ratio)"""
    r, i, j = ratios(core.gauge(v))
    lo, hi = r.min(), r.max()
    return ([(int(a), int(b)) for r_, a, b in zip(r, i, j) if r_ <= lo * (1 + TOL)],
            [(int(a), int(b)) for r_, a, b in zip(r, i, j) if r_ >= hi * (1 - TOL)])


def oriented(v):
    """endpoint gauge, arc drawn above its chord"""
    v = core.gauge(v)
    return np.conj(v) if v.imag.mean() < 0 else v


def draw(ax, v, pairs, colour):
    """v must already be in the drawing orientation"""
    ax.plot(v.real, v.imag, "-", lw=0.9, color="0.72", zorder=1)
    ax.plot([0, 1], [0, 0], ":", lw=0.8, color="0.55", zorder=0)
    for a, b in pairs:
        ax.plot(v[[a, b]].real, v[[a, b]].imag, "-", lw=1.7, color=colour,
                solid_capstyle="round", zorder=2)
    ax.plot(v.real, v.imag, "o", ms=3.2, color="0.15", mec="w", mew=0.7, zorder=3)
    for k, z in enumerate(v):
        ax.annotate(str(k), (z.real, z.imag), textcoords="offset points",
                    xytext=(3.6, 3.2), fontsize=7, zorder=4)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])


def main():
    c = 3.0 ** -ALPHA
    a, b = (1 - c) / 2, np.sqrt(c ** 2 - ((1 - c) / 2) ** 2)
    configs = [
        ("$m=3$", np.array([0.0, 0.5 + 1j * np.sqrt(3 ** -1 - 0.25), 1.0]),
         "every ratio equal"),
        ("$m=4$", np.array([0.0, a + 1j * b, 1 - a + 1j * b, 1.0]), "four chords $=1/3$"),
        ("$m=5$", np.load(ROOT / "data/m5_best.npy"), "four steps $=1/3$"),
        ("$m=6$", np.load(ROOT / "data/m6_best.npy"), "five steps equal"),
    ]
    rows = ("tight minima", "tight maxima")
    fig, axes = plt.subplots(2, 4, figsize=(6.6, 2.45))
    for col, (name, v, gloss) in enumerate(configs):
        tmin, tmax = tight(v)
        w = oriented(v)
        xlo, xhi = w.real.min() - 0.10, w.real.max() + 0.16
        ylo, yhi = min(w.imag.min(), 0.0) - 0.10, w.imag.max() + 0.16
        D = ratios(core.gauge(v))[0]
        D = D.max() / D.min()
        for row, (pairs, colour, label) in enumerate(((tmin, CMIN, rows[0]),
                                                      (tmax, CMAX, rows[1]))):
            ax = axes[row, col]
            draw(ax, w, pairs, colour)
            ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
            ax.annotate(f"{len(pairs)} pairs" + (" (same three)" if len(pairs) == 3 and
                                                 row == 1 and col == 0 else ""),
                        (0.5, -0.03), xycoords="axes fraction", ha="center", va="top",
                        fontsize=7.2, color=colour)
            if row == 0:
                ax.set_title(f"{name}: {gloss}\n$D={D:.6f}$", fontsize=8.4)
            if col == 0:
                ax.set_ylabel(label, fontsize=8.6)
    handles = [Line2D([], [], color=CMIN, lw=1.9, label="pair at the minimum ratio"),
               Line2D([], [], color=CMAX, lw=1.9, label="pair at the maximum ratio"),
               Line2D([], [], color="0.72", lw=1.0, label="the configuration")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.6, frameon=False,
               bbox_to_anchor=(0.5, 0.002))
    for name, v, _ in configs:
        r = ratios(core.gauge(v))[0]
        print(f"  {name}: D = {r.max() / r.min():.9f}  "
              f"tight {len(tight(v)[0])}+{len(tight(v)[1])}")
    fig.subplots_adjust(left=0.052, right=0.996, top=0.79, bottom=0.175,
                        wspace=0.05, hspace=0.34)
    out = ROOT / "figures/small-m-optima.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
