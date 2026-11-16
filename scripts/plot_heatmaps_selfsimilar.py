"""Where the distortion lives, for the self-similar arcs and for the ladder.

For a parametrised arc with nodes (x_i, v_i) plot

    r(x_i,x_j) = |v_i - v_j| / |x_i - x_j|^alpha ,     alpha = log3/log4,

on I x I.  The self-similar arcs are drawn on *their own* parameter grids, which are not uniform:
a cell is narrow where the substitution puts many nodes in a short parameter interval.  Using a
uniform grid instead would require interpolating the polyline, which creates collinear nodes whose
r tends to zero and would drive the apparent minimum down.

Each panel is on its own logarithmic scale with a diverging map: blue at the minimum, white at the
geometric mean of that panel's range, red at the maximum, so the colour change falls at the
logarithmic centre of the range.  The geometric-mean level is also drawn as a contour.  Nothing is
clipped.  The diagonal is masked.

    python3 scripts/plot_heatmaps_selfsimilar.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import explicit_curves as ec          # noqa: E402
import align_core as core             # noqa: E402

ALPHA = ec.ALPHA


def field(x, v):
    d = np.abs(v[:, None] - v[None, :])
    g = np.abs(x[:, None] - x[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(g > 0, d / g ** ALPHA, np.nan)


def panel(ax, x, v, title):
    r = field(x, v)
    lo, hi = np.nanmin(r), np.nanmax(r)
    ax.pcolormesh(x, x, r, cmap="RdBu_r", norm=LogNorm(vmin=lo, vmax=hi),
                  shading="nearest", rasterized=True)
    ax.contour(x, x, np.nan_to_num(r, nan=np.nanmean(r)), levels=[np.sqrt(lo * hi)],
               colors="k", linewidths=0.35, alpha=0.6)
    ax.plot([0, 1], [0, 1], "-", lw=0.5, color="0.3")
    ax.set_title(f"{title}\n$r\\in[{lo:.3f},{hi:.3f}]$, $D={np.nanmax(r)/np.nanmin(r):.4f}$",
                 fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.tick_params(labelsize=8)
    ax.set_xlabel("$x$", fontsize=9); ax.set_ylabel("$y$", fontsize=9)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 5.6))
    z = ec.koch_nodes(5); z = (z - z[0]) / (z[-1] - z[0])
    panel(axes[0, 0], np.linspace(0, 1, len(z)), z, "Koch arc, level 5")
    # the two proved extremal pairs of the Koch arc (and the mirror images)
    for (a, b), mk in (((5 / 16, 11 / 16), "o"), ((63 / 256, 321 / 1024), "s")):
        axes[0, 0].plot([a, 1 - a, a, 1 - a], [b, 1 - b, 1 - b, b], mk, ms=5,
                        mfc="none", mec="k", mew=0.9)

    p, v, _ = ec.two_piece(t=0.563, level=9)
    o = np.argsort(p)
    panel(axes[0, 1], p[o], v[o] / (v[o][-1] - v[o][0]), "two-piece reversible, $t=0.563$, level 9")

    p, v, _ = ec.three_piece(level=5)
    o = np.argsort(p)
    panel(axes[1, 0], p[o], v[o] / (v[o][-1] - v[o][0]), "three-piece reversible, level 5")

    panel(axes[1, 1], np.linspace(0, 1, 800), core.gauge(np.load(ROOT / "results/m800_converged.npy")),
          "S-branch ladder, $m=800$")

    fig.tight_layout(h_pad=1.6, w_pad=1.0)
    out = ROOT / "figures/distortion-heatmaps.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
