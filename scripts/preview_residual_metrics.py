"""Compare the two residual metrics on one interval, as pictures.

Left: the gap at equal parameter, which is affine between knots and therefore
V-shaped.  Right: the distance from each curve to the piecewise-linear
interpolant of the other, which does not care where the two curves cross.

    python3 scripts/preview_residual_metrics.py            # [259,511] of m=800
    python3 scripts/preview_residual_metrics.py 455 689
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core  # noqa: E402

SRC = ROOT / "data/inputs/s_m800_plain.npy"


def main():
    a = int(sys.argv[1]) if len(sys.argv) > 2 else 259
    b = int(sys.argv[2]) if len(sys.argv) > 2 else 511
    v = np.load(SRC)
    W = core.gauge(v)
    rep = core.interval_report(v, a, b)
    kind = rep["best"]
    t = np.linspace(0, 1, 2001)
    w = core.block(v, a, b, t)
    g = core.whole_at(W, t, kind)
    gap = np.abs(w - g)

    tv = np.arange(len(v)) / (len(v) - 1)
    gv = core.whole_at(W, tv, kind)
    bt = np.arange(b - a + 1) / (b - a)
    bv = core.block(v, a, b, bt)
    fwd = core.dist_to_polyline(w.real, w.imag, gv.real, gv.imag)
    bwd = core.dist_to_polyline(g.real, g.imag, bv.real, bv.imag)

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 3.4), sharey=True)
    ax[0].plot(t, gap, color="#c8342a", lw=1.1)
    ax[0].set_title(f"gap at equal parameter: sup {gap.max():.4f}", fontsize=9)
    ax[0].set_xlabel("t", fontsize=8)
    ax[1].plot(t, fwd, color="#2c9c4a", lw=1.2, label=f"block to curve, sup {fwd.max():.4f}")
    ax[1].plot(t, bwd, color="#2c9c4a", lw=0.9, ls="--", label=f"curve to block, sup {bwd.max():.4f}")
    ax[1].set_title("distance to the other curve's interpolant", fontsize=9)
    ax[1].set_xlabel("t", fontsize=8)
    ax[1].legend(fontsize=7)
    for a_ in ax:
        a_.tick_params(labelsize=7)
        a_.set_ylim(0, max(gap.max(), fwd.max(), bwd.max()) * 1.08)
    fig.suptitle(f"interval [{a},{b}] of the m=800 curve, orientation {kind}", fontsize=9, y=1.0)
    fig.tight_layout()
    out = ROOT / "figures/residual-metrics-preview.png"
    fig.savefig(out, dpi=110, facecolor="white", bbox_inches="tight", pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")
    print(f"gap sup {gap.max():.5f} | distance block->curve {fwd.max():.5f} "
          f"curve->block {bwd.max():.5f}")


if __name__ == "__main__":
    main()
