"""Static preview of what the alignment page shows, drawn with matplotlib.

No browser runs here, so this is the only rendering of the app's geometry that
can be looked at: the left panel with the interval, the endpoint-gauged overlay,
and the pointwise residual.  Numbers come from ``scripts/align_core.py``.

    python3 scripts/preview_align_app.py            # default interval [259,511]
    python3 scripts/preview_align_app.py 455 689
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
    v = core.gauge(np.load(SRC))
    rep = core.interval_report(v, a, b)
    # same parameters as the page: both knot sets plus a fine grid, so the plotted
    # maximum is the exact sup
    t = np.unique(np.concatenate([core.knots(len(v), b - a), np.linspace(0, 1, 1201)]))
    W = core.gauge(v)
    g = core.whole_at(W, t, rep["best"])
    w = core.block(v, a, b, t)
    res = np.abs(w - g)

    fig = plt.figure(figsize=(11.6, 4.5))
    ax1 = fig.add_axes([0.04, 0.08, 0.30, 0.82])
    ax2 = fig.add_axes([0.39, 0.08, 0.28, 0.82])
    ax3 = fig.add_axes([0.72, 0.08, 0.26, 0.82])

    ax1.plot(v.real, v.imag, "-", color="0.78", lw=1.0)
    seg = v[a:b + 1]
    ax1.plot(seg.real, seg.imag, "-", color="#e08a2e", lw=2.0)
    ax1.plot([seg[0].real, seg[-1].real], [seg[0].imag, seg[-1].imag], ":", color="#c8342a", lw=1.0)
    ax1.plot(seg[[0, -1]].real, seg[[0, -1]].imag, "o", color="#c8342a", ms=5, mec="w")
    ax1.set_aspect("equal")
    ax1.set_title(f"whole curve, interval [{a},{b}]", fontsize=9)
    ax1.tick_params(labelsize=7)

    ax2.plot(g.real, g.imag, "-", color="0.72", lw=1.2, label="whole")
    ax2.plot(w.real, w.imag, "-", color="#1f6fd0", lw=1.8, label="shortening")
    ax2.set_aspect("equal")
    ax2.set_title(f"overlay, {rep['best']}", fontsize=9)
    ax2.legend(fontsize=7, loc="lower right")
    ax2.tick_params(labelsize=7)

    ax3.plot(t, res, color="#c8342a", lw=1.2)
    ax3.axhline(rep["sup"], color="0.6", ls=":", lw=0.9)
    ax3.set_title("residual |w - gW| against t", fontsize=9)
    ax3.tick_params(labelsize=7)
    ax3.set_ylim(0, max(rep["sup"] * 1.15, 1e-6))

    fig.suptitle(
        f"sup {rep['sup']:.4f}   rms {rep['rms']:.4f}   free fit {rep['free_sup']:.4f}   "
        f"D(whole) {rep['D_whole']:.9f}   D(interval) {rep['D_block']:.9f}   "
        f"r/R {rep['r_over_R']:.4f}   percentile {100 * rep['r_percentile']:.1f}%",
        fontsize=9, y=0.99)
    out = ROOT / "figures/align-app-preview.png"
    fig.savefig(out, dpi=110, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)}  ({out.stat().st_size/1024:.0f} kB)")
    print(f"sup={rep['sup']:.6f} rms={rep['rms']:.6f} best={rep['best']} "
          f"D_block={rep['D_block']:.9f} pct={100*rep['r_percentile']:.1f}")


if __name__ == "__main__":
    main()
