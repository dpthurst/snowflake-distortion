"""Draw the best partial overlaps found by search_overlaps.py.

Each panel compares a block with an arc after both are normalised to unit
diameter and the fitted rotation is applied — the same scale-free comparison
the search scores, so the visual agreement matches the number.

    python3 scripts/plot_overlaps.py [results/overlaps_m1600.json]
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core  # noqa: E402

CURVE = ROOT / ".ultramath/smwhdqmht/m1/s_m1600_plain.npy"


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/overlaps_m1600.json"
    hits = json.loads(src.read_text())[:4]
    V = core.gauge(np.load(CURVE))
    n = len(V) - 1
    fig, ax = plt.subplots(1, len(hits), figsize=(3.2 * len(hits), 3.4))
    if len(hits) == 1:
        ax = [ax]
    for k, (r, a_) in enumerate(zip(hits, ax)):
        a, L, c, K = r["a"], r["L"], r["c"], r["K"]
        t = np.linspace(0, 1, 801)
        w = core._interp(V, (a + t * L) / n)
        q = core._interp(V, (c + t * K) / n)
        wc, qc = w - w.mean(), q - q.mean()
        dw, dq = np.abs(wc).max(), np.abs(qc).max()
        s = np.vdot(wc / dw, qc / dq) / np.vdot(wc / dw, wc / dw).real
        # both curves normalised by their own diameter; the arc is then turned
        # back by the fitted rotation so that the two can be compared directly.
        # (An earlier version of this figure plotted s*block, i.e. the block
        # against its own rotation, which of course agreed.)
        arc = qc / dq / s
        shown = float(np.abs(arc - wc / dw).max())
        a_.plot(wc.real / dw, wc.imag / dw, "-", color="#1f6fd0", lw=1.6, label="block")
        a_.plot(arc.real, arc.imag, "--", color="#c8342a", lw=1.0,
                label="arc, turned back")
        a_.set_aspect("equal"); a_.set_xticks([]); a_.set_yticks([])
        a_.set_title(f"[{a:.0f},{a + L:.0f}] vs [{c:.0f},{c + K:.0f}]\n"
                     f"resid {r['residual']:.3f} (drawn {shown:.3f}), size {r['size']:.2f}, "
                     f"block->arc rotation {r['argC']:+.0f}°", fontsize=8)
        if k == 0:
            a_.legend(fontsize=7, loc="lower right")
    fig.suptitle("partial overlaps in the m=1600 curve: each shape normalised to unit "
                 "diameter, the arc turned back by the fitted rotation", fontsize=9, y=1.0)
    fig.tight_layout()
    out = ROOT / "figures/overlaps-m1600.png"
    fig.savefig(out, dpi=120, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size / 1024:.0f} kB)")
    for r in hits:
        print(f"  [{r['a']:.0f},{r['a'] + r['L']:.0f}] vs [{r['c']:.0f},{r['c'] + r['K']:.0f}]: "
              f"resid {r['residual']:.4f} size {r['size']:.3f} arg {r['argC']:+.1f}")


if __name__ == "__main__":
    main()
