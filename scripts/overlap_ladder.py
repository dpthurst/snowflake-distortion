"""One alignment, followed up the m ladder.

An overlap found at one m is written in "m=800 index units"; at each rung the
indices are rescaled by (m-1)/799, the shapes compared scale-free, and the pair
drawn: the block in blue, the arc normalised by its own diameter and turned back
by the fitted rotation, dashed.

The residual is reported twice: on the whole piece (401 samples) and on its
coarse shape only (24 samples). If the second stays put while the first grows,
what worsens is fine detail; if both grow, the coarse shape itself is changing.

    python3 scripts/overlap_ladder.py
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
import search_overlaps as so  # noqa: E402

PATHS = {200: ROOT / ".scratch/descended_m200.npy",
         400: ROOT / ".ultramath/smwhdqmht/m1/s_m400.npy",
         800: ROOT / ".ultramath/smwhdqmht/m1/s_m800_plain.npy",
         1600: ROOT / ".ultramath/smwhdqmht/m1/s_m1600_plain.npy"}
LADDER = (200, 400, 800, 1600)

# alignments in m=800 index units: (a, L, c, K)
ALIGNMENTS = [((120, 89, 24, 76), "one of the best found at m=800"),
              ((533, 61, 648, 30), "a short, high-quality one")]


def main():
    V = {m: core.gauge(np.load(p)) for m, p in PATHS.items()}
    for idx, (hit, label) in enumerate(ALIGNMENTS, start=1):
        a, L, c, K = hit
        fig, axes = plt.subplots(1, 4, figsize=(12.4, 3.5))
        print(f"\n{label}: block [{a},{a+L}] vs arc [{c},{c+K}] in m=800 units")
        print(f"{'m':>5s} {'index range':>20s} {'arc range':>20s} {'resid 401':>10s} "
              f"{'resid 24':>9s}")
        for ax, m in zip(axes, LADDER):
            r = (m - 1) / 799.0
            aa, LL, cc, KK = a * r, L * r, c * r, K * r
            n = len(V[m]) - 1
            t = np.linspace(0, 1, 801)
            w = so.sample(V[m], aa + t * LL)
            q = so.sample(V[m], cc + t * KK)
            wc, qc = w - w.mean(), q - q.mean()
            dw, dq = np.abs(wc).max(), np.abs(qc).max()
            s = np.vdot(wc / dw, qc / dq) / np.vdot(wc / dw, wc / dw).real
            arc = qc / dq / s
            r_fine = float(np.abs(arc - wc / dw).max())
            tt = np.linspace(0, 1, 24)
            w2 = so.sample(V[m], aa + tt * LL)
            q2 = so.sample(V[m], cc + tt * KK)
            r_coarse = so.shape_score(w2, q2)[1]
            print(f"{m:5d} [{aa:7.1f},{aa+LL:7.1f}] [{cc:7.1f},{cc+KK:7.1f}] "
                  f"{r_fine:10.3f} {r_coarse:9.3f}")
            ax.plot(wc.real / dw, wc.imag / dw, "-", color="#1f6fd0", lw=1.5,
                    label="block" if m == LADDER[0] else None)
            ax.plot(arc.real, arc.imag, "--", color="#c8342a", lw=1.0,
                    label="arc, turned back" if m == LADDER[0] else None)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"m={m}\nresid {r_fine:.3f}, coarse {r_coarse:.3f}", fontsize=8)
            if m == LADDER[0]:
                ax.legend(fontsize=7, loc="lower right")
        fig.suptitle(f"one alignment followed up the ladder — {label}; "
                     f"block [{a},{a+L}] vs arc [{c},{c+K}] at m=800", fontsize=9, y=1.02)
        fig.tight_layout()
        out = ROOT / f"figures/overlap-ladder-{idx}.png"
        fig.savefig(out, dpi=120, facecolor="white", bbox_inches="tight",
                    pil_kwargs={"optimize": True})
        print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
