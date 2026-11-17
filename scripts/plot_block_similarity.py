"""Parts of the m=1600 curve that resemble other parts, possibly at another scale.

Reuses the search of scripts/search_overlaps.py: for index intervals I=[a,a+L] and J=[c,c+K] it
asks whether V(c+Kt) ~ s V(a+Lt) + d for a complex s, scoring scale free (both blocks centred and
divided by their own diameter, only a rotation fitted).  Here the best few hits at *different size
ratios* are drawn: the block solid, the arc dashed after being normalised and turned back.

    python3 scripts/plot_block_similarity.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core       # noqa: E402
import search_overlaps as so    # noqa: E402

SRC = ROOT / ".ultramath/smwhdqmht/m1/s_m1600_plain.npy"


def main():
    V = core.gauge(np.load(SRC))
    m = len(V)
    ts = np.linspace(0, 1, 24)
    lengths = [int(f * m) for f in (0.10, 0.15, 0.22, 0.30)]
    ratios = [0.5, 0.7, 1.0, 1.4, 2.0]
    hits = so.coarse(V, lengths, ratios, ts, 16, 16, 0.32)
    print(f"{len(hits)} coarse hits")
    ref = [so.refine(V, h) for h in sorted(hits)[:400]]
    # drop the trivial correspondence: same place and same size
    ref = [r for r in ref
           if r["residual"] < 0.25
           and not (abs(r["size"] - 1) < 0.15 and abs(r["argC"]) < 10.0)]
    ref.sort(key=lambda r: r["residual"])
    print(f"{len(ref)} non-trivial refined hits below 0.25; best {ref[0]['residual']:.4f}")

    # best hit per bin of size ratio, and no two sharing index range
    bins = [(0.3, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 1.4), (1.4, 1.8), (1.8, 2.5), (2.5, 4.0)]
    chosen, taken = [], []
    for lo, hi in bins:
        for r in ref:
            if not (lo <= r["size"] < hi):
                continue
            iv = (r["a"], r["a"] + r["L"])
            if any(not (iv[1] < t0 or iv[0] > t1) for t0, t1 in taken):
                continue
            taken.append(iv)
            chosen.append(r)
            break
        if len(chosen) == 6:
            break
    print("chosen (residual, size ratio, rotation, block, arc):")
    for r in chosen:
        print(f"  resid {r['residual']:.4f}  size {r['size']:.2f}  rot {r['argC']:+.1f} deg  "
              f"block [{r['a']:.0f},{r['a']+r['L']:.0f}] arc [{r['c']:.0f},{r['c']+r['K']:.0f}]")

    t = np.linspace(0, 1, 401)
    fig, axes = plt.subplots(2, 5, figsize=(6.9, 3.9))
    for ax, r in zip(axes[0], chosen):
        w = so.sample(V, r["a"] + t * r["L"])
        q = so.sample(V, r["c"] + t * r["K"])
        s, resid, size = so.shape_score(w, q)
        wn = (w - w.mean()) / np.abs(w - w.mean()).max()
        qn = (q - q.mean()) / np.abs(q - q.mean()).max() / s
        ax.plot(wn.real, wn.imag, "-", lw=1.5, color="#1f6fd0")
        ax.plot(qn.real, qn.imag, "--", lw=1.1, color="#c8342a")
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{resid:.3f}, size {size:.2f}", fontsize=7)
        # where the two pieces sit on the whole curve
        bx = axes[1][list(axes[0]).index(ax)]
        bx.plot(V.real, V.imag, "-", lw=0.5, color="0.6")
        ia, ib = int(round(r["a"])), int(round(r["a"] + r["L"]))
        ic, id_ = int(round(r["c"])), int(round(r["c"] + r["K"]))
        bx.plot(V[ia:ib].real, V[ia:ib].imag, "-", lw=1.8, color="#1f6fd0")
        bx.plot(V[ic:id_].real, V[ic:id_].imag, "-", lw=1.8, color="#c8342a")
        bx.set_aspect("equal"); bx.set_xticks([]); bx.set_yticks([])
    for ax in axes[1][len(chosen):]:
        ax.axis("off")
    axes[0][0].set_ylabel("shapes", fontsize=7)
    axes[1][0].set_ylabel("on the curve", fontsize=7)
    fig.suptitle("parts of the $m=1600$ curve that resemble other parts\n"
                 "top: the block (blue) against the matching arc (red, normalised and turned back); "
                 "bottom: where the two sit on the whole curve", fontsize=7.5)
    fig.tight_layout(h_pad=0.6)
    out = ROOT / "figures/block-similarity.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
