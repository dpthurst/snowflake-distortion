"""The optimal configurations at m = 3, 4, 5, 6.

m=3 has D = 1 exactly (three points can carry all three ratios equal); m=4 is the closed-form
optimum of the note; m=5 and m=6 are the certified ones.  Drawn in the endpoint gauge (v_0 = 0, v_{m-1} = 1)
with the chord dotted, the ratios printed underneath, and the nodes marked.

    python3 scripts/plot_small_m.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core      # noqa: E402

ALPHA = core.ALPHA


def ratios(v):
    m = len(v)
    i, j = np.triu_indices(m, 1)
    return np.abs(v[i] - v[j]) / np.abs(i - j) ** ALPHA


def panel(ax, v, D, title):
    v = core.gauge(v)
    ax.plot(v.real, v.imag, "-", lw=1.3, color="#1f3f6f", zorder=2)
    ax.plot([0, 1], [0, 0], ":", lw=0.9, color="0.55", zorder=1)
    ax.plot(v.real, v.imag, "o", ms=5, color="#c8342a", mec="w", zorder=3)
    for k, z in enumerate(v):
        ax.annotate(str(k), (z.real, z.imag), textcoords="offset points",
                    xytext=(4, 4), fontsize=7)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.margins(0.22)
    ax.set_title(f"{title}\n$D={D:.6f}$", fontsize=8)


def main():
    fig, axes = plt.subplots(1, 4, figsize=(6.8, 2.1))
    # m = 3: all three ratios equal 2^{-alpha}; an isoceles triangle
    v3 = np.array([0.0, 0.5 + 1j * np.sqrt(3 ** -1 - 0.25), 1.0])
    panel(axes[0], v3, 1.0, "$m=3$: every ratio equal")
    # m = 4: the closed form of the note
    c = 3.0 ** -ALPHA
    a, b = (1 - c) / 2, np.sqrt(c ** 2 - ((1 - c) / 2) ** 2)
    v4 = np.array([0.0, a + 1j * b, 1 - a + 1j * b, 1.0])
    panel(axes[1], v4, ratios(v4).max() / ratios(v4).min(), "$m=4$: four chords $=1/3$")
    # m = 5: the certified optimum
    v5 = np.conj(np.load(ROOT / "results/m5_best.npz")["v"])   # draw it above the chord
    panel(axes[2], v5, ratios(v5).max() / ratios(v5).min(), "$m=5$: four steps $=1/3$")
    v6 = np.conj(np.load(ROOT / "results/m6_best.npz")["v"])
    panel(axes[3], v6, ratios(v6).max() / ratios(v6).min(), "$m=6$")
    for ax, v in zip(axes, (v3, v4, v5, v6)):
        r = ratios(core.gauge(v))
        print("  ratios:", np.round(r, 6), " D =", np.round(r.max() / r.min(), 9))
    fig.tight_layout()
    out = ROOT / "figures/small-m-optima.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
