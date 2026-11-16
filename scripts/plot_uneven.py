"""Figure for the uneven-spacing experiments."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core                       # noqa: E402

ALPHA = np.log(3) / np.log(4)
COL = {1: "#1f3f6f", 2: "#3a6ea5", 4: "#2c9c4a", 8: "#e08a2e", 16: "#c8342a", 32: "#8b2fa0"}


def load(rho, m):
    z = np.load(ROOT / f"results/uneven_bounded_{rho:g}_m{m}.npz")
    return core.gauge(z["v"]), z["x"]


def main():
    rows = json.loads((ROOT / "results/uneven-bounded.json").read_text())
    fig, (ax, bz, cz) = plt.subplots(1, 3, figsize=(6.8, 2.3),
                                 gridspec_kw=dict(width_ratios=[1.0, 1.15, 0.95]))

    # (a) D against m for each gap-ratio bound
    for rho in (1, 2, 4, 8, 16, 32):
        ms = [r["m"] for r in rows if r["par"] == rho]
        Ds = [r["D"] for r in rows if r["par"] == rho]
        ax.plot(ms, Ds, "o-", ms=2.6, lw=0.9, color=COL[rho],
                    label=f"gap ratio $\\leq{rho}$" + (" (uniform)" if rho == 1 else ""))
    ax.set_xticks([8, 12, 16, 24, 32])
    ax.set_xticklabels(["8", "12", "16", "24", "32"], fontsize=6)
    ax.set_xlim(6, 34)
    ax.set_xlabel("$m$")
    ax.set_ylabel("best $D_m$")
    ax.grid(alpha=0.25)
    ax.set_title("uneven spacing lowers the rungs", fontsize=6.8)
    ax.legend(fontsize=5.4, loc="lower right", handlelength=1.6, labelspacing=0.3)

    # (b) the m=16 shapes for three bindings
    for rho in (1, 4, 32):
        v, x = load(rho, 16)
        bz.plot(v.real, v.imag, "-", lw=0.9, color=COL[rho], label=f"$\\rho={rho}$")
        bz.plot(v.real, v.imag, ".", ms=1.8, color=COL[rho])
    bz.set_aspect("equal"); bz.set_xticks([]); bz.set_yticks([])
    bz.set_title("$m=16$: same profile, uneven sampling", fontsize=6.8)
    bz.legend(fontsize=5.4, loc="lower left")

    # (c) where the points go: cumulative parameter against index
    for rho in (1, 4, 32):
        v, x = load(rho, 32)
        cz.plot(np.arange(32) / 31, x, "o-", ms=2.0, lw=0.8, color=COL[rho],
                label=f"$\\rho={rho}$")
    cz.plot([0, 1], [0, 1], ":", lw=0.8, color="0.5")
    cz.set_xlabel("index fraction $i/(m-1)$")
    cz.set_ylabel("parameter $x_i$")
    cz.grid(alpha=0.25)
    cz.set_title("$m=32$: optimal spacings", fontsize=6.8)
    cz.legend(fontsize=5.4, loc="upper left")

    fig.tight_layout()
    out = ROOT / "figures/uneven-spacing.png"
    fig.savefig(out, dpi=150, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
