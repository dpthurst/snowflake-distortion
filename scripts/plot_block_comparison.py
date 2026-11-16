"""Does a contiguous block of a configuration see the same distortion as the whole?

For a block of nodes [i0, i1] with its own parameter values x, the restricted distortion is the
max/min of r_ij = |v_i - v_j| / |x_j - x_i|^alpha over pairs *inside* the block.  If the block is a
copy of the whole -- which is what self-similarity gives at the substitution scales -- this equals
the global D exactly.  For the numerical ladder it is a property, not a theorem.

    python3 scripts/plot_block_comparison.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import explicit_curves as ec          # noqa: E402
import align_core as core             # noqa: E402

ALPHA = ec.ALPHA


def restricted(x, v, i0, i1):
    s = slice(i0, i1 + 1)
    d = np.abs(v[s][:, None] - v[s][None, :])
    g = np.abs(x[s][:, None] - x[s][None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(g > 0, d / g ** ALPHA, np.nan)
    return np.nanmax(r) / np.nanmin(r)


def main():
    # the ladder configuration
    v = core.gauge(np.load(ROOT / "results/m800_converged.npy"))
    m = len(v)
    x = np.linspace(0, 1, m)
    D = restricted(x, v, 0, m - 1)
    lens = [10, 20, 40, 80, 160, 320, 640]
    rows = []
    for L in lens:
        devs = [abs(restricted(x, v, i, i + L - 1) / D - 1.0) for i in range(0, m - L + 1)]
        devs = np.array(devs)
        rows.append((L, devs.max(), np.median(devs)))
        print(f"  m=800, blocks of {L:4d}: {len(devs):4d} blocks, "
              f"max |D_block/D-1| = {devs.max():.3e}, median {np.median(devs):.3e}", flush=True)

    # the self-similar arcs: aligned blocks are copies of the whole
    z = ec.koch_nodes(6); z = (z - z[0]) / (z[-1] - z[0])
    xk = np.linspace(0, 1, len(z))
    Dk = restricted(xk, z, 0, len(z) - 1)
    print(f"\n  Koch level 6: D = {Dk:.9f}")
    for j in (5, 3):
        L = 4 ** j + 1
        devs = [abs(restricted(xk, z, i, i + L - 1) / Dk - 1.0)
                for i in range(0, len(z) - L + 1, L - 1)]
        Db = restricted(xk, z, 0, L - 1)
        print(f"    blocks of {L:4d} nodes = the level-{j} pieces ({len(devs)} of them): "
              f"D_block = {Db:.12f}, max |D_block/D-1| = {max(devs):.3e}")
    p, v2, _ = ec.two_piece(t=0.563, level=10); o = np.argsort(p)
    p, v2 = p[o], v2[o]
    D2 = restricted(p, v2, 0, len(v2) - 1)
    n1 = int(np.searchsorted(p, 0.563))
    d1 = abs(restricted(p, v2, 0, n1 - 1) / D2 - 1.0)
    d2 = abs(restricted(p, v2, n1, len(v2) - 1) / D2 - 1.0)
    print(f"  two-piece level 10: D = {D2:.9f}; the two pieces as blocks: "
          f"deviations {d1:.2e}, {d2:.2e}")

    fig, ax = plt.subplots(figsize=(4.7, 2.9))
    L = np.array([r[0] for r in rows], float)
    ax.semilogy(L, [r[1] for r in rows], "o-", ms=4, lw=1.0, color="#1f3f6f", label="max over blocks")
    ax.semilogy(L, [r[2] for r in rows], "s--", ms=3.6, lw=0.9, color="#c8342a", label="median")
    ax.axhline(1e-15, color="#2c9c4a", lw=1.0, ls=":",
               label="self-similar arc, aligned blocks (exact)")
    ax.set_xlabel("block length $L$ (of $m=800$)")
    ax.set_ylabel("$|D_{\\rm block}/D-1|$")
    ax.set_title("blocks of the $m=800$ configuration", fontsize=9)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    out = ROOT / "figures/block-comparison.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"\nwrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
