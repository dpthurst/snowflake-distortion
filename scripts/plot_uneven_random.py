"""Figures for the uneven/random-spacing experiments (preliminary).

    figures/uneven-random-shapes.png      random fixed spacings at m=16,32,64 against m=800
    figures/uneven-random-diagnostics.png D against m, the m=64 draws overlaid, and distances
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
import align_core as core                       # noqa: E402

COL = ["#1f3f6f", "#c8342a", "#2c9c4a"]


def gauge(v):
    return core.gauge(np.asarray(v, complex))


def dense(v, N=4000):
    t = np.linspace(0, 1, N)
    return np.interp(t, np.linspace(0, 1, len(v)), v.real) + \
        1j * np.interp(t, np.linspace(0, 1, len(v)), v.imag)


def orient(v, ref):
    """pick v or its conjugate, whichever is closer to the reference curve"""
    Rd = dense(ref)
    d = lambda s: np.abs(dense(s)[:, None] - Rd[None, :]).min(axis=1).max()
    return v if d(v) <= d(np.conj(v)) else np.conj(v)


def dist_to(v, ref):
    Rd = dense(ref)
    d = lambda s: np.abs(dense(s)[:, None] - Rd[None, :]).min(axis=1).max()
    return min(d(v), d(np.conj(v)))


def main():
    ref = gauge(np.load(ROOT / "results/m800_converged.npy"))
    rows = json.loads((ROOT / "results/uneven-random.json").read_text())
    ms = sorted({r["m"] for r in rows})

    # ---- row of shapes
    fig, ax = plt.subplots(1, 4, figsize=(6.9, 1.55))
    for a, m in zip(ax, [16, 32, 64]):
        v = orient(gauge(np.load(ROOT / f"results/uneven_random_0_m{m}.npz")["v"]), ref)
        a.plot(dense(ref).real, dense(ref).imag, "-", lw=1.4, color="0.8")
        a.plot(v.real, v.imag, "o-", ms=2.2, lw=0.7, color=COL[0])
        a.set_title(f"random spacing, $m={m}$", fontsize=6.5)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    ax[3].plot(ref.real, ref.imag, "-", lw=0.9, color="0.35")
    ax[3].set_title("$m=800$ uniform", fontsize=6.5)
    ax[3].set_aspect("equal"); ax[3].set_xticks([]); ax[3].set_yticks([])
    fig.tight_layout()
    out = ROOT / "figures/uneven-random-shapes.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")

    # ---- diagnostics
    uni = {r["m"]: r["D"] for r in json.loads((ROOT / "results/uneven-bounded.json").read_text())
           if r["par"] == 1}
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(6.8, 2.3),
                                     gridspec_kw=dict(width_ratios=[1, 1.05, 0.95]))
    a1.plot(sorted(uni), [uni[m] for m in sorted(uni)], "o-", color="0.35", ms=2.6, lw=0.9,
            label="uniform")
    for k in range(3):
        xs = [r["m"] for r in rows if int(r["par"]) == k]
        ys = [r["D"] for r in rows if int(r["par"]) == k]
        a1.plot(xs, ys, "s--", color=COL[k], ms=2.8, lw=0.8, label=f"draw {k}")
    a1.set_xlabel("$m$"); a1.set_ylabel("best $D_m$"); a1.grid(alpha=0.25)
    a1.set_title("random spacings give lower $D$", fontsize=6.8)
    a1.legend(fontsize=5.2, loc="lower right")

    a2.plot(dense(ref).real, dense(ref).imag, "-", lw=1.6, color="0.7", label="$m=800$ uniform")
    for k in range(3):
        v = orient(gauge(np.load(ROOT / f"results/uneven_random_{k}_m64.npz")["v"]), ref)
        a2.plot(v.real, v.imag, "o-", ms=2.0, lw=0.7, color=COL[k], label=f"draw {k}, $m=64$")
    a2.set_aspect("equal"); a2.set_xticks([]); a2.set_yticks([])
    a2.set_title("three draws at $m=64$", fontsize=6.8)
    a2.legend(fontsize=5.0, loc="upper left")

    for k in range(3):
        ds = [dist_to(gauge(np.load(ROOT / f"results/uneven_random_{k}_m{m}.npz")["v"]), ref)
              for m in ms]
        a3.plot(ms, ds, "s--", color=COL[k], ms=2.6, lw=0.8, label=f"draw {k}")
    um = sorted(uni)[-5:]
    a3.plot(um, [dist_to(gauge(np.load(ROOT / f"results/uneven_bounded_1_m{m}.npz")["v"]), ref)
                 for m in um], "o-", color="0.35", ms=2.6, lw=0.9, label="uniform")
    a3.set_xlabel("$m$"); a3.set_ylabel("distance to $m=800$")
    a3.grid(alpha=0.25)
    a3.set_title("same curve, same distances", fontsize=6.8)
    a3.legend(fontsize=5.0)
    fig.tight_layout()
    out = ROOT / "figures/uneven-random-diagnostics.png"
    fig.savefig(out, dpi=150, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
