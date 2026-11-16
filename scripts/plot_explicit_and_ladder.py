"""Figures for the note on improving the Koch extremal.

  figures/explicit-curves.png   Koch, the two-piece and the three-piece arc
  figures/ladder-m.png          best known positions at m = 100 .. 1600
  figures/ladder-D.png          D_m against m, with the explicit bounds
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core            # noqa: E402
import explicit_curves as ec         # noqa: E402

KOCK_D = np.sqrt(633) * 23 ** -ec.ALPHA

LADDER = [
    (100, ROOT / ".scratch/final_m100.npz", "npz"),
    (200, ROOT / ".scratch/descended_m200.npy", "npy"),
    (400, ROOT / "results/m400_converged.npy", "npy"),
    (800, ROOT / "results/m800_converged.npy", "npy"),
    (1600, ROOT / ".ultramath/smwhdqmht/m1/s_m1600_plain.npy", "npy"),
]


def load(path, kind):
    if kind == "npz":
        z = np.load(path)
        return core.gauge(z["W"][int(np.argmin(z["Ds"]))])
    return core.gauge(np.load(path))


def fig_curves():
    p2, z2, _ = ec.two_piece(level=10)
    p3, z3, _ = ec.three_piece(level=7)
    zk = ec.koch_nodes(5)
    fig, ax = plt.subplots(1, 3, figsize=(12.6, 4.1))
    panels = [(zk, "Koch arc", f"$D=2.096805$  (sharp)"),
              (z2, "two-piece reversible", f"$D=2.04204$  (level 16)"),
              (z3, "three-piece reversible", f"$D=2.017935$  (levels 9–11)")]
    for a, (z, name, sub) in zip(ax, panels):
        w = (z - z[0]) / (z[-1] - z[0])
        if w.imag.mean() < 0:
            w = np.conj(w)
        a.plot(w.real, w.imag, "-", color="#1f3f6f", lw=1.1)
        a.plot(w[[0, -1]].real, w[[0, -1]].imag, "o", ms=4, color="#c8342a", mec="w")
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_title(f"{name}\n{sub}", fontsize=9)
    fig.suptitle("explicit arcs with distortion below the Koch constant", fontsize=10, y=1.0)
    fig.tight_layout()
    out = ROOT / "figures/explicit-curves.png"
    fig.savefig(out, dpi=130, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


def fig_ladder():
    fig, ax = plt.subplots(1, len(LADDER), figsize=(3.0 * len(LADDER), 3.3))
    rows = []
    for a, (m, path, kind) in zip(ax, LADDER):
        v = load(path, kind)
        D = core.distortion(v)
        rows.append((m, D))
        a.plot(v.real, v.imag, "-", color="#1f3f6f", lw=0.9)
        a.plot(v[[0, -1]].real, v[[0, -1]].imag, "o", ms=3.5, color="#c8342a", mec="w")
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_title(f"$m={m}$\n$D_m={D:.6f}$", fontsize=8.5)
    fig.suptitle("best known positions, endpoint-gauged, for large $m$", fontsize=10, y=1.02)
    fig.tight_layout()
    out = ROOT / "figures/ladder-m.png"
    fig.savefig(out, dpi=130, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")
    for m, D in rows:
        print(f"  m={m:5d}  D={D:.9f}")


def fig_ladder_D():
    ms = np.array([40, 100, 200, 400, 800, 1600], float)
    Ds = np.array([1.3523099782, 1.4362283999, 1.5031685289,
                   1.5417920000, 1.5900261898, 1.6355308997])
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogx(ms, Ds, "o-", color="#1f3f6f", label="best known $D_m$")
    ax.axhline(KOCK_D, color="#c8342a", lw=1.2, ls="--",
               label=f"Koch $D={KOCK_D:.6f}$ (sharp)")
    ax.axhline(2.04204, color="#e08a2e", lw=1.2, ls=":",
               label="two-piece arc, $D=2.04204$")
    ax.axhline(2.017935, color="#2c9c4a", lw=1.2, ls="-.",
               label="three-piece arc, $D=2.017935$")
    ax.set_xlabel("$m$"); ax.set_ylabel("$D$")
    ax.set_ylim(1.3, 2.15)
    ax.grid(alpha=0.25); ax.legend(fontsize=8)
    ax.set_title("the finite-$m$ ladder against the explicit arcs", fontsize=9.5)
    fig.tight_layout()
    out = ROOT / "figures/ladder-D.png"
    fig.savefig(out, dpi=140, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    fig_curves()
    fig_ladder()
    fig_ladder_D()
