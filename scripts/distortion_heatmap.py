"""Heatmaps of the node-level distortion density on I x I.

For a configuration of m nodes with x_i = i/(m-1) and chord 1,
    r(x_i, x_j) = |v_i - v_j| / |x_i - x_j|^alpha ,   alpha = log3/log4,
so max r / min r is the discrete distortion D of the project.

Two configurations:

* the standard Koch arc at level 5 (1025 nodes), the case whose structure is
  known exactly: the four sub-arcs are similar copies of the whole with ratio
  1/3 = 4^-alpha, so r is exactly invariant under the level-n rescaling;
* the stored m=800 S-branch configuration.

The continuum r along the Koch arc is *not* bounded below: inside a straight
piece F is affine, so |F(x)-F(y)|/|x-y|^alpha ~ |x-y|^{1-alpha} -> 0. The
bounded object is the node-level matrix used here, whose ratio is the D of the
finite configurations and stabilises at 2.0968 as the level grows.

    python3 scripts/distortion_heatmap.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core  # noqa: E402

ALPHA = core.ALPHA
SRC = ROOT / ".ultramath/smwhdqmht/m1/s_m800_plain.npy"


def koch_nodes(level):
    z = np.array([0 + 0j, 1 + 0j])
    for _ in range(level):
        out = [z[0]]
        for A, B in zip(z[:-1], z[1:]):
            d = (B - A) / 3
            out += [A + d, A + d + d * np.exp(1j * np.pi / 3), A + 2 * d, B]
        z = np.array(out)
    return z


def density(v):
    x = np.arange(len(v)) / (len(v) - 1)
    d = np.abs(v[:, None] - v[None, :])
    g = np.abs(x[:, None] - x[None, :]) ** ALPHA
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(g > 0, d / g, np.nan), x


def stats(r, x, tol=0.01):
    lo, hi = np.nanmin(r), np.nanmax(r)
    am = np.unravel_index(np.nanargmax(r), r.shape)
    an = np.unravel_index(np.nanargmin(r), r.shape)
    return dict(lo=float(lo), hi=float(hi), D=float(hi / lo),
                amax=(x[am[1]], x[am[0]]), amin=(x[an[1]], x[an[0]]),
                near_hi=np.argwhere(r >= hi * (1 - tol)),
                near_lo=np.argwhere(r <= lo * (1 + tol)), tol=tol)


def main():
    data = [("Koch arc, level 5 (1025 nodes)", *density(koch_nodes(5))),
            ("m=800 S-branch (stored nodes)", *density(core.gauge(np.load(SRC))))]
    st = [stats(r, x) for _, r, x in data]
    for (name, _, _), s in zip(data, st):
        print(f"{name}: r in [{s['lo']:.6f}, {s['hi']:.6f}]  D={s['D']:.9f}")
        print(f"   argmax at x,y = {s['amax']}   argmin at x,y = {s['amin']}")
        print(f"   pairs within 1% of max {len(s['near_hi'])}, of min {len(s['near_lo'])}")

    vmin = min(s["lo"] for s in st)
    vmax = max(s["hi"] for s in st)
    # r = 1 is the alpha-power law itself, so centre the diverging scale there;
    # the limits are exactly the extremes over both panels, so nothing is clipped
    for (name, r, _), s in zip(data, st):
        assert np.nanmin(r) >= vmin and np.nanmax(r) <= vmax
    assert vmin < 1.0 < vmax, (vmin, vmax)
    norm = Normalize(vmin=vmin, vmax=vmax)
    print(f"shared linear scale [{vmin:.6f}, {vmax:.6f}]: nothing falls outside it, "
          f"so nothing is clipped")
    fig, ax = plt.subplots(2, 2, figsize=(11.6, 10.8))

    for col, ((name, r, x), s) in enumerate(zip(data, st)):
        a0 = ax[0, col]
        im = a0.imshow(r, origin="lower", extent=[0, 1, 0, 1], cmap="YlOrRd",
                       norm=norm, interpolation="nearest")
        # the level set r = 1 is the alpha-power law itself
        a0.contour(np.linspace(0, 1, r.shape[1]), np.linspace(0, 1, r.shape[0]),
                   np.nan_to_num(r, nan=1.0), levels=[1.0], colors="k", linewidths=0.5)
        a0.set_title(f"{name}\n"
                     f"r in [{s['lo']:.4f}, {s['hi']:.4f}],  D = {s['D']:.4f}", fontsize=9)
        a0.set_xlabel("x", fontsize=8); a0.set_ylabel("y", fontsize=8)
        a0.tick_params(labelsize=7)
        cb = fig.colorbar(im, ax=a0, fraction=0.046, pad=0.02,
                          ticks=[vmin, 0.8, 1.0, vmax])
        cb.set_label("r(x,y);  contour and tick at the power law r = 1", fontsize=7)
        cb.ax.tick_params(labelsize=7)
        cb.ax.axhline(1.0, color="k", lw=0.8)
        a0.plot([0, 1], [0, 1], color="0.35", lw=0.4, alpha=0.7)

        a1 = ax[1, col]
        a1.imshow(r, origin="lower", extent=[0, 1, 0, 1], cmap="Greys",
                  interpolation="nearest", alpha=0.22)
        a1.plot(x[s["near_hi"][:, 1]], x[s["near_hi"][:, 0]], ".", ms=0.8,
                color="#c8342a", alpha=0.6, label=f"within 1% of max: {len(s['near_hi'])}")
        a1.plot(x[s["near_lo"][:, 1]], x[s["near_lo"][:, 0]], ".", ms=0.8,
                color="#1f6fd0", alpha=0.6, label=f"within 1% of min: {len(s['near_lo'])}")
        for pos, colr, lab in ((s["amax"], "#c8342a", "max"), (s["amin"], "#1f6fd0", "min")):
            a1.plot([pos[0]], [pos[1]], "x", color=colr, ms=8, mew=1.5)
            a1.annotate(f"{lab} ({pos[0]:.3f}, {pos[1]:.3f})", (pos[0], pos[1]),
                        textcoords="offset points", xytext=(6, 4), fontsize=7, color=colr)
        a1.set_xlim(0, 1); a1.set_ylim(0, 1)
        a1.set_xlabel("x", fontsize=8); a1.set_ylabel("y", fontsize=8)
        a1.tick_params(labelsize=7)
        a1.legend(fontsize=7, loc="lower right", markerscale=8, framealpha=0.85)
        a1.set_title(f"where the extremes sit (tolerance {int(100 * s['tol'])}%)", fontsize=9)

    fig.suptitle(r"distortion density $r(x,y)=|v_i-v_j|\,/\,|x_i-x_j|^{\alpha}$ on $I\times I$,"
                 r"  $\alpha=\log 3/\log 4$, both curves of chord 1", fontsize=10, y=0.996)
    fig.tight_layout()
    out = ROOT / "figures/distortion-heatmap.png"
    fig.savefig(out, dpi=130, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
