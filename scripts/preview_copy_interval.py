"""Picture one copy-mode alignment: a block, and a similarity copy of an arc.

    python3 scripts/preview_copy_interval.py CURVE.npy A B S0 S1 [ABS_C ARG_DEG] OUT.png

CURVE is the whole; the block is [A,B] of the same curve unless a second path is
given, and the copy covers the index range [S0, S1] (both in node units). With no
|c| and angle the similarity is fitted; the translation is always the
least-squares one, as in the page's sliders.
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


def main():
    path = Path(sys.argv[1])
    a, b, s0i, s1i = (float(x) for x in sys.argv[2:6])
    a, b = int(a), int(b)
    if len(sys.argv) >= 9 and sys.argv[6] not in ("fit", "-"):
        abs_c, arg_c = float(sys.argv[6]), float(sys.argv[7])
        out = Path(sys.argv[8])
        c = abs_c * np.exp(1j * np.deg2rad(arg_c))
        fitted = False
    else:
        c = None
        out = Path(sys.argv[-1])
        fitted = True

    W = core.gauge(np.load(path))
    m = len(W)
    s0, s1 = s0i / (m - 1), s1i / (m - 1)
    t = np.linspace(0, 1, 2001)
    w = core.block(W, a, b, t)
    q = core._interp(W, s0 + s1 * t)
    if c is None:
        qc, wc = q - q.mean(), w - w.mean()
        c = np.vdot(qc, wc) / np.vdot(qc, qc).real
    d = w.mean() - c * q.mean()
    res = np.abs(w - (c * q + d))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    k = np.linspace(s0 * (m - 1), (s0 + s1) * (m - 1), 500).astype(int).clip(0, m - 1)
    ax[0].plot(W.real, W.imag, "-", color="0.82", lw=0.8)
    ax[0].plot(W[a:b + 1].real, W[a:b + 1].imag, "-", color="#e08a2e", lw=2.2,
               label=f"block [{a},{b}]")
    ax[0].plot(W[k].real, W[k].imag, "-", color="#2c9c4a", lw=2.2,
               label=f"arc [{s0i:.0f},{s0i + s1i:.0f}]")
    ax[0].set_aspect("equal"); ax[0].legend(fontsize=7)
    ax[0].set_title("whole curve, block and covered arc", fontsize=9)
    ax[1].plot((c * q + d).real, (c * q + d).imag, "-", color="#2c9c4a", lw=1.8, label="copy")
    ax[1].plot(w.real, w.imag, "-", color="#1f6fd0", lw=1.4, label="block")
    ax[1].set_aspect("equal"); ax[1].legend(fontsize=7)
    ax[1].set_title(f"overlay: |c| {abs(c):.4f}, arg {np.angle(c, deg=True):+.2f} deg", fontsize=9)
    for x in ax:
        x.tick_params(labelsize=7)
    fig.suptitle(f"sup {res.max():.5f}   rms {np.sqrt((res ** 2).mean()):.5f}   "
                 f"d = {d.real:+.6f}{d.imag:+.6f}i   "
                 f"(|c| and angle {'fitted' if fitted else 'given'})", fontsize=9, y=0.99)
    fig.tight_layout()
    fig.savefig(out, dpi=120, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out} ({out.stat().st_size/1024:.0f} kB)")
    print(f"|c|={abs(c):.6f} arg={np.angle(c, deg=True):+.3f} d={d.real:+.6f}{d.imag:+.6f}i "
          f"sup={res.max():.6f} rms={np.sqrt((res**2).mean()):.6f}")


if __name__ == "__main__":
    main()
