"""Static picture of the copy mode: a small-m curve fitted into a piece of the
m=800 curve by a similarity and an affine reparametrisation.

Mirrors what the page draws (no browser here), so the geometry can be looked at.

    python3 scripts/preview_align_copy.py
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

W = np.load(ROOT / "data/inputs/s_m800_plain.npy")


def main():
    S = np.load(ROOT / "data/inputs/final_m40.npz")["W"][0]
    a, b = 0, len(S) - 1
    fit = core.fit_copy_multires(W, a, b, S)
    s0, s1 = fit["s0"], fit["s1"]
    m = len(W)
    Wg = core.gauge(W)
    t = np.linspace(0, 1, 1201)
    w = core.block(S, a, b, t)
    q = core._interp(Wg, s0 + s1 * t)
    copy = fit["c"] * q + fit["d"]
    res = np.abs(w - copy)
    nat = core.profile_at(Wg, t, w, a / (m - 1), (b - a) / (m - 1))

    fig = plt.figure(figsize=(11.6, 4.6))
    ax1 = fig.add_axes([0.03, 0.09, 0.30, 0.78])
    ax2 = fig.add_axes([0.36, 0.09, 0.28, 0.78])
    ax3 = fig.add_axes([0.69, 0.09, 0.29, 0.78])

    ax1.plot(W.real, W.imag, "-", color="0.8", lw=0.9)
    lo, hi = min(s0, s0 + s1) * (m - 1), max(s0, s0 + s1) * (m - 1)
    k = np.linspace(lo, hi, 400).astype(int).clip(0, m - 1)
    ax1.plot(W[k].real, W[k].imag, "-", color="#2c9c4a", lw=2.4,
             label=f"covered arc  s=[{s0:.3f},{s0 + s1:.3f}]")
    ax1.set_aspect("equal"); ax1.legend(fontsize=7, loc="lower right")
    ax1.set_title("m=800 curve and the piece the copy covers", fontsize=9)
    ax1.tick_params(labelsize=7)

    ax2.plot(copy.real, copy.imag, "-", color="#2c9c4a", lw=2.0, label="copy of the whole")
    ax2.plot(w.real, w.imag, "-", color="#1f6fd0", lw=1.3, label="m=40 curve")
    ax2.set_aspect("equal"); ax2.legend(fontsize=7, loc="lower right")
    ax2.set_title(f"overlay, |c|={abs(fit['c']):.3f}, arg c={np.angle(fit['c'], deg=True):.1f}°", fontsize=9)
    ax2.tick_params(labelsize=7)

    q_nat = core._interp(Wg, a / (m - 1) + (b - a) / (m - 1) * t)
    res_nat = np.abs(w - (nat["c"] * q_nat + nat["d"]))
    ax3.plot(t, res, color="#c8342a", lw=1.2, label=f"fitted copy: sup {fit['sup']:.4f}")
    ax3.plot(t, res_nat, color="0.6", lw=0.9, ls="--",
             label=f"index correspondence: sup {nat['sup']:.4f}")
    ax3.set_xlabel("t on the m=40 curve", fontsize=8)
    ax3.set_ylabel("|w - copy|", fontsize=8)
    ax3.legend(fontsize=7); ax3.tick_params(labelsize=7)
    ax3.set_title("residual", fontsize=9)

    fig.suptitle(f"m=40 optimum fitted into the m=800 curve: sup {fit['sup']:.4f} "
                 f"(index correspondence {nat['sup']:.4f}), |c|={abs(fit['c']):.2f}", fontsize=9, y=0.99)
    out = ROOT / "figures/align-app-copy-preview.png"
    fig.savefig(out, dpi=110, facecolor="white", bbox_inches="tight", pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")
    print(f"s0={s0:.4f} s1={s1:.4f} |c|={abs(fit['c']):.4f} sup={fit['sup']:.5f} "
          f"index sup={nat['sup']:.5f}")


if __name__ == "__main__":
    main()
