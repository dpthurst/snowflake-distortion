"""The first few stages of the two explicit constructions.

Each construction is a substitution: the whole arc is the union of N copies of
itself, copy i scaled by |b_i|, rotated, reflected if eps_i, and traversed
backwards if sigma_i = -1.  Level k applies the substitution k times, giving a
polyline with N^k segments, and consecutive levels refine each other: the
level-(k-1) vertices are vertices of level k.

    python3 scripts/plot_construction_stages.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import explicit_curves as ec   # noqa: E402

COL = ["#1f3f6f", "#c8342a", "#2c9c4a"]


def two_piece_maps(t=0.563):
    """sigma, eps, b, C, w for the two-piece reversible arc of the note"""
    r0, r1 = t ** ec.ALPHA, (1 - t) ** ec.ALPHA
    cosphi = (r0 ** 2 + 1 - r1 ** 2) / (2 * r0)
    b0 = r0 * np.exp(-1j * np.arccos(np.clip(cosphi, -1, 1)))
    b1 = b0 - 1
    sigma, eps = [+1, -1], [1, 0]
    b = [b0, b1]
    A = np.zeros(3, dtype=complex)
    for i in range(2):
        A[i + 1] = A[i] + sigma[i] * b[i]
    C = [A[i] if sigma[i] == 1 else A[i] - b[i] for i in range(2)]
    return sigma, eps, b, C, np.array([t, 1 - t])


def three_piece_maps():
    """sigma, eps, b, C, w for the closed-form three-piece reversible arc"""
    w = np.array([np.sqrt(2) - 1, 3 - 2 * np.sqrt(2), np.sqrt(2) - 1])
    sigma, eps = [+1, -1, -1], [1, 0, 0]
    r = w ** ec.ALPHA
    b = np.zeros(3, dtype=complex)
    b[1] = -r[1]
    re0 = (1 - r[1]) / 2
    b[0] = re0 + 1j * np.sqrt(max(r[0] ** 2 - re0 ** 2, 0.0))
    b[2] = -np.conj(b[0])
    A = np.zeros(4, dtype=complex)
    for i in range(3):
        A[i + 1] = A[i] + sigma[i] * b[i]
    C = [A[i] if sigma[i] == 1 else A[i] - b[i] for i in range(3)]
    return sigma, eps, b, C, w


def chain(maps, level):
    """the level-k polyline in increasing parameter.

    Returns (P, Z, cuts): the parameter and position of every vertex, and the
    index at which each first-level piece starts (N+1 of them).  Piece i covers
    the parameter interval [a_i, a_i + w_i], a the running sum of the weights.
    """
    sigma, eps, b, C, w = maps
    N = len(b)
    if level == 0:
        return np.array([0.0, 1.0]), np.array([0 + 0j, 1 + 0j]), [0, 1]
    P0, Z0, _ = chain(maps, level - 1)
    a = np.concatenate(([0.0], np.cumsum(w)))
    Ps, Zs, cuts = [], [], [0]
    for i in range(N):
        if sigma[i] == 1:
            P, Z = a[i] + w[i] * P0, Z0
        else:                                   # reversed: the copy is run backwards
            P, Z = a[i] + w[i] * (1 - P0[::-1]), Z0[::-1]
        Zs.append(C[i] + b[i] * (np.conj(Z) if eps[i] else Z))
        Ps.append(P)
        cuts.append(cuts[-1] + len(P))
    P, Z = np.concatenate(Ps), np.concatenate(Zs)
    keep = np.concatenate(([True], np.diff(P) > 0))   # drop the repeated joint
    starts = [int(np.searchsorted(P[keep], a[i] - 1e-15)) for i in range(N)]
    return P[keep], Z[keep], starts + [int(keep.sum()) - 1]


def distortion(P, Z):
    i, j = np.triu_indices(len(Z), 1)
    r = np.abs(Z[i] - Z[j]) / np.abs(P[i] - P[j]) ** ec.ALPHA
    return float(r.max() / r.min()), float(r.min()), float(r.max())


def check_against_explicit_curves(maps, level, fun):
    """the curves must be the ones the note's formulas produce"""
    P, Z, _ = chain(maps, level)
    p, z, _ = fun(level=level)
    d = np.abs(z[:, None] - Z[None, :]).min(axis=1).max()
    assert d < 1e-9, f"level {level}: {d}"


def draw(ax, maps, level, flip, finest):
    P, Z, cuts = chain(maps, level)
    Z = flip(Z)
    for i in range(len(cuts) - 1):
        s = slice(cuts[i], cuts[i + 1] + 1)
        ax.plot(Z[s].real, Z[s].imag, "-", lw=0.3 if finest else 0.9,
                color=COL[i % 3], solid_capstyle="round")
    if level > 1:                       # the level-1 pattern, for orientation
        Z1 = flip(chain(maps, 1)[1])
        ax.plot(Z1.real, Z1.imag, ":", lw=0.7, color="0.6", zorder=0)
    if not finest:
        ax.plot(Z.real, Z.imag, "o", ms=1.5, color="0.2", zorder=3)
    if level == 1:
        for i in range(len(cuts) - 1):
            A, B = Z[cuts[i]], Z[cuts[i + 1]]
            m, d = (A + B) / 2, 0.14 * (B - A) / abs(B - A)
            ax.annotate("", xy=((m + d).real, (m + d).imag),
                        xytext=((m - d).real, (m - d).imag), zorder=4,
                        arrowprops=dict(arrowstyle="-|>", color=COL[i % 3],
                                        lw=1.0, mutation_scale=9))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.margins(0.10)


def main():
    two, three = two_piece_maps(), three_piece_maps()
    rows = [(two, "two-piece\nreversible", 9, ec.two_piece),
            (three, "three-piece\nreversible", 6, ec.three_piece)]
    fig, ax = plt.subplots(2, 4, figsize=(6.9, 3.55))
    for r, (maps, name, deep, fun) in enumerate(rows):
        N = len(maps[2])
        for k in (1, 2, 3, deep):
            check_against_explicit_curves(maps, k, fun)
        flip = (lambda z: np.conj(z)) if chain(maps, 1)[1].imag.mean() < 0 \
            else (lambda z: z)
        Ds = []
        for c, k in enumerate((1, 2, 3)):
            Ds.append(distortion(*chain(maps, k)[:2]))
            draw(ax[r, c], maps, k, flip, finest=False)
            ax[r, c].set_title(f"level {k}   ${N}^{k}$ segments\n"
                               f"$D={Ds[-1][0]:.5f}$", fontsize=6.5)
        draw(ax[r, 3], maps, deep, flip, finest=True)
        D = distortion(*chain(maps, deep)[:2])
        ax[r, 3].set_title(f"level {deep}   ${N}^{{{deep}}}$ segments\n"
                           f"$D={D[0]:.5f}$", fontsize=6.5)
        ax[r, 0].set_ylabel(name, fontsize=6.5)
        print(f"  {name.split(chr(10))[0]:11s} " +
              "  ".join(f"k={k}: {distortion(*chain(maps, k)[:2])[0]:.9f}"
                        for k in (1, 2, 3, 4, 6, deep)))
    fig.suptitle("the substitution rule: every segment is replaced by a copy of the whole, "
                 "scaled, rotated, reflected and/or reversed\n"
                 "(colour: the first-level piece;  arrows: its traversal direction;  "
                 "dotted grey: the level-1 pattern)", fontsize=6.8, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = ROOT / "figures/construction-stages.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
