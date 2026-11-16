"""The stress measures at m = 200, 400, 800, canonically selected, and their trend.

`slp2.stress` returns *a* feasible point of the multiplier LP, which at these m is a vertex of a
badly degenerate polytope (thousands of tight pairs) and therefore not a meaningful approximation
of anything.  The canonical choice used here is the maximum-entropy point of the feasible set:

    max  -sum z log z    s.t.  M z = b,  z >= 0,

with z = (lam, mu), M = [[1,0],[0,1],[A+, -A-]], b = (1,1,0).  Its dual is smooth and unconstrained,

    min_nu  sum_i exp((M^T nu)_i - 1) - nu.b,

solved by L-BFGS in dimension 2 + nfree.  Unlike a vertex, this selection does not depend on the
order the pairs happen to be listed in.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".ultramath/smwgx79ya/m8"))
sys.path.insert(0, str(ROOT / "scripts"))
from core2 import Problem          # noqa: E402
from slp2 import grad_rows         # noqa: E402
import align_core as core          # noqa: E402

CONFIGS = [(200, "results/m200_S.npy"), (400, "results/m400_converged.npy"),
           (800, "results/m800_converged.npy")]


def tight(P, v, tol=1e-8):
    A, L = grad_rows(P, v)
    mx, mn = L.max(), L.min()
    imax = np.where(L > mx - tol)[0]
    imin = np.where(L < mn + tol)[0]
    return A.tocsr(), imax, imin


def maxent(A, imax, imin, nfree, iters=60, tol=1e-13):
    """maximum-entropy point of the multiplier polytope, by Newton on the (smooth) dual"""
    Am = A[imax].T.tocsr()             # nfree x na
    Bm = A[imin].T.tocsr()
    na, nb = len(imax), len(imin)
    I = sparse.identity
    M = sparse.vstack([
        sparse.hstack([sparse.csr_matrix(np.ones((1, na))), sparse.csr_matrix((1, nb))]),
        sparse.hstack([sparse.csr_matrix((1, na)), sparse.csr_matrix(np.ones((1, nb)))]),
        sparse.hstack([Am, -Bm]),
    ]).tocsr()
    b = np.concatenate([[1.0, 1.0], np.zeros(nfree)])
    nu = np.zeros(2 + nfree)
    for it in range(iters):
        t = M.T @ nu
        z = np.exp(np.clip(t - 1.0, -700, 700))
        val = z.sum() - b @ nu
        grad = M @ z - b
        if np.abs(grad).max() < tol:
            break
        H = (M @ sparse.diags(z) @ M.T).toarray()
        step = np.linalg.solve(H + 1e-14 * np.eye(len(b)), -grad)
        # damped Newton: keep the dual value decreasing
        t0, s = val, 1.0
        for _ in range(40):
            zt = np.exp(np.clip(M.T @ (nu + s * step) - 1.0, -700, 700))
            if zt.sum() - b @ (nu + s * step) < t0 - 1e-4 * s * (grad @ step):
                break
            s *= 0.5
        nu = nu + s * step
    t = M.T @ nu
    z = np.exp(np.clip(t - 1.0, -700, 700))
    return z[:na], z[na:], dict(iters=it + 1, resid=float(np.abs(M @ z - b).max()))


def main():
    out = {}
    for m, path in CONFIGS:
        v = core.gauge(np.load(ROOT / path))
        P = Problem(m, gauge="v2")
        A, imax, imin = tight(P, v)
        lam, mu, info = maxent(A, imax, imin, P.nfree)
        print(f"m={m:4d}: tight {len(imax)}+{len(imin)}  entropy solution: "
              f"constraint residual {info['resid']:.2e}  "
              f"mass on {int((lam > 1e-9).sum())}+{int((mu > 1e-9).sum())} pairs  "
              f"[{info['iters']} Newton steps]", flush=True)
        out[m] = dict(iu_max=P.iu[imax], ju_max=P.ju[imax], lam=lam,
                      iu_min=P.iu[imin], ju_min=P.ju[imin], mu=mu,
                      n_max=len(imax), n_min=len(imin))
    np.savez(ROOT / "results/stress-maxent.npz",
             **{f"{k}_{n}": a for k, d in out.items() for n, a in d.items()})
    print("wrote results/stress-maxent.npz")
    figure(out)


def marginal(m, ii, jj, w):
    """mass on each node, from a measure on pairs"""
    mu = np.zeros(m)
    np.add.at(mu, ii, w)
    np.add.at(mu, jj, w)
    return mu / mu.sum()


def figure(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(5.6, 4.2))
    ax, bz, cz, dz = axes.ravel()
    # (a) where the stress lives, m = 800
    d = out[800]; m = 800
    for ii, jj, w, col, lab in ((d["iu_max"], d["ju_max"], d["lam"], "#1f3f6f", "max pairs"),
                                (d["iu_min"], d["ju_min"], d["mu"], "#c8342a", "min pairs")):
        k = w > 1e-6
        ax.scatter(ii[k] / (m - 1), jj[k] / (m - 1), s=90 * w[k] ** 0.4 + 0.3,
                   color=col, alpha=0.7, linewidth=0, label=lab)
    ax.plot([0, 1], [0, 1], ":", lw=0.6, color="0.6")
    ax.set_xlabel("$i/m$"); ax.set_ylabel("$j/m$")
    ax.set_title("stress at $m=800$", fontsize=7)
    ax.legend(fontsize=5.4, loc="upper left")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    cols = {200: "#9ecae1", 400: "#4292c6", 800: "#08519c"}
    # (b) marginal of the min measure on the parameter
    for mm in (200, 400, 800):
        dd = out[mm]
        mu = marginal(mm, dd["iu_min"], dd["ju_min"], dd["mu"])
        bz.plot(np.arange(mm) / (mm - 1), mu * mm, "-", lw=0.7, color=cols[mm],
                label=f"$m={mm}$")
    bz.set_xlabel("parameter"); bz.set_ylabel("marginal density")
    bz.set_title("marginal of the min measure", fontsize=7)
    bz.legend(fontsize=5.4)
    # (c) the min measure in lattice units
    for mm in (200, 400, 800):
        dd = out[mm]
        sep = dd["ju_min"] - dd["iu_min"]
        h, e = np.histogram(sep, bins=np.arange(0.5, 13.5), weights=dd["mu"])
        cz.plot(0.5 * (e[1:] + e[:-1]), h, "o-", ms=2.2, lw=0.7, color=cols[mm],
                label=f"$m={mm}$")
    cz.set_xlabel("separation $j-i$ (lattice spacings)")
    cz.set_ylabel("mass")
    cz.set_title("min measure is local", fontsize=7)
    cz.legend(fontsize=5.4)
    # (d) the max measure in ordinary units
    for mm in (200, 400, 800):
        dd = out[mm]
        sep = (dd["ju_max"] - dd["iu_max"]) / (mm - 1)
        h, e = np.histogram(sep, bins=10, range=(0, 1), weights=dd["lam"])
        dz.plot(0.5 * (e[1:] + e[:-1]), h, "o-", ms=2.2, lw=0.7, color=cols[mm],
                label=f"$m={mm}$")
    dz.set_xlabel("separation $|i-j|/m$")
    dz.set_ylabel("mass")
    dz.set_title("max measure is not yet settled", fontsize=7)
    dz.legend(fontsize=5.4)
    for a in (bz, cz, dz):
        a.grid(alpha=0.25)
    fig.tight_layout()
    outp = ROOT / "figures/stress-measures.png"
    fig.savefig(outp, dpi=160, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {outp.relative_to(ROOT)} ({outp.stat().st_size/1024:.0f} kB)")

    # convergence diagnostic on the marginals
    print("\nmarginal convergence (L1 distance of 20-bin histograms)")
    prev = None
    for mm in (200, 400, 800):
        d = out[mm]
        mu = marginal(mm, d["iu_min"], d["ju_min"], d["mu"])
        x = np.arange(mm) / (mm - 1)
        h, _ = np.histogram(x, bins=20, range=(0, 1), weights=mu)
        if prev is not None:
            print(f"  m={prev} -> {mm}: L1 = {np.abs(h - prev_h).sum():.4f}")
        prev, prev_h = mm, h
    for mm in (200, 400, 800):
        d = out[mm]
        s = np.abs(d["ju_min"] - d["iu_min"]) / (mm - 1)
        print(f"  m={mm}: tight {d['n_max']}+{d['n_min']}, "
              f"min-pair separations: median {np.median(s):.3f}, "
              f"90th pct {np.percentile(s, 90):.3f}")


if __name__ == "__main__":
    main()
