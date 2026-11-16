"""How far is a configuration from the stress (KKT) balance, in units of the typical force?

The certificate asks for weights lam on the tight-max pairs and mu on the tight-min pairs with
   sum_a lam_a A_a = sum_b mu_b A_b,   sum lam = sum mu = 1,   lam, mu >= 0,
where A_a is the gradient of log r_a w.r.t. the free coordinates.  `slp2.stress` tests
feasibility exactly.  When it fails we want the size of the failure, so here we minimise the
L1 norm of the residual and divide by the mean row norm of A -- a dimensionless number.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".ultramath/smwgx79ya/m8"))
sys.path.insert(0, str(ROOT / "scripts"))
from core2 import Problem            # noqa: E402
from slp2 import grad_rows           # noqa: E402


def residual(P, v, tol=1e-8):
    A, L = grad_rows(P, v)
    A = A.tocsr()
    mx, mn = L.max(), L.min()
    imax = np.where(L > mx - tol)[0]
    imin = np.where(L < mn + tol)[0]
    Am, Bm = A[imax], A[imin]
    nfree = P.nfree
    na, nb = len(imax), len(imin)
    # variables: lam(na), mu(nb), s+(nfree), s-(nfree)
    nv = na + nb + 2 * nfree
    c = np.concatenate([np.zeros(na + nb), np.ones(2 * nfree)])
    Aeq = sparse.vstack([
        sparse.hstack([sparse.csr_matrix(np.ones((1, na))), sparse.csr_matrix((1, nb)),
                       sparse.csr_matrix((1, 2 * nfree))]),
        sparse.hstack([sparse.csr_matrix((1, na)), sparse.csr_matrix(np.ones((1, nb))),
                       sparse.csr_matrix((1, 2 * nfree))]),
        sparse.hstack([Am.T, -Bm.T, -sparse.identity(nfree), sparse.identity(nfree)]),
    ]).tocsr()
    beq = np.concatenate([[1.0, 1.0], np.zeros(nfree)])
    res = linprog(c, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    if res.status != 0:
        return None, na, nb
    s = res.x[na + nb:]
    resid = s[:nfree] + s[nfree:]
    rows = np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
    scale = rows[imax].mean()
    return float(np.sum(resid) / (2 * scale)), na, nb


def main():
    import explicit_curves as ec
    import align_core as core
    print("stress residual (L1 residual / mean row norm); 0 means a certificate exists")
    print(f"{'configuration':>28} {'m':>6} {'D':>10} {'residual':>10} {'#max':>6} {'#min':>6}")
    for k in (2, 3, 4):
        z = ec.koch_nodes(k)
        v = (z - z[0]) / (z[1] - z[0])
        P = Problem(len(v), gauge="v2")
        r, na, nb = residual(P, v)
        print(f"{'Koch level ' + str(k):>28} {len(v):6d} {core.distortion(v):10.6f} "
              f"{r:10.4f} {na:6d} {nb:6d}")
    for name, path in (("S branch m=200 (certified)", "results/m200_S.npy"),
                       ("S branch m=400 (certified)", "results/m400_converged.npy"),
                       ("S branch m=1600 (not cert.)", ".ultramath/smwhdqmht/m1/s_m1600_plain.npy")):
        p = ROOT / path
        if not p.exists():
            continue
        v = core.gauge(np.load(p))
        m = len(v)
        if m > 800:
            print(f"{name:>28} {m:6d} {core.distortion(v):10.6f}     (LP too large)")
            continue
        P = Problem(m, gauge="v2")
        r, na, nb = residual(P, v)
        print(f"{name:>28} {m:6d} {core.distortion(v):10.6f} {r:10.5f} {na:6d} {nb:6d}")


if __name__ == "__main__":
    main()
