"""SLP for the exact minimax with row-normalised gradients.

slp2.slp uses the raw gradient rows A_i = d log r_i / d(positions), whose norms are 1/chord and so
spread over a factor equal to the chord ratio -- 10^3 at i.i.d. spacings.  The trust region then
means something different from pair to pair.  Here each row is divided by its norm, which makes the
trust region uniform in "relative change per pair", and the step is then converted back.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
from slp2 import grad_rows            # noqa: E402


def slp_scaled(P, v, rho0=0.02, iters=300, verbose=False):
    x = P.pack(v)
    v = P.unpack(x)
    r = P.logrange(v)
    rho = rho0
    n = P.nfree
    for it in range(iters):
        A, L = grad_rows(P, v)
        A = A.tocsr()
        nrm = np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
        nrm[nrm == 0] = 1.0
        S = A.multiply(1.0 / nrm[:, None]).tocsr()          # normalised rows
        # minimise t  s.t.  S d + t >= -L,  -S d + t >= L,  |d| <= rho,  t >= 0
        npair = len(L)
        I = sparse.identity(n)
        A_ub = sparse.vstack([
            sparse.hstack([S, -sparse.csr_matrix(np.ones((npair, 1)))]),
            sparse.hstack([-S, -sparse.csr_matrix(np.ones((npair, 1)))]),
            sparse.hstack([I, sparse.csr_matrix((n, 1))]),
            sparse.hstack([-I, sparse.csr_matrix((n, 1))]),
        ]).tocsr()
        b_ub = np.concatenate([-L, L, rho * np.ones(n), rho * np.ones(n)])
        c = np.zeros(n + 1); c[-1] = 1.0
        res = linprog(c, A_ub=A_ub, b_ub=b_ub,
                      bounds=[(None, None)] * n + [(0, None)], method="highs")
        if res.status != 0:
            rho *= 0.5
            if rho < 1e-12:
                break
            continue
        d = res.x[:n]
        step, improved = 1.0, False
        for _ in range(20):
            vn = P.unpack(x + step * d)
            rn = P.logrange(vn)
            if rn < r - 1e-16:
                improved = True
                break
            step *= 0.5
        if not improved:
            rho *= 0.5
            if rho < 1e-12:
                break
            continue
        x, v, r = x + step * d, vn, rn
        rho = min(rho * 1.5, 0.2)
        if verbose and it % 20 == 0:
            print(f"    it {it}: logrange {r:.12f} rho {rho:.2e}", flush=True)
    return v, r
