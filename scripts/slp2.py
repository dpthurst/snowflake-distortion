"""SLP for the exact minimax, gauge-agnostic (works with core2.Problem)."""
import numpy as np
from scipy.optimize import linprog
from scipy import sparse
from core2 import ALPHA


def grad_rows(P, v):
    dv, d2 = P._geom(v)
    L = 0.5 * np.log(d2)[P.iu, P.ju] - ALPHA * P.loggap
    G = dv / d2
    g = G[P.iu, P.ju]
    n = len(P.free_idx)
    rows = np.arange(P.npair)
    data, ii, jj = [], [], []
    for off, comp in ((0, g.real), (n, g.imag)):
        for a, sgn in ((P.iu, 1.0), (P.ju, -1.0)):
            pa = P.pos[a]
            keep = pa >= 0
            data.append(sgn * comp[keep]); ii.append(rows[keep]); jj.append(off + pa[keep])
    data = np.concatenate(data); ii = np.concatenate(ii); jj = np.concatenate(jj)
    return sparse.coo_matrix((data, (ii, jj)), shape=(P.npair, P.nfree)).tocsr(), L


def slp(P, v, rho0=0.01, iters=200, verbose=0):
    x = P.pack(v); v = P.unpack(x)
    r = P.logrange(v); rho = rho0
    nfree, npair = P.nfree, P.npair
    c = np.concatenate([np.zeros(nfree), [1.0, 1.0]])
    for it in range(iters):
        A, L = grad_rows(P, v)
        em = sparse.csr_matrix(-np.ones((npair, 1))); Z = sparse.csr_matrix((npair, 1))
        A_ub = sparse.vstack([sparse.hstack([A, em, Z]),
                              sparse.hstack([-A, Z, em])]).tocsr()
        b_ub = np.concatenate([-L, L])
        bounds = [(-rho, rho)] * nfree + [(0, None), (0, None)]
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if res.status != 0:
            rho *= 0.5
            if rho < 1e-13: break
            continue
        delta = res.x[:nfree]
        step = 1.0; improved = False
        for _ in range(14):
            y = x + step * delta; vn = P.unpack(y); rn = P.logrange(vn)
            if rn < r - 1e-15:
                improved = True; break
            step *= 0.5
        if not improved:
            rho *= 0.5
            if rho < 1e-13: break
            continue
        x, v, r = y, vn, rn
        rho = min(rho * 2.0, 0.2)
        if verbose and it % 10 == 0: print(f"   slp it={it} logrange={r:.12f} rho={rho:.2e}")
    return v, r


def stress(P, v, tol=1e-8):
    A, L = grad_rows(P, v)
    mx, mn = L.max(), L.min()
    imax = np.where(L > mx - tol)[0]; imin = np.where(L < mn + tol)[0]
    Am = A[imax].toarray(); Bm = A[imin].toarray()
    na, nb = len(imax), len(imin); nv = na + nb
    Aeq = np.zeros((2 + P.nfree, nv)); beq = np.zeros(2 + P.nfree)
    Aeq[0, :na] = 1; beq[0] = 1; Aeq[1, na:] = 1; beq[1] = 1
    Aeq[2:, :na] = Am.T; Aeq[2:, na:] = -Bm.T
    res = linprog(np.zeros(nv), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    ok = res.status == 0
    return dict(status=res.status, imax=imax, imin=imin,
                lam=res.x[:na] if ok else None, mu=res.x[na:] if ok else None,
                pairs_max=list(zip(P.iu[imax], P.ju[imax])),
                pairs_min=list(zip(P.iu[imin], P.ju[imin])))


def tight_sets(P, v, tol=1e-8):
    A, L = grad_rows(P, v)
    mx, mn = L.max(), L.min()
    return A, L, np.where(L > mx - tol)[0], np.where(L < mn + tol)[0]


def cone_dim(P, v, tol=1e-8, K=30, seed=0):
    """dimension of the first-order flat cone {d: g_a.d<=0 (max), g_b.d>=0 (min)}"""
    A, L, ia, ib = tight_sets(P, v, tol)
    Ad = A.toarray(); n = P.nfree
    A_ub = np.vstack([Ad[ia], -Ad[ib]]); b_ub = np.zeros(len(ia) + len(ib))
    rng = np.random.default_rng(seed); rays = []
    for k in range(K):
        c = rng.normal(size=n)
        res = linprog(-c, A_ub=A_ub, b_ub=b_ub, bounds=[(-1, 1)] * n, method="highs")
        if res.status == 0 and np.abs(res.x).max() > 1e-9:
            rays.append(res.x / np.abs(res.x).max())
    if not rays: return 0, np.zeros((0, n))
    R = np.array(rays); s = np.linalg.svd(R, compute_uv=False)
    return int((s > 1e-7 * s[0]).sum()), R


def multiplier_rank(P, v, tol=1e-8):
    A, L, ia, ib = tight_sets(P, v, tol)
    st = stress(P, v, tol)
    if st['lam'] is None: return None
    Ad = A.toarray()
    sup = np.concatenate([ia[st['lam'] > 1e-9], ib[st['mu'] > 1e-9]])
    s = np.linalg.svd(Ad[sup], compute_uv=False)
    return dict(na=len(ia), nb=len(ib), support=len(sup), rank=int((s > 1e-8).sum()),
                nfree=P.nfree, lam=st['lam'], mu=st['mu'], ia=ia, ib=ib)


def face_spectrum(P, v, tol=1e-8, max_lp=400):
    """min over faces d_j=+-1 of h(d); h(d)=max_A g.d - min_B g.d.  Sorted."""
    A, L, ia, ib = tight_sets(P, v, tol)
    Ad = A.toarray(); n = P.nfree
    na, nb = len(ia), len(ib)
    A_ub = np.zeros((na + nb, n + 2))
    A_ub[:na, :n] = Ad[ia]; A_ub[:na, n] = -1
    A_ub[na:, :n] = -Ad[ib]; A_ub[na:, n + 1] = -1
    bb = np.zeros(na + nb); c = np.r_[np.zeros(n), 1.0, 1.0]
    out = []
    idxs = range(n) if n <= max_lp // 2 else np.linspace(0, n - 1, max_lp // 2).astype(int)
    for j in idxs:
        for sgn in (1.0, -1.0):
            bounds = [(-1, 1)] * n + [(0, None), (0, None)]
            bounds[j] = (sgn, sgn)
            res = linprog(c, A_ub=A_ub, b_ub=bb, bounds=bounds, method="highs")
            if res.status == 0:
                d = res.x[:n]
                out.append((res.fun, d / np.abs(d).max()))
    out.sort(key=lambda z: z[0])
    return out
