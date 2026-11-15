"""Reference implementation of the geometry used by the alignment app.

``scripts/test_align_app.py`` runs the JavaScript core embedded in the app and
compares it with these functions, so the two must stay in step.  Everything
works on a curve ``v`` of ``m`` complex nodes with uniform index spacing; ``t``
is the normalised index parameter.

Conventions, as in ``results/``:

* endpoint gauge ``(v - v[0]) / (v[-1] - v[0])``
* block ``[a,b]``: ``w(t) = (V(a + t (b-a)) - v[a]) / (v[b] - v[a])``
* variants ``id``, ``conj``, ``rev``, ``revconj`` with ``rev(t) = 1 - W(1-t)``
* residual: ``min over variants of sup_t |w(t) - g(t)|``.  Both curves are
  gauged to ``0 -> 1``, so the unit is the endpoint chord.  The sup is exact on
  the union of the two knot sets, where the difference is affine in ``t``.
* RMS: over 4001 uniform ``t`` (stated because it is a sampling choice).
* free fit: complex least squares mapping the block onto the whole.
"""

import numpy as np

ALPHA = np.log(3.0) / np.log(4.0)
VARIANTS = ("id", "conj", "rev", "revconj")
N_RMS = 4001


def gauge(v):
    v = np.asarray(v, dtype=complex)
    return (v - v[0]) / (v[-1] - v[0])


def _interp(arr, t):
    arr = np.asarray(arr, dtype=complex)
    x = np.asarray(t, dtype=float) * (len(arr) - 1)
    k = np.clip(np.floor(x).astype(int), 0, len(arr) - 2)
    f = x - k
    return arr[k] * (1 - f) + arr[k + 1] * f


def block(v, a, b, t):
    v = np.asarray(v, dtype=complex)
    idx = a + np.asarray(t, dtype=float) * (b - a)
    return (_interp(v, idx / (len(v) - 1)) - v[a]) / (v[b] - v[a])


def whole_at(W, t, kind):
    """Evaluate the gauged whole ``W`` in one of the four variants."""
    t = np.asarray(t, dtype=float)
    if kind == "id":
        return _interp(W, t)
    if kind == "conj":
        return np.conj(_interp(W, t))
    if kind == "rev":
        return 1 - _interp(W, 1 - t)
    if kind == "revconj":
        return 1 - np.conj(_interp(W, 1 - t))
    raise ValueError(kind)


def knots(m, L):
    return np.unique(np.concatenate([np.arange(L + 1) / L, np.arange(m) / (m - 1)]))


def residuals(whole, a, b, source=None):
    """Per-variant (sup, rms) residuals for the block ``[a,b]``."""
    v = np.asarray(whole, dtype=complex)
    src = v if source is None else np.asarray(source, dtype=complex)
    m, L = len(v), b - a
    W = gauge(v)
    t = knots(m, L)
    w = block(src, a, b, t)
    td = np.linspace(0.0, 1.0, N_RMS)
    wd = block(src, a, b, td)
    out = {}
    for kind in VARIANTS:
        e = np.abs(w - whole_at(W, t, kind))
        ed = np.abs(wd - whole_at(W, td, kind))
        out[kind] = dict(sup=float(e.max()), rms=float(np.sqrt(np.mean(ed ** 2))))
    return out


def procrustes(whole, a, b, kind="id", source=None):
    """Free complex similarity mapping the block onto the whole variant."""
    v = np.asarray(whole, dtype=complex)
    src = v if source is None else np.asarray(source, dtype=complex)
    t = np.linspace(0.0, 1.0, N_RMS)
    z = block(src, a, b, t)
    u = whole_at(gauge(v), t, kind)
    zc, uc = z - z.mean(), u - u.mean()
    c = np.vdot(zc, uc) / np.vdot(zc, zc).real
    return complex(c), float(np.abs(c * zc - uc).max())


def distortion(v):
    v = gauge(v)
    n = len(v)
    i, j = np.triu_indices(n, 1)
    r = np.abs(v[i] - v[j]) / (j - i) ** ALPHA
    return float(r.max() / r.min())


def ladder(v, L):
    """``r_k = |v[k+L] - v[k]| / L^alpha`` for every start index."""
    v = gauge(v)
    return np.abs(v[L:] - v[:-L]) / L ** ALPHA


def interval_report(whole, a, b, source=None):
    v = np.asarray(whole, dtype=complex)
    src = v if source is None else np.asarray(source, dtype=complex)
    m, L = len(v), b - a
    res = residuals(v, a, b, src)
    best = min(VARIANTS, key=lambda k: res[k]["sup"])
    c, free_sup = procrustes(v, a, b, best, src)
    seg = src[a:b + 1]
    dia = float(np.abs(seg[:, None] - seg[None, :]).max())
    r_ab = abs(src[b] - src[a]) / L ** ALPHA
    R_full = abs(v[-1] - v[0]) / (m - 1) ** ALPHA
    rep = dict(best=best, sup=res[best]["sup"], rms=res[best]["rms"],
               sup_all={k: res[k]["sup"] for k in VARIANTS},
               free_sup=free_sup, free_abs_c=abs(c), free_arg_c=float(np.angle(c)),
               D_whole=distortion(v), D_block=distortion(seg),
               s_block=dia / abs(src[b] - src[a]),
               r_over_R=float(r_ab / R_full))
    if source is None:
        lad = ladder(v, L)
        others = lad[np.arange(lad.size) != a]
        rep["r_percentile"] = float((others < r_ab).mean())
    return rep


# ---------------------------------------------------------------- copy mode
# A similarity image of the whole with an affine reparametrisation:
#   w(t) ~ c W(s0 + s1 t) + d,   c, d in C,  [s0, s0+s1] inside [0, 1].

def sample_block(v, a, b, n=301):
    t = np.linspace(0.0, 1.0, n)
    return t, block(v, a, b, t)


def copy_residual(whole, a, b, source, s0, s1, c, d, t=None):
    """sup and rms of w(t) - (c W(s0+s1 t) + d) at the given parameters."""
    if t is None:
        t = np.linspace(0.0, 1.0, 301)
    w = block(source, a, b, t)
    q = _interp(gauge(whole), s0 + s1 * t)
    r = np.abs(w - (c * q + d))
    return float(r.max()), float(np.sqrt(np.mean(r ** 2)))


def fit_similarity(whole, a, b, source, s0, s1, t=None):
    """Least-squares complex c, d for a fixed parameter correspondence."""
    if t is None:
        t = np.linspace(0.0, 1.0, 301)
    w = block(source, a, b, t)
    q = _interp(gauge(whole), s0 + s1 * t)
    qc, wc = q - q.mean(), w - w.mean()
    c = np.vdot(qc, wc) / np.vdot(qc, qc).real
    return complex(c), complex(w.mean() - c * q.mean())


def fit_copy_reference(whole, a, b, source, power=12.0, n=301):
    """Reference six-parameter fit with scipy's Nelder-Mead on the soft maximum."""
    from scipy.optimize import minimize

    t, w = sample_block(source, a, b, n)
    W = gauge(whole)
    m = len(np.asarray(whole))
    nat = (a / (m - 1), (b - a) / (m - 1))
    c0, d0 = fit_similarity(whole, a, b, source, *nat, t=t)

    def soft(p):
        q = _interp(W, p[4] + p[5] * t)
        r = np.abs(w - ((p[0] + 1j * p[1]) * q + (p[2] + 1j * p[3])))
        pen = 0.0
        lo, hi = min(p[4], p[4] + p[5]), max(p[4], p[4] + p[5])
        if lo < -1e-12:
            pen += 1e3 * (-lo)
        if hi > 1 + 1e-12:
            pen += 1e3 * (hi - 1)
        if abs(p[5]) < 1e-3:
            pen += 1e3 * (1e-3 - abs(p[5]))
        return (np.mean(r ** power)) ** (1 / power) + pen

    x0 = [c0.real, c0.imag, d0.real, d0.imag, nat[0], nat[1]]
    res = minimize(soft, x0, method="Nelder-Mead",
                   options=dict(maxiter=4000, xatol=1e-10, fatol=1e-12, adaptive=True))
    p = res.x
    sup, rms = copy_residual(whole, a, b, source, p[4], p[5],
                             complex(p[0], p[1]), complex(p[2], p[3]),
                             t=np.linspace(0, 1, 1201))
    nat_sup, _ = copy_residual(whole, a, b, source, nat[0], nat[1], c0, d0,
                               t=np.linspace(0, 1, 1201))
    return dict(params=p, sup=sup, rms=rms, natural_sup=nat_sup, nat=nat)


def profile_at(Wg, t, w, s0, s1, power=12.0):
    """Least-squares (c,d) at a fixed correspondence, with sup and soft-max."""
    q = _interp(Wg, s0 + s1 * t)
    qc, wc = q - q.mean(), w - w.mean()
    den = np.vdot(qc, qc).real
    c = np.vdot(qc, wc) / den
    d = w.mean() - c * q.mean()
    r = np.abs(w - (c * q + d))
    return dict(c=c, d=d, sup=float(r.max()),
                soft=float((np.mean(r ** power)) ** (1 / power)))


def fit_copy_multires(whole, a, b, source):
    """The app's copy-mode fit: coarse grid, then two refinements."""
    from scipy.optimize import minimize

    Wg = gauge(whole)
    m = len(np.asarray(whole))

    def penal(s0, s1):
        p = 0.0
        if s0 < 0:
            p += 1e3 * (-s0)
        if s0 + s1 > 1:
            p += 1e3 * (s0 + s1 - 1)
        if abs(s1) < 0.02:
            p += 1e3 * (0.02 - abs(s1))
        return p

    def stage(n, init=None):
        t = np.linspace(0.0, 1.0, n)
        w = block(source, a, b, t)

        def obj(p):
            s1 = min(p[1], 1 - p[0])
            return profile_at(Wg, t, w, p[0], s1)["soft"] + penal(p[0], s1)

        if init is None:
            best = (np.inf, 0.0, 0.25)
            for s0 in np.linspace(0, 0.95, 61):
                for s1 in np.linspace(1 / 24, 1.0, 24):
                    if s0 + s1 > 1:
                        continue
                    v = profile_at(Wg, t, w, s0, s1)["sup"]
                    if v < best[0]:
                        best = (v, s0, s1)
            init = [best[1], best[2]]
        res = minimize(obj, init, method="Nelder-Mead",
                       options=dict(maxiter=800, xatol=1e-10, fatol=1e-12, adaptive=True))
        p = res.x
        return [p[0], min(p[1], 1 - p[0])], t, w

    p, t1, w1 = stage(101)
    p, t2, w2 = stage(301, p)
    s0, s1 = p
    pr = profile_at(Wg, t2, w2, s0, s1)
    nat = profile_at(Wg, t2, w2, a / (m - 1), (b - a) / (m - 1))
    return dict(s0=s0, s1=s1, c=pr["c"], d=pr["d"], sup=pr["sup"],
                natural_sup=nat["sup"], natural_c=nat["c"], natural_d=nat["d"])


# ---------------------------------------------------------------- set distance
def _dist_point_segment(px, py, ax, ay, bx, by):
    vx, vy = bx - ax, by - ay
    wx, wy = px - ax, py - ay
    c1 = vx * wx + vy * wy
    c2 = vx * vx + vy * vy
    t = np.where(c2 > 0, np.clip(c1 / np.where(c2 > 0, c2, 1.0), 0.0, 1.0), 0.0)
    return np.hypot(px - (ax + t * vx), py - (ay + t * vy))


def dist_to_polyline(px, py, PX, PY):
    """Distance from points (px,py) to the polyline through (PX,PY); brute force."""
    px = np.atleast_1d(px)[:, None]
    py = np.atleast_1d(py)[:, None]
    out = np.full(px.shape[0], np.inf)
    for i in range(len(PX) - 1):
        d = _dist_point_segment(px[:, 0], py[:, 0], PX[i], PY[i], PX[i + 1], PY[i + 1])
        out = np.minimum(out, d)
    return out


def set_distance_sup(whole, a, b, source=None, kind="id", n=601):
    """Symmetric supremum of the distance between the two polylines.

    Returns (block_to_curve, curve_to_block, max).
    """
    v = np.asarray(whole, dtype=complex)
    src = v if source is None else np.asarray(source, dtype=complex)
    W = gauge(v)
    t = np.linspace(0.0, 1.0, n)
    w = block(src, a, b, t)
    g = whole_at(W, t, kind)
    tv = np.arange(len(v)) / (len(v) - 1.0)
    gv = whole_at(W, tv, kind)
    bt = np.arange(b - a + 1) / (b - a)
    bv = block(src, a, b, bt)
    fwd = dist_to_polyline(w.real, w.imag, gv.real, gv.imag)
    bwd = dist_to_polyline(g.real, g.imag, bv.real, bv.imag)
    return float(fwd.max()), float(bwd.max()), float(max(fwd.max(), bwd.max()))


def heat_stats(v):
    """min, max, geometric mean and arg-extrema of r(x_i,x_j), x_i = i/(m-1)."""
    v = np.asarray(v, dtype=complex)
    n = len(v)
    x = np.arange(n) / (n - 1)
    i, j = np.triu_indices(n, 1)
    r = np.abs(v[i] - v[j]) / np.abs(x[i] - x[j]) ** ALPHA
    imin = int(r.argmin())
    imax = int(r.argmax())
    lo, hi = float(r[imin]), float(r[imax])
    return dict(lo=lo, hi=hi, gmean=float(np.sqrt(lo * hi)), D=hi / lo,
                imin=[int(i[imin]), int(j[imin])], imax=[int(i[imax]), int(j[imax])])
