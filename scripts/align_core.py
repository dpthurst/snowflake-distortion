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
