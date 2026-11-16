"""Is the Koch curve stationary for the distortion functional?

Take the four-piece self-similar family with parameter weights w = (w1,w2,w2,w4), sum 1, and
step lengths r_i = w_i^alpha (which is what makes the sidelengths equal to the parameter
lengths to the power alpha, i.e. the curve lives at the right exponent).  Closure with the
Koch directions (0, theta, -theta, 0) forces cos(theta) = (1 - r1 - r4)/(r2 + r3).  At
w = (1/4,1/4,1/4,1/4) this is exactly the Koch curve (theta = 60 degrees, r = 1/3).

So D(w1,w4) is a two-parameter family of admissible competitors through the Koch curve, and
the sign of its gradient at (1/4,1/4) decides whether the Koch curve can be stationary.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ALPHA = np.log(3) / np.log(4)


def curve(w, level):
    """(parameters, positions) of the level-k polyline for weights w and exponent ALPHA"""
    w = np.asarray(w, float)
    r = w ** ALPHA
    cosa = (1.0 - r[0] - r[3]) / (r[1] + r[2])
    if abs(cosa) > 1:
        return None, None
    th = np.arccos(cosa)
    steps = r * np.array([1.0, np.exp(1j * th), np.exp(-1j * th), 1.0])
    A = np.concatenate([[0.0], np.cumsum(steps)])
    if abs(A[-1] - 1.0) > 1e-12:
        return None, None
    a = np.concatenate([[0.0], np.cumsum(w)])
    P = np.array([0.0, 1.0])
    Z = np.array([0.0 + 0.0j, 1.0 + 0.0j])
    for _ in range(level):
        Ps, Zs = [], []
        for i in range(4):
            Ps.append(a[i] + w[i] * P)
            Zs.append(A[i] + steps[i] * Z)
        P = np.concatenate(Ps)
        Z = np.concatenate(Zs)
    return P, Z


def distortion(P, Z, chunk=512):
    n = len(Z)
    lo, hi = np.inf, -np.inf
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        d = np.abs(Z[s:e, None] - Z[None, :])
        g = np.abs(P[s:e, None] - P[None, :])
        mask = g > 0
        r = d[mask] / g[mask] ** ALPHA
        lo = min(lo, r.min()); hi = max(hi, r.max())
    return hi / lo


def main():
    level = 6
    base = np.array([0.25, 0.25, 0.25, 0.25])
    print(f"level {level}: {4**level + 1} nodes")
    for eps in (0.0, 0.01, -0.01):
        row = []
        for delta in (0.0, 0.01, -0.01):
            w = base.copy()
            w[0] += eps; w[3] += delta; w[1] -= 0.5 * (eps + delta); w[2] -= 0.5 * (eps + delta)
            P, Z = curve(w, level)
            if P is None:
                row.append("   n/a ")
                continue
            row.append(f"{distortion(P, Z):.7f}")
        print(f"  dw1={eps:+.2f}: " + "  ".join(row))
    P, Z = curve(base, level)
    print(f"\nKoch itself at level {level}: D = {distortion(P, Z):.9f}   (level-limit 2.096805161)")


if __name__ == "__main__":
    main()
