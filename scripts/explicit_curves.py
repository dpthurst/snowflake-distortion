"""The two explicit self-similar arcs that beat the Koch distortion.

Koch: D = sqrt(633) * 23^-alpha = 2.096805160875 (sharp, two-sided).

Two-piece reversible arc (from the earlier swarm, audited): for x <= t,
gamma(x) = b0 * conj(gamma(x/t)); for x >= t, gamma(x) = 1 + b1*gamma((1-x)/(1-t)),
with |b0| = t^alpha, |b1| = (1-t)^alpha and b0 - b1 = 1; t = 0.563 gives the
champion, D = 2.04204 (level 16 in the report).

Three-piece reversible arc (closed form, no free parameters): symmetric weights
w = (sqrt2-1, 3-2 sqrt2, sqrt2-1), signs sigma=(+1,-1,-1), reflections
eps=(1,0,0) as in the source note, D = 2.017935389433 (levels 9-11).

Both are numerical constructions: the quoted D is a deep-level measurement, not a
certificate.
"""

import numpy as np

ALPHA = np.log(3) / np.log(4)


def two_piece(t=0.563, level=12):
    """nodes (parameter, position) of the two-piece reversible arc"""
    r0, r1 = t ** ALPHA, (1 - t) ** ALPHA
    # b0 - b1 = 1 with |b0| = r0, |b1| = r1: solve for the argument of b0
    # b0 = r0 e^{i phi}, b1 = b0 - 1, |b1| = r1  ->  |r0 e^{i phi} - 1|^2 = r1^2
    cosphi = (r0 ** 2 + 1 - r1 ** 2) / (2 * r0)
    phi = np.arccos(np.clip(cosphi, -1, 1))
    b0 = r0 * np.exp(-1j * phi)
    b1 = b0 - 1
    pts = [(0.0, 0j)]
    params = [0.0]
    pos = [0j]
    for _ in range(level):
        new_p, new_z = [], []
        for p, z in zip(params, pos):
            new_p.append(t * p);                     new_z.append(b0 * np.conj(z))
            new_p.append(1 - (1 - t) * p);           new_z.append(1 + b1 * z)
        params, pos = new_p, new_z
    order = np.argsort(params)
    return np.array(params)[order], np.array(pos)[order], dict(b0=b0, b1=b1, t=t, phi=phi)


def three_piece(level=8):
    """nodes of the closed-form three-piece reversible arc"""
    w = np.array([np.sqrt(2) - 1, 3 - 2 * np.sqrt(2), np.sqrt(2) - 1])
    sigma = np.array([+1, -1, -1])
    eps = np.array([1, 0, 0])
    r = w ** ALPHA
    b = np.zeros(3, dtype=complex)
    b[1] = -r[1]
    re0 = (1 - r[1]) / 2
    im0 = np.sqrt(max(r[0] ** 2 - re0 ** 2, 0.0))
    b[0] = re0 + 1j * im0
    b[2] = -np.conj(b[0])
    A = np.zeros(4, dtype=complex)
    for i in range(3):
        A[i + 1] = A[i] + sigma[i] * b[i]
    C = np.array([A[i] if sigma[i] == 1 else A[i] - b[i] for i in range(3)])
    a = np.array([0.0, w[0], w[0] + w[1]])

    def phi(i, s):
        return a[i] + w[i] * s if sigma[i] == 1 else a[i] + w[i] * (1 - s)

    def S(i, z):
        return C[i] + b[i] * (np.conj(z) if eps[i] else z)

    params, pos = [np.array([0.0])], [np.array([0j])]
    for _ in range(level):
        new_p, new_z = [], []
        for i in range(3):
            new_p.append(phi(i, params[-1]))
            new_z.append(S(i, pos[-1]))
        params = [np.concatenate(new_p)]
        pos = [np.concatenate(new_z)]
    p, z = params[0], pos[0]
    order = np.argsort(p)
    return p[order], z[order], dict(w=w, sigma=sigma, eps=eps, b=b, A=A)


def distortion_of_nodes(z, p):
    n = len(z)
    i, j = np.triu_indices(n, 1)
    gap = np.abs(p[i] - p[j])
    keep = gap > 0
    r = np.abs(z[i] - z[j])[keep] / gap[keep] ** ALPHA
    return float(r.max() / r.min()), float(r.min()), float(r.max())


def koch_nodes(level):
    z = np.array([0 + 0j, 1 + 0j])
    for _ in range(level):
        out = [z[0]]
        for A, B in zip(z[:-1], z[1:]):
            d = (B - A) / 3
            out += [A + d, A + d + d * np.exp(1j * np.pi / 3), A + 2 * d, B]
        z = np.array(out)
    return z


if __name__ == "__main__":
    print(f"alpha = {ALPHA:.10f}, Koch D = {np.sqrt(633)*23**-ALPHA:.12f}")
    for K in (4, 6, 8, 10):
        p, z, _ = two_piece(level=K)
        D, lo, hi = distortion_of_nodes(z, p)
        print(f"  two-piece   K={K:2d} m={len(z):6d}: D={D:.9f}  r in [{lo:.6f},{hi:.6f}]")
    for K in (4, 5, 6, 7):
        p, z, _ = three_piece(level=K)
        D, lo, hi = distortion_of_nodes(z, p)
        print(f"  three-piece K={K:2d} m={len(z):6d}: D={D:.12f}  r in [{lo:.6f},{hi:.6f}]")
    for K in (4, 5):
        z = koch_nodes(K)
        D, lo, hi = distortion_of_nodes(z, np.arange(len(z)) / (4 ** K))
        print(f"  Koch        K={K:2d} m={len(z):6d}: D={D:.12f}")
