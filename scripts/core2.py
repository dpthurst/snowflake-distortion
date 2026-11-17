"""Problem class with selectable gauge.

gauge='v2' : v_1=0, v_2=1   (optimisation gauge, anchors the first step)
gauge='vm' : v_1=0, v_m=1   (MEASUREMENT gauge: the free-coordinate sup-norm
                             equals the normalised shape sup-distance, so local
                             stiffness computed there is in shape units)
"""
import numpy as np
from scipy.optimize import minimize

ALPHA = np.log(3.0) / np.log(4.0)


class Problem:
    def __init__(self, m, gauge="v2"):
        self.m = m
        self.gauge = gauge
        iu, ju = np.triu_indices(m, 1)
        self.iu, self.ju = iu, ju
        self.loggap = np.log((ju - iu).astype(float))
        self.npair = len(iu)
        if gauge == "v2":
            self.free_idx = np.arange(2, m)
        elif gauge == "vm":
            self.free_idx = np.arange(1, m - 1)
        else:
            raise ValueError(gauge)
        self.nfree = 2 * len(self.free_idx)
        self.pos = -np.ones(m, int)
        self.pos[self.free_idx] = np.arange(len(self.free_idx))

    def pack(self, v):
        z = np.asarray(v)[self.free_idx]
        return np.concatenate([z.real, z.imag])

    def unpack(self, x):
        m = self.m
        n = len(self.free_idx)
        v = np.zeros(m, dtype=complex)
        if self.gauge == "v2":
            v[1] = 1.0
        else:
            v[m - 1] = 1.0
        v[self.free_idx] = x[:n] + 1j * x[n:]
        return v

    def _geom(self, v):
        dv = v[:, None] - v[None, :]
        d2 = dv.real ** 2 + dv.imag ** 2
        np.fill_diagonal(d2, 1.0)
        return dv, d2

    def Lvec(self, v):
        _, d2 = self._geom(v)
        return 0.5 * np.log(d2)[self.iu, self.ju] - ALPHA * self.loggap

    def D(self, v):
        L = self.Lvec(v)
        return float(np.exp(L.max() - L.min()))

    def logrange(self, v):
        L = self.Lvec(v)
        return float(L.max() - L.min())

    def fg(self, x, beta):
        m = self.m
        v = self.unpack(x)
        dv, d2 = self._geom(v)
        L = 0.5 * np.log(d2)[self.iu, self.ju] - ALPHA * self.loggap
        a = beta * (L - L.max()); ea = np.exp(a); Z1 = ea.sum()
        b = -beta * (L - L.min()); eb = np.exp(b); Z2 = eb.sum()
        F = (beta * L.max() + np.log(Z1)) / beta + (-beta * L.min() + np.log(Z2)) / beta
        w = ea / Z1 - eb / Z2
        W = np.zeros((m, m))
        W[self.iu, self.ju] = w
        W[self.ju, self.iu] = w
        g = (W * (dv / d2)).sum(axis=1)
        gf = g[self.free_idx]
        return F, np.concatenate([gf.real, gf.imag])

    def Fbeta(self, x, beta):
        return self.fg(x, beta)[0]


def optimise(prob, v0, betas, maxiter=400, gtol=1e-11, verbose=False):
    x = prob.pack(v0)
    for beta in betas:
        res = minimize(prob.fg, x, args=(beta,), jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "gtol": gtol, "ftol": 1e-16,
                                "maxcor": 30})
        x = res.x
        if verbose:
            print(f"  beta={beta:9.3g} F={res.fun: .8f} D={prob.D(prob.unpack(x)):.8f}")
    return prob.unpack(x)


SLOW = [0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
MED = [3, 10, 30, 100, 300, 1000, 3000, 10000]
FAST = [30, 100, 300, 1000, 3000, 10000]


def norm_shape(v):
    v = np.asarray(v, dtype=complex)
    return (v - v[0]) / (v[-1] - v[0])


def orbit(v):
    v = np.asarray(v, dtype=complex)
    return [norm_shape(v), norm_shape(np.conj(v)),
            norm_shape(v[::-1]), norm_shape(np.conj(v[::-1]))]


def supdist(A, B):
    a = norm_shape(A)
    return min(float(np.max(np.abs(a - b))) for b in orbit(B))


def supdist_variants(A, B):
    a = norm_shape(A)
    return [float(np.max(np.abs(a - b))) for b in orbit(B)]


def canonical(v):
    cands = orbit(v)
    k = int(np.argmax([np.mean(c.imag) for c in cands]))
    return cands[k]


def save_cfg(path, v):
    v = np.asarray(v, dtype=complex)
    np.savetxt(path, np.column_stack([np.arange(1, len(v) + 1), v.real, v.imag]),
               fmt="%d %.17g %.17g", header="index re im")


def load_cfg(path):
    a = np.loadtxt(path)
    return a[:, 1] + 1j * a[:, 2]
