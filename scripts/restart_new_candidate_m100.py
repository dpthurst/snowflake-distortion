"""Restart the 100-point discrete optimizer from the three-piece IFS curve.

Samples t_j=j/(m-1), NOT equally spaced indices of the IFS grid.
The saved curve coefficient file supplies b, weights and traversal flags.
Run from the project root with python3 -u scripts/restart_new_candidate_m100.py.
"""
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
from core2 import Problem, norm_shape, supdist
from slp2 import slp, stress, grad_rows

ALPHA = np.log(3.0) / np.log(4.0)
DATA = np.load(ROOT / 'data/inputs/coords_3piece_sym.npz')
w, b, signs, refl = DATA['w'], DATA['b'], DATA['signs'], DATA['refl']
a = np.r_[0.0, np.cumsum(w)[:-1]]
A = np.r_[0j, np.cumsum(signs * b)[:-1]]
C = np.where(signs > 0, A, A - b)


def gamma(x, depth=64):
    """Evaluate the IFS curve by composing maps along the parameter address."""
    scalar = np.asarray(x).ndim == 0
    x = np.atleast_1d(np.asarray(x, dtype=float)).copy()
    c0, c1 = np.zeros(len(x), complex), np.ones(len(x), complex)
    parity = np.zeros(len(x), bool)
    for _ in range(depth):
        i = (x > w[0]).astype(int) + (x > w[0] + w[1]).astype(int)
        s = np.where(signs[i] > 0, (x - a[i]) / w[i],
                     (a[i] + w[i] - x) / w[i])
        x = np.clip(s, 0.0, 1.0)
        bb, cc, ee = b[i], C[i], refl[i]
        n0 = c0 + c1 * np.where(parity, np.conj(cc), cc)
        n1 = c1 * np.where(parity, np.conj(bb), bb)
        c0, c1, parity = n0, n1, parity ^ ee
    return c0[0] if scalar else c0


def stats(v):
    i, j = np.triu_indices(len(v), 1)
    r = np.abs(v[i] - v[j]) / (j - i) ** ALPHA
    return dict(D=float(r.max() / r.min()), rmin=float(r.min()), rmax=float(r.max()),
                argmin=[int(i[r.argmin()]), int(j[r.argmin()])],
                argmax=[int(i[r.argmax()]), int(j[r.argmax()])],
                minimum_step=float(np.min(np.abs(np.diff(v)))),
                distinct_vertices=len(np.unique(v)))


def optimise_one(label, initial, seed):
    P = Problem(len(initial), gauge='v2')
    # P.unpack anchors v[0]=0 and v[1]=1, not the final endpoint.
    initial = (initial - initial[0]) / (initial[1] - initial[0])
    x = P.pack(initial)
    assert stats(P.unpack(x))['distinct_vertices'] == len(initial)
    best = stats(initial)['D']
    best_shape = norm_shape(initial)
    best_x = x.copy()
    print(label, 'start', stats(norm_shape(initial)), flush=True)
    betas = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000,
             10000, 30000, 100000, 300000, 1000000]
    for beta in betas:
        for _ in range(3):
            res = minimize(P.fg, x, args=(beta,), jac=True, method='L-BFGS-B',
                           options=dict(maxiter=350, gtol=1e-15, ftol=1e-16,
                                        maxcor=30))
            x = res.x
            current_shape = norm_shape(P.unpack(x))
            current = stats(current_shape)['D']
            if np.isfinite(current) and current < best:
                best, best_shape, best_x = current, current_shape, x.copy()
        print(label, 'after beta', beta, 'raw D', current, 'best valid D', best,
              flush=True)
        if not np.isfinite(stats(current_shape)['D']):
            print(label, 'singular annealing output: revert to best-valid checkpoint',
                  flush=True)
            x = best_x.copy()
            break
    annealed = best_shape
    annealed_for_slp = (annealed - annealed[0]) / (annealed[1] - annealed[0])
    assert stats(annealed_for_slp)['distinct_vertices'] == len(annealed)
    polished, _ = slp(P, annealed_for_slp, rho0=0.005, iters=400)
    polished = norm_shape(polished)
    st = stress(P, polished, tol=1e-8)
    A, L = grad_rows(P, polished)
    residual = np.full(P.nfree, np.nan)
    if st['status'] == 0:
        ak = A[st['imax']].toarray()
        bk = A[st['imin']].toarray()
        residual = st['lam'] @ ak - st['mu'] @ bk
    record = dict(label=label, seed=seed, annealed=stats(annealed),
                  polished=stats(polished),
                  tight_max=len(st['pairs_max']), tight_min=len(st['pairs_min']),
                  stress_status=int(st['status']),
                  multiplier_residual_inf=float(np.max(np.abs(residual))),
                  pairs_max=[list(map(int, p)) for p in st['pairs_max']],
                  pairs_min=[list(map(int, p)) for p in st['pairs_min']])
    print(label, 'FINAL', json.dumps(record), flush=True)
    return norm_shape(initial), annealed, polished, st, record


def main():
    m = 100
    t = np.linspace(0.0, 1.0, m)
    v0 = gamma(t)
    assert len(np.unique(v0)) == m
    print('UNIFORM t start:', stats(v0), flush=True)
    rng = np.random.default_rng(271828)
    runs = [('unperturbed', v0, -1)]
    for seed in (41, 73):
        rng = np.random.default_rng(seed)
        vv = v0.copy()
        vv[1:-1] += 1e-5 * np.median(np.abs(np.diff(v0))) * (
            rng.normal(size=m - 2) + 1j * rng.normal(size=m - 2))
        runs.append((f'perturb_1e-5_seed_{seed}', vv, seed))
    arrays, records = {}, []
    for label, initial, seed in runs:
        initial, annealed, polished, st, record = optimise_one(label, initial, seed)
        arrays[label + '_initial'] = initial
        arrays[label + '_annealed'] = annealed
        arrays[label + '_polished'] = polished
        if st['status'] == 0:
            arrays[label + '_multipliers_max'] = st['lam']
            arrays[label + '_multipliers_min'] = st['mu']
        records.append(record)
        np.savez(ROOT / 'results/new_candidate_uniform_m100_restarts.npz',
                 parameters=t, **arrays)
        (ROOT / 'results/new_candidate_uniform_m100_restarts.json').write_text(
            json.dumps(records, indent=2) + '\n')
        time.sleep(0.05)
    base = ROOT / 'data/inputs/v2_100best.npy'
    if base.exists():
        vv = norm_shape(np.load(base))
        print('previously saved best D', stats(vv), flush=True)
        for rec in records:
            label = rec['label']
            print('normalized shape distance', label, 'to old best:',
                  supdist(arrays[label + '_polished'], vv), flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
