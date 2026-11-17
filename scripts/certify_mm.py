"""Interval branch-and-bound for the m-point global distortion, for m = 4 and 5.

Generalises scripts/certify_m4.py: the free coordinates are the m-2 interior points (2(m-2) real
numbers, so 4 dimensions at m=4 and 6 at m=5).  A box is discarded when a rigorous lower bound for
D on it exceeds the target; the a priori box comes from D <= Dbar with Dbar an achieved value.

    python3 scripts/certify_mm.py --m 5 --eps 1e-2 --Dbar 1.093091082
"""
import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 50
ALPHA = mp.log(3) / mp.log(4)


def outward(x, down=True, n=4):
    y = float(x)
    for _ in range(n):
        y = np.nextafter(y, -np.inf if down else np.inf)
    return y


def gap_intervals(m):
    """|i-j|^alpha as outward-rounded intervals, one per gap 1..m-1"""
    out = [None]
    for g in range(1, m):
        v = mp.power(g, ALPHA)
        out.append((outward(v, True, 6), outward(v, False, 6)))
    return out


def inflate(lo, hi, n=3):
    for _ in range(n):
        lo = np.nextafter(lo, -np.inf)
        hi = np.nextafter(hi, np.inf)
    return lo, hi


def abs_iv(relo, rehi, imlo, imhi):
    def sq(a, b):
        lo = np.where(a > 0, a * a, np.where(b < 0, b * b, 0.0))
        hi = np.maximum(a * a, b * b)
        return inflate(lo, hi, 1)
    s1, S1 = sq(relo, rehi)
    s2, S2 = sq(imlo, imhi)
    L, H = inflate(s1 + s2, S1 + S2, 1)
    return (np.nextafter(np.sqrt(np.maximum(L, 0.0)), -np.inf),
            np.nextafter(np.sqrt(np.maximum(H, 0.0)), np.inf))


def ratios(boxes, m, G):
    """lower/upper bounds of every r_ij for each box (boxes: (N, 2(m-2), 2))"""
    nb = len(boxes)
    lo_list, hi_list = [], []
    # positions: point 0 is 0, point m-1 is 1, the rest are free
    for i in range(m):
        for j in range(i + 1, m):
            if i == 0 and j == m - 1:
                lo = np.full(nb, 1.0 / G[m - 1][1])
                hi = np.full(nb, 1.0 / G[m - 1][0])
                lo_list.append(lo); hi_list.append(hi); continue
            ri = None if i == 0 else (i - 1)
            rj = None if j == m - 1 else (j - 1)
            if ri is None:       # point i is the origin
                relo = -boxes[:, 2 * rj, 1]
                rehi = -boxes[:, 2 * rj, 0]
                imlo = -boxes[:, 2 * rj + 1, 1]
                imhi = -boxes[:, 2 * rj + 1, 0]
            elif rj is None:     # point j is 1
                relo = boxes[:, 2 * ri, 0] - 1.0
                rehi = boxes[:, 2 * ri, 1] - 1.0
                imlo = boxes[:, 2 * ri + 1, 0]
                imhi = boxes[:, 2 * ri + 1, 1]
            else:
                relo = boxes[:, 2 * ri, 0] - boxes[:, 2 * rj, 1]
                rehi = boxes[:, 2 * ri, 1] - boxes[:, 2 * rj, 0]
                imlo = boxes[:, 2 * ri + 1, 0] - boxes[:, 2 * rj + 1, 1]
                imhi = boxes[:, 2 * ri + 1, 1] - boxes[:, 2 * rj + 1, 0]
                relo, rehi = inflate(relo, rehi, 1)
                imlo, imhi = inflate(imlo, imhi, 1)
            lo, hi = abs_iv(np.asarray(relo), np.asarray(rehi), np.asarray(imlo), np.asarray(imhi))
            g = j - i
            lo_list.append(lo / G[g][1]); hi_list.append(hi / G[g][0])
    lo = np.array(lo_list); hi = np.array(hi_list)
    return inflate(lo, hi, 2)


def D_lower(boxes, m, G):
    lo, hi = ratios(boxes, m, G)
    best = np.zeros(lo.shape[1])
    for i in range(len(lo)):
        for j in range(len(lo)):
            if i != j:
                best = np.maximum(best, lo[i] / hi[j])
    return best


def initial_box(m, Dbar, G):
    nb = 2 * (m - 2)
    lo = np.full(nb, -np.inf); hi = np.full(nb, np.inf)
    # |v_i - v_j| <= Dbar * g_ij / 3   (since min r >= r(0,m-1)/Dbar = 1/(3 Dbar))
    # r(0,m-1) = 1/g_{m-1} is the reference pair, so every competitor has
    # |v_i - v_j| <= Dbar * g_ij / g_{m-1}
    ref = G[m - 1][0]
    for i in range(1, m - 1):
        lim = Dbar * G[i][1] / ref
        lo[2 * (i - 1)], hi[2 * (i - 1)] = -lim, lim
        lo[2 * (i - 1) + 1], hi[2 * (i - 1) + 1] = -lim, lim
        lim_r = Dbar * G[m - 1 - i][1] / ref     # distance to the endpoint 1
        lo[2 * (i - 1)] = max(lo[2 * (i - 1)], 1.0 - lim_r)
        hi[2 * (i - 1)] = min(hi[2 * (i - 1)], 1.0 + lim_r)
    assert np.all(lo < hi), (lo, hi)
    return np.stack([lo, hi], axis=1)[None, :, :]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--eps", type=float, default=1e-2)
    ap.add_argument("--Dbar", type=float, default=None)
    ap.add_argument("--maxboxes", type=int, default=3_000_000)
    args = ap.parse_args()
    m = args.m
    G = gap_intervals(m)
    if args.Dbar is None:
        v = np.load(ROOT / f"data/m{m}_best.npy")
        i, j = np.triu_indices(m, 1)
        r = np.abs(v[i] - v[j]) / np.abs(i - j) ** float(ALPHA)
        args.Dbar = float(r.max() / r.min())
    Dbar = args.Dbar
    L = Dbar - args.eps
    boxes = initial_box(m, Dbar, G)
    print(f"m={m}, Dbar={Dbar:.9f}, target L={L:.9f}, initial box "
          f"{np.round(boxes[0], 3).tolist()}", flush=True)
    t0 = time.time(); total = 1; level = 0
    while True:
        lb = D_lower(boxes, m, G)
        alive = lb < L
        na = int(alive.sum())
        if level % 2 == 0 or na == 0:
            print(f"  level {level:2d}: {len(boxes):9d} boxes, {na:8d} below target "
                  f"[{time.time()-t0:.1f}s]", flush=True)
        if na == 0:
            print(f"CERTIFIED: D_{m} >= {L:.9f}")
            break
        w = (boxes[alive][:, :, 1] - boxes[alive][:, :, 0]).max()
        if w < 1e-12:
            print("refinement limit; certified only", float(lb[alive].min())); break
        b = boxes[alive]
        width = b[:, :, 1] - b[:, :, 0]
        k = np.argmax(width, axis=1)
        mid = 0.5 * (b[np.arange(len(b)), k, 0] + b[np.arange(len(b)), k, 1])
        left = b.copy(); right = b.copy()
        left[np.arange(len(b)), k, 1] = mid
        right[np.arange(len(b)), k, 0] = mid
        boxes = np.concatenate([left, right])
        total += len(boxes); level += 1
        if total > args.maxboxes:
            print("box budget exceeded; certified only", float(lb[alive].min())); break


if __name__ == "__main__":
    main()
