"""Certify a lower bound for D_4, hence for D_infinity.

The m=4 problem is four real parameters after the similarity gauge v0=0, v3=1:
v1, v2 in C.  Its six ratios are

    r01=|v1|,  r12=|v1-v2|,  r23=|v2-1|,  r02=|v2|/sqrt3,  r13=|v1-1|/sqrt3,  r03=3^-alpha,

because |i-j|^alpha equals 1, sqrt3, 3^alpha for gaps 1, 2, 3 (2^alpha = sqrt3 and 4^alpha = 3).

The configuration  v1 = a+bi,  v2 = 1-a+bi,  a=(1-c)/2,  b=sqrt(c^2-a^2),  c=3^-alpha,
has r01=r12=r23=r03=c and r02=r13=sqrt(c(1+c))/sqrt3, so it attains
K = sqrt((1+c)/(3c)) = 1.0627659230763...

This script proves D_4 > K - eps for a chosen eps, by interval branch-and-bound over the
box that can contain any configuration with D <= K (outside it D > K by the a priori bound
below).  All arithmetic on boxes is outward rounded: every float operation is followed by
an inflation of a few ulps, so the computed D-lower-bound is a true lower bound.

    python3 scripts/certify_m4.py [--eps 1e-4] [--out results/m4-global.json]
"""

import argparse
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 50
ALPHA_MP = mp.log(3) / mp.log(4)
C_MP = mp.power(3, -ALPHA_MP)                      # 3^-alpha
K_MP = mp.sqrt((1 + C_MP) / (3 * C_MP))            # attained value


def outward(x, down=True, ulps=4):
    y = float(x)
    for _ in range(ulps):
        y = np.nextafter(y, -np.inf if down else np.inf)
    return y


# interval endpoints of the two transcendental constants, outward rounded
C_IN = (outward(C_MP, True, 8), outward(C_MP, False, 8))          # 3^-alpha
S3_IN = (outward(mp.sqrt(3), True, 4), outward(mp.sqrt(3), False, 4))
INV_C_IN = (outward(1 / C_MP, True, 8), outward(1 / C_MP, False, 8))   # 3^alpha


def inflate(lo, hi, ulps=4):
    return (np.nextafter(lo, -np.inf), np.nextafter(hi, np.inf)) if ulps == 1 else \
           (np.nextafter(np.nextafter(np.nextafter(lo, -np.inf), -np.inf), -np.inf),
            np.nextafter(np.nextafter(np.nextafter(hi, np.inf), np.inf), np.inf))


def sq_interval(lo, hi):
    """[lo,hi]^2 as an interval"""
    lo2 = np.where(lo > 0, lo * lo, np.where(hi < 0, hi * hi, 0.0))
    hi2 = np.maximum(lo * lo, hi * hi)
    return inflate(lo2, hi2)


def abs_interval(relo, rehi, imlo, imhi):
    """|x| for x in a rectangle, lower and upper bounds (outward)"""
    slo, shi = sq_interval(relo, rehi)
    tlo, thi = sq_interval(imlo, imhi)
    lo_in, hi_in = inflate(slo + tlo, shi + thi)
    return np.nextafter(np.sqrt(np.maximum(lo_in, 0.0)), -np.inf), \
        np.nextafter(np.sqrt(np.maximum(hi_in, 0.0)), np.inf)


def ratios(boxes):
    """boxes: (N,4,2); returns (6,N) lower and upper bounds of r01,r12,r23,r02,r13,r03"""
    x1l, x1h = boxes[:, 0, 0], boxes[:, 0, 1]
    y1l, y1h = boxes[:, 1, 0], boxes[:, 1, 1]
    x2l, x2h = boxes[:, 2, 0], boxes[:, 2, 1]
    y2l, y2h = boxes[:, 3, 0], boxes[:, 3, 1]
    ones = np.ones_like(x1l)
    out_lo, out_hi = [], []
    # (0,1): re = x1, im = y1
    lo, hi = abs_interval(x1l, x1h, y1l, y1h)
    out_lo.append(lo); out_hi.append(hi)
    # (1,2): re = x1 - x2, im = y1 - y2
    a, b = inflate(x1l - x2h, x1h - x2l)
    c, d = inflate(y1l - y2h, y1h - y2l)
    lo, hi = abs_interval(a, b, c, d)
    out_lo.append(lo); out_hi.append(hi)
    # (2,3): re = x2 - 1, im = y2
    a, b = inflate(x2l - ones, x2h - ones)
    lo, hi = abs_interval(a, b, y2l, y2h)
    out_lo.append(lo); out_hi.append(hi)
    # (0,2): re = x2, im = y2, divided by sqrt3
    lo, hi = abs_interval(x2l, x2h, y2l, y2h)
    out_lo.append(lo / S3_IN[1])
    out_hi.append(hi / S3_IN[0])
    # (1,3): re = x1 - 1, im = y1, divided by sqrt3
    a, b = inflate(x1l - ones, x1h - ones)
    lo, hi = abs_interval(a, b, y1l, y1h)
    out_lo.append(lo / S3_IN[1])
    out_hi.append(hi / S3_IN[0])
    # (0,3): 3^-alpha
    out_lo.append(np.full_like(x1l, C_IN[0]))
    out_hi.append(np.full_like(x1l, C_IN[1]))
    lo = np.array(out_lo); hi = np.array(out_hi)
    return inflate(lo, hi)


def D_lower(boxes):
    """lower bound for max/min over each box, as max over ordered pairs lo_i / hi_j"""
    lo, hi = ratios(boxes)
    n = lo.shape[1]
    best = np.zeros(n)
    for i in range(6):
        for j in range(6):
            if i != j:
                est = lo[i] / hi[j]
                best = np.maximum(best, est)
    return best


def initial_box():
    """any configuration with D <= K has r's in [c/K, K c]"""
    cK_lo = float(C_IN[0]) / float(K_MP)
    cK_hi = float(C_IN[1]) * float(K_MP)
    s3cK_hi = float(S3_IN[1]) * cK_hi
    x1 = (max(-cK_hi, 1.0 - s3cK_hi), cK_hi)
    y1 = (-cK_hi, cK_hi)
    x2 = (1.0 - cK_hi, s3cK_hi)
    y2 = (-cK_hi, cK_hi)
    return np.array([[x1, y1, x2, y2]], dtype=float)      # shape (1,4,2)


def subdivide(boxes):
    """bisect every box along its widest coordinate"""
    width = boxes[:, :, 1] - boxes[:, :, 0]
    k = np.argmax(width, axis=1)
    lo = boxes[:, :, 0].copy(); hi = boxes[:, :, 1].copy()
    mid = 0.5 * (lo[np.arange(len(boxes)), k] + hi[np.arange(len(boxes)), k])
    lo2 = lo.copy(); hi2 = hi.copy()
    hi[np.arange(len(boxes)), k] = mid
    lo2[np.arange(len(boxes)), k] = mid
    left = np.stack([lo, hi], axis=2)
    right = np.stack([lo2, hi2], axis=2)
    return np.concatenate([left, right], axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--maxboxes", type=int, default=4_000_000)
    ap.add_argument("--out", default=str(ROOT / "results/m4-global.json"))
    args = ap.parse_args()

    K = float(K_MP)
    L = K - args.eps
    t0 = time.time()
    boxes = initial_box()
    print(f"K = {K:.15f}, target L = {L:.15f}, initial box {boxes[0].tolist()}", flush=True)
    total = 1
    level = 0
    while True:
        lb = D_lower(boxes)
        alive = lb < L
        n_alive = int(alive.sum())
        print(f"  level {level:2d}: {len(boxes):9d} boxes, {n_alive:8d} below L, "
              f"max lower bound {lb.max():.12f}  [{time.time()-t0:.1f}s]", flush=True)
        if n_alive == 0:
            print(f"\nCERTIFIED: D_4 >= {L:.15f}  (K = {K:.15f})", flush=True)
            break
        width = (boxes[alive][:, :, 1] - boxes[alive][:, :, 0]).max()
        if width < 1e-12:
            print("refinement limit reached; certified only D_4 >= "
                  f"{lb[alive].min():.15f}", flush=True)
            L = float(lb[alive].min())
            break
        boxes = subdivide(boxes[alive])
        total += len(boxes)
        level += 1
        if total > args.maxboxes:
            print("box budget exceeded; certified only D_4 >= "
                  f"{lb[alive].min():.15f}", flush=True)
            L = float(lb[alive].min())
            break
    summary = dict(eps=args.eps, K=K, L_certified=L, boxes=total, levels=level,
                   seconds=time.time() - t0, c=float(C_MP))
    Path(args.out).write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
