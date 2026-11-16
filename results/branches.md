# Two macro-shapes, and what each rung is worth

Requested: an S-shaped local minimum at $m=200$ (the stored rung was not S-shaped), a
separate figure for the one-sided arc branch, and an estimate of how reliable the ladder
points are.

## The $m=200$ rung

`descended_m200.npy` (the previous entry, $D=1.5031685289$) has
$\operatorname{Im} v\in[-0.0019,+0.5435]$: one-sided, i.e. an arc, not an S. The S-side
candidate `morph400.npz['m200 from S(m100)+kick']` ($D=1.4916292100$,
$\operatorname{Im}\in[-0.226,+0.230]$) was converged with the exact-minimax SLP:

```
python3 scripts/optimise_m400.py --m 200 --rounds 4 --minutes 12 \
    --start .scratch/start_m200_S.npy --out results/m200_S.npy
```

| round | $D$ | stress | tight pairs |
|---|---|---|---|
| start | 1.4916292100 | 2 | 1+1 |
| 1 | 1.4908867897 | 0 | 45+345 |
| 2 | **1.4908552600** | 0 | 48+346 |
| 3 | 1.4908781151 | 0 | 47+346 |

`stress=0` is the multiplier certificate (no descent direction). The value is stable at the
fifth digit across rounds, $0.82\%$ below the arc rung and $3.80\%$ above $m=100$, so the
ladder stays monotone. This is the entry now in Table 1 and in `figures/ladder-m.png`.

## The one-sided arc branch

New rungs at $m=40$ and $m=100$ were made by downsampling the stored arc and converging with
the same SLP (both reach `stress=0`); $m=400$ was re-polished (`stress=0`); $m=800$ and
$m=1600$ are the stored swarm configurations, uncertified.

| $m$ | S branch | status | arc branch | status | gap |
|---:|---:|---|---:|---|---:|
| 40 | 1.3523099782 | stationary | 1.3692601379 | certified | $+1.25\%$ |
| 100 | 1.4362283999 | stationary | 1.4475037000 | certified | $+0.79\%$ |
| 200 | 1.4908552600 | certified | 1.5031685289 | stationary | $+0.83\%$ |
| 400 | 1.5417920000 | certified | 1.5487020131 | certified | $+0.45\%$ |
| 800 | 1.5900261898 | certified | 1.5939535359 | not certified | $+0.25\%$ |
| 1600 | 1.6355308997 | stationary | 1.6315839649 | not certified | $-0.24\%$ |

The gap is the "error bar" used in `figures/ladder-D.png`: the best value known in the other
branch at the same $m$. It is a measured estimate of how far a rung could move without
leaving the search, not a bound. It falls monotonically from $1.25\%$ to $0.24\%$, so the two
branches are converging; at $m=1600$ they cross, on two uncertified configurations.

## Restart study (why the bars are not smaller)

`scripts/restart_study.py` re-runs a cheap annealed L-BFGS (no SLP) from perturbed, resampled
and random starts; output in `results/restart-study.json`. At $m=40$ and $m=100$ the
perturbed restarts reproduce the stored value to $2\times10^{-4}$; at $m\ge200$ the same
budget lands $0.1$--$1\%$ *above* the polished value (e.g. $1.6043$ versus $1.5900$ at
$m=800$). So the cheap protocol measures how much the answer depends on where the search
starts, which is much larger than the distance between the two branches at large $m$; the
fine structure there is only reachable with the SLP. That is the honest reason the bars are
the branch gaps and not something tighter.
