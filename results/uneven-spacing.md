# Uneven parameter spacings at finite $m$

Asked: do finite-$m$ optima with uneven spacings converge to the same shapes as the uniform
ladder?  Scripts: `scripts/optimise_uneven.py` (free and softmax spacings),
`scripts/uneven_spacing.py` (bounded gap ratio and fixed profiles), `scripts/plot_uneven.py`.
Data: `results/uneven-*.json`, `results/uneven_*.npz`.

## Free spacings are degenerate

With $0=x_0<\dots<x_{m-1}=1$ completely free the infimum of $D$ is $1$ for every $m\ge3$:

* three points already give $D=1$ exactly (sides $c,c,2^\alpha c$);
* splitting each of them into near-duplicates --- $\varepsilon$ apart in parameter, $c\,\varepsilon^\alpha$
  apart in space --- keeps every ratio at $1$.

The optimiser finds exactly this: at $m=6$, $D=1.0008$ with gaps spanning $5\times10^{-5}$ to
$0.56$ and the six points in three pairs. So the free problem approximates no curve and says nothing
about $D_\infty$; it cannot be used to study shapes.

## Bounded unevenness: well posed, and faster

Impose $\max\Delta x/\min\Delta x\le\rho$ ($\rho=1$ is the uniform ladder). Comparable gaps keep the
parameters asymptotically uniform, so $D_m(\rho)\to D_\infty$ for each fixed $\rho$, and every value
is still $\le D_\infty$.

```
          m=8      m=12     m=16     m=24     m=32
 rho=1   1.15661  1.22013  1.24535  1.29708  1.33192
 rho=2   1.15012  1.20140  1.23488  1.28427  1.31850
 rho=4   1.11664  1.18004  1.21941  1.27031  1.30524
 rho=8   1.09427  1.15733  1.18935  1.25117  1.28663
 rho=16  1.08978  1.13516  1.17524  1.23142  1.27044
 rho=32  1.06277  1.09517  1.16348  1.20199  1.26423
```

Unevenness buys a few per cent at fixed $m$, and it is worth having: the uniform ladder needs
$m=100$ to reach $1.436$, while $\rho=32$ gets $1.26$ by $m=32$.

## The shapes: the same S, reached sooner

* $m\ge12$: every optimum is the same S profile as the uniform ladder. Endpoint-gauged Hausdorff
  distance between the $\rho=1$ and $\rho=32$ curves: $0.171$ at $m=16$, $0.045$ at $m=32$. So the
  uneven problem converges to the same curve, just faster.
* $m=8$: with $\rho\ge4$ the optimum leaves the S branch and becomes one-sided --- the $m=4$ pattern
  with points inserted; at $\rho=32$ it reproduces $D=1.0628$, the certified $D_4$.
* Where the points go is not a smooth grading. At $\rho=32$ the gaps sit at the bound in clusters:
  at $m=32$ they alternate between $0.0026$ and $0.082$ --- near-duplicate groups separated by jumps,
  which is what a self-similar structure wants.

## Smooth profiles are much worse

Fixing $x_i=(i/(m-1))^\gamma$ and optimising only the positions:

```
          m=12     m=16     m=24     m=32
 gamma=1  1.22013  1.25629  1.29768  1.33242   (uniform)
 gamma=0.5  5.13709  6.01134  7.19555  8.59356
 gamma=2   11.16436 14.41858 20.42381 25.99210
```

So it is not the size of the unevenness that matters but its structure: a smooth profile destroys the
matching of chords, while optimised spacings place near-duplicate points where the curve needs them.

## Caveats

These are local optima from a handful of starts. At fixed $\rho$ the values are not monotone in $m$
($\rho=32$: $1.0628$ at $m=8$, $1.0952$ at $m=12$), which cannot happen for exact optima and is the
signature of under-convergence at $m\ge12$. Seeding $m=12$ and $m=16$ from the $m=8$ optimum improves
them slightly ($1.1352\to1.1223$, $1.1635\to1.1516$) but not to $1.0628$.
