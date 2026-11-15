# Converging a configuration, and what it does to the overlaps

The question this was meant to settle: when the overlap search finds poorer
matches at larger $m$, is that because the finer configurations are less
converged, or because finer detail is genuinely different?

## The run

`scripts/optimise_m400.py --m 400` (the script takes `--m`), starting from `s_m400.npy`
($D=1.5422713964$), five rounds of $\beta$-annealed L-BFGS followed by the
exact-minimax SLP (200 iterations, $\rho_0=0.01$). Three minutes in total.

```
start:  D=1.5422713964  stress=2   (1+1 tight pairs)
round 0: D=1.5437942979 stress=2   (1+1)      <- L-BFGS only
round 0: D=1.5418525759 stress=0   (80+698)   <- after SLP
round 1: D=1.5418180845 stress=0   (73+707)
round 2: D=1.5418046959 stress=0   (71+707)
round 3: D=1.5417920000 stress=0   (70+706)
round 4: D=1.5418644386 stress=0   (73+700)
```

`stress` is the multiplier certificate: it solves for weights on the tight max
and min pairs, so $0$ means a first-order KKT point of the minimax. The stored
point has one tight max pair and one tight min pair and **no** certificate — the
singleton-pair contradiction recorded earlier. The converged point has ~70 tight
max pairs and ~700 tight min pairs and a certificate.

So: the stored $m=400$ configuration really was not stationary, by the only
test we have.

## What changed, and what did not

```
D             1.5422713964  ->  1.5417920000     (-0.031%)
shape         0.0077 of the diameter apart (size ratio 1.002, rotation -0.1 deg)
```

Converging moved the curve by less than one per cent of its diameter. The
landscape is very flat here: a $0.031\%$ gain in $D$ for a $0.77\%$ change of
shape.

The overlap search, run on both with the same relative grid:

```
           best    #<0.12   #<0.15   #<0.20
stored     0.0481    176      787     4578
converged  0.0566    185      799     4591
```

Essentially unchanged — the counts move by a few per cent, the best residual
shifts within the scatter of the grid.

## The same at $m=800$

The run takes 19 minutes at $m=800$ (1.2 s per SLP iteration):

```
start:  D=1.5910519607  stress=2   (1+1 tight pairs)
round 0: D=1.5924480930 stress=2   (1+1)        <- L-BFGS only
round 0: D=1.5900725043 stress=0   (143+1403)
round 1: D=1.5900396023 stress=0   (148+1409)
round 2: D=1.5900382880 stress=0   (148+1410)
round 3: D=1.5900336157 stress=0   (142+1407)
round 4: D=1.5900261898 stress=0   (147+1408)
```

$D$ improves by $0.065\%$, the shape moves by $0.0107$ of the diameter (size
ratio 1.004, rotation $-0.4^\circ$), and the overlaps are again unchanged:

```
           best    #<0.12   #<0.15   #<0.20
stored     0.0785     70      418     3386
converged  0.0800     69      412     3426
```

Comparing the two *certified* configurations at the same relative grid, the
thinning is still there: $787 \to 412$ below $0.15$ from $m=400$ to $m=800$.

## What this says

At both $m=400$ and $m=800$ the deficiency in optimisation is **not** what the
overlaps were measuring. In each case the stored point was technically
non-stationary — one tight max and one tight min pair, no multipliers — but
already within about $1\%$ of a certified stationary point, and the overlaps
survive that move intact. So the thinning of overlaps as $m$ grows is not
explained by the optimiser having stopped early.

That leaves the resolution explanation: finer configurations carry structure at
scales the coarse ones do not have, and the coarse-shape comparison is degraded
by it. The earlier observation that the *coarse* residual also grows with $m$
suggests something beyond sampling, but it is not optimisation.

## Caveats

* One starting point per $m$. The converged points are certified KKT points, not
  global optima, and each sits in the basin of the stored configuration it came
  from.
* The certificate is first-order only: no second-order or global claim.
* The overlap comparison is on a coarse grid; small changes are below its
  resolution.
* The certificate is first-order; the stored points may still differ from the
  converged ones in ways the overlap grid cannot see.

## Reproduce

```bash
python3 scripts/optimise_m400.py --m 400 --rounds 5 --minutes 20
python3 scripts/optimise_m400.py --m 800 --rounds 5 --minutes 40
```
Writes `results/m{400,800}_converged.npy` and `.json`.
