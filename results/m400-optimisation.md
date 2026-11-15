# Converging the $m=400$ configuration, and what it did to the overlaps

The question this was meant to settle: when the overlap search finds poorer
matches at larger $m$, is that because the finer configurations are less
converged, or because finer detail is genuinely different?

## The run

`scripts/optimise_m400.py`, starting from `s_m400.npy`
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

## What this says

At $m=400$, the deficiency in optimisation is **not** what the overlaps were
measuring. The stored point was technically non-stationary but already within
$0.8\%$ of a certified stationary point, and the overlaps survive that move
intact. So the thinning of overlaps as $m$ grows is not explained by the
optimiser having stopped early — at this $m$ at least, and the $m=400$ point is
the one with the *largest* documented defect.

That leaves the resolution explanation: finer configurations carry structure at
scales the coarse ones do not have, and the coarse-shape comparison is degraded
by it. The earlier observation that the *coarse* residual also grows with $m$
suggests something beyond sampling, but it is not optimisation.

## Caveats

* One $m$, one starting point. The converged point is a certified KKT point, not
  a global optimum, and it sits in the basin of `s_m400.npy`.
* The certificate is first-order only: no second-order or global claim.
* The overlap comparison is on a coarse grid; small changes are below its
  resolution.
* $m=800$ would test the same question closer to where the loss appears, at
  roughly five times the cost (about 15–20 min).

## Reproduce

```bash
python3 scripts/optimise_m400.py --rounds 5 --minutes 20
```
Writes `results/m400_converged.npy` and `results/m400_converged.json`.
