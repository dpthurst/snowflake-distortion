# Partial self-overlaps across the ladder

The search of `results/overlaps.md` repeated at several $m$, plus two cross-$m$
tests. All numbers are the scale-free shape residual defined there (both
intervals centred and normalised by their own diameter, rotation fitted, affine
parameter); a residual of $0.10$ is a $10\%$ mismatch of the diameter.

## Within one curve, same relative grid

Lengths and start steps are held fixed as fractions of $m$, so every curve is
searched with the same relative configurations (about $10^6$ each). Controls
have the curve's own step lengths and turning angles in shuffled order.

```
    m     curve    best  #<0.12  #<0.15  #<0.20
  200    stored  0.0559     142     664    3442
  200  shuffled  0.1324       0       3      50
  400    stored  0.0481     176     787    4578
  400  shuffled  0.1389       0       6     127
  800    stored  0.0785      70     418    3386
  800  shuffled  0.1418       0       1      41
 1600    stored  0.0862      35     271    2574
 1600  shuffled  0.1639       0       0      11
```

## Same node count

The coarse stage samples $24$ points whatever the interval length, so a finer
curve is relatively under-sampled. To remove that, every curve is decimated to
$200$ nodes and searched on one identical grid:

```
    m    best  #<0.12  #<0.15  #<0.20
  200  0.0709      97     550    3008
  400  0.0741      52     298    1899
  800  0.0772      30     212    1534
 1600  0.0677      18     151    1293
```

(Decimating a fine optimum is not the same object as a genuine optimum at the
coarse resolution — its chords cut the fine wiggles — so this is a control for
the sampling, not a substitute for the ladder.)

## Do overlaps found at one $m$ survive at the next?

Eight distinct overlaps at $m=800$, mapped to $m=1600$ by the index ratio
$1599/799$ and rescored there:

```
resid at m=800   resid at m=1600   on shuffled m=1600   percentile among 300 random placements
0.0926           0.1252            0.911                0.0%
0.1053           0.1380            0.814                0.0%
0.1220           0.1463            0.549                0.0%
0.1148           0.1261            0.733                0.0%
0.0962           0.1171            0.823                0.0%
0.1187           0.1576            0.800                0.0%
0.1048           0.1384            0.911                0.0%
0.0888           0.1142            0.600                0.0%
```

The alignments are inherited: they sit at the same place at doubled indices and
are far better there than a shuffled curve, though the residual grows by about
a third.

## Cross-$m$

Block from one curve, arc from another:

```
pair        best (away from   #<0.15 away    corresponding    shuffled arcs
            the correspondence)   from it     piece scores        best
800 -> 1600      0.0879          326           0.0166            0.183
400 -> 800       0.0635          616           0.0296            0.144
800 -> 400       0.0794          572           0.0300            0.128
```

Two things are visible. The two samplings of the same piece of curve agree to
$0.017$–$0.030$ — that is the floor of this metric across $m$, not zero, so
values near $0.03$ are not evidence of anything. And away from the
correspondence the curves still match at $0.06$–$0.09$, an order of magnitude
better than against shuffled arcs: the self-overlap structure is *shared*
between neighbouring rungs of the ladder.

## What this says about a limiting object

* Self-overlaps are real at every $m$: hundreds below $0.15$, against $0$–$6$
  for the shuffled controls.
* Their **number falls** as $m$ grows — $787 \to 418 \to 271$ below $0.15$ on
  the relative grid, and $550 \to 298 \to 212 \to 151$ at fixed node count —
  and the **best residual does not improve** ($0.048$ at $m=400$, $0.086$ at
  $m=1600$).
* The overlaps that exist are **inherited**: the same alignments reappear at
  doubled indices and are shared across neighbouring curves.

So the coarse shapes of the finite optima resemble themselves, and the
resemblance is stable enough to survive refinement, but it *thins out* rather
than sharpening. A family converging to a self-similar limit curve would be
expected to hold or accumulate such coincidences as detail is added; here they
dissolve. That is evidence against a strongly self-similar $m\to\infty$ object
at these relative scales, and equally against reading these coincidences as
structure of a limiting curve rather than of the smoother coarse members of the
ladder.

The residual level is also worth keeping in view: $0.05$–$0.09$ is one to two
orders of magnitude above the $0.003$–$0.02$ at which exact copies sit in this
project. None of these is a near-copy.

## Reproduce

```bash
python3 scripts/overlaps_vs_m.py        # writes results/overlaps-vs-m.json and .out
```
