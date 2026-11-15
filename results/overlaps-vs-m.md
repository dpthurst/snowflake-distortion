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

## Two confounds in the ladder, found afterwards

**The ladder is not one family.** Checking the macro-morphology (whether the
curve crosses its own endpoint chord):

```
curve              min Im w   max Im w   morphology
m=100 final         -0.2049    +0.2338    S, crosses its chord
m=200 descended     -0.0019    +0.5435    arc, one side
m=400 s_m400        -0.2444    +0.2268    S
m=800 s_m800_plain  -0.2470    +0.2281    S
m=1600 s_m1600_plain -0.2479   +0.2302    S
```

The stored $m=200$ curve is an *arc*; the rest of the S-branch rungs were seeded
from the S shape at $m=100$ and carried upward by upsampling. So the $m=200$ row
of the tables above compares different morphologies, not different resolutions.
The trend from $m=400$ to $1600$ is within one branch and stands; the $m=200$
entry should not be read as the bottom of it.

**The chain is inherited by construction.** `morph.py` builds the S shape at
$m=400$ from the $m=100$ S shape; `m800.py` starts from that and upsamples to
$m=800$; `m1600.py` starts from `s_m800_plain` and upsamples to $m=1600$, each
followed by a local descent. The persistence of an alignment from $m=400$ to
$m=1600$ is therefore partly inherited rather than independently found at each
rung.

**Where the optima actually sit.** At $m=400$ four stored configurations give
$D=1.5423,\ 1.5488,\ 1.5491,\ 1.5526$ — a spread of $0.67\%$; at $m=800$ three
give $1.5911,\ 1.5940,\ 1.5946$; at $m=1600$ two give $1.6355,\ 1.6316$. And the
KKT check recorded in `limiting_curve_directions.md` gives relative residuals
$1.07$ and $1.04$ for the $m=400$ and $m=800$ points, i.e. they are not
stationary. These configurations are accordingly not samples of "the optimum" at
their resolution, and their fine structure in particular is whatever the
upsampling and local descent produced.

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

## One alignment followed up the ladder

`figures/overlap-ladder-1.png` and `-2.png`: the same pair of intervals, with
indices rescaled by $(m-1)/799$, at $m=200,400,800,1600$; the block in blue and
the arc normalised and turned back by the fitted rotation.

```
block [120,209] vs arc [24,100] (m=800 units)     block [533,594] vs arc [648,678]
   m     resid(401)  resid(24 samples)                resid(401)  resid(24)
  200      0.242        0.223                          0.310       0.267
  400      0.109        0.101                          0.086       0.067
  800      0.102        0.078                          0.092       0.088
 1600      0.119        0.106                          0.113       0.104
```

Two readings:

* At every rung the residual falls when only the coarse shape is sampled — fine
  detail costs something.
* But the coarse residual also *grows* from $m=400$ to $m=1600$
  ($0.067\to0.088\to0.104$ for the second alignment). So "more detail is
  visible" is not the whole story: the coarse correspondence itself degrades.
* At $m=200$ the alignment is simply absent, which is explained by the
  morphology confound above rather than by resolution.
