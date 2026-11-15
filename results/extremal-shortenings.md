# Can long extremal chords suggest self-shortenings?

**Qualified yes, as a cheap candidate generator; no, as a criterion.** On the
saved `s_m800_plain` configuration, approximate long-index **maximisers** of
$$r_{ij}=\frac{|v_i-v_j|}{|i-j|^\alpha},\qquad \alpha=\frac{\log3}{\log4}$$
supply some of the better candidate shortenings. Approximate **minimisers**
perform badly. The best match found is still not a close geometrical copy.
This is a test on one finite configuration, not a continuum self-similarity
theorem or a global rule for all configurations.

## What extremality determines, and what it does not

Put $m=799$ and normalise $V_i=(v_i-v_0)/(v_{799}-v_0)$. A similarity taking the
whole arc to the restriction $[a,b]$, if one exists, is forced by the endpoints:
its linear complex multiplier must be
$$\lambda=V_b-V_a.$$
This gives both a proposed scale and rotation. It imposes no condition on the
interior. A real self-copy must further satisfy
$$\gamma\!\left(\frac{(1-t)a+tb}{m}\right)-\gamma\!\left(\frac am\right)
 =\lambda\,\gamma(t)$$
for every $t$ (up to parameter reversal or conjugation), not just at $t=0,1$.
Thus one extremal scalar $r_{ab}$ cannot certify such a relation.

I discard short **index gaps**, not short Euclidean chord lengths. A chord of
short Euclidean length spanning a long index interval measures folding and is
not automatically discretisation noise.

## Test

For each candidate $[a,b]$, endpoint-normalise the block. Compare it to the full
saved 800-node configuration at matching parameter fractions, allowing identity,
conjugation, index reversal, and both together. Since both endpoint chords are
normalised to $1$, there is no additional unconstrained fitted rotation or
scale. The reported residual is the maximum pointwise disagreement divided by
the **block's own diameter**.

The full reference is sampled by linear-in-index interpolation for this shape
comparison only. This is not offered as a finite-distortion continuous curve.
For a control, the same comparison returns supremum error
$1.11\times10^{-16}$ for an exact Koch quarter-copy.

I used both $b-a\ge0.2(799)$ and $b-a\ge0.5(799)$. Besides the most extreme
chords, I tested all chords in the near-extremal bands
$$r_{ab}\ge r_{\max}/1.02\quad\text{and}\quad r_{ab}\le1.02r_{\min},$$
respectively. The table gives the best supremum match within each band. For
each selected block, every interval of **exactly the same index gap** was also
tested with the same score and all four orientation variants.

| index-gap cutoff | kind | near-extremal chords | best interval | variant | sup residual | RMS residual | sup rank among intervals of that length | best same-length residual |
|---|---|---:|---|---|---:|---:|---:|---:|
| $0.2m$ | max | 3587 | $(277,479)$ | reversal | $0.141146$ | $0.075870$ | 1 / 598 | $0.141146$ |
| $0.2m$ | min | 102 | $(346,720)$ | reversal + conjugation | $0.539112$ | $0.350895$ | 298 / 426 | $0.164211$ |
| $0.5m$ | max | 878 | $(112,669)$ | identity | $0.144740$ | $0.076216$ | 18 / 243 | $0.132376$ |
| $0.5m$ | min | 129 | $(6,418)$ | conjugation | $0.572507$ | $0.362555$ | 285 / 388 | $0.131670$ |

The single largest-ratio chord in the $0.2m$ band is $(455,689)$. Its residual
is $0.314118$, rank $111$ out of $566$ same-length intervals. The largest-ratio
chord with gap at least $0.5m$ is $(100,666)$, with residual $0.162350$, rank
$60$ out of $234$. Looking at a small near-maximising band therefore helps more
than blindly taking the top chord, but none of these establishes a close copy.

The selected blocks have discrete distortion very close to that of the whole
configuration, $D\approx1.59105196$; even the poorly matching near-minimal
blocks do. Consequently $D$ preservation must not be used as a shortcut for the
shape test. The previous point-symmetry test on $[25,739]$ answers yet another
question and does not imply self-shortening either.

## A necessity counterexample

Even exact self-shortening does **not** imply extremality. The standard Koch arc
satisfies $F(t/4)=F(t)/3$, so its first quarter is an exact self-copy. Its
parameter-normalised endpoint ratio is
$$\frac{|F(1/4)-F(0)|}{(1/4)^\alpha}=1.$$
But the sharp global ratio range is
$$[2\,6^{-\alpha},\,2\sqrt{211}\,69^{-\alpha}]
 \approx[0.4834590784,\;1.0137194906],$$
so this self-shortening is neither a maximiser nor a minimiser. Hence a search
only among extremal chords can also miss genuine examples.

## Practical use

A reasonable cheap procedure is:

1. truncate by index gap, not Euclidean chord length;
2. collect approximate long-index **maximisers**, not just the top maximiser;
3. use their endpoints to propose a similarity scale and rotation;
4. test the entire normalized restriction against the whole arc, including
   reversal/conjugation, and compare against arbitrary intervals of the same
   length;
5. treat a small residual as evidence worth repeating across $m$, not as a
   continuum construction.

The numerical evidence supports that heuristic weakly for maxima, and not for
minima on this configuration. It does not support asserting a rigid
self-similar construction. A broader cross-$m$ test would be needed before
calling the weak positive signal structural.

## Reproduction and checks

- `scripts/test_extremal_shortenings.py`: all rankings, fits, and same-length
  comparisons; writes the full JSON and a short summary JSON.
- `results/extremal_shortenings.json` and
  `results/extremal-shortenings-summary.json`.
- `scripts/calibrate_extremal_shortenings.py`: exact Koch quarter-copy control;
  output in `results/extremal-shortenings-calibration.out`.
- The best $[277,479]$ match was independently recomputed without importing the
  checker; all four variant errors agreed to floating-point precision.

## Follow-up

Independent exhaustive and cross-resolution checks refine this conclusion.
Approximate long-index maximisers are a useful empirical ranking screen, but
no tested nontrivial block comes close to the genuine-copy controls. Extremality
is neither necessary nor sufficient. All $205{,}120$ blocks of index gap at least
$160$ have restricted $D$ within relative error $2.8613\times10^{-6}$ of the
whole $D$, so this preservation is qualitatively vacuous for distinguishing
shape. See [the verified follow-up](extremal-shortenings-followup.md), including
the exact scope, corrected correlations and reproduction commands.
