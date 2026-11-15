# Follow-up: what the exhaustive checks establish

## Bottom line

Approximate long-index **maximisers** of $r_{ij}$ are a useful **screening and
ranking device**, not a self-similarity criterion. On `s_m800_plain.npy`, they
often rank better than arbitrary intervals of the same index length. But the
best tested nontrivial restrictions are still far from being geometric copies
of the whole, and approximate minimisers perform particularly badly.

“Near” is relative to the fixed-index piecewise-linear shape-discrepancy metrics
used below. These experiments do not establish a continuum theorem, select a
unique limiting curve, or rule out another correspondence or configuration.

## Why endpoint extremality alone is not evidence

Suppose a restriction has index length $\tau=L/799$ and is similar to the whole
via a complex dilation $\lambda$. Both curves are endpoint-normalised before
comparing them. The **interior** points must still satisfy the similarity, not
just the endpoints. The chord ratio satisfies
$$r_{a,a+L}=|\lambda|\tau^{-\alpha}R_{\rm full},$$
where $R_{\rm full}$ is the full-arc endpoint ratio. Neither the dilation nor
this scale factor is fixed by geometric resemblance alone. Requiring scale
factor one means imposing an additional absolute-$r$-profile-preservation
condition; it is not part of the definition of a self-shortening.

An exact Koch quarter-copy already demonstrates this: it is a genuine
self-shortening but its ratio is $1$, while the full sharp Koch ratio range is
$[0.4834590784,1.0137194906]$. Thus an exact self-shortening need not be a
global maximiser or minimiser. A straight segment sampled on 800 nodes also
shows that $799=17\cdot47$ is only a node-commensurability obstruction, not an
obstruction to a geometric self-copy in its interpolant.

## What the exhaustive scans found

There are $205{,}120$ contiguous blocks with index gap $L\ge160$, including
the whole block. One scan used the endpoint-gauged maximum discrepancy at the
block's actual nodes, minimising over conjugation, reversal, and both. A second
computed the supremum between the two piecewise-linear interpolants exactly
from their union of index knots. That is an exact maximum of a piecewise-affine
complex difference, not a dense-grid approximation.

* Among all $124{,}920$ such blocks with $L/799\le1/2$, the smallest node-based
  discrepancy was $0.1021$, at $[259,511]$. That interval is **not** an
  extremal chord.
* With $L\le600$, the best exact-polyline discrepancy was $0.0866895$, at
  $[165,623]$, again a mid-$r$ interval rather than the literal maximum.
* Genuine-copy controls attain residuals of about $0.003$–$0.02$, depending on
  the construction and resolution. The gap between those controls and the
  measured candidates is therefore material, although their residuals are not
  universal resolution floors.
* Cross-resolution checks at the tested $1/2$ and $1/4$ index scales likewise
  found no close copy. Their best residuals are about $0.116$–$0.165$, roughly
  $9$–$37$ times the corresponding same-branch residuals. These are empirical
  thresholds for those scales and metrics, not a nonexistence theorem.

An independent literal knot-union calculation gave $E=0.3156823538$ for
$[455,689]$, $E=0.0866895052$ for $[165,623]$, and $E=0.1327098944$ for
$[25,739]$. The last value is the earlier symmetry-check interval: its point
symmetry does not make it an approximate self-copy.

## Is the ranking still useful?

At fixed index length, larger $r$ tends to have smaller shape error in this
file. For example, node-based Spearman correlations are approximately
$-0.57$, $-0.76$, $-0.68$, $-0.80$, $-0.90$, $-0.83$, and $-0.77$ at
$L=160,200,252,320,400,500,600$. This is a real ranking signal, not a
certificate. It weakens as the block approaches the whole ($+0.02$ at $L=700$),
so pooled correlations without length controls should not be used.

A corrected cross-resolution Procrustes computation (after initially
conjugating the fitted complex multiplier incorrectly) gave consistent rather
than inverted correlations: $-0.896$, $-0.894$, and $-0.622$ on its three
scans. An exact rotated-copy test of that formula returns zero residual.

The fair practical recommendation is therefore to use approximate long-index
maxima to **propose** endpoint pairs, compare them with arbitrary intervals of
the same index length, and test the entire restriction. Do not use near-minima
here as an analogous selection rule, and do not stop the search after selecting
the top chord. If nothing scores near a calibrated positive control, the
screen has not found a geometric shortening.

## Restricted distortion really is uninformative here

For every one of the $205{,}120$ blocks with $L\ge160$, exact pairwise chord
extrema were queried using per-row prefix maxima/minima and a one-dimensional
sparse range table. Literal random-submatrix checks passed for both extrema.

The full value is $D=1.5910519607023$. Across all tested blocks:
$$\max|D_{\rm block}-D|=4.55237604546\times10^{-6},\qquad
\max|D_{\rm block}/D-1|=2.86123656417\times10^{-6}.$$
Only $50{,}616$ blocks ($24.68\%$) match $D$ to relative tolerance $10^{-9}$.
Thus it is **approximate**, not an identity; in particular, none of the $640$
blocks of length $160$ matches that tightly. But $D$ is far too close to
constant to discriminate shape: even controls with clearly non-similar
macroscopic profiles have essentially the same value.

This strengthens, rather than replaces, the earlier concrete $[25,739]$ check:
that particular equality is valid, but it is not evidence of self-similar
geometry or an optimal shortening.

## Review and reproducibility

The four-member exploratory run ended with **two accepted solutions**:
[m1's exhaustive-polyline report](.ultramath/smwhr9283/m1/report.md) and
[m3's revised synthetic-control report](.ultramath/smwhr9283/m3/report.md).
The other directions stopped without an accepted final report. Their corrected
numerical findings were checked but are not treated as additional accepted
solutions. The introductory $2\%$-band experiment and its positive Koch
calibration remain in [the original note](extremal-shortenings.md).

Reproduction:

```bash
python3 scripts/test_extremal_shortenings.py
python3 scripts/calibrate_extremal_shortenings.py
python3 scripts/check_block_distortion.py
```

The last script writes [the exhaustive $D$ data](block-distortion-exhaustive.json).
The main accepted scans and their scripts are in `.ultramath/smwhr9283/m1/` and
`.ultramath/smwhr9283/m3/`. The fixed-index PL interpolants in these tests are
comparison devices, not continuum curves whose snowflake distortion is finite.
