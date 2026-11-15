# Shortening similarity and symmetry: what was actually checked

We compared the raw saved node sets, not an interpolant fitted as a new curve.
The tested candidates were
`.ultramath/smwhdqmht/m1/s_m1600_plain.npy`,
`s_m1600_kick.npy`, `arc_m1600.npy`, and the earlier exact point-symmetric candidates
`.ultramath/smwhdqmht/m2/w_smoothS_t2_200.npy` and `w_smoothS_t2_400.npy`.

For each contiguous interval of relative length $1/2,1/4,1/8,1/16,1/32$, we
reparameterised it by its endpoint chord, centred it, and fitted a complex
translation-dilation. Conjugation and index reversal were tested as four variants
($0$ identity, $1$ conjugation, $2$ reversal, $3$ both). The comparison starts with
$64$ nodes per interval and is repeated with $512$ nodes (capped at the number of
original node records). Every error below is the **maximum pointwise residual
normalised by the shortening's own diameter**. All-pairs $D$ values were also
recomputed from the raw coordinates with ordinary index gaps $|i-j|^\alpha$,
$\alpha=\log3/\log4$.

## Best fitted shortenings

| candidate | interval(s), approximately | variant | max residual at 64 samples | max residual at 512 samples |
|---|---|---|---:|---:|
| S plain | $[0.624765,0.687305]$ | identity | 0.12638 | 0.13402 |
| one-sided arc | $[0.624765,0.749844]$ | identity | 0.15596 | 0.15834 |
| exact point-symmetric $m=400$ | $[0.343358,0.375940]$ and $[0.624060,0.656642]$ | conjugation + reversal | 0.09992 | 0.10294 |

These are not near-zero residuals and should not be read as proving
self-similarity. A coarser local structural analogy is possible, but there is no
support here for pruning to a rigid self-similar construction. Other interval
positions/depths were not exhaustively proved non-self-similar; these were the
best of all tested contiguous node blocks, not cherry-picked examples.

## Exact versus approximate point symmetry

Our reported inversion test is
$$e_{\rm inv}=\sup_i\left|\frac{w_i+w_{m-1-i}-1}{\operatorname{diam}(w)}\right|,$$
so zero means exact point symmetry under $z\mapsto1-z$ together with index reversal.
A different mirror operation, $w_i\mapsto1-\overline{w_{m-1-i}}$, was also measured
separately; it is not the same test and must not be confused with parameter reversal.

| configuration | inversion supremum error | positive peak | negative peak | peak/trough indices |
|---|---:|---:|---:|---:|
| S plain, $m=1600$ | $0.143105$ | $0.230180$ | $-0.247871$ | $0.26079$ / $0.68480$ |
| S kick, $m=1600$ | $0.143368$ | $0.230767$ | $-0.247984$ | $0.26079$ / $0.68480$ |
| arc, $m=1600$ | $0.985821$ | $0.657157$ | $0$ (no negative lobe) | $0.43965$ / — |
| exact symmetric, $m=200$ | $0$ | $\pm0.238845$ | $\pm0.238845$ | $0.27638$ / $0.72362$ |
| exact symmetric, $m=400$ | $1.1\times10^{-16}$ | $\pm0.248537$ | $\pm0.248537$ | $0.26817$ / $0.73183$ |

The exact symmetric $m=400$ configuration therefore already is an example of the
requested **perfect** pair of lobes: every corresponding vertex pair sums to $1$.
It is not just approximately symmetric in lobe height. Its exact all-pairs
$D=1.543697915377456$ is only about $0.09\%$ above the newer S-branch $m=400$
value $1.5422713964167$. At $m=200$, the exact symmetric value $1.497314357482232$
is about $0.38\%$ above the newer S-kick value $1.4916292099925$. These are feasible
achieved configurations, not certified stationary points or global optima.

## Relation to restriction optimality

The best similar-subarc direction is not the same as the restriction direction.
If $\gamma$ is globally optimal and denoted $\gamma_{a,b}$ after endpoint
normalisation, then every restriction has the same optimal distortion. There is no
additional implication that a subarc must be an approximate **geometric copy** of
the whole curve. Our data illustrate the distinction: some blocks have very similar
$D$ to the whole configuration but substantially different fitted shape. Evaluating
both a restriction's $D$ and its profile is therefore preferable to assuming that
restriction constancy obliges self-similarity.

## Reproducibility

- `scripts/analyze_shortenings.py`: all fits, deviations, and the 512-node checks.
- `results/shortening-similarity.json`: per-interval residuals and fitted coefficients.
- `results/shortening-similarity.out`: complete tables over all the tested intervals.

No new continuous curve is claimed to be certified by this finite-$m$ comparison.
The practical symmetric alternative is the exact point-symmetric $m=400$ seed above;
it is a candidate for a higher-$m$ continuation in the same invariant subspace, not
a substitute for a stationarity or hierarchy certificate.
