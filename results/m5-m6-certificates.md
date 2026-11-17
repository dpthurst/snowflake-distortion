# Certified windows for $D_5$ and $D_6$, and the improved lower bound

The $m=4$ certificate of `results/m4-global.md` generalises: `scripts/certify_mm.py --m 5 --eps ...`
does the same interval branch-and-bound in $2(m-2)$ dimensions (six at $m=5$, eight at $m=6$). The a
priori box comes from $D\le\bar D$: every ratio lies in $[r_{\min},\bar D r_{\min}]$ with
$r_{\min}=1/|i-j|^{\alpha}$ at the reference pair $(0,m-1)$ after the gauge $v_0=0$, $v_{m-1}=1$.
Boxes are discarded when $\max_{i\ne j}\underline r_i/\overline r_j$ exceeds the target.

| $m$ | certified window | width | boxes | time |
|---|---|---|---|---|
| 4 | $[1.0627659230747,\ 1.0627659230764]$ | $1.7\times10^{-12}$ | 54 163 | 0.07 s |
| 5 | $[1.093090076,\ 1.093091076]$ | $10^{-6}$ | 80 | 0.7 s |
| 6 | $[1.127430051,\ 1.128130050]$ | $7.0\times10^{-4}$ | $1.15\times10^{7}$ | 11 s |
| $\ge7$ | fails: $1.5\times10^6$ boxes at $10^{-2}$, bound still below 1.01 | | | |

Upper ends are achieved configurations: at $m=5$ the optimum found here
($D=1.0930910764117$, better than the $1.093091572$ of the earlier census) and at $m=6$ one with
$D=1.1281300496046$. Hence

$$D_\infty \ge 1.127430051,$$

up from $1.0627659230$, and from $1.119534712$ as this table stood before the $m=6$ search.

## The $m=6$ rung, re-done

The earlier $m=6$ row was $[1.119534712,\ 1.129534712]$ at width $10^{-2}$: the upper end was the
best configuration then in hand ($1.129534712$, while a 40-start census had already reported
$1.128130873$) and the lower end was whatever the branch-and-bound could reach with that $\bar D$.

`scripts/optimise_m6.py` re-searches the eight-dimensional problem: $20\,000$ random starts,
annealed L-BFGS on the smooth surrogate of $\log D$, then three rounds of sequential LP on the exact
range. Every polished start lands on the same configuration,
$$D_6 = 1.1281300496046,$$
with the next distinct basin $0.55\%$ away at $1.134347$. So the value is almost certainly the
optimum as well as an achieved value, though nothing here proves that --- and the certificate does
not need it: it needs only *an* upper bound.

The configuration is rigid, like the $m=5$ one. Writing $c$ for the common minimum ratio,
$$|v_i-v_{i+1}| = c\ \ (i=0,\dots,4), \qquad |v_5-v_0| = 5^{\alpha}c,$$
so the five steps and the end-to-end chord are tight minima (six pairs), and the tight maxima are
the five pairs at spans two and three, $(0,3),(1,3),(1,4),(2,5),(3,5)$. Eleven tight pairs against
eight free coordinates: a degenerate KKT point, with `slp2.stress` returning status $0$ and an
explicit multiplier vector on those pairs, so no descent direction leaves it.

Lowering $\bar D$ from $1.129534712$ to $1.1281300506$ shrinks the a priori box, and the same code
then closes the window to $7\times10^{-4}$ instead of $10^{-2}$:

```
python3 scripts/certify_mm.py --m 6 --eps 7e-4 --Dbar 1.1281300506 --maxboxes 12000000
CERTIFIED: D_6 >= 1.127430051   (11497959 boxes, 11.3s)
```

The bound is first order in the box width, so this is near the end of what the method gives at
$m=6$: $\epsilon=6\times10^{-4}$ already exhausts $1.2\times10^7$ boxes, and even at
$\epsilon=7\times10^{-4}$ the default budget of $3\times10^6$ is not enough --- the run above needs
`--maxboxes 12000000`, which is why the command carries it. Controls run alongside:
at $m=5$ the known setting ($\epsilon=10^{-6}$, $\bar D=1.093091076$) still certifies, while a
target that is false ($D_6\ge1.12815$ with $\bar D=1.1282$) is not certified. The full log is
`results/m6-certificate.log`, the search record `results/m6-search.json`.

## Why $m=5$ is cheaper than $m=6$

The $m=5$ optimum is rigid: four consecutive chords exactly $1/3$ (the same value as the end-to-end
chord, since $4^\alpha=3$) and the three chord-of-two steps equal, so the a priori box is already a
small neighbourhood and eighty boxes close the window to $10^{-6}$. At $m=6$, even with the rigid
configuration above, the first-order bound resolves only to $7\times10^{-4}$: it is linear in the box
width, so each extra decimal digit in $\bar D-\epsilon$ costs a factor of roughly a thousand boxes.
At $m=7$ (ten coordinates) it fails at $10^{-2}$.

Validation: the interval bound was tested against $25$--$30$ sample points in $300$ random boxes per
dimension (worst margin $10^{-6}$), the a priori box was checked to contain the known optimum, and
at a degenerate box centred on the optimum the bound returns $D$ to machine precision.

## What would improve it

A second-order (centred/mean-value) bound would cut the box count by orders of magnitude and might
bring $m=7,8$ into range; that is the natural next step if the lower bound is to be pushed toward
the ladder at all.
