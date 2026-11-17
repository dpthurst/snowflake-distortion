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
| 6 | $[1.119534712,\ 1.129534712]$ | $10^{-2}$ | 268 | 1.8 s |
| $\ge7$ | fails: $1.5\times10^6$ boxes at $10^{-2}$, bound still below 1.01 | | | |

Upper ends are achieved configurations: at $m=5$ the optimum found here
($D=1.0930910764117$, better than the $1.093091572$ of the earlier census) and at $m=6$ one with
$D=1.129534712$ (the census reports $1.128130873$, $0.12\%$ better, which would tighten the upper
end but not the lower). Hence

$$D_\infty \ge 1.119534712,$$

up from $1.0627659230$.

## Why $m=5$ is cheap and $m=6$ is not

The $m=5$ optimum is rigid: four consecutive chords exactly $1/3$ (the same value as the end-to-end
chord, since $4^\alpha=3$) and the three chord-of-two steps equal, so the a priori box is already a
small neighbourhood and eighty boxes close the window to $10^{-6}$. At $m=6$ the box is still
tractable but the first-order interval bound cannot resolve below $10^{-2}$ without exhausting the
budget; at $m=7$ (ten coordinates) it fails at $10^{-2}$.

Validation: the interval bound was tested against $25$--$30$ sample points in $300$ random boxes per
dimension (worst margin $10^{-6}$), the a priori box was checked to contain the known optimum, and
at a degenerate box centred on the optimum the bound returns $D$ to machine precision.

## What would improve it

A second-order (centred/mean-value) bound would cut the box count by orders of magnitude and might
bring $m=7,8$ into range; that is the natural next step if the lower bound is to be pushed toward
the ladder at all.
