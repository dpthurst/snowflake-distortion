# Direct symmetry check of one shortened S configuration

## Scope

This is a small reproducible computation on the saved **finite configuration**
`.ultramath/smwhdqmht/m1/s_m800_plain.npy`, which contains $800$ complex nodes.
It is not a search result: only the full configuration and the single interval
$25\le i\le739$ are evaluated. In particular, this note does **not** establish
that the interval is optimal, that its two lobes are exactly symmetric, or that
a particular continuous limiting curve exists or is unique.

The earlier broad exploratory sweep is intentionally not retained in the checking
script. All quantities below are recomputed with direct array operations on the
saved coordinates, without its sparse range-maximum routine.

## Definitions and interval convention

For an arc with node coordinates $v_i$, use the endpoint gauge
$$w_i=\frac{v_i-v_{\rm first}}{v_{\rm last}-v_{\rm first}},$$
so the endpoints are $0$ and $1$. For $N$ nodes let
$$e_i^{\rm pt}=|w_i+w_{N-1-i}-1|,\qquad
 e_i^{\rm mir}=|w_i-(1-\overline{w_{N-1-i}})|.$$
The reported supremum and RMS errors are respectively
$\max_i e_i/\operatorname{diam}(w)$ and
$\sqrt{N^{-1}\sum_i(e_i)^2}/\operatorname{diam}(w)$. The first tests invariance
under point inversion $z\mapsto1-z$ with index reversal; the second is a distinct
mirror test and is not a test of similarity between the two lobes.

With $\alpha=\log3/\log4$, the discrete distortion is
$$D=\frac{\max_{i<j}|w_i-w_j|/(j-i)^\alpha}
{\min_{i<j}|w_i-w_j|/(j-i)^\alpha}.$$
Index gaps $j-i$ retain their spacing on the original $800$-node grid;
endpoint normalisation does not change $D$.

The selected block includes its two endpoints, so it has $715$ nodes. Its
parameter interval is $[25/799,739/799]\approx[0.0312891,0.9249061]$, and its
relative parameter length is $714/799\approx0.893617$.

## Numerical result

| quantity | full $800$ nodes | shortened block |
|---|---:|---:|
| point-inversion supremum | $0.1431998031$ | $0.0379849072$ |
| point-inversion RMS | $0.0663846085$ | $0.0237535320$ |
| mirror supremum | $0.4473882296$ | $0.4299759323$ |
| mirror RMS | $0.3017448096$ | $0.2680665197$ |
| discrete $D$ | $1.591051960702300$ | $1.591051960702299$ |

The block's diameter in its own endpoint gauge is $1.0006878424$; the full
configuration's is $1.0078468725$. The largest point-inversion residual in the
block occurs at original node $267$. Thus point symmetry improves substantially
under this particular shortening, but a $3.80\%$-of-diameter worst residual is
not exact point symmetry or a perfect equal-lobe construction.

The restricted $D$ agrees with the full value to the displayed floating-point
precision for a concrete reason: both extremal pairs remain in the block. The
reported minimising pair is $(509,599)$ and the reported maximising pair is
$(455,689)$, in original node indices. This finite-$D$ preservation is not
evidence of an optimal restriction, self-similarity, or a distinguished
continuum limit.

## Reproduction

Run `python3 scripts/check_shortening_symmetry.py`. It writes the same values to
`results/shortening_symmetry.json`; the captured run is
`results/shortening_symmetry.out`.
