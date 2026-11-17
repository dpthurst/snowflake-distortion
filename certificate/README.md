# Enclosures behind the Koch constants

The companion note `notes/koch-extremals.pdf` proves

    inf rho = 2*6^(-alpha)   and   sup rho = 2*sqrt(211)*69^(-alpha),   alpha = log3/log4,

with both constants sharp. Its §Certification cites the four documents here. Each one is an
unedited transcript of the working note that produced that part of the computation; they were
written in parallel by independent lines of attack, so their notation is their own
(`omega = e^{i pi/3}`, `lambda = 1+omega`, `D` a normalised numerator, `R` the full ratio).

| file | what it establishes |
|---|---|
| `lower-bound-enclosure.md` | the lower inequality: exact reduction to the corner square, the renormalisation that makes the pinch an exact self-copy of the whole problem, the equality families, and the branch-and-bound runs (`513` boxes at gap $\ge 1/2$ with no residual; `17372` boxes at gap $\ge 2^{-10}$ with relative slack $10^{-4}$) |
| `maximum-enclosure.md` | the upper inequality: `Phi(a,b) <= 2 sqrt(211) 69^(-alpha)` on the fundamental domain, by exact `Z[omega][1/3]` arithmetic plus directed-interval signs on the residual cells |
| `maximum-enclosure-output.txt` | the run that closes it: `109620` cells examined, `156` residuals, all `156` closed, least margin `2.53e-7` |
| `minimum-local-checks.md` | an independent pass at the minimum: a fully rational per-cell certificate on the four cells around the diagonal extremum, plus a signed-derivative check of the maximum on `2` of its `8` residual cells (a check on the argument above, not a substitute for it) |

The scripts named inside these documents (`code/...`, `max/...`) were
scratch files of the analysis and are not part of this repository; the numbers they produced are
quoted in the documents and reproduced in the output file. The certificate code that *is* in the
repository, `scripts/certify_mm.py`, is a separate implementation: it is what produces the certified
finite-`m` bounds of the main note, not the two constants here.
