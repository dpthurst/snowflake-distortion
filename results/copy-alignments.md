# Copy mode: a small-$m$ curve fitted into a piece of a larger one

Exploratory, and the caveats at the end matter as much as the numbers.
Computed with `scripts/align_core.py` (`fit_copy_multires`); the same
computation is in the page's copy mode, checked against these functions.

## What is fitted

For the block $w$ of the source curve and the whole $W$,
$$w(t)\approx c\,W(\sigma_0+\sigma_1 t)+d,$$
with complex $c,d$ and real $\sigma_0,\sigma_1$ keeping
$[\sigma_0,\sigma_0+\sigma_1]\subset[0,1]$. That is a similarity image of the
whole placed over the block, with its parameter allowed to move affinely. Six
real parameters, found by a coarse grid over the parameter range followed by
two Nelder–Mead refinements (101 then 301 samples), translation always re-fitted
in closed form.

Residuals are supremum norms in units of the block's endpoint chord, sampled on
1201 parameters.

## Checked against known answers

- An exact similarity image of a window of the $m=800$ curve (constructed so
  that its nodes sit on whole nodes) is recovered with
  $\sup=2.0\times10^{-10}$ and $\sigma$ to four decimals.
- The same construction with every other whole node dropped is chord-limited:
  the fit returns the chord floor, $\sup=5\times10^{-3}$.
- Breaking the affine parameter in the fit loses the exact recovery, so the
  test can fail.

## Measured, $m=800$ as the whole

The block is the small curve in its entirety; $\sigma_1(m-1)$ is how many nodes
of the $m=800$ curve the fitted copy covers.

```
source    sup     index-sup   |c|     s0      s1      covered nodes
m=40     0.0744   0.308       5.18    0.306   0.120    96
m=100    0.0633   0.243       0.99    0.000   0.984   786
m=200    0.0993   0.328       5.49    0.315   0.160   128
```

"index-sup" is the same fit restricted to the index correspondence
($\sigma$ fixed to the block's own range); the affine parameter helps by a
factor of three to four.

Shape-free control: random-turn curves with the same step-length multiset as
the source, twelve draws each. Sup medians 0.49, 0.40, 0.49 (minima 0.13, 0.13,
0.22). So the measured curves sit below every draw for $m=40$ and $m=200$, and
below the median for $m=100$ but above its minimum.

## Caveats

- Six free parameters against 40–200 nodes is a lot of freedom; the control is
  a weak one and is not a calibrated null for smooth curves.
- The $m=40$ and $m=200$ fits use a 5× blowup of a short arc (96–128 nodes), and
  a short arc of a large-$m$ curve is comparatively smooth. Some of the
  agreement may be smoothness rather than shape.
- The $m=100$ fit is a different regime: $|c|\approx1$ and $\sigma$ almost all
  of $[0,1]$, i.e. the two curves agree at the same scale over the whole range,
  which is what convergence of the ladder would produce.
- Both curves are optima of the same functional, so resemblance is expected.
  Nothing here says which mechanism is at work, and no continuum statement
  follows.

## Reproduce

```bash
python3 scripts/preview_align_copy.py    # the m=40 picture
python3 scripts/test_align_app.py        # fit checked against numpy
```

## The fitted scale is the reciprocal of the chord ratio

Both curves are gauged to unit endpoint chord, so for a block that is a piece
of the whole the exact fit has $|c| = 1/\rho$ with
$\rho=\lvert v_b-v_a\rvert/\lvert v_{m-1}-v_0\rvert$ the chord ratio of the
piece. Verified on seven intervals of the $m=800$ curve (residuals at
$10^{-12}$):

```
interval     tau      rho     1/rho    |c| fit    rho/tau^alpha
[0,799]     1.000    1.0000   1.0000   1.0000    1.0000
[0,399]     0.499    0.4777   2.0932   2.0932    0.8283
[0,199]     0.249    0.3525   2.8368   2.8368    1.0607
[259,511]   0.315    0.3722   2.6865   2.6865    0.9289
[455,689]   0.293    0.4284   2.3341   2.3341    1.1338
```

So a piece of the same curve always fits with $|c|>1$: in this gauge the
"shrunk copy" appears as $\rho<1$, and the content of the scale is $\rho$
against $\tau^\alpha$, not $|c|$ itself. The page now prints $\rho$,
$\tau^\alpha$ and their ratio in copy mode.

Across curves $|c|<1$ does occur, but only marginally: the $m=400$, $m=800$ and
$m=1600$ curves as the whole against the $m=100$ curve as the block give
$|c|=0.995$, $0.995$, $0.994$.

Full-curve blocks, fitted $|c|$ (whole in rows, block source in columns):

```
          m40     m100     m200     m400     m800    m1600
m40      1.000    1.959    1.859    1.942    1.939    1.937
m100    14.428    1.000    5.135    1.026    1.025    1.025
m200     4.474    2.022    1.000    5.349    5.344    5.345
m400    13.073    0.995    5.414    1.000    1.000    1.000
m800     5.179    0.995    5.485    1.000    1.000    1.000
m1600    2.469    0.994    5.461    1.000    1.000    1.000
```
