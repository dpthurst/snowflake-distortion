# Where the distortion lives on $I\times I$

For a configuration with $m$ nodes, $x_i=i/(m-1)$ and chord 1, plot
$$r(x_i,x_j)=\frac{|v_i-v_j|}{|x_i-x_j|^{\alpha}},\qquad \alpha=\frac{\log3}{\log4},$$
whose max/min ratio is the project's discrete distortion $D$. Two
configurations, both at chord 1 so the fields are comparable:

| | nodes | $r$ range | $D$ |
|---|---|---|---|
| Koch arc, level 5 | 1025 | $[0.483459,\ 1.013719]$ | $2.096805161$ |
| $m=800$ S-branch | 800 | $[0.712611,\ 1.133801]$ | $1.591051961$ |

Figure: `figures/distortion-heatmap.png`, with the pairs within $1\%$ of each
extreme marked underneath.

Colour convention: a shared linear scale whose limits are **exactly** the global
minimum and maximum over both panels, $[0.483459,\ 1.133801]$, so no value is
clipped (the script asserts it). Pale is the minimum, dark red the maximum; the
black contour is the level set $r=1$, the $\alpha$-power law itself. The
diagonal $x=y$, where $r$ is undefined, is masked and drawn as a thin line.

An earlier rendering of the same field used `magma` on a log scale, whose low
end is black: the pixels for $r\approx0.48$–$0.55$ then came out black and read
as missing data. They were the minima. On that scale those values occupy only
the first $15\%$ of the map, which is why the low band looked like holes.

## Why the node level, not the continuum

Along the Koch arc the ratio is **not** bounded below: inside a straight piece
$F$ is affine, so $|F(x)-F(y)|/|x-y|^\alpha \asymp |x-y|^{1-\alpha}\to0$, and
the continuum infimum is $0$. The bounded object, and the one whose ratio is the
$D$ of the finite configurations, is the node-level matrix at a fixed level;
$D$ there is $2.0896$, $2.0934$, $2.0968$ at levels $3,4,5$, stabilising at the
known $2.0968$.

## What the Koch heatmap shows, and why it is the reference

The four sub-arcs are similar copies of the whole with ratio $1/3=4^{-\alpha}$,
so the field is exactly invariant under the 4-adic rescaling:
$$r(4i,4j)=r(i,j)\quad\text{to }1.4\times10^{-14}.$$
That is the block structure visible in the picture — the pattern in the four
sub-squares repeats the whole field, three levels deep. Nothing is
concentrated at any particular scale, which is exactly what a self-similar
candidate should look like.

## What the $m=800$ heatmap shows instead

The field is flatter ($D$ smaller by a quarter) and its extremes sit at
different places than the Koch's:

```
                      near-min pairs        near-max pairs
Koch, level 5         168  (336 symmetric)  988  (1976 symmetric)
   median index gap    6                    10
   gaps <= 5           0 of 168             340 of 988
   largest gap         385                  640

m=800                 1866 (3732)           1704 (3408)
   median index gap    3                    214
   gaps <= 5           1169 of 1866          24 of 1704
   largest gap         269                  573
```

Read off from that:

* In the Koch configuration the small-$r$ pairs are spread across scales —
  none at gap $\le5$, the minimiser at gap $12$, the near-min pairs reaching
  gap $385$. Exact self-similarity keeps all scales at the same ratio.
* In the $m=800$ configuration the near-minimum is dominated by the shortest
  gaps: $1169$ of $1866$ near-min pairs have index gap $\le5$. That is the
  discretisation floor, not structure, and it is the finite-$m$ analogue of
  the Koch's straight-piece degeneration.
* The near-maximum is the opposite: median index gap $214$, only $24$ pairs at
  gap $\le5$. The distortion is carried by long-range pairs, in a coarse block
  pattern that echoes the Koch's without being scale-invariant.

Restricting to long chords changes $D$ only slightly, because the minimiser is
already long-range: at gap $\ge 0.2n$, $D=2.0896$ for Koch (the level-4 value,
i.e. the same structure one level coarser) and $D=1.5893$ for $m=800$ against
$1.5911$ unrestricted.

## Reproduce

```bash
python3 scripts/distortion_heatmap.py
```
