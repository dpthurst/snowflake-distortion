# $D_4$ exactly, and the first lower bound for $D_\infty$

## The value

After the similarity gauge $v_0=0$, $v_3=1$ the six ratios are
$r_{01}=|v_1|$, $r_{12}=|v_1-v_2|$, $r_{23}=|v_2-1|$, $r_{02}=|v_2|/\sqrt3$, $r_{13}=|v_1-1|/\sqrt3$,
$r_{03}=3^{-\alpha}$, because $2^\alpha=\sqrt3$ and $4^\alpha=3$. With

$$c=3^{-\alpha}=0.418687843545147,\qquad a=\frac{1-c}{2}=0.290656078227427,\qquad b=\sqrt{c^2-a^2}=0.301361169565585,$$

the configuration $v_1=a+bi$, $v_2=1-a+bi$ has

$$r_{01}=r_{12}=r_{23}=r_{03}=c,\qquad r_{02}=r_{13}=\frac{\sqrt{c(1+c)}}{\sqrt3},$$

so it attains

$$D_4\le K:=\sqrt{\frac{1+c}{3c}}=1.06276592307635403\ldots$$

(the identity $(1-a)^2+b^2=(1+c)^2/4+c^2-(1-c)^2/4=c(1+c)$ is exact, so all six ratios are exact
in $\mathbb{Q}(c)$ apart from the square roots).

## The certificate

`scripts/certify_m4.py` proves $D_4\ge K-1.7\times10^{-12}$ by interval branch-and-bound:

* **a priori box.** If $D\le K$ then every ratio lies in $[c/K,Kc]$; hence
  $|v_1|\le Kc$, $|v_2|\le\sqrt3Kc$, $|v_1-1|\le\sqrt3Kc$, $|v_2-1|\le Kc$, which confines
  $(v_1,v_2)$ to $[0.2293,0.4450]\times[-0.4450,0.4450]\times[0.5550,0.7707]\times[-0.4450,0.4450]$.
  Outside the box $D>K$, so a bound inside it is a bound everywhere.
* **interval bound.** For a box, each $r_{ij}$ is enclosed by $[\underline r_{ij},\overline r_{ij}]$
  (rectangle arithmetic, outward rounded by a few ulps at every operation), and
  $D(v)\ge\max_{i\ne j}\underline r_i/\overline r_j$ for every $v$ in the box, since
  $D\ge r_i/r_j$ pairwise.
* **refinement.** Boxes whose lower bound exceeds the target are discarded, the rest are bisected
  along the widest coordinate.

Run at target $K-10^{-12}$: 54 163 boxes, depth 156, 0.07 s, certified lower end
$1.0627659230747$. So

$$1.0627659230747\ \le\ D_4\ \le\ 1.0627659230764,$$

a window of $1.7\times10^{-12}$, and therefore, since $D_m\le D_\infty$,

$$\boxed{D_\infty\ge 1.0627659230747.}$$

This is the first explicit lower bound in the workspace. The upper bound is still Koch's
$2.096805160875$, so $1.0628\le D_\infty\le2.0968$; the ladder suggests the truth is near $1.64$
or below, but nothing proves it.

## Checks

* The closed-form configuration reproduces $K$ to $10^{-42}$ in 40-digit arithmetic, and its
  minimum ratios are $c=0.418687843545146687$, its maximum $\sqrt{c(1+c)}/\sqrt3=0.444967172526105914$.
* The interval bound was tested against 16 000 sample points in 400 random boxes (scales
  $10^0$ to $10^{-8}$): the sampled $D$ never fell below the box bound.
* 200 000 random configurations in $[-1,2]^2$ for $(v_1,v_2)$: none with $D\le K$ outside the
  a priori box.
* 300 random-start minimisations all land at $K$ or above.

## What it does not give

It sits 0.57 below the $m=1600$ value, which is exactly the gap the note is about; but it is
a *lower* bound, which is what the ladder cannot give. The "floors" direction stops here:
$m=4$ is four-dimensional after the gauge, $m=6$ is eight-dimensional, and the same box
count would grow by roughly $(\text{volume ratio})^2$ --- out of reach.
