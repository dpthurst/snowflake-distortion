# The $\alpha=\tfrac12$ endpoint is excluded: primary-source check

The attached paper is exactly the relevant one:

A. S. Besicovitch and I. J. Schoenberg,
**On Jordan arcs and Lipschitz classes of functions defined on them**,
*Acta Mathematica* **106** (1961), 113–136. It is kept as `refs/besicovitch-schoenberg-1961.pdf`
in this repository (the user's attachment).

## Why it settles the endpoint

Their Theorem 6 (p. 117) reads:

> If $J$ is a plane Jordan arc and the function $G(v)$ is uniformly of the Lipschitz class
> $x^2$ along $J$, then $G(v)$ is necessarily a constant.

Suppose $\gamma$ were a finite-distortion embedding of the $1/2$-snowflaked interval $[0,1]$ in
$\mathbb R^2$, i.e. there exist $0<c\le C<\infty$ with
$$c|x-y|^{1/2}\le|\gamma(x)-\gamma(y)|\le C|x-y|^{1/2}.$$
Its image $J=\gamma([0,1])$ is a Jordan arc in the plane. Write $G=\gamma^{-1}$ on $J$.
For $v=\gamma(x)$ and $v'=\gamma(y)$, the lower inequality gives
$$|G(v)-G(v')|=|x-y|\le c^{-2}|\gamma(x)-\gamma(y)|^2=c^{-2}|v-v'|^2.$$
Hence $G$ has Lipschitz class $x^2$ uniformly along $J$, so Theorem 6 forces it to be constant.
But a constant inverse implies that $J$ has only one distinct point, contrary to $x\ne y$ in $[0,1]$.

Therefore **there is no finite-distortion embedding of $([0,1],|x-y|^{1/2})$ in $\mathbb R^2$**.
Equivalently, writing $D_\infty$ for the infimum of distortion over injective arcs of $[0,1]$, we have
$D_\infty(\alpha=1/2)=+\infty$. This also rules out a bounded infinite version obtained by repeatedly
sampling/rescaling compact subintervals: such an arc, restricted to a compact interval, would possess
the inverse used above.

So the earlier statement in this project that $\alpha=1/2$ was open was wrong. The ladder at that
exponent should be read as evidence diverging to an infinite continuum value. This does not say an
injective Osgood arc of Hausdorff dimension 2 does not exist (those exist); it says there is no
bi-Hölder $1/2$ parametrisation of $[0,1]$ by one.

The whole proof is the inverse-class observation above,应 along the Jordan arc rather than trying
to apply the forward Holder inequalities. The paper itself also shows the contrast at $\alpha<1$:
it constructs nonconstant Lipschitz classes along suitably chosen Jordan arcs in the plane, but
planes cannot support an inverse of class $x^2$. That is why the endpoint is really singular.

The change is recorded in `notes/improving-koch-extremals.tex` and `README.md`, with the proof
sketched above and the primary journal/page reference. These files and the reference PDF are copied
locally only; they were not pushed.
