# Lean fragment: the value at the minimising configuration

`KochSnowflake.lean` formalises the exact value

```lean
theorem KochSnowflake.rho_value_at_minimiser :
    KochSnowflake.kochRatio 5 11 2 = 2 * (6 : ℝ) ^ (-KochSnowflake.alpha)
```

Here `kochRatio 5 11 2` is the ratio at the code pair `5/16`, `11/16`, and
`alpha = log 3 / log 4`. In the notation of `notes/koch-extremals.tex`, this is

\[
  \rho(5/16,11/16)=2\cdot 6^{-\alpha}.
\]

Compile from this directory with the toolchain and Mathlib on `LEAN_PATH`:

```sh
lean KochSnowflake.lean
```

The file has no `sorry`, `admit`, or custom axioms. `#print axioms
KochSnowflake.rho_value_at_minimiser` reports only `propext`,
`Classical.choice`, and `Quot.sound`.

## Scope

This is deliberately only the endpoint computation, not the global statement

\[
  \forall x\ne y,\qquad \rho(x,y)\ge 2\cdot 6^{-\alpha}.
\]

The global lower bound uses the renormalisation
\(\Phi(a/4,b/4)=\Phi(a,b)\), a finite compact-domain branch-and-bound, exact
Eisenstein arithmetic, and rational enclosures of \(\alpha\)-powers. Formalising
that certificate requires an interval-arithmetic layer for `Real.rpow` and a
machine-checkable enumeration; it is not contained in this file. The Lean
statement should therefore not be read as a proof that `2 * 6^(-alpha)` is the
global infimum.

## Source correspondence

* `kochMap_matches_source` checks that the four maps of `kochMap` are exactly
  \(z/3,\ \omega z/3+1/3,\ \bar\omega z/3+(1+\omega)/3,\ z/3+2/3\).
* `kochPoint_five` and `kochPoint_eleven` compute the two code points as
  \((3+\omega)/9\) and \((5+\omega)/9\).
* `kochPoint_five_sub_eleven` gives the chord \(-2/9\).
* `rpow_four` and `rpow_two` record \(4^\alpha=3\) and \(2^\alpha=\sqrt3\);
  `rho_value_algebra` performs the scalar simplification.
