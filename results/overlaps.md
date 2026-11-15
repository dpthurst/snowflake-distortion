# Partial self-overlaps of the $m=1600$ curve

Prompted by one alignment found by eye: block $[898,1215]$ against the arc
$[1068.1,1455.1]$, scale $3.26$, rotation $-24.6°$. This is a systematic search
for such pairs, with a control.

## What is being asked

For index intervals $I=[a,a+L]$ and $J=[c,c+K]$, is
$$V(c+Kt)\ \approx\ s\,V(a+Lt)+d,\qquad t\in[0,1],$$
with $s,d\in\mathbb{C}$? That is exactly the page's copy-mode model: a
similarity carrying the block onto the arc, with the parameter taken affinely.

The score is **scale free**: both sampled curves are centred and divided by
their own diameter, and only a rotation is fitted. A residual of $0.10$ means
"the shapes agree to $10\%$ of the diameter"; the size ratio is reported
separately. (Scoring the raw least-squares fit instead lets a small arc match
almost anything, because its errors are never amplified — a first version of
this search did that and produced spurious $3\%$ hits.)

Four parameters $(a,L,c,K)$ make an exhaustive search hopeless: a coarse grid
over $L\in\{120,\dots,640\}$, $K/L\in\{0.5,\dots,2.2\}$, block start and arc
start (steps of $12$ nodes, $24$ sample points) is followed by a local refit of
all four at $201$ points. Configurations whose index map is close to the
identity on the overlap are discarded.

## Control

The point of the control is that a smooth curve has more self-similarity than a
random one, so the question is whether the *order* of the steps matters. Two
controls, both with the $m=1600$ curve's own step lengths:

* **shuffled turns** — the same steps and the same multiset of turning angles,
  in a random order: same step-length and turn statistics, no sequence;
* **random turns** — i.i.d. uniform turning angles.

```
curve                best resid   #<0.12   #<0.15   (all length ratios)
as stored                 0.0981      34      261
shuffled turns 1          0.1639       0       0
shuffled turns 2          0.1448       0       1
shuffled turns 3          0.1688       0       0
random turns 1                -       0       0
random turns 2            0.1810       0       0
random turns 3            0.1834       0       0
```

Restricting to intervals of comparable length ($0.85\le K/L\le1.2$), where the
two pieces carry the same fine structure:

```
curve                best resid   #<0.15   #<0.20
as stored                 0.1071      85     1427
shuffled turns 1          0.1292      10       54
shuffled turns 2          0.1260       7       87
shuffled turns 3          0.1164       9       56
```

So the stored configuration has roughly ten times as many comparable-length
partial overlaps below $0.15$ as curves with the same step lengths and the same
distribution of turns. The overlaps are a property of the *order* of the steps,
not of smoothness statistics.

## The best few

```
resid   block            arc              size   angle
0.103   [ 947,1070]      [ 768, 824]      0.53   -66°
0.104   [1394,1577]      [  12, 138]      0.75    +3°
0.105   [ 240, 420]      [  48, 201]      0.90   +26°
0.108   [ 947,1127]      [ 768, 859]      0.59   -77°
0.110   [ 876,1204]      [ 732, 892]      0.55   -80°
0.111   [ 238, 490]      [  48, 254]      0.88   +31°
```

Figure: `figures/overlaps-m1600.png` (shapes normalised, fitted rotation
applied). Several of these are the same coincidence found at slightly different
bounds: the region around indices $730$–$1200$ is where most of them live.

Note the column `size`: the strongest hits mostly have $K/L\approx0.5$ with
$|s|\approx 0.5$—$0.6$, i.e. shorter arcs matching longer blocks at
approximately the $\alpha$-power size ratio ($0.5^{\alpha}=0.58$). Those are
matches of *coarse* shape: with half the index length the arc has fewer fine
wiggles, so only the large-scale shape can agree.

## How far each alignment reaches

Extending both intervals by $\delta$ nodes (keeping the affine map) and
re-measuring:

```
block            arc              δ=0   δ=10  δ=20  δ=40  δ=60  δ=80 δ=120
[ 947,1070]  [ 768, 824]        0.104  0.220 0.229 0.269 0.230 0.174 0.149
[1394,1577]  [  12, 138]        0.104  0.115    -     -     -     -     -
[ 240, 420]  [  48, 201]        0.112  0.251 0.390 0.411    -     -     -
[ 947,1127]  [ 768, 859]        0.111  0.105 0.120 0.154 0.156 0.255 0.184
```

Typically the agreement is destroyed by ten more nodes on each side — the
residual roughly doubles. One survivor holds together out to $120$ nodes. This
is the quantitative form of "it lines up, and diverges when I move a little
farther along".

## Caveats

* $0.10$ of the diameter is not a close copy: exact or near-exact copies in this
  project sit at $0.003$–$0.02$. These are coarse shape coincidences.
* The counts depend on the grid (steps of $8$–$12$ nodes) and on the $24$-point
  coarse stage; small or very short overlaps are missed, and the thresholds are
  arbitrary. The real-versus-control contrast is the robust part, not the counts.
* Only $m=1600$ was searched.
* The correspondence is affine. A non-affine one would find more (and would be a
  different question).
* "Scale free shape distance after fitting a rotation" is one metric among
  several; it ignores how the pieces sit in the plane.

## Reproduce

```bash
python3 scripts/search_overlaps.py --out results/overlaps_m1600.json
python3 scripts/search_overlaps.py --shuffle-turns 1     # control
python3 scripts/search_overlaps.py --surrogate 1         # control
python3 scripts/plot_overlaps.py
```
