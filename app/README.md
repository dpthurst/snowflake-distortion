# Alignment explorer

`app/snowflake-align.html` is a single self-contained page: open it in a browser,
no server and no network needed.

It lines up a saved configuration (the small-$m$ ladder and the arc and
three-piece branches) with a contiguous shortening of itself or of another
curve, in two modes:

* **index correspondence** — overlay at the same parameter, endpoint-gauged,
  with the four symmetries available; residual in units of the endpoint chord;
* **residual metric** — either the gap at equal parameter (affine between
  knots, hence V-shaped where the curves cross), or the distance from each
  curve to the piecewise-linear interpolant of the other, taken in both
  directions and reported as a symmetric maximum;
* **distortion heatmap** — a third panel draws $r(x,y)$ on $I\times I$ for the
  whole curve: blue at the minimum, red at the maximum, white at their
  geometric mean (a log scale, so the colour change sits at the logarithmic
  centre). The minimum and maximum pairs are marked, as are crosses at the
  current chord. It is computed in the browser once per curve and cached;
* **copy mode** — fit ``w(t) ~ c W(s0 + s1 t) + d``: a similarity image of the
  whole with a free affine parameter. Six free parameters, found by a coarse
  grid over the parameter range followed by two local refinements; the sliders
  set scale, rotation, s0 and s1 and re-fit the translation.

With the block taken from the same curve the index correspondence matches it
exactly, so copy mode is informative only when the two curves differ.

Changing the interval or either curve leaves the copy where you put it: the
placement is independent of which block is being matched, and only *refit all
six* recomputes it. The same holds for the parameter-range sliders — shortening
the copy's arc keeps its scale and position. The scale and rotation sliders do
re-anchor the translation, so that a change of size keeps the copy over the
block rather than swinging it away.

In copy mode the green copy can be dragged directly on the overlay: with the
drag control on *moves the copy* a drag translates it; set it to *rotates and
scales it* (or hold Shift) and the same drag applies $k=z_1/z_0$ about the
copy's centroid, so the gesture rotates and scales together. The overlay frame
is fixed by the whole curve and the block, so the view does not chase the copy
while it moves. The sliders remain, and *refit all six* re-runs the search.

The page is generated from `scripts/align_app_template.html` plus the saved
arrays; rebuild it after changing either:

    python3 scripts/build_align_app.py

Checks (see `scripts/`):

    python3 scripts/test_align_app.py    # browser geometry against numpy, under node
    ~/node/bin/node scripts/test_align_ui.js   # wiring in a DOM stub, not a GUI test
    python3 scripts/preview_align_app.py      # static PNG, index mode
    python3 scripts/preview_align_copy.py     # static PNG, copy mode
    python3 scripts/preview_residual_metrics.py  # the two metrics side by side

Coordinates are embedded at full precision. Rounding them to 7 decimals moved
the $m=800$ arg-min from the gap-90 pair $(509,599)$ to a gap-1 pair, because
the shortest chords are sensitive at the $10^{-5}$ level; a test now compares
the embedded arrays with the source `.npy` files bit for bit.

There is no browser on this machine, so drawing, dragging and layout have not
been tested; only the geometry, the heatmap statistics and the event wiring
have.
