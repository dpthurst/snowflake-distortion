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
* **copy mode** — fit ``w(t) ~ c W(s0 + s1 t) + d``: a similarity image of the
  whole with a free affine parameter. Six free parameters, found by a coarse
  grid over the parameter range followed by two local refinements; the sliders
  set scale, rotation, s0 and s1 and re-fit the translation.

With the block taken from the same curve the index correspondence matches it
exactly, so copy mode is informative only when the two curves differ.

The page is generated from `scripts/align_app_template.html` plus the saved
arrays; rebuild it after changing either:

    python3 scripts/build_align_app.py

Checks (see `scripts/`):

    python3 scripts/test_align_app.py    # browser geometry against numpy, under node
    ~/node/bin/node scripts/test_align_ui.js   # wiring in a DOM stub, not a GUI test
    python3 scripts/preview_align_app.py      # static PNG, index mode
    python3 scripts/preview_align_copy.py     # static PNG, copy mode
    python3 scripts/preview_residual_metrics.py  # the two metrics side by side

There is no browser on this machine, so drawing, dragging and layout have not
been tested; only the geometry and the event wiring have.
