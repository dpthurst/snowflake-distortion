# Alignment explorer

`app/snowflake-align.html` is a single self-contained page: open it in a browser,
no server and no network needed.

It lines up a saved configuration (the small-$m$ ladder and the arc and
three-piece branches) with a contiguous shortening of itself or of another
curve. The overlay is endpoint-gauged, with the four symmetries available and
the residual in units of the endpoint chord.

The page is generated from `scripts/align_app_template.html` plus the saved
arrays; rebuild it after changing either:

    python3 scripts/build_align_app.py

Checks (see `scripts/`):

    python3 scripts/test_align_app.py    # browser geometry against numpy, under node
    ~/node/bin/node scripts/test_align_ui.js   # wiring in a DOM stub, not a GUI test
    python3 scripts/preview_align_app.py # static PNG of the same view

There is no browser on this machine, so drawing, dragging and layout have not
been tested; only the geometry and the event wiring have.
