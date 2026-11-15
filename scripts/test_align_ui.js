/* Wiring test for app/snowflake-align.html, without a browser.
 *
 * There is no browser on this machine, so this is not a GUI test: it builds a
 * minimal DOM stub, runs the page's own scripts in it, and checks that the page
 * loads, that every element id the code asks for exists in the HTML, that the
 * readout appears, and that changing the interval, the curve or a preset
 * updates the readout without throwing.  Canvas drawing is stubbed out.
 *
 * Run:  ~/node/bin/node scripts/test_align_ui.js [path-to-html]
 */
const fs = require('fs');
const vm = require('vm');
const path = require('path');

const file = process.argv[2] || path.join(__dirname, '..', 'app', 'snowflake-align.html');
const html = fs.readFileSync(file, 'utf8');

/* ids that actually appear in the static markup */
const ids = new Set([...html.matchAll(/id="([^"]+)"/g)].map(m => m[1]));

function canvasCtx() {
  const base = {
    createImageData: (w, h) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
    putImageData: () => undefined,
    getImageData: (x, y, w, h) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
    drawImage: () => undefined,
    save: () => undefined, restore: () => undefined,
    beginPath: () => undefined, moveTo: () => undefined, lineTo: () => undefined,
    arc: () => undefined, fill: () => undefined, stroke: () => undefined,
    clearRect: () => undefined, fillRect: () => undefined, fillText: () => undefined,
    setLineDash: () => undefined,
  };
  return new Proxy(base, {
    get: (t, k) => (k in t ? t[k] : () => undefined),
    set: (t, k, v) => { t[k] = v; return true; },
  });
}

function el(tag = 'div', id = '') {
  const listeners = {};
  const e = {
    tagName: tag.toUpperCase(), id, value: '', textContent: '', innerHTML: '',
    children: [], style: {}, width: 440, height: 440, max: '0',
    appendChild(c) { this.children.push(c); return c; },
    addEventListener(t, f) { (listeners[t] = listeners[t] || []).push(f); },
    dispatch(t, ev) { (listeners[t] || []).forEach(f => f(ev)); },
    getContext() { return canvasCtx(); },
    getBoundingClientRect() { return { left: 0, top: 0, width: this.width, height: this.height }; },
    setPointerCapture() {},
  };
  return e;
}

const nodes = new Map();
for (const id of ids) nodes.set('#' + id, el('div', id));

const document = {
  querySelector(sel) { return nodes.get(sel) || null; },
  createElement(tag) { return el(tag); },
};
const window = { addEventListener(t, f) { (window._l = window._l || {})[t] = f; } };

const ctx = vm.createContext({ document, window, console, Math, Float64Array, Float32Array, Array, Object, Number, JSON, Uint8ClampedArray });

const blocks = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);
if (blocks.length !== 3) { console.error(`expected 3 script blocks, found ${blocks.length}`); process.exit(1); }
for (const b of blocks) vm.runInContext(b, ctx, { filename: 'inline.js' });

const readout = () => nodes.get('#readout').textContent;
const problems = [];
function check(name, cond, detail) {
  if (cond) { console.log(`  ok   ${name}`); }
  else { console.log(`  FAIL ${name}${detail ? ': ' + detail : ''}`); problems.push(name); }
}

console.log(`ids in the page: ${[...ids].join(', ')}`);
check('page runs and writes a readout', readout().includes('interval'), readout().slice(0, 80));
check('curve pickers were populated', nodes.get('#whole').children.length > 1,
      `${nodes.get('#whole').children.length} options`);

const app = window.__app;
check('app handle exposed', !!app);

/* put the page on the m=800 curve, as a user would */
const wholeSel = nodes.get('#whole'), sourceSel = nodes.get('#source');
wholeSel.value = 'S-branch m=800';
wholeSel.onchange({ target: wholeSel });
sourceSel.value = 'S-branch m=800';
sourceSel.onchange({ target: sourceSel });
check('default curve is the m=800 branch', readout().includes('D whole'), readout().slice(0, 60));

/* moving the interval changes the numbers */
const before = readout();
app.setInterval(259, 511);
const after = readout();
check('setInterval updates the readout', before !== after);
check('the readout names the new interval', after.includes('[259, 511]'), after.split('\n')[0]);

/* the sup printed for a known interval matches the numpy value in the reports */
const rep = app.report('S-branch m=800', 'S-branch m=800', 259, 511);
check('sup for [259,511] matches the report', Math.abs(rep.sup - 0.102082610) < 5e-8,
      `${rep.sup}`);

/* a preset moves the interval and switches the curve if needed */
const preset = nodes.get('#preset');
preset.value = '0';
preset.onchange({ target: preset });
check('preset applied', readout().includes('[259, 511]'), readout().split('\n')[0]);

/* switching the whole curve does not throw and updates the label */
const whole = nodes.get('#whole');
whole.value = 'S-branch m=400';
whole.onchange({ target: whole });
check('switching curve updates the readout', readout().includes('D whole'),
      readout().slice(0, 60));

/* variant switch */
const variant = nodes.get('#variant');
variant.value = 'rev';
variant.onchange({ target: variant });
check('variant switch keeps a readout', readout().includes('reversal') || readout().includes('rev'),
      readout().split('\n')[1]);

/* distortion heatmap panel */
const heatTitle = readout();   // just to force a draw before reading
nodes.get('#an').value = '0'; nodes.get('#an').oninput({ target: nodes.get('#an') });
nodes.get('#bn').value = '799'; nodes.get('#bn').oninput({ target: nodes.get('#bn') });
const ht = nodes.get('#heatTitle').textContent;
check('heatmap panel titles the current curve', /r\(x,y\) for/.test(ht) && /D =/.test(ht), ht);
const beforeCurve = nodes.get('#heatTitle').textContent;
const wholeH = nodes.get('#whole');
wholeH.value = 'S-branch m=100'; wholeH.onchange({ target: wholeH });
check('heatmap follows the curve switch', nodes.get('#heatTitle').textContent !== beforeCurve &&
      /m=100/.test(nodes.get('#heatTitle').textContent), nodes.get('#heatTitle').textContent);
wholeH.value = 'S-branch m=800'; wholeH.onchange({ target: wholeH });

/* residual metric switch */
const metric = nodes.get('#metric');
const beforeMetric = readout();
metric.value = 'setdist';
metric.onchange({ target: metric });
check('distance metric changes the readout', readout() !== beforeMetric &&
      /distance sup/.test(readout()), readout().split('\n')[3]);
metric.value = 'pointwise';
metric.onchange({ target: metric });
check('switching back restores the gap readout', /residual sup/.test(readout()));

/* copy mode: similarity + affine parameter */
const mode = nodes.get('#mode');
mode.value = 'copy';
let copyThrew = null;
try { mode.onchange({ target: mode }); } catch (e) { copyThrew = e; }
check('switching to copy mode does not throw', !copyThrew, copyThrew && copyThrew.message);
const txt = readout();
check('copy mode reports a similarity', /\|c\|/.test(txt) && /parameter s/.test(txt),
      txt.split('\n').slice(3, 6).join(' | '));
check('copy-mode sliders were initialised', nodes.get('#lscale').value !== '' &&
      nodes.get('#s0').value !== '', `lscale=${nodes.get('#lscale').value} s0=${nodes.get('#s0').value}`);
const beforeScale = readout();
const ls = nodes.get('#lscale');
ls.value = String(Number(ls.value) + 0.2);
ls.oninput({ target: ls });
check('moving the scale changes the residual', readout() !== beforeScale);
let refitThrew = null;
try { nodes.get('#refit').onclick({}); } catch (e) { refitThrew = e; }
check('refit does not throw', !refitThrew, refitThrew && refitThrew.message);

/* dragging the copy on the overlay */
const app2 = window.__app;
const srcSel = nodes.get('#source');
srcSel.value = 'S-branch m=40'; srcSel.onchange({ target: srcSel });
const cvR = nodes.get('#cvRight');
function ptr(x, y, extra) { return Object.assign({ clientX: x, clientY: y, pointerId: 3,
                                                  preventDefault() {} }, extra || {}); }
const c0 = { ...app2.state.copy };
cvR.dispatch('pointerdown', ptr(200, 200));
cvR.dispatch('pointermove', ptr(240, 210));
cvR.dispatch('pointerup', ptr(240, 210));
const c1 = { ...app2.state.copy };
check('dragging translates the copy', Math.hypot(c1.dx - c0.dx, c1.dy - c0.dy) > 1e-6 &&
      Math.hypot(c1.cx - c0.cx, c1.cy - c0.cy) < 1e-12,
      `d moved by ${Math.hypot(c1.dx - c0.dx, c1.dy - c0.dy).toExponential(2)}, ` +
      `|c| ${Math.hypot(c0.cx, c0.cy).toFixed(4)} -> ${Math.hypot(c1.cx, c1.cy).toFixed(4)}`);
cvR.dispatch('pointerdown', ptr(200, 200));
cvR.dispatch('pointermove', ptr(240, 210));
cvR.dispatch('pointerup', ptr(240, 210));
const c2 = { ...app2.state.copy };
check('the frame stays put while dragging (same pixel move, same shift)',
      Math.abs((c2.dx - c1.dx) - (c1.dx - c0.dx)) < 1e-9 &&
      Math.abs((c2.dy - c1.dy) - (c1.dy - c0.dy)) < 1e-9,
      `${(c1.dx - c0.dx).toFixed(6)} vs ${(c2.dx - c1.dx).toFixed(6)}`);

/* transform: check the algebra exactly, not just that something changed */
const view = app2.view(), pivot = app2.pivot();
const toData = (px, py) => [view.cx + (px - view.w / 2) / view.s,
                            view.cy - (py - view.h / 2) / view.s];
const p0d = toData(260, 240), p1d = toData(300, 260);
const z0 = [p0d[0] - pivot[0], p0d[1] - pivot[1]];
const z1 = [p1d[0] - pivot[0], p1d[1] - pivot[1]];
const kd = z0[0] * z0[0] + z0[1] * z0[1];
const kx = (z1[0] * z0[0] + z1[1] * z0[1]) / kd, ky = (z1[1] * z0[0] - z1[0] * z0[1]) / kd;
const cPre = { ...app2.state.copy };
const dm = nodes.get('#dragmode');
dm.value = 'transform';
const scaleBefore = Math.hypot(app2.state.copy.cx, app2.state.copy.cy);
cvR.dispatch('pointerdown', ptr(260, 240));
cvR.dispatch('pointermove', ptr(300, 260));
cvR.dispatch('pointerup', ptr(300, 260));
const scaleAfter = Math.hypot(app2.state.copy.cx, app2.state.copy.cy);
check('transform mode changes the scale', Math.abs(scaleAfter - scaleBefore) > 1e-6,
      `|c| ${scaleBefore.toFixed(4)} -> ${scaleAfter.toFixed(4)}`);
check('transform mode keeps the copy finite', Number.isFinite(app2.state.copy.dx) &&
      scaleAfter > 1e-3 && scaleAfter < 1e3, `|c| ${scaleAfter}`);
{
  const cNow = app2.state.copy;
  const expC = [cPre.cx * kx - cPre.cy * ky, cPre.cx * ky + cPre.cy * kx];
  const dx = cPre.dx - pivot[0], dy = cPre.dy - pivot[1];
  const expD = [pivot[0] + (kx * dx - ky * dy), pivot[1] + (kx * dy + ky * dx)];
  const errC = Math.hypot(cNow.cx - expC[0], cNow.cy - expC[1]);
  const errD = Math.hypot(cNow.dx - expD[0], cNow.dy - expD[1]);
  check('the transform is k = z1/z0 about the copy centroid', errC < 1e-9 && errD < 1e-9,
        `errors ${errC.toExponential(2)}, ${errD.toExponential(2)}`);
}

/* pointer drag on the left canvas */
const cv = nodes.get('#cvLeft');
let threw = null;
try {
  cv.dispatch('pointerdown', { clientX: 100, clientY: 100, pointerId: 1 });
  cv.dispatch('pointermove', { clientX: 120, clientY: 110, pointerId: 1 });
  cv.dispatch('pointerup', { pointerId: 1 });
} catch (e) { threw = e; }
check('pointer drag does not throw', !threw, threw && threw.message);

/* keyboard nudge */
let kthrew = null;
try { window._l.keydown({ key: 'ArrowRight', shiftKey: true, preventDefault() {} }); }
catch (e) { kthrew = e; }
check('arrow key does not throw', !kthrew, kthrew && kthrew.message);

if (problems.length) { console.error(`${problems.length} problem(s)`); process.exit(1); }
console.log('wiring test passed (no browser: drawing itself is untested)');
