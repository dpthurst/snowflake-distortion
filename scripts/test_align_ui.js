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
  return new Proxy({}, {
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

const ctx = vm.createContext({ document, window, console, Math, Float64Array, Array, Object, Number, JSON });

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
