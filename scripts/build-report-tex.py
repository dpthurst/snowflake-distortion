"""Build a report: reports/<name>.md -> .tex -> .pdf.

Usage (from the project root):
    python3 scripts/build-report-tex.py                  # expansion-prior-work
    python3 scripts/build-report-tex.py diameter-growth-2cat

The Markdown is the master copy; the .tex is generated and the .pdf compiled.
Run this after every edit to a report's Markdown so that its PDF never lags
behind.  LaTeX intermediates go to build/report-tex/<name>/.
"""
import io, os, re, shutil, subprocess, sys

NAME = sys.argv[1] if len(sys.argv) > 1 else 'expansion-prior-work'
TITLE = {
    'expansion-prior-work': 'Related prior work for Expansion',
    'diameter-growth-2cat': 'Diameter growth for the 2Cat endomorphism',
}.get(NAME, NAME.replace('-', ' '))

SRC = 'reports/%s.md' % NAME
MD = 'build/report-tex/%s-tex.md' % NAME
RAW = 'build/report-tex/%s-raw.tex' % NAME
OUT = 'reports/%s.tex' % NAME
PDF = 'reports/%s.pdf' % NAME
BUILD = 'build/report-tex/%s' % NAME

# ---- 1. preprocess Markdown for pandoc -------------------------------------
lines = io.open(SRC, encoding='utf-8').read().split('\n')
assert lines[0].startswith('# '), lines[0]
lines = lines[1:]
lines = [l for l in lines if not l.strip().startswith('*Prepared by Wren')]

out = []
for l in lines:
    if l.startswith('### '):
        out.append('## ' + re.sub(r'^###\s+\d+(\.\d+)*\.?\s+', '', l))
    elif l.startswith('## '):
        out.append('# ' + re.sub(r'^##\s+\d+(\.\d+)*\.?\s+', '', l))
    else:
        out.append(l)

res = []
for i, l in enumerate(out):
    res.append(l)
    if (re.match(r'^\*\*[^*].*\*\*$', l.strip())
            and i + 1 < len(out) and re.match(r'^\s*[-*] ', out[i + 1])):
        res.append('')
io.open(MD, 'w', encoding='utf-8').write('\n'.join(res))

# ---- 2. pandoc -------------------------------------------------------------
subprocess.run(['pandoc', MD, '-f', 'markdown', '-t', 'latex', '--standalone',
                '--ascii', '--number-sections',
                '--metadata', 'title=%s' % TITLE,
                '--metadata', 'author=Wren, an Open Math OS Agent', '-o', RAW], check=True)

# ---- 3. post-process -------------------------------------------------------
tex = io.open(RAW, encoding='utf-8').read()

def sub_once(old, new):
    global tex
    n = tex.count(old)
    assert n == 1, (n, old[:60])
    tex = tex.replace(old, new)

PREAMBLE = ('\\usepackage[margin=1in]{geometry}\n'
            '\\allowdisplaybreaks\n'
            '\\setlength{\\emergencystretch}{3em}')
anchor = '\\usepackage{longtable,booktabs,array}'
if anchor in tex:                       # pandoc emits this only when the
    sub_once(anchor, anchor + '\n' + PREAMBLE)   # source has a table
else:
    sub_once('\\usepackage{amsmath,amssymb}',
             '\\usepackage{amsmath,amssymb}\n' + PREAMBLE)
sub_once('\\maketitle', '\\maketitle\n\\tableofcontents\n\\newpage')
n = tex.count('A\\(^\\#\\)')
tex = tex.replace('A\\(^\\#\\)', 'A$^{\\#}$')
tex = tex.replace('\\date{}', '\\date{Compiled \\today}')
nsec = tex.count('\u00a7')
tex = tex.replace('\u00a7', 'Section ')
tex = tex.replace('\u2248', '$\\sim$')
# pandoc renders a backticked URL as \texttt{...}, which cannot break across
# lines; use \url instead (the previous build did the same).
def _url(m):
    inner = m.group(1)
    inner = inner.replace('\\textasciitilde{}', '~')
    assert '{' not in inner and '}' not in inner, inner
    return '\\url{' + inner + '}'
tex, nurl = re.subn(r'\\texttt\{(https?://(?:[^{}]|\\textasciitilde\{\})*)\}', _url, tex)
assert all(ord(c) < 128 for c in tex), 'non-ASCII left in the .tex'

io.open(OUT, 'w', encoding='utf-8').write(tex)
print('wrote %s (%d lines), patched %d A-sharp, %d section signs, %d urls'
      % (OUT, tex.count('\n') + 1, n, nsec, nurl))

# ---- 4. compile the PDF ----------------------------------------------------
os.makedirs(BUILD, exist_ok=True)
for attempt in range(4):
    r = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                        '-output-directory', BUILD, OUT],
                       capture_output=True, text=True)
    log = os.path.join(BUILD, '%s.log' % NAME)
    tail = io.open(log, encoding='utf-8', errors='replace').read() if os.path.exists(log) else ''
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-4000:])
        sys.exit('pdflatex failed (pass %d); see %s' % (attempt + 1, log))
    if 'Rerun to get' not in tail:
        break
else:
    sys.exit('pdflatex still asking for a rerun after 4 passes')

built = os.path.join(BUILD, '%s.pdf' % NAME)
shutil.copyfile(built, PDF)
over = tail.lower().count('overfull')
info = subprocess.run(['pdfinfo', PDF], capture_output=True, text=True).stdout
pages = re.search(r'^Pages:\s+(\d+)', info, re.M)
print('wrote %s (%s pages, %d overfull boxes)'
      % (PDF, pages.group(1) if pages else '?', over))
