# Imported from the Expansion project

On 2026-11-15 this repository was populated by copying the working tree of the
project titled **Expansion** (project id `p-bcf364e629b24a6d9bd061bbabc49e45`,
`repo/`), because the new project's root was a fresh empty git repository and
the old tree was not shared. The git history came with it: 64 commits, HEAD
`65f742f`.

The import is deliberately narrowed to the snowflake-distortion material. The
Expansion paper's own sources and literature thread stayed behind. Nothing was
deleted from the source; the old project still holds every item listed below.

## Kept

| path | contents |
|---|---|
| `.ultramath/` | all 38 swarm reports and coordinate arrays, six runs |
| `.scratch/` | snowflake optimisation code and logs from the earlier sessions |
| `results/`, `scripts/` | the distortion reports, checks and plotting scripts |
| `figures/` | Koch, finite-$m$ cloud, IFS and tie figures |
| `notes/` | `koch-renormalisation.md`, `koch-extremals.tex`/`.pdf` |
| `lean/` | `KochSnowflake.lean`, `Test.lean` |
| `refs/` | nine conformal-dimension and self-similar-set PDFs only |
| `.skill/` | project skills |
| `stationary-cloud-m100.md` | the $m=100$ stationary cloud note |

## Left behind in the Expansion project

- The paper sources (`Expansion.tex`, `preamble.tex`, `intro.tex`,
  `alternatives.tex`, `clippings.tex`, `comb-models.tex`, `covers.tex`,
  `torus.tex`, `visual-metrics.tex`, `min-growth.tex`, `Makefile`,
  `Expansion.xmpdata`) and `min-growth-fix.patch`.
- `reports/expansion-prior-work.{md,tex,pdf}` and `handoffs/prior-work.md`,
  the Expansion literature memo and its handoff.
- `attachments/`, and the expansion-dynamics PDFs in `refs/` (expansivity,
  expanding maps, Lyapunov exponents, Lebesgue numbers, chain recurrence,
  expansion rates, Hruška, Milnor, Kifer, Schreiber, Dedekind).

Those PDFs are re-downloadable; `.gitignore` already treats `refs/*.pdf` as
disposable. The removal of the three tracked files is recorded in the import
commit.

## Checks

- `git fsck` reports no corruption (only dangling trees, the blobs of the
  Expansion sources that were staged but never committed in the source repo).
- `git log` reproduces the source history: `65f742f`, `4381587`, `b1d5925`, ...
- 38 swarm reports present; `.ultramath` 122 MB, `.scratch` 35 MB,
  `figures` 5.1 MB, 202 MB total.
- Worktree state matches the source apart from the pruned items: 41 untracked
  entries, the same figures, notes, Lean and skill directories.
