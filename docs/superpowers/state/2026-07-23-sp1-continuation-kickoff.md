# SP1 Continuation Kickoff — Stock/TEA Kernel (mid-execution handoff)

I'm continuing **SP1 (exact stock model + TEA kernel)** of the certified
adaptive-clearing programme in **`/Users/jelle/Code/CADCAM/compas_cgal_prs`**,
branch **`jf/adaptive-clearing-sp1`**. Commit `d01879f`
(`feat(engagement): typed engagement audit replay (Task 6b)`) is the last
hard checkpoint; **a merge-repair agent was IN FLIGHT at handoff** — collect
or complete it first (§2).

## 1. What I want you to do

Resume subagent-driven execution of the SP1 plan: use
`superpowers:subagent-driven-development` (fresh Opus subagent per task —
NEVER spawn Fable agents; parent-session rules in `~/.claude/CLAUDE.md`
apply, including the ETH/DeepMind bar and the choice-menu rule). The
progress ledger is the single source of truth for what is DONE — do not
re-dispatch completed tasks.

## 2. Read first, in dependency order

1. `.superpowers/sdd/progress.md` — the execution ledger: Tasks 1–7 + 6a/6b
   CLOSED with review verdicts, all adjudications, the open Minor-findings
   triage list, and the merge-hole record. Trust it over recollection.
2. `docs/plans/2026-07-23-sp1-stock-tea-kernel.md` — the SP1 plan WITH its
   amendment blocks (Tasks 4/5/6/7 amendments supersede their older text).
   Remaining: Task 8 (baseline), Task 9 (wall-clock gate + CHANGELOG).
3. `docs/plans/2026-07-23-adaptive-clearing-design.md` — the programme spec
   (§2 v1 claim, §3 stock model incl. conservative under-covering doctrine).
4. `CLAUDE.md` — repo law: "Exact-Kernel Discipline", "Numeric comparison is
   the exact-kernel idiom", build/test commands. Full doctrine:
   `docs/exactness.md` (mkdocs Design Notes page) + verbatim reference
   `docs/dev/exact-cgal-idioms.md`.
5. `.superpowers/sdd/merge-repair-report.md` — IF PRESENT, the in-flight
   repair's report. If absent, see the collection protocol below.
6. `.superpowers/sdd/task-N-brief.md` / `task-N-report.md` — per-task record
   (N ∈ 1..7, 6a, 6b; plus `review-*.diff` packages).

**Merge-repair collection protocol.** The repair closes a certified
soundness hole: `max_run_tea` jumps discontinuously at run-MERGES
(O(1) vs the lemma's O(√d)), so `certify_segment_tea` was unsound for merge
configurations. Approved repair (user-sanctioned, option A): **gap-closure
pessimism** — `engagement_at` gains `gap_close_ratio` (chord form
`4·sin²(γ/2)`, default 0.0 = old semantics); at each station, void gaps with
span ≤ γ are absorbed (a chain of runs+gaps is a SINGLE arc → existing exact
chord/orientation predicates apply — no angle sums); the certifier passes
`γ_guard = 2·GROWTH(hs)` per level; cap decision on PESSIMISTIC runs,
reported totals stay TRUE measures; lemma comment rewritten with the
contradiction argument (a gap that closes within the step was already ≤
γ_guard at the station). Required evidence tests (a) raw-measure merge jump
EXCEEDS raw lemma (inverted assertion, documents the hole), (b) pessimism
unit test, (c) `cap_certified=False` across the merge at cap < merged span,
(d) pessimistic-measure growth stays ≤ raw lemma, (e) grazing-incidence
docstring softening + observed-headroom pin in
`tests/test_engagement_oracle.py`. Check `git log --oneline -3` for
`fix(engagement): gap-closure pessimism closes certificate merge hole`:
- **Committed** → generate a review package
  (`<skill-dir>/scripts/review-package <parent> <sha>`) and dispatch a
  scoped Opus reviewer verifying the fixpoint closure, the rewritten lemma
  argument, and tests a–d; then proceed to Task 8.
- **Uncommitted / partial** (at handoff: `src/engagement_2.{h,cpp}` modified,
  tests untouched) → dispatch a fresh Opus implementer with the spec above
  (full detail also in the ledger), starting from the tree; TDD; single
  commit; then the scoped review.
- **If test (a) fails to demonstrate or (d) fails to pass** → STOP, report
  the measurements to the user; do not tune geometry or widen bounds.

## 3. Boundary contract (what is already landed and callable)

```python
# compas_cgal._stock_2  (built extension; sources src/stock_2.*, src/engagement_2.*)
Stock2(boundary: float64[N,3], holes: list[float64[M,3]])
Stock2.contains(x, y) -> bool;  Stock2.is_empty() -> bool
Stock2.subtract_capsule(x0, y0, x1, y1, radius) -> None       # certified under-covering disk chain
Stock2.subtract_disk(cx, cy, radius) -> None                   # exact
Stock2.subtract_arc_sweep(cx, cy, sx, sy, ex, ey, cw, tool_radius) -> None
engagement_at(stock, cx, cy, tool_radius, cap_chord_ratio) -> (total_tea, max_run_tea, exceeded)
    # cap_chord_ratio = 4*sin^2(cap/2) in (0, 4]; `exceeded` is the EXACT verdict
    # AFTER the repair: + gap_close_ratio: float = 0.0 (pessimistic merge closing)
certify_segment_tea(stock, x0, y0, x1, y1, tool_radius, cap_radians) -> (max_tea, cap_certified, stations)
    # 0 < cap <= pi validated at the binding; guarded exact thresholds per level
_sign_mixed_radical(a, b, c, d, alpha, beta) -> int            # test-only exact primitive

# compas_cgal.stock:  Stock(polygon, holes=None) — typed wrapper, .raw -> Stock2
# compas_cgal.engagement:  audit_toolpath_engagement(polygon, result, tool_diameter, tea_cap, holes=None)
#   -> EngagementReport(tool_diameter, tea_cap, operations: list[OperationEngagement],
#                       max_tea, cap_violations, engaged_ops)   # frozen dataclasses
```

Suite at `d01879f`: **154 passed** (45 toolpath untouchable + stock/
engagement/oracle/audit). All decisions in Epeck-land are exact predicates;
`to_double`/`atan2` are reporting-only.

## 4. Deliverables remaining (plan Tasks 8–9 + wrap)

1. Merge-repair collection + scoped review (§2).
2. **Task 8** — `scripts/engagement_baseline.py` +
   `docs/benchmarks/engagement_baseline.{md,json}` + one subprocess test in
   `tests/test_engagement_audit.py`. MUST run AFTER the repair (baseline
   reflects the SOUND certifier). Full-toolpath audits are expensive
   (Task 6b measured >120 s for one full square audit) — run with generous
   timeouts, RECORD wall-clock per pocket in the report (BLUF first;
   "This baseline is the reference SP2 must beat.").
3. **Task 9** — `tests/test_stock_gate.py` wall-clock exit gate
   (≤ 10 s/pocket commit+certify throughput; an honest FAIL is a result to
   report, not tune) + CHANGELOG entry (plan Task 9 has the text).
4. **Whole-branch final review** — Opus subagent per
   `superpowers:requesting-code-review`, package
   `review-package $(git merge-base jf/toolpath-redesign HEAD) HEAD`. Hand
   it the ledger's Minor-findings list for triage. The acceptance bar is the
   user's, verbatim: *"if its not ETH deepmind grade engineering it doesnt
   go in our repo."*
5. Then `superpowers:finishing-a-development-branch`. NEVER merge/rebase/
   push without the user's explicit instruction.

## 5. Exit gates (plan §8, quoted)

- "(a) TEA kernel agrees with Monte-Carlo oracle within sampling bounds on
  the benchmark suite" — DONE (Task 7, 8 seeds).
- "(b) engagement audit of the existing trochoid generator produces a
  truthful per-op TEA report … the number SP2 must beat" — audit DONE
  (6b); the BASELINE REPORT is Task 8.
- "(c) measured wall-clock ≤ 10 s per benchmark pocket for commit+certify
  throughput — a gate, not a guess; if missed, the fallback stays within
  CGAL (coarser exact station grids), never an approximate clipper" —
  Task 9. If the gate slips, write
  `docs/superpowers/state/sp1-gate-c-analysis.md` BEFORE any fallback work
  and present options to the user (choice menu pre-filtered to the bar).

## 6. Decoupling rule

Everything remaining runs against committed interfaces — no dependency on
uncommitted state except the merge-repair tree (§2). Task 8/9 test
`compas_cgal.engagement`/`_stock_2` only; the 45 `tests/test_toolpath.py`
tests are a frozen contract (NEVER modify to make anything pass).

## 7. Load-bearing conventions

- Exact-Kernel Discipline (CLAUDE.md): no epsilon/deflation/ulp handling in
  ANY decision path; `CGAL::sign`/`CGAL::compare` at the number-type level;
  cross-root `Sqrt_extension` ARITHMETIC is UB (comparisons fine); doubles
  cross into FT exactly ONCE at documented boundaries; analytic lemma
  bounds only for refinement, with integer safety factors + stated safe
  direction.
- One agent per file set — a directive to a stood-down agent is a
  RE-ACTIVATION (this bit us once: see ledger COLLISION entry). Reviews are
  read-only and may pipeline; implementers serialize on shared files.
- Commits: author+committer Jelle Feringa <jelleferinga@gmail.com>, concise
  conventional-commit messages, NO attribution lines. Never commit without
  standing instruction (plan tasks carry commit steps = standing).
- `docs/examples/example_isolines.py` is a PRE-EXISTING user modification:
  never stage, never revert.
- Linear history; never rewrite committed history.
- pytest ALWAYS `-n auto`. Ruff before commits (Python).

## 8. Environment booby traps

- Python: `PY=/Users/jelle/mambaforge/envs/cgal-dev/bin/python` (py312).
  The env's editable compas_cgal install points at the CLEAN sibling
  checkout — ALWAYS run with `PYTHONPATH=src` from this repo, or you'll
  silently import the wrong package.
- Rebuild C++ (no pip involved):
  `cmake -S . -B build/cp312-abi3-macosx_15_0_arm64 && cmake --build build/cp312-abi3-macosx_15_0_arm64 --target _stock_2 -j 8`
  (`--target _toolpath` for the trochoid module). `src/compas_cgal/*.so`
  are symlinks into `build/` (gitignored).
- Test: `PYTHONPATH=src $PY -m pytest tests/ -n auto -q`.
- nanobind: NB_STATIC — one module per registry; `_stock_2` is built from
  TWO sources (stock_2.cpp + engagement_2.cpp, CMakeLists ~line 242). NEVER
  give NB_MODULE functions eager `std::vector<...>` DEFAULT args
  (std::bad_cast at import).
- The 10 pytest warnings are pre-existing compas `LooseVersion`
  deprecations — not ours, don't chase them.
- mkdocs deps (mkdocs-material, markdown-callouts, mkdocstrings[python])
  are installed in cgal-dev; `$PY -m mkdocs build` validates docs.

## 9. Suggested file structure (remaining work)

```
compas_cgal_prs/
├── src/engagement_2.{h,cpp}      # MODIFY (merge-repair: gap_close_ratio)
├── tests/test_stock.py           # MODIFY (repair tests b, c)
├── tests/test_engagement_oracle.py  # MODIFY (repair tests a, d + grazing/headroom)
├── scripts/engagement_baseline.py   # CREATE (Task 8)
├── docs/benchmarks/engagement_baseline.{md,json}  # CREATE (Task 8, generated+committed)
├── tests/test_engagement_audit.py   # MODIFY (Task 8 subprocess test)
├── tests/test_stock_gate.py      # CREATE (Task 9)
└── CHANGELOG.md                  # MODIFY (Task 9)
```

## 10. Risk + budget

- Gate (c) is the SP1 risk: exact-boolean cost grows with accumulated
  subtractions (Task 6b: full square audit >120 s). The gate measures
  commit+certify throughput (60 passes), not full-audit cost — keep the two
  separate in reporting. Fallback direction is pre-decided in the spec:
  coarser exact station grids, never an approximate clipper.
- Merge-repair evidence tests may legitimately FAIL (that's data, not a
  bug) — escalate with measurements, never tune.

## 11. Reviews & quality machinery

- Skill scripts:
  `/Users/jelle/.claude/plugins/cache/claude-plugins-official/superpowers/6.1.1/skills/subagent-driven-development/scripts/{task-brief,review-package}`.
- Reviewer dispatches: read-only, diff-file based, named verifications,
  no pre-judging findings; Important findings → fix subagent → re-review.
- Ledger discipline: append one line per event to
  `.superpowers/sdd/progress.md`; it is the recovery map.

## 12. Start

```bash
cd /Users/jelle/Code/CADCAM/compas_cgal_prs && git log --oneline -3 && git status -s && cat .superpowers/sdd/progress.md | tail -20
```

Then apply §2's merge-repair collection protocol.

End of prompt.
