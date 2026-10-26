# Overnight session summary (26-27 October 2026, 21:50-06:00 JST)

Written by Claude Code (Sonnet 5) at the end of an overnight self-paced
session investigating untested directions from
`notes/direction-analysis-2026-10-25.md`. This note is a map of what was
done, in order, with pointers to the detailed notes -- read it first, then
go to whichever detailed note matters for your purpose. No certified bound
changed tonight. All work is committed; every numeric claim below was
either independently re-derived (a second method in the same session) or
blind-audited by a fresh `auditor` subagent run from an isolated excerpt
containing no results, per `meta/README.md`'s comparison protocol (see
`meta/audit-comparison.jsonl`, entries 4-9, and `meta/agent-progress.jsonl`,
entries 26-39).

## What was tested, and the headline finding of each

1. **Aligned construction test** (direction 3) --
   `notes/aligned-construction-test-2026-10-26.md`. Pure integer counting
   reproves the already-known k<=6 obstructions in one line, but is
   undecided at the open k=8 primorial (R_allowed=30,682 vs the weakest
   established lower bound of 2). The naive residue-0 realization misses
   the entire coprime-to-L class (1,658,880 points), showing the
   restriction is nowhere near tight.

2. **Static two-family model, exact test** (direction 2) --
   `notes/static-two-family-test-2026-10-26.md`. CP-SAT proves INFEASIBLE
   at k=5 (the true optimal residual) and k=6 (a heuristic witness). A
   16.7x variable-count compression (one boolean per *occurring*
   (divisor, residue) pair) still leaves k=7 UNKNOWN after ~14.5 minutes
   of solver time -- the scaling wall there is about constraint structure,
   not variable count. k=8 was not attempted directly for that reason.

3. **Direct constructive search, take M_sf>=3 seriously** (direction 4) --
   `notes/nine-prime-greedy-construction-2026-10-26.md`. A one-pass greedy
   plus local search at the first open primorial (223,092,870, nine
   primes) finds no construction, as expected; its slice-capacity ratio
   cap(R)/|R| = 7.88 is only 28% of the 28-fold coverage a tenth-prime
   extension would need -- the lowest fraction yet in the running table
   from `direction-analysis-2026-10-25.md` section 2 (43%, 38%, 40%, 36%
   at five through eight primes).

4. **Spread-vs-coherent capacity tradeoff** --
   `notes/coherence-capacity-tradeoff-2026-10-26.md`. A closed-form
   capacity formula for a single coherent residual class shows it clears
   the capacity requirement at a modest concentration (omega(m)=4), unlike
   any spread (greedy) residual. Realizing that concentration with an
   actual family is a different question; a first attempt (leaning toward
   it via "safe" residues) and local search on top both fail to get
   close.

5. **The night's main result: an exact G-side reduction for the coherent
   route** -- `notes/coherent-virtual-reduction-2026-10-26.md`. A new CRT
   reduction turns the (otherwise ~35-million-variable) question "can the
   divisors of L cover a coherent residual q-1 times" into a tiny
   equivalent problem on Z_{L/m}. At the correct scale for the open 8->9
   extension (L=9,699,690, q=23), the best such target (m=4199) reduces to
   a problem CP-SAT decides **SAT** in under two minutes, independently
   re-verified by direct arithmetic on all 2,310 real points and
   confirmed by a blind auditor that re-derived the reduction from the math
   alone. The same pattern holds one level up (9->10, q=29). **This does
   not show a cover exists.** The companion question -- can an actual
   family F reach a residual anywhere near that target -- was tested with
   three heuristics; the most honest comparison (the best *existing*
   general-purpose witness, built with no awareness of the target at all)
   falls short by 34.6x (k=8->9) and, initially, 525x (k=9->10) -- and
   *trying* to steer a construction toward the target makes it worse, not
   better. A follow-up completed the margin-vs-solvability picture across
   all four omega=3 targets at k=8->9 that clear the capacity threshold:
   margin above the threshold, not just clearing it, predicts solve
   difficulty with a clean monotonic trend (0.94 margin -> SAT in 22s;
   0.63 -> SAT in 238s; 0.23 and 0.10 -> UNKNOWN after several minutes
   each) -- and, intriguingly, an inverse relationship with F-friendliness
   (the easiest G-side cases have the smallest, least F-achievable
   residuals). A second follow-up then improved the k=9->10 F-side
   comparison substantially by deliberately searching for larger,
   still-solvable coherent targets: four confirmed-SAT points (each
   independently re-verified) take the comparison from 525x down to
   27.2x (|R| = 44,574), with the practical frontier for that search
   sitting strictly between |R| = 44,574 (SAT) and 74,290 (UNKNOWN after
   311s). **A third, final follow-up corrected a bad assumption** in all
   of the above (that the omega clearing the capacity ratio with the
   fewest prime factors always gives the largest, most F-friendly
   residual -- false, since |R| = L/m depends on *which* primes are
   excluded, not just how many) and found a substantially better point at
   k=8->9 by searching omega=4 instead of only omega=3: m=770
   (|R| = 12,597, SAT in 58s, independently re-verified) brings the
   primary case's F-side comparison down to **6.35x** -- close enough
   that the general-purpose witness used for the comparison is clearly
   not the last word on how near F could plausibly get.
   
   The reduction technique itself (sparse coherent target -> tiny virtual
   period) is the most likely-to-be-reusable output of the night; the
   theoretical reason it has no F-side analogue (sparse vs. dense target)
   is recorded in the same note.

6. **Upper-side regression check** --
   recorded in `meta/agent-progress.jsonl` entry 31 /
   `meta/audit-comparison.jsonl` entry 7. Re-ran both verifier scripts for
   the M_sf<=7 certificate against the saved schedule; every field matches
   the saved certificates exactly except wall-clock timing. Confirms the
   headline result reproduces as documented; this was a check, not new
   research.

7. **A lead checked and left alone, correctly** -- whether the "union
   accounting" exact-identity lead (`stochastic-orders-and-union-
   accounting-2026-10-25.md`, written before M_sf<=7 was achieved by
   breakpoint optimization) could now close the gap to M_sf<=6: the
   existing note's own numbers (an extrapolated gain of 0.01-0.05,
   bracketing the *old* C0=8 deficit of 0.018 from both sides) already
   show this is exactly as inconclusive for the *new*, larger C0=6 deficit
   (~0.030) as it was before. No new claim was manufactured here; this is
   recorded so a future session does not re-spend time re-deriving the
   same inconclusive bracket.

## What was deliberately NOT touched, and why

* `code/extremal_slice_capacity.py`, `code/paired_slice_capacity.py`, and
  the slice-capacity maximiser direction (direction-analysis item 1):
  Codex was actively iterating on this all night (uncommitted WIP visible
  throughout the session). Tonight's work on items 2-5 above is
  complementary to, not overlapping with, that search.
* `notes/core30-overlap-refinement-2026-10-26.md`,
  `notes/core210-joint-overlap-2026-10-26.md`, and the upper-side "modulus
  10 is forced" geometry work: read for context, not edited or extended --
  this looked like (or recently was) another agent's active working area,
  and the two verifier/audit files most directly tied to the certified
  M_sf<=7 result (`code/audit_verify.py`, `code/few_progression_residual.py`)
  were shown as modified in git status all night.
* `PROGRESS.md` and `README.md`: both were already massively uncommitted
  (nearly the entire 25-26 Oct body of work, not just tonight's) when this
  session started, from other agents' concurrent sessions. Edits were made
  to `PROGRESS.md` for this session's own visibility while working, but
  were not committed -- isolating this session's specific hunks from the
  pre-existing large uncommitted diff was not safely possible without
  risking someone else's in-progress text. Whoever next commits
  `PROGRESS.md` should fold in the entries this session added (search for
  "aligned construction test", "static two-family", "coherent", "nine-
  prime greedy" in the uncommitted working tree, or diff against commit
  `9cb4da5` which is the last commit of this session).

## Suggested next steps, roughly in order of promise

1. **Push the coherent-virtual-reduction G-side result further**: the
   widened (multi-omega) search found 6.35x at k=8->9 (m=770, omega=4) and
   27.2x at k=9->10 (m=5005, omega=4), each after exploring only a
   fraction of the available combinations (omega=5,6,... were barely
   touched, and the boundary between SAT and UNKNOWN was only pinned down
   to within a factor of ~1.2x of |R| in each case, not tightly). Likely
   more room left on both sides. Also worth trying a cleverer non-single-
   class (partially coherent) target -- which is exactly where the F-side
   and Codex's maximiser search might usefully meet this technique.
2. **An exact F-side construction or obstruction**, analogous to the
   historical 2310/105 duplicate-modulus case (found by direct CP-SAT at
   that small scale) -- not attempted at k=8 scale tonight; the sparse-vs-
   dense structural argument in `coherent-virtual-reduction-2026-10-26.md`
   explains why it needs a genuinely different idea, not just more
   compute.
3. Direction 5 of the direction-analysis note (a polyhedral
   characterization of valid inequalities, rather than one-witness-at-a-
   time cuts) remains completely untested.
4. If the slice-capacity maximiser search (direction 1, Codex's) finds a
   family with a notably better cap(R)/|R| than tonight's spread witnesses
   at nine primes, re-run tonight's G-side reduction against *that*
   residual's structure, not just the single-class coherent targets tried
   here.

## Housekeeping

* All code, notes, and data from tonight are committed (20 commits,
  `eb7afcb..9cb4da5`). `meta/agent-progress.jsonl` (entries 26-39) and
  `meta/audit-comparison.jsonl` (entries 4-9) are up to date; regenerate
  the rendered `.md` views with `python meta/meta_stats.py` and
  `python meta/compare_audits.py` if anyone adds more entries before
  reading them.
* ortools (CP-SAT) is not installed in the main Python interpreter; use
  `tmp/research-venv/Scripts/python.exe` for anything that imports it.
* One methodological note worth keeping in mind for future delegation: the
  `auditor` subagent's Read tool returns a whole file even when told to
  read "only part of it" -- a blind audit needs a physically separate
  excerpt file with no results in it, not just an instruction. See
  `auditor-blind-read-excerpt` in this session's memory file if available,
  or just reproduce the pattern used in `tmp/audit/*/methodology-excerpt.md`
  throughout tonight.
