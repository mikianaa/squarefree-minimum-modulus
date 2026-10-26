# The static two-family model (direction 2) fails at k=5, on the actual optimal witness

Written 26 October 2026 (Claude Code, Sonnet 5). Direction 2 of
`notes/direction-analysis-2026-10-25.md` is the general (non-aligned) model
that `notes/aligned-construction-test-2026-10-26.md` only tested a special
case of (G-residues forced to 0). This note tests the general model
directly, exactly, at a small scale where it is fully decidable, using
*actual* optimal data rather than counting.

## The model and the test

For primes dividing L with next prime q: a hypothetical cover with lcm Lq
needs a q-free family F over Z_L (distinct moduli >= 3) with residual R,
plus a single family G -- one residue a_d per divisor d > 1 of L (any
residue, not just 0) -- such that every point of R is covered by at least
q-1 of the resulting progressions mod d. (WLOG G may use every divisor of L,
since more capacity can only help; F's own moduli d and G's actual moduli
qd never collide since q does not divide L.)

This is a small, exact decision problem once R is fixed: one boolean
"which residue" choice per divisor (sum to 1), and one linear constraint per
point of R (sum over divisors of the indicator that its chosen residue
matches that point's residue mod d, >= q-1). `code/static_two_family_test.py`
builds and solves this with CP-SAT (OR-Tools, from `tmp/research-venv`,
where it is installed; not in the main interpreter).

## Result at L=2310, q=13 (k=5, five primes)

Used the actual CP-SAT-certified optimal witness
`data/min-uncovered-lcm2310-feasible121.json` (u(2310)=121, the true
minimum, not a lower bound): F uses all 30 divisors >= 3 of 2310, leaving
exactly 121 explicit uncovered residues. G has 31 divisors (2 through 2310)
available and needs every one of those 121 points covered q-1=12 times.

**INFEASIBLE** (solved in well under the 120s time limit; no partial/unknown
status). The solver's own witness format gives no certificate to replay
independently of trusting the solver's INFEASIBLE claim, but the companion
FEASIBLE direction (see below) is verified independently inside the script
by recomputing every point's multiplicity from the saved assignment without
trusting the solver -- this is the standard asymmetry of SAT/UNSAT results
(a SAT witness is independently checkable; an UNSAT proof, without the proof
log, is not, here).

This reproves -- via a genuinely different mechanism than
`code/slice_count.py`'s pigeonhole counting -- that this specific optimal
family's residual cannot be slice-extended to 30030, which was already
known (the direction-analysis note's own table: counting alone already
settles k=5). The point of running it anyway was to measure how much slack
the fuller (non-aligned) model buys over the aligned one, on real data, not
just bounds. Here: none. Both say impossible, on this witness.

Since |R| = 121 here is the true minimum for this period, and smaller R
only makes G's job easier (same total G-capacity spread over fewer required
points), the other four CP-SAT-optimal anchor patterns in
`data/min-uncovered-2310-five-patterns-optimal.json` (uncovered 123, 124,
134, 141) are expected to fail at least as badly and were not separately
run.

## Why this does not scale to k=8 (9699690, q=23) directly

The same exact-CP-SAT encoding needs sum_{d|L, d>1} d boolean variables (one
per divisor, per its own residues) -- about 35 million at L=9699690 (sigma(L)
minus 1 and 2), far beyond what CP-SAT can even build in memory, let alone
solve, in one evening. The k=5 case above only worked because L=2310 is
small (sigma(2310)-3 = 6,909 variables). Resolving the open k=8 case this
way needs either a genuinely different (polyhedral/profile-compressed)
encoding -- exactly direction 5 of the direction-analysis note -- or a
heuristic witness-construction attempt instead of an exact solver, as in
`notes/nine-prime-greedy-construction-2026-10-26.md`'s direct-cover
experiment (which took the complementary approach of attacking the full
period directly rather than this slice/G-family decomposition, since the
decomposition's own bookkeeping does not fit in an exact solver at this
scale).

## Scope

* Not a new certified bound; this reproves an already-known fact about
  2310 -> 30030 by a different method, as a calibration of the model before
  attempting (and failing, for scale reasons) to run it at k=8.
* Confirms the fuller, non-aligned static two-family model buys nothing over
  the aligned special case on the one piece of real data tested here.
* No change to `primes_needed` (9) or the certified bounds.

## Reproduction

```
tmp/research-venv/Scripts/python.exe code/static_two_family_test.py \
    --source data/min-uncovered-lcm2310-feasible121.json --q 13 --output data/static-two-family-2310-13.json
```

ortools is not installed in the main interpreter; use the project's own
`tmp/research-venv`. Runs in a few seconds.
