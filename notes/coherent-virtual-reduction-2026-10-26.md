# The G-side of the coherent-residual route is exactly solvable at k=8 -> 9 -- and SAT

Written 27 October 2026 (Claude Code, Sonnet 5). This is the most significant
single finding of tonight's session. **It is not a new covering system and
does not change any certified bound.** Read the "Scope" section before
drawing any conclusion from the headline result.

## Setting the scale right

L = 9699690 = 2*3*5*7*11*13*17*19 (eight primes), q = 23 (the ninth prime,
correctly *excluded* from L -- the open question is whether a cover with
lcm L*q = 223092870 exists with least modulus 3). This is the first
untested primorial. (An earlier draft of this investigation mistakenly ran
some checks with q=23 against L=223092870, which already includes 23 --
a nonsense combination, since q must be a prime *not* already in L.
Those results are not reported.)

The static two-family model (direction 2 of
`notes/direction-analysis-2026-10-25.md`, already tested exactly at k=5,6,7
in `notes/static-two-family-test-2026-10-26.md`) splits into two questions
for a candidate residual R:

1. **F-side**: does some actual distinct-moduli family over Z_L leave
   residual R (or a subset of it)?
2. **G-side**: can the divisors of L (d > 1, each usable independently any
   number of times across *different* real divisor values) assign residues
   so that every point of R is hit at least q-1 = 22 times?

Both are needed for R to yield a slice-extended cover. Tonight's other
notes (`aligned-construction-test`, `static-two-family-test`,
`coherence-capacity-tradeoff`) mostly bounded or tested the F-side, or
tested both together on *spread* (greedy) residuals, where the G-side
always failed badly (ratio ~7-8 against a requirement of 22-28). This note
isolates the G-side on the best possible *coherent* target and solves it
exactly.

## The coherent target and its exact reduction

For divisor m of L, R = {x in [0, L) : m | x}, |R| = L/m. Scanning every
divisor of L = 9699690 by omega(m) (number of prime factors), the smallest
omega(m) whose best capacity ratio cap(R)/|R| clears the requirement 22 is
**omega(m) = 3**, at **m = 4199 = 13*17*19**, with ratio 22.94 and
**|R| = 2310** (the complementary primes 2,3,5,7,11 multiply to exactly
2310). Every omega(m) <= 2 option has ratio < 12 (omega=1) or < 12.89
(omega=2) -- all well short of 22. So m = 4199 is the *most achievable*
(largest |R|, i.e. least restrictive) coherent target that even has a
chance on the G-side.

The G-side question for this R -- can the ~254 divisors of L assign
residues so every one of the 2310 multiples of 4199 gets >= 22 matching
divisors? -- naively needs about sum_{d|L,d>1} d ~ 35 million boolean
variables, intractable directly (the same wall found for the k=7 case in
`static-two-family-test-2026-10-26.md`).

**Key fact, proved and checked**: writing x = m*k for k in [0, L/m), a real
divisor d (with g = gcd(d,m)) can only match points of R via a residue a
that is a multiple of g (by the CRT solvability criterion for simultaneous
congruences), and whenever it does, the matched points correspond to a
*single* residue class of k modulo d' := d/g (since m/g is invertible
mod d', as gcd(d', m) = 1 by construction). So every real divisor acts
exactly like a smaller "virtual" divisor d' of L/m on the k-domain.
Consequences:

* Divisors of m itself (d' = 1) match *every* k unconditionally, for free.
  There are 7 of them (d > 1, all odd since m is odd: 13,17,19,221,247,323,4199).
* Every other real divisor d of L maps to some d' > 1 dividing L/m = 2310,
  and distinct real divisors of L mapping to the *same* d' act as
  independent "copies" of that virtual divisor, since they are different
  integers and each gets its own residue choice. The number of real
  divisors mapping to a given d' > 1 is exactly 2^omega(m) = 8 (one per
  divisor of m used to build back up to a real d = d' * s, s | m),
  independent of which d'.

So the G-side question reduces *exactly* to: on Z_2310 (31 divisors d' > 1,
8 independent copies each, 55,288 boolean variables after the same
occurring-residue compression used for k=7), does every k get at least
22 - 7 = 15 hits? This is now the same scale as the original k=5 problem,
not k=8, and CP-SAT decides it directly.

`code/coherent_virtual_reduction.py` implements this reduction and the
reduced CP-SAT solve, and -- critically -- translates any SAT solution's
virtual assignment back into an explicit list of 255 real (divisor,
residue) pairs on Z_9699690 and re-verifies the original claim (every real
multiple of 4199 gets >= 22 matching real divisors) using ordinary integer
arithmetic on the real domain, not trusting the reduction's own algebra.
Validated first at small scale (L = 2310, m = 11 and m = 6, matching the
direct unreduced CP-SAT model in `static_two_family_test.py` exactly on
both INFEASIBLE verdicts; a trivial m = 1155 case where the free baseline
alone suffices confirmed SAT handling).

## Result

**SAT**, in 118.5s (55,288 compressed boolean variables), with the
translated-back real system of 255 (divisor, residue) pairs independently
re-verified by direct integer arithmetic on all 2310 real points of
Z_9699690 that are multiples of 4199: `data/coherent-virtual-reduction-4199-23.json`.

So: **if** some distinct-moduli family F over Z_9699690 (moduli >= 3)
could leave its uncovered residual contained in {x : 4199 | x}, **then**
the G-side (slice-coverage) requirement for extending to a cover with lcm
223092870 is achievable, by an explicit, verified construction. This is a
new exact result (not previously computed in this repository) and is the
tightest G-side witness found so far at this scale by a wide margin (ratio
22.94 against a requirement of 22, versus ~7-8 for every spread witness
tested elsewhere tonight).

## The F-side: heuristically, nowhere close

Does an actual family F achieve a residual anywhere near {4199 | x}
(2310 points, 0.024% of Z_9699690)? `code/safe_residue_construction.py`
(restricting every candidate modulus to residues that provably can't touch
a multiple of 4199, per `coherence-capacity-tradeoff-2026-10-26.md`) gives
an actual family whose uncovered set is **4,751,373** points -- 2,056 times
the coherent target. A further local-search refinement
(`code/safe_residue_local_search.py`, 10 sweeps) converges after a single
sweep with **zero improvement**: the safe-residue greedy is already stuck
at a strict local optimum, and local search cannot move it at all here
(contrast with the m=390, k=9 case, where local search bought a modest 5%).

This is much weaker evidence than the G-side's exact result: it only rules
out one specific (and, by construction, quite conservative) heuristic
family-building strategy, the same way the earlier coherence-capacity note
found for a different, larger coherent target at k=9. It does not attempt,
and does not rule out, a genuinely different exact construction for the
F-side, analogous to the one that worked for the historical 2310/105
duplicate case (found by direct CP-SAT at that small scale, not a greedy
heuristic). Whether an analogous exact F-side reduction exists for k=8 is
not attempted here and is the natural next step.

## Same pattern one level up: 9 -> 10 (q = 29)

Repeating the whole exercise at L = 223092870 (nine primes), q = 29: the
smallest omega(m) whose capacity ratio clears 28 is **omega(m) = 4**, at
**m = 96577 = 13*17*19*23** (ratio 46.88, |R| = 2310 -- the same size by
coincidence, since the complementary primes are again 2,3,5,7,11). The
reduced G-side model (free = 15, need_extra = 13, 31 virtual divisors of
2310, 16 copies each, 110,576 compressed variables) is **SAT in 22.1s**,
with the translated-back 511-divisor real system re-verified directly on
all 2310 real multiples of 96577 in Z_223092870.

F-side, same story but worse: the safe-residue construction for m = 96577
leaves **89,929,666** uncovered -- 38,931 times the coherent target, an
order of magnitude further off than the m=4199 case at k=8->9 (2,056x).
So the pattern holds, and gets worse, one level up: **the G-side capacity
requirement is consistently achievable (exactly) for the smallest-omega
coherent target at both tested extensions, while the F-side gets further
and further from achievable as the target gets smaller/more concentrated.**
This is a specific, nontrivial reason the coherent-residual route looks
like a dead end that is different from (and sharper than) "not enough raw
capacity": the capacity is there; getting a real family anywhere near the
needed residual is what fails, and fails harder as primes are added.

## Scope (read this before citing the headline result)

* This does **not** show a 9-prime covering system exists, and does not
  change `primes_needed` (9), the certified bracket on u(9699690), or any
  other certified bound in `PROGRESS.md`.
* It solves, exactly, one half (the G-side) of one specific necessary
  condition (coherent residual = {4199 | x}) for one specific route to
  such a cover. The other half (F actually achieving that residual) is
  only tested by a conservative heuristic here, and that heuristic fails
  badly -- which is suggestive but not a proof of F-side infeasibility.
* The exact reduction technique (reducing a sparse coherent-target
  G-coverage problem to a tiny virtual period) is new to this repository
  and may generalize usefully to other coherent targets or other k, which
  is the main reason this result is worth recording even though it settles
  nothing about the global conjecture.

## Reproduction

```
tmp/research-venv/Scripts/python.exe code/coherent_virtual_reduction.py \
    --m 4199 --q 23 --primes 2 3 5 7 11 13 17 19 --time 120 \
    --output data/coherent-virtual-reduction-4199-23.json

python code/safe_residue_construction.py --m 4199 --primes 2 3 5 7 11 13 17 19 \
    --output data/safe-residue-4199-at-9699690.json
python code/safe_residue_local_search.py --init data/safe-residue-4199-at-9699690.json \
    --sweeps 10 --output data/safe-residue-4199-local-search.json

# the 9 -> 10 extension (q=29), same pattern:
tmp/research-venv/Scripts/python.exe code/coherent_virtual_reduction.py \
    --m 96577 --q 29 --primes 2 3 5 7 11 13 17 19 23 --time 180 \
    --output data/coherent-virtual-reduction-96577-29.json
python code/safe_residue_construction.py --m 96577 --output data/safe-residue-96577.json
```

ortools is only in `tmp/research-venv`, not the main interpreter. The G-side
solve takes about two minutes; the F-side scripts run in a few seconds at
this smaller L.
