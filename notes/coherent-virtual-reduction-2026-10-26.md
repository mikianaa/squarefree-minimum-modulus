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

## Not every passing target is equally easy -- margin matters, with a clean trend

Four omega(m)=3 combinations clear ratio >= 22 at k=8->9: m=4199 (ratio
22.94, margin 0.94), 3553 (22.63, margin 0.63), 2717 (22.23, margin 0.23),
and 2431 (22.10, margin 0.10). All four reduced models are essentially the
same size (55-92 thousand compressed variables). Testing all four gives a
clean monotonic trend:

| m | ratio | margin above 22 | \|R\| | status | time to resolve |
| --- | --- | --- | --- | --- | --- |
| 4199 | 22.94 | 0.94 | 2310 | **SAT** | 22.1s |
| 3553 | 22.63 | 0.63 | 2730 | **SAT** | 238.2s |
| 2717 | 22.23 | 0.23 | 3570 | UNKNOWN | (400s, not resolved) |
| 2431 | 22.10 | 0.10 | 3990 | UNKNOWN | (~15 min total, not resolved) |

The capacity ratio is a necessary, not sufficient, condition, and clearing
it by more headroom makes the actual constraint-satisfaction problem
reliably easier to decide, not just "sometimes": the ordering by margin
exactly matches the ordering by solve difficulty across all four points
tested, with a roughly 10x jump in solve time between the two SAT cases
and no resolution at all (within several minutes) for the two closer
calls. Note the inverse relationship with F-friendliness: the two easiest
G-side cases (4199, 3553) have the *smallest* |R| (hardest for F to reach,
per the section below), while the two hardest-to-decide G-side cases
(2717, 2431) have the *largest* |R| (easiest for F, per their general-
witness T-shortfall ratios of 22.4x and 20.0x respectively -- vs 29.3x for
3553 and 34.6x for 4199 -- both still far short but consistently better
than the easy-to-decide G-side cases). This is a genuine tension in the
overall route, not just a numerical curiosity: whichever end of this
spectrum is more promising for one side of the problem is less promising
for the other.

## A much better k=8->9 point: omega=4 can beat omega=3, down to 6.35x

The margin-matters section above only searched omega(m)=3 (the *minimum*
omega clearing the ratio at this scale), on the assumption that higher
omega always means smaller |R|. That assumption is wrong: |R| = L/m
depends on *which* primes are excluded, not just how many, and some
omega=4 choices exclude a product of primes larger than any omega=3
choice's excluded set. Searching omega in {3,4,5} together for the
largest |R| that still clears ratio 22:

| m | omega | ratio | margin | Lm = \|R\| | compressed vars | status | time | general-witness F-ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1155 (3*5*7*11) | 4 | 27.81 | 5.81 | 8,398 | 241,904 | **SAT** | 42.3s | **9.52x** |
| 770 (2*5*7*11) | 4 | 24.61 | 2.61 | 12,597 | 322,544 | **SAT** | 58.1s | **6.35x** |
| 714 (2*3*7*17) | 4 | 22.74 | 0.74 | 13,585 | 322,544 | UNKNOWN | (251s) | (5.89x, if it held) |
| 462 (2*3*7*11) | 4 | 22.05 | 0.05 | 20,995 | -- | not attempted | -- | (3.81x, if it held) |

Both confirmed points were independently re-verified on the real domain
(255 divisor/residue pairs each), exactly as for every other SAT result
tonight. **6.35x is the best confirmed F-side comparison found for the
primary open case (k=8->9) tonight** -- a large improvement over the
34.6x/29.3x figures found before this search was widened past omega=3,
and over an order of magnitude better than the first (omega=3-only)
pass. The same margin-vs-virtual-period-size pattern holds: m=714, with a
smaller margin (0.74) and a larger virtual period (13,585) than m=770,
did not resolve in comparable time. m=462, with an even larger |R| but a
razor-thin margin (0.05), was not attempted -- expected to be harder
still, consistent with the pattern throughout this note.

Trying the same trick at k=9->10 (q=29): m=1155 (=3*5*7*11, omega=4,
margin 1.06, |R|=193,154 -- which would give F-ratio 6.28x, as good as
the k=8->9 result, if it held) is **UNKNOWN after 291s** (5.8 million
compressed variables -- far more than any case tried tonight). Not
pursued further given the time budget; the confirmed best for k=9->10
remains 27.2x (m=5005, section above). The same omega=4-vs-omega=3
oversight likely applies there too (the search above only tried omega=4
options with margin high enough relative to their own |R|, not
specifically the largest-|R| omega=4 options), so there is likely a
better confirmable point at k=9->10 as well; it was not found tonight.

This changes the "next steps" framing in
`overnight-session-summary-2026-10-27.md`: the k=8->9 omega=3 search was
exhausted, but omega=4 (and untried omega=5) was not, and turned out to
be where the best point actually is. The same widened search should be
tried at k=9->10 too (only omega=4 was searched there, section below);
there may be an even better point at omega=5 or 6.

## The F-side: three heuristics tried, all far short -- and biasing toward m makes it WORSE

Does an actual family F achieve a residual anywhere near {4199 | x}
(2310 points, 0.024% of Z_9699690)? Three heuristics, in increasing order
of actual performance:

1. **Safe-residue construction** (`code/safe_residue_construction.py`,
   restricting every candidate modulus to residues that provably can't
   touch a multiple of 4199): uncovered = **4,751,373** (2,056x the
   target). Local search (`code/safe_residue_local_search.py`, 10 sweeps)
   converges after one sweep with **zero improvement**.
2. **Direct T-targeting** (`code/cover_complement_of_coherent.py`): drop
   the "safe" restriction entirely -- it is not actually needed, since
   accidentally covering some multiples of m is harmless, not forbidden;
   all that matters is covering every point that is *not* a multiple of m.
   One-pass greedy, scored only on genuine non-multiples of m: uncovered
   (non-multiples of m still missing) = **123,687** (53.6x the target), a
   38x improvement in method over (1). Local search
   (`code/cover_complement_local_search.py`) converges after 2 sweeps to
   **122,653** (53.1x), almost no further gain. Total residual size (the
   123,687-ish non-m points plus the handful of m-multiples still missing)
   is 122,767 -- cap(R)/|R| = 6.40, slightly *worse* than the untargeted
   general-purpose witness below.
3. **The best already-existing general-purpose witness**
   (`data/min-uncovered-9699690-witness.json`, |R| = 79,979, built by this
   project's ordinary local search with *no* awareness of m = 4199 at
   all): checking directly, **0 of its 79,979 uncovered points happen to
   be multiples of 4199**, so its T-shortfall is exactly 79,979 -- **34.6x
   the target**, better than either method that specifically tried to
   steer toward m.

So biasing the construction toward the coherent target, by either method
tried, makes things *worse* than simply running the best known
general-purpose covering search and accepting whatever (here: zero)
overlap with the target occurs by chance. The honest headline number for
"how far is F from the coherent target" is **34.6x** for m=4199 itself,
from a construction that was not even trying to get there, and the two
deliberate methods both did worse. (Using the best *confirmed-SAT*
coherent target overall instead of m=4199 specifically -- m=770, the
omega=4 point found later in this note -- the same general witness gives
|R|/|R_target| = 79979/12597 = **6.35x**, substantially better, for
exactly the same reason: this witness's residual also happens to have
zero multiples of 770.) That is itself informative: it suggests the
coherent target is not merely hard to *reach deliberately*, it may not sit
anywhere near where general-purpose good covers naturally land either --
though 6.35x is close enough that this specific general-purpose witness is
clearly not the last word on how close F could plausibly get.

This is still much weaker evidence than the G-side's exact result: it only
rules out these three heuristic family-building strategies, not every
possible construction, and does not rule out a genuinely different exact
construction for the F-side, analogous to the one that worked for the
historical 2310/105 duplicate case (found by direct CP-SAT at that small
scale, not a greedy heuristic). Whether an analogous exact F-side
reduction exists for k=8 is not attempted here and is the natural next
step.

## Why the G-side trick does not have an F-side analogue

Worth recording why no amount of cleverness is likely to find a G-side-style
exact reduction for the F-side question. The G-side reduction worked
because its target R = {m | x} is *sparse* (|R| = L/m, a small fraction of
L): every real divisor's contribution to covering R collapses to a single
residue class on a *much smaller* virtual domain, because R itself only
occupies a thin slice of each divisor's progression.

The F-side target T = Z_L \\ {m | x} is the *complement* -- dense, nearly
all of Z_L. Writing Z_L = Z_m x Z_n (n = L/m) by CRT, a divisor d = d_m*d_n
(d_m | m, d_n | n) with residue (a_m, a_n) covers the *cylinder*
{i ≡ a_m (mod d_m)} x {k ≡ a_n (mod d_n)} -- a 1/d_m fraction of the
i-coordinate and 1/d_n of the k-coordinate, simultaneously. Small d_m
(hence large 1/d_m) is what makes a divisor efficient at covering T broadly,
exactly as in ordinary unbiased covering; there is no sparse structure to
exploit, so the natural effective strategy is just "cover efficiently",
not "cover this specific small complement cleverly" -- consistent with the
empirical finding just above that trying to steer toward the target
backfires. This is a structural (sparse-vs-dense), not just empirical,
reason the F-side resisted the G-side's trick tonight.

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
leaves **89,929,666** uncovered -- 38,931 times the coherent target. But
as above, that specific heuristic is not the best available comparison:
the existing general-purpose nine-prime witness used throughout
`nine-prime-greedy-construction-2026-10-26.md`
(`data/local-search-nine-prime.json`, |R| = 1,213,060, built with no
awareness of m = 96577) has **0 overlap** with multiples of 96577 too, so
its T-shortfall is the full 1,213,060 -- **525x the target**. Smaller than
the safe-residue number, but still two orders of magnitude off, and still
worse (as a *ratio*) than the k=8->9 case's best comparison figure
(34.6x). So the pattern holds, and gets worse, one level up either way it
is measured: **the G-side capacity requirement is consistently achievable
(exactly) for the smallest-omega coherent target at both tested
extensions, while the F-side gets further and further from achievable as
the target gets smaller/more concentrated -- and in both cases, deliberately
steering toward the target is no better (here, much worse) than just
running the best known unaware general-purpose search.**
This is a specific, nontrivial reason the coherent-residual route looks
like a dead end that is different from (and sharper than) "not enough raw
capacity": the capacity is there; getting a real family anywhere near the
needed residual is what fails, and fails harder as primes are added.

## Pushing for a better F-side comparison at 9 -> 10: four confirmed-SAT points, down to 27x

The single m = 96577 target above was the *best-margin* omega=4 option,
but far from the most F-friendly (|R| = 2310 is tiny, giving the weak
525x comparison). There are 98 omega=4 combinations clearing ratio >= 28
at this scale; searching among the larger-|R| ones for a combination of
decent margin and a virtual period small enough to actually resolve:

| m | margin above 28 | Lm = \|R\| | compressed vars | status | time | general-witness F-ratio |
| --- | --- | --- | --- | --- | --- | --- |
| 96577 | 18.88 | 2,310 | 110,576 | SAT | 22.1s | 525x |
| 37145 | 13.97 | 6,006 | 258,032 | **SAT** | 46.1s | **202.0x** |
| 10659 | 7.99 | 20,930 | 774,128 | **SAT** | 158.1s | **58.0x** |
| 5005 | 8.22 | 44,574 | 1,658,864 | **SAT** | 271.8s | **27.2x** |
| 3003 | 4.49 | 74,290 | 2,488,304 | UNKNOWN | (311s, not resolved) | (16.3x, if it held) |
| 2145 | 2.90 | 104,006 | 3,317,744 | UNKNOWN | (246s, not resolved) | (11.7x, if it held) |

Every SAT row's real system (511 divisor/residue pairs) was re-verified
directly on all real multiples of m in Z_223092870, same as the headline
case. The pattern confirms and refines the margin-matters finding from the
k=8->9 section: it is not just margin above the ratio threshold that
matters, but margin *relative to the size of the resulting virtual
problem* -- m=2145 has a smaller margin (2.90) than every SAT case here,
and also by far the largest virtual period (104,006 vs <=44,574), and it
is the one that did not resolve. The practical frontier for this specific
exact method, tonight, sits between |R| = 44,574 (SAT, m=5005) and 74,290
(UNKNOWN after 311s, m=3003): the F-side
comparison ratio for the best *confirmed* target has improved from 525x to
**27.2x** by choosing a better coherent target, without weakening the
G-side guarantee at all (every SAT row is an exact, re-verified result,
not a heuristic). Whether an even better point exists past the 2145
frontier (with more solver time, or a smarter CP-SAT configuration) is not
resolved here.

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

# the better (unrestricted) F-side probe at k=8->9:
python code/cover_complement_of_coherent.py --m 4199 --primes 2 3 5 7 11 13 17 19 \
    --output data/cover-complement-4199.json
python code/cover_complement_local_search.py --init data/cover-complement-4199.json \
    --sweeps 10 --output data/cover-complement-4199-local-search.json
```

ortools is only in `tmp/research-venv`, not the main interpreter. The G-side
solve takes about two minutes; the F-side scripts run in a few seconds at
this smaller L.
