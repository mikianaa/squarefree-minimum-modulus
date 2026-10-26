# A direct greedy construction attempt at the first open primorial (nine primes)

Written 26 October 2026 (Claude Code, Sonnet 5), direction 4 of
`notes/direction-analysis-2026-10-25.md`: "take the possibility M_sf >= 3
seriously and search constructively... with heuristics, not exact solvers."
This is a heuristic experiment, not a proof in either direction. It adds no
certified bound.

## What was tried

L = 223092870 = 2*3*5*7*11*13*17*19*23 (nine primes) is the first primorial
not excluded by the eight-prime obstruction (every counterexample omitting
modulus 2 needs >= 9 distinct prime factors). A direct cover of Z_L using
distinct divisors >= 3 of L, if one existed, would be an actual squarefree
covering system with least modulus 3 -- i.e. M_sf >= 3, contradicting the
usual conjecture M_sf = 2. No such cover is claimed or expected here; this
is a cheap sanity probe of how close a naive constructive heuristic gets.

`code/greedy_cover_nine_primes.py`: a one-pass greedy. Process the 510
divisors >= 3 of L in a fixed order; for each, pick the residue class that
covers the most still-uncovered points (via one `numpy` reshape + column-sum
per divisor -- O(L) time independent of d, so no per-residue lookup table is
built, unlike `code/improve_cover_witness.py`, which is not memory-tractable
at this L: summing divisors of L gives about 8.4*10^8, far too many explicit
residue entries). No backtracking or subsequent local search.

Three orders tested, same divisor set, each ~63 seconds on this machine:

| order | final uncovered | density |
| --- | --- | --- |
| ascending (prime-count, then value) | 1,325,141 | 0.594% |
| random (seed 7) | 5,237,242 | 2.35% |
| descending | 9,687,783 | 4.34% |

Ascending order (small, high-value divisors first) dominates, as expected
from standard covering-system heuristics; it is the only one carried further.
Self-check: the saved `(modulus, residue)` system was independently rebuilt
into a fresh coverage array and reproduced the same uncovered count exactly
(`data/greedy-nine-prime-smalld.json`).

## Result

The best heuristic witness leaves **1,325,141** of 223,092,870 points
uncovered (density 0.594%). This is an upper witness only (`u(223092870) <=
1325141`, almost certainly far from tight -- the one-pass greedy has no
chance to fix early bad commitments, unlike the tuned local search already
used at L = 9699690). For calibration, the SAME crude method at L = 9699690
(eight primes) gives 93,123 uncovered (0.96%), slightly worse than the
already-published tuned witness there (79,979, 0.82%). So this new number
likely overstates the true minimum at nine primes by a similar margin, but
the order of magnitude and density trend are informative.

Despite being cruder, the nine-prime density (0.594%) is *lower* than the
tuned eight-prime witness's density (0.82%) -- consistent with, and
continuing, the steady decrease already recorded in section 1 of the
direction-analysis note (10.5% -> 5.2% -> 1.7-2.6% -> 0.1-1.5% -> <=0.82% for
k = 4..8). This on its own says nothing new about feasibility; smaller
residual density was already the expected, unsurprising trend.

## The informative part: cap(R)/|R| continues to fall behind

Direction 4 asks specifically to track cap(R)/|R| (the slice-capacity ratio
of section 2 of the direction-analysis note) for an *actual* constructed
family, not just its residual size. For the best witness here (R = the
1,325,141 uncovered points), cap(R) := sum over divisors d > 1 of L of the
largest single-residue intersection with R, computed exactly by the same
reshape-and-max trick applied to the fixed R (511 divisors, ~63 seconds):

    cap(R) = 10,307,316,  cap(R)/|R| = 7.778

The next prime after the nine is q = 29, requiring every point of R to have
>= q-1 = 28 pairwise-coherent divisors through it for a slice extension to
L*q = 6,469,693,230 (ten primes) to even be worth attempting. Extending the
note's own table by this one row:

| L (primes) | q | \|R\| | cap(R)/\|R\| | required (q-1) | shortfall | ratio to requirement |
| --- | --- | --- | --- | --- | --- | --- |
| 2310 (5) | 13 | 121 | 5.17 | 12 | 6.8 | 0.43 |
| 30030 (6) | 17 | 792 | 6.04 | 16 | 10.0 | 0.38 |
| 510510 (7) | 19 | 7436 | 7.15 | 18 | 10.9 | 0.40 |
| 9699690 (8) | 23 | 79979 | 7.87 | 22 | 14.1 | 0.36 |
| **223092870 (9)** | **29** | **1325141** | **7.78** | **28** | **20.2** | **0.28** |

cap(R)/|R| has essentially *stopped growing* between eight and nine primes
(7.87 -> 7.78, actually slightly down, likely because this witness is
cruder than the tuned 8-prime one rather than a real reversal), while the
requirement keeps climbing with the next prime (22 -> 28). The achieved
fraction of the requirement drops to 0.28 -- almost exactly "a third of the
requirement" at nine primes, which direction 4 flagged in advance as "the
strongest available evidence for the conjecture" if it held at nine and ten
primes. It held at nine. This is one data point from one crude heuristic,
not a proof; a properly tuned witness (local search, as in
`improve_cover_witness.py`, scaled up) could move cap(R)/|R| up somewhat, as
the gap between the crude 8-prime number (this method, 0.96% density) and
the tuned one (0.82%) shows is possible -- but the gap would have to widen
by an order of magnitude to approach 28, not shrink by the ~15% that tuning
bought at eight primes.

## Follow-up: tuned local search barely moves the ratio

The note above flagged a tuned local search as the natural next step and
predicted the gap "would have to widen by an order of magnitude... not
shrink by the ~15% that tuning bought at eight primes" to threaten the
conclusion. `code/local_search_nine_primes.py` implements this properly at
the same L = 223092870 scale: a coordinate-descent (Gauss-Seidel) sweep over
all 510 moduli, each revisited and allowed to pick a better residue given
every other modulus's current choice, using an int16 coverage-count array
(not the one-pass greedy's one-shot commitment) so each sweep can only keep
or reduce the uncovered count.

Starting from the one-pass witness above, five sweeps (~93s each) converge
to **1,213,060** uncovered -- an 8.5% reduction from 1,325,141, consistent
with the 15% tuning bought at eight primes. Recomputing the capacity ratio
on this improved residual:

    cap(R) = 9,556,220,  cap(R)/|R| = 7.878,  fraction of requirement (28) = 0.2814

Essentially unchanged from the one-pass witness (7.778 -> 7.878, 0.278 ->
0.281). Tuning shrank the residual by 8.5% but barely moved the capacity
ratio at all. This is the strongest evidence this note has for its own
conclusion: the 28% figure is not an artifact of the crude one-pass method,
it survives an actual (if still modest) local-search refinement at the same
scale.

Reproduction:

```
python code/local_search_nine_primes.py --init data/greedy-nine-prime-smalld.json --sweeps 6 --output data/local-search-nine-prime.json
```

~8 minutes total (5 sweeps to convergence), numpy required.

## Scope and honesty about limits

* This is a heuristic upper witness, not a lower bound and not a claim that
  no better construction exists. It does not touch `primes_needed` (still 9)
  or the certified bracket on u(223092870) (there is none yet; this note
  supplies the first recorded one: u(223092870) <= 1325141).
* The one-pass greedy is deliberately cheap (no backtracking); a tuned local
  search at this scale is architecturally possible (the O(L) reshape trick
  used here could replace the per-residue dictionary in
  `improve_cover_witness.py` to make it memory-tractable at L =
  223092870) but was not attempted tonight.
* The ten-primes extension (L*29 = 6,469,693,230) was not attempted directly;
  only the capacity ratio of a nine-prime witness against the q=29
  requirement was computed, which is direction 4's own suggested cheap
  proxy, not a full search.

## Reproduction

```
python code/greedy_cover_nine_primes.py --order smalld --output data/greedy-nine-prime-smalld.json
```

About 63 seconds, numpy required, ~223 MB peak for the coverage array. The
cap(R)/|R| figure was computed by a short follow-up script reusing the saved
system (not committed separately; see `data/greedy-nine-prime-smalld-cap-check.json`
for the numeric result).
