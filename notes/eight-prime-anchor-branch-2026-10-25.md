# An entire normalized eight-prime branch is excluded

Continuation: branch 001 has also now been independently excluded; see
`eight-prime-anchor001-2026-10-25.md`. The remaining eight-prime low-anchor
branches are 010 and 011. The 001 screening status below is historical.

## Audited theorem and limits

Let N=9699690=2*3*5*7*11*13*17*19. Consider a family with one class for
each divisor of N greater than 2, in the usual normalized prime/twice-prime
form: the class p has coordinate label 0, and the class 2p has coordinate
label 1 in parity b_p. Missing classes can first be added, and redundant
twice-prime classes can be moved, only increasing coverage.

**If b_3=b_5=b_7, the family leaves at least 6696 residues uncovered.**
After interchanging the two parity labels, it suffices to take
(b_3,b_5,b_7)=(0,0,0). All sixteen choices of b_11,b_13,b_17,b_19 are
included in the bound. Unlike the earlier fixed-core result, this theorem
covers every choice of the other classes in the 210-core in this branch.

The global conjecture is not settled. The low-anchor patterns 001, 010,
and 011 remain open in the eight-prime problem, apart from previously
excluded specific cores. Arbitrarily larger prime supports also remain
unresolved. The global prime-support lower bound is still eight.

## Why the duals can be transferred between cores

The aggregate conditional-tail proof in
`notes/conditional-tail-obstruction-2026-10-25.md` has two ingredients:

1. Pointwise inequalities depending only on the four higher primes, their
   anchor parities, and local plane/tail hit counts.
2. Upper bounds on the concentration of those hit counts on a core set.

Only the second ingredient depends on the core. Thus the audited local
dual inequalities can be reused after replacing every concentration
right-hand side by another valid upper bound. The old subset table for
the 22-point example is not assumed valid for another core.

For the 000 branch the anchor boxes have dimensions (1,3,5) and (2,4,6).
For a subset with parity sizes (m0,m1), let C(m0,m1) be the sharp
arbitrary-subset bound for the sum of the 14 core projection capacities.
Coordinate compression computes C by assigning the seven even shapes to
the two layers, sorting the induced point weights, and taking the best
sum of the m0 and m1 largest weights over all 128 assignments. The
independent audit reconstructs these weights from actual integers modulo
210, not from the NumPy producer.

Adding core shapes 1 and 2 gives the valid tail-group upper bound

    C_16(m0,m1) = C(m0,m1) + m0 + m1 + max(m0,m1).

For each threshold pair, use the maximum of C_j(m0,m1)-t0*m0-t1*m1 over
the appropriate count rectangle. Multiplying these upper bounds by the
nonpositive dual multipliers preserves the lower bound for the uncovered
count.

## Reducing larger core residuals to nine count pairs

Every 210-core residual in the 000 branch has parity counts that
componentwise dominate one of

```
(0,32), (1,30), (2,29), (3,28), (4,27),
(5,25), (6,24), (7,23), (8,22).
```

The independent auditor enumerates every nonempty ordinary-residue choice
of the odd pair classes 15,21,35 and their even doubles 30,42,70, without
label symmetry. It then permits the final classes 105 and 210 to delete
up to one point in each parity plus another point in one parity. The
resulting optimistic count frontier is exactly the list above.

If an actual residual has larger counts, select an arbitrary subset of
it with one of these exact count pairs. A cover would also cover that
subset. More generally, the original family's number of uncovered points
is at least the number uncovered on this subset. The concentration bound
C applies to every such subset of the anchor boxes.

This subset argument avoids assuming that the final dual expression is
monotone in the two counts. It also avoids the feasibility-only normalized
profile dominance used in earlier work.

All 9*16 combinations have strictly positive integer lower bounds. The
smallest is 6696, at core counts (7,23) and higher pattern 1101. These are
certificate bounds, not claimed optima.

## Independent verification

`code/audit_eight_prime_anchor_branch.py` uses only the standard library.
It checks:

* 15819570 ordinary pair-union/final-even-line combinations, giving the
  nine-pair core frontier;
* the compressed capacity table via explicit integer point weights and
  all 128 even-shape allocations;
* all local states for the sixteen reused aggregate duals, using Python
  integer arithmetic;
* the new concentration constants and all 144 final lower bounds.

The audit passed in about 20 seconds. It does not depend on the new orbit
reduction or the provisional screening described below. It also does not
depend on optimization solver status or floating-point feasibility.

Reproduce it with:

```
python code/audit_eight_prime_anchor_branch.py
```

The producer is `code/eight_prime_anchor_branch.py`; the data and report
are `data/eight-prime-anchor000-certificate.json` and
`data/audit-eight-prime-anchor000.json`. This is not proof-assistant
verification, and no novelty claim is made.

## Work toward the remaining branches

`code/core_orbit_reduction.py` gives an exact representation of the
anchor-preserving coordinate permutations: 120 permutations of the free
7-labels, followed by sorting the three free 5-label rows. The group has
720 elements. Direct coordinate-permutation tests validate the encoding;
an exhaustive comparison confirms that canonical-even and full-even pair
enumerations have identical orbit sets. Inclusion reduction then keeps
representatives whose residuals contain no relabelled smaller retained
residual.

| Low pattern | Raw canonical pair masks | Relabelling orbits | Retained inclusion representatives |
| --- | ---: | ---: | ---: |
| 000 | 61604 | 50145 | 5330 |
| 001 | 45415 | 35994 | 3174 |
| 010 | 44836 | 35400 | 2494 |
| 011 | 39174 | 29997 | 1911 |

These are residuals before 105 and 210. The orbit equality audit is in
`data/audit-core-even-orbits.json`; the reduction catalogue is
`data/core-orbit-reduction.json`. The inclusion reduction has regression
checks, but a complete independent witness audit has not yet been added.

The portable screening implementation additionally bounds H by joint
parity slice counts for each modulus, intersected with the compression
bound. It treats the last two core classes by selecting the optimistic
smaller subsets described above. Its output is in
`data/portable-tail-core-screen.json`. Those batch screening results are
research guidance, not a substitute for the independent whole-branch
audit above.

In branch 001, the original duals leave 320 representative/high-pattern
combinations unresolved. Splitting the two optimistic deletion choices
and merging identical capacity inputs leaves 167 distinct optimization
inputs, saved in `data/anchor001-unresolved-capacity-inputs.json`.

An abstract aggregate model now accepts these count/capacity inputs
directly. Eight newly generated duals reduce the unresolved input count
from 167 to 41. They are saved with the flag `independent_audit: false` in
`data/anchor001-refined-dual-pool.json`; their independent arithmetic audit
and completeness audit are still required. No exclusion of the whole
001 branch is claimed.

The 47-test regression suite passed. A subsequently added abstract-model
regression and the three adjacent conditional-tail tests also passed.

Next work should continue the remaining 001 capacity inputs, then audit
both the new duals and coverage of the original core cases. The 000 branch
does not need to be searched again.
