# Seven-prime investigation and the exact minimum residual at 2310

Subsequent result: the seven-prime case is now excluded by the audited plane-load
argument in [the seven-prime obstruction note](seven-prime-obstruction-2026-10-25.md).
The unsuccessful relaxations below are retained as research history.

The global minimum-modulus-2 problem remains open. The preceding solver-free
six-prime obstruction is unchanged. This investigation strengthens a finite
residual bound and tests whether that improvement can exclude 510510.

## Minimum uncovered count at 2310

Write u(2310) for the minimum number of uncovered residues when selecting one
class for each divisor d>=3 of 2310. The exact parity-box CP-SAT model gives
u(2310)=121. The complete case coverage is:

| Anchor pattern b_3,b_5,b_7,b_11 | Lower bound | Evidence |
| --- | ---: | --- |
| 0000 | 156 | Audited integer capacity enumeration |
| 0001 | 138 | Audited integer capacity enumeration |
| 0010 | 126 | Audited two-parity frontier enumeration |
| 0011 | 141, attained | CP-SAT OPTIMAL, 52.2 s |
| 0100 | 134, attained | CP-SAT OPTIMAL, 43.1 s |
| 0101 | 124, attained | CP-SAT OPTIMAL, 32.3 s |
| 0110 | 123, attained | CP-SAT OPTIMAL, 57.1 s |
| 0111 | 121, attained | CP-SAT OPTIMAL, 43.1 s |

All five solver witnesses were checked directly on all 2310 integer residues.
The lower-bound conclusions for these five rows rely on OR-Tools' OPTIMAL status
and the normalization/encoding argument; no independently checked proof log for
the 121 bound has been obtained. This is a solver-certified result, unlike the
fully integer-enumerated six-prime obstruction.

Saved results: `data/min-uncovered-2310-five-patterns-optimal.json`. The other
three cases are in `data/audit-parity-capacity-30030.json`. Reproduce a solver
case with, for example:

    python code/parity_cover.py 2310 --pattern 0111 --minimize --residual-lb 90 --time 60 --workers 8 --output out.jsonl

The input lower bounds used were already certified by the integer enumeration.
They were 109, 103, 90, 98, and 90 in the five solver cases listed above.

## What the stronger bound buys

Intersecting the certified two-parity frontier at pattern 0111 with total
residual >=121 and propagating through 13 gives a lower bound 282 at pattern
01111 modulo 30030, compared with 74 without the new input. Its 261 minimal
count pairs all still pass the capacity tests for adjoining 17, for either
parity of modulus 34. This is inconclusive, not evidence of a covering.

Reproduce this refinement with:

    python code/propagate_solver_bound.py --pattern 0111 --new-parity 1 --output out.json

For comparison, the simple summed-capacity test would exclude both extensions
at this pattern if a total residual lower bound of 786 could be established.
That number is a sufficient threshold for this specific pattern and test,
not a universal necessary threshold or a claim about the true minimum.

## Linked projection models

The new `core_projection.py` keeps the actual CRT positions in the primes beyond
a 210 core. Each class is allowed an optimistic capacity in the core residual.
This retains interactions between the high-prime slices that the previous
step-by-step count recurrence discards.

For a fixed core residual count r, its per-cell residual lower bound is

    max(0, r - sum_i min(r,C_i)) = max(0, r - sum_i C_i).

Consequently, when no absolute intermediate bounds are imposed, this residual
is nondecreasing in r, so minimization can restrict the core counts to a proved
componentwise lower frontier. The model never assumes those lower pairs are
actually realizable.

An independently checked feasible relaxation at pattern 01111 modulo 30030
has 268 uncovered units, while the solver lower bound after 60 s was only 66.
Thus this relaxation alone cannot prove the sufficient threshold 786.
The direct 510510 feasibility probe at pattern 011111 and core counts (8,14)
returned UNKNOWN after 60 s; that single probe proves no seven-prime case.
Using the exact core residue sets extracted from the optimal 2310 witness,
with counts (10,12), also returned UNKNOWN after 60 s. This last test uses
the actual per-modulus intersection capacities of that one core family; its
scope does not extend to other core families.

`layered_core_projection.py` adds residual-count variables at each prefix,
constraining every child count between zero and its parent and below by parent
minus new class capacities. Certified absolute prefix bounds can constrain their
sums. Because these absolute bounds cannot simply be scaled down, the root
counts in this version range over every pair dominating a core frontier point,
not just the minimal frontier itself.

With prefix bounds 121 and 282, the layered model has a checked feasible
relaxation with 390 uncovered units at 30030. Its lower bound remained 282.
The corresponding 510510 feasibility probe returned UNKNOWN after 60 s.
Neither relaxation's feasible assignments are actual congruence coverings.

All saved feasible projection assignments were verified independently against
every aggregate cell and every prefix inequality. Small positive and negative
controls are in `code/test_projection_models.py`.

## Independent SAT and proof-log checks

`sat_uncovered_audit.py` independently constructs the fixed-pattern problem with
ordinary integer residues and a cardinality network, instead of CRT box points
and CP-SAT optimization. At 210, pattern 011, it refuted the possibility of at
most 21 uncovered residues. A separate Minisat unit-propagation check accepted
all 104976 Glucose proof additions and obtained the final contradiction.
The large compressed proof remains in the ignored temporary directory; the
solver-free base certificate remains the preferred permanent proof artifact.

The analogous 2310, pattern 0111, <=120-uncovered audit was inconclusive:
Glucose returned UNKNOWN after 120 s and CaDiCaL after a 60 s conflict-budget
loop. These timeouts neither validate nor refute the CP-SAT optimum.

## Next mathematical issue

The next investigation is recorded in `core-geometry-2026-10-25.md`: all four
210 core patterns now have audited dominating capacity catalogues. Their
weighted total-capacity tests still leave every seven-prime anchor pattern
open. A verified partial family at 30030 leaves 792 residues uncovered, while
the lower bound at its pattern remains 282.

Size-only capacities still permit core intersection patterns that cannot coexist
in a single residual set. The next refinement should track the residual's
intersections with small core congruences, or the three surviving cells modulo 6,
and preserve their shared residue choices. Merely extending the current count
recurrence or repeating these timeout probes is not a demonstrated route to the
global statement.
