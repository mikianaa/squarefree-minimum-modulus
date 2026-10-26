# Three-curve search and missing pair compatibility

The global conjecture and the fixed nine-prime core/pattern remain unresolved.
The new facts are an audited three-curve candidate and an exact explanation
of one compatibility condition missing from the fractional case relaxation.

## Combining curves

`code/optimize_threshold_support.py` now optionally uses a small basis of
reward curves. At each group/point, each basis curve is normalized by its
own total mass, and the nonnegative basis coefficients sum to at most one.
The original curve's total mass caps the combined reward. This allows new
shapes without increasing the point/group reward budget.

The three-curve run uses the original curve and the two free-threshold
candidates from `notes/free-threshold-support-2026-10-25.md`. It retains the
40,197 actual assignments, the 16,240 fractional scenarios, and the 41,375
known valid local count states. Three complete pricing iterations finished;
each closed the constraints of those sampled pools before full local pricing.

Exact catalogue excesses were -2065872.863901, -1619166.412346, and finally
-1109526.907486. The last candidate's unrestricted affine excess is
-3714315.319351. Both remain negative. Exact full-domain pricing still finds
weak local states outside the active seed; 2,895 states are now retained.

The independent audit checks all 9,453,375 plane tuples, 349 exceptional
tuples, and 45 integer tail DPs for the final candidate:

* `data/nine-prime-three-curve-support.json`.
* `data/audit-nine-prime-three-curve-support.json`.

Closing sampled constraints is not a claim that the entire local-state or
group-assignment problem has been optimized. Improving group bounds alone
cannot make these fixed curves work, because they already fail against the
actual finite catalogue. Further weight changes and local-state constraints
remain necessary.

## A specific missing compatibility condition

For N core resources with moduli d, two points x,y satisfy

    count(x) + count(y) <= N + #{d : x == y modulo d}.

Every resource contributes at most one hit unless its chosen class can hit
both points; only the listed congruent resources can contribute two. Hence
events count(x)>=t and count(y)>=u cannot both occur when t+u exceeds this
capacity. Their probabilities must then sum to at most one. This is the
existing pair-threshold inequality in `code/point_pair_capacity.py`, not a
new general covering-system theorem.

The fractional barrier witness from the fixed-curve investigation violates
one such inequality *after averaging its group scenarios*. In group
{11,17}, the core resources are all 16 divisors of 210. The two residual
points are 67 and 164. Their difference 97 is coprime to 210, so only d=1
can hit both, and the pair capacity is 17. Therefore

    P(count(67)>=16) + P(count(164)>=2) <= 1.

The saved fractional mixture instead gives exactly

    1 + 4961018949 / 50000000000 = 1.09922037898.

Thus it cannot represent actual group assignments and does not survive the
pair-compatible relaxation. This does not invalidate the earlier theorem
about the weaker, explicitly stated relaxation. It identifies information
that the weaker relaxation discarded.

## Conditional pair constraints

Inside a residue case, let h_x,h_y be the fixed hits, r the number of
remaining resources, and s the number of those remaining resources for
which x and y have the same residue. Then the sharper case capacity is
h_x+h_y+r+s. The boundary threshold pairs t+u=h_x+h_y+r+s+1 suffice because
higher thresholds only decrease probabilities.

The diagnostic found 37 individual scenarios violating global pair bounds
and 41 violating conditional pair bounds. These are counts of scenarios
with a recorded worst violation, not counts of all inequalities. The single
averaged violation above and all recorded individual violations are checked
by the independent standard-library auditor:

* `data/nine-prime-case-barrier-pair-check.json`.
* `data/audit-nine-prime-case-barrier-pair-check.json`.
* `code/audit_case_pair_check.py`.

Separately, checking the existing 327 certified group-event cuts found two
individual violating states but no averaged violation. The diagnostic is
`data/nine-prime-case-barrier-event-cut-check.json`; the prior full audit of
those cuts remains `data/audit-nine-prime-transferred-group-primal.json`.

## Pair-aware exploratory pricing

`code/fractional_group_scenarios.py` accepts `--pairs`. Its LP adds violated
conditional pair inequalities and resolves until all pair inequalities have
been checked numerically. A small exhaustive test reduces a fractional
score from 3 to the actual maximum 2, with fixed and unit resources tested.
These LP outputs are still exploratory; they are not rational upper-bound
certificates.

Under the final three-curve weights, repricing the 40 distinct cases behind
the 41 recorded violations improves 38 numerical relaxation optima. The
largest score decrease is approximately 21490.464624, and 1,416 pair
constraints were added in total. The selected cases and distributions are in
`data/nine-prime-pair-compatible-case-probe.json`. This does not cover every
case and proves no group-wide maximum.

All 150 tests pass in `data/tests-three-curve-and-pair-support.log`.

## Continuation

Preserve the three-curve basis, its 2,895-state seed, and the 41,375-state
local separation pool. Before using the stronger pair relaxation in a
proof, extend the affine dual certificates with nonnegative pair-cut
multipliers and independently verify their rational inequalities and pair
capacities. Filter or regenerate old fractional scenarios that violate pair
compatibility; merely adding new valid scenarios leaves the old invalid
search constraints in place. Then reoptimize and perform exact local
pricing. A different feasible fractional mixture may survive, so eliminating
this particular witness is not an optimality or impossibility result for
the stronger model.

Continuation: the pair-aware affine certificates and their independent
auditor are now implemented. They verify all 40 selected cases and a full
group-{11,17} case tree. The search pool has been repaired and reoptimized;
the new candidate is still insufficient. See
`notes/pair-affine-certificates-2026-10-25.md`.
