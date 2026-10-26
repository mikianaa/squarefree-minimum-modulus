# Rational pair-aware affine certificates

The global covering conjecture remains unresolved. This continuation makes
the pair-compatible case bounds independently checkable and repairs the
fractional search pool. It does not turn the fixed nine-prime case into a
covering obstruction.

## Certificate inequality

After fixing some resource residues, remove their hits from every event
threshold and collect already satisfied events into an integer constant.
Let r be the number of remaining resources and b=1 if the unit resource is
remaining, otherwise b=0. At each point x, the remaining count k lies in
{b,...,r}, and Q_x(k) is its shifted nonnegative threshold reward.

For each pair cut c on points x,y and thresholds t,u, verify directly that

    t + u > r + #{d remaining : x == y modulo d}.

Then I(count(x)>=t)+I(count(y)>=u) <= 1 for every actual assignment.
Choose nonnegative pair multipliers mu_c and nonnegative slopes a_x, and
intercepts q_x satisfying, for every integer k in the count domain,

    Q_x(k) <= q_x + a_x(k-b)
                 + sum_(c incident to x) mu_c I(k >= threshold_(c,x)).

Summing over points and resources gives the upper bound

    sum_x q_x
      + sum_(d remaining, d != 1) max_a sum_(x == a mod d) a_x
      + sum_c mu_c.

Add the already satisfied integer constant and round down. The final reward
is integer-valued, so this rounding is valid. A separate pointwise reward
bound can also be taken if smaller.

Numerical LP duals only suggest a_x and mu_c. The generator rounds them to
rational values and reconstructs each intercept as the exact maximum needed
to satisfy all its count inequalities. It recomputes all resource maxima and
the final floor exactly. If the new certificate does not improve the old
affine bound, it retains the old bound.

## Independent verification

* `code/pair_event_majorant.py` generates rational certificates.
* `code/audit_pair_event_majorant.py` independently checks every pair
  incompatibility, count inequality, resource capacity, and final bound.
  This checker uses only the standard library and imports no generator or
  optimization solver.
* `code/branched_event_bound.py --pairs` uses these bounds at case-tree nodes.
  The existing independent tree checker verifies the full case partition and
  inherited bounds, dispatching pair certificates to the new checker.

The primal LP's redundant upper bounds on probability variables were removed
so their dual values do not obscure the mean-resource and pair multipliers.
The mass equations still imply that probabilities are at most one. This
change does not alter the relaxation's feasible set.

Tests compare the bounds with exhaustive ordinary assignments on small
instances, including fixed cases, a unit resource, and fully fixed cases.
They reject false pair thresholds, negative multipliers, incorrect
intercepts, and changed fixed-event constants.

## Verified improvements

For the 40 selected cases from the preceding numerical probe, 38 bounds
strictly improve after rational reconstruction. The independent audit checks
353 positive pair multipliers. The largest bound decrease is exactly
21490464624 in integer certificate units, or 21490.464624 after dividing by
the event scale 1000000.

* `data/nine-prime-pair-case-certificates.json`.
* `data/audit-nine-prime-pair-case-certificates.json`.

These 40 cases do not form a complete partition of every group. A separate
complete depth-six tree was therefore tested on group {11,17}, axis 6, under
the weights in `data/nine-prime-three-curve-support.json`:

* Ordinary affine tree: bound 121950699883, 106 nodes.
* Pair-aware tree: bound 119163880384, 116 nodes, 1,253 pair multipliers.

Both bounds cover every actual residue assignment in that group and pass
independent audit. Their difference is 2786819499, or 2786.819499 in normalized
units. The adaptive tree partitions differ; the controlled 40-case comparison
above isolates improvements within fixed cases.

The complete-tree audit is `data/audit-nine-prime-pair-tree-axis6.json`, with
trees in `data/nine-prime-plain-tree-axis6-depth6.json` and
`data/nine-prime-pair-tree-axis6-depth6.json`. These numerical values belong
to the named old weights and must not be reused for newly optimized weights.

## Repairing the search pool

Simply adding pair-compatible states would leave incompatible old states in
the optimization constraints. `code/repair_pair_scenarios.py` instead replaces
each offending count distribution with a mixture of itself and its minimum
count. Resource residue probabilities remain unchanged, and expected count
demand only decreases.

If the largest pair-event probability sum is L>1, scale all nonminimum
count masses by slightly less than 1/L. Any potentially violated pair has
both thresholds strictly above their minimum counts, so both probabilities
scale by this factor. Pairs with an event already forced at the minimum have
their other event impossible and cannot violate the inequality.

Of the 16,240 old scenarios, 2,793 required repair. The largest original pair
sum was 16/9; the smallest retained nonminimum mass is approximately 9/16.
The replacement pool is `data/nine-prime-pair-repaired-scenarios.json`.
Its checks are numerical and it is used only to guide optimization. Final
proofs rely on independent rational upper bounds, not these scenario values.
The pool remains a sample, not a complete description of the pair relaxation.

All 156 tests pass in `data/tests-pair-affine-certificates.log`.

## Continuation

The completed repaired-pool run saved
`data/nine-prime-three-curve-pair-repaired-support.json`. The selected
candidate has independently verified catalogue excess -546600.730794 and
ordinary unrestricted affine excess -2560420.711351. These bounds use its
new weights; the complete pair-aware group-6 bounds above use the older
three-curve weights. The new candidate still fails against the actual pool,
so pair-aware group-bound refinement alone cannot make these fixed weights
an obstruction.

The audit `data/audit-nine-prime-three-curve-pair-repaired-support.json`
checks all 9,453,375 plane tuples, 612 exceptional cases, and 90 integer DPs.
Three pricing passes retained 3,381 local support states. The last sampled
master was not fully closed when its search budget ended; the exact audit
is independent of that numerical stopping condition. Candidate selection
uses sampled-scenario margin, which need not select the iteration with the
best actual-catalogue margin.

The three-curve search has been rerun with the repaired scenarios while
retaining the 41,375 known nonpositive local states and all 40,197 actual
assignments. Exact local pricing remains essential: sampled local constraints
can be closed while undiscovered integer states still destroy the margin.
Preserve every newly found support state. If this bottleneck persists,
generate additional valid states from diverse plane tuples by solving the
fractional tail problem and rounding tail counts upward; these provide cheap
extra constraints without replacing exact final local pricing.

Further pair-aware case pricing can enlarge the scenario pool. Any positive
finite-pool or sampled-scenario margin would still require newly computed
complete group bounds for its exact weights, all core/pattern cases, and an
argument addressing arbitrary squarefree supports before the global goal
could be considered settled.
