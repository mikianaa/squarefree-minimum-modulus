# Mutually exclusive event families strengthen the nine-prime relaxation

The squarefree modulus-2 conjecture remains unresolved. The new result is
an independently checked obstruction to all 329 fractional states used in
the preceding three-direction barrier witness, followed by stronger exact
upper bounds for conditional group scores. It does not yet exclude a
covering, a whole nine-prime core, or all possible barrier mixtures.

## The missing condition

For a group with fixed core residues, let S be the remaining core resources
and let K_x count those resources whose chosen residue contains core point x.
Two events E(x,t) = {K_x >= t} and E(y,u) are incompatible whenever

    t + u > sum_{d in S} max_a (1[x = a mod d] + 1[y = a mod d]).

The previous relaxation imposed P(E)+P(F)<=1 for each incompatible pair.
But if a whole family of events is pairwise incompatible, the stronger
valid inequality is

    sum_{E in the family} P(E) <= 1.

The preceding triple screen tested whether three events could all hold,
using a bound of two on the sum of their probabilities. It did not test
the bound of one that follows when every pair is incompatible. Thus its
negative screening result did not cover these clique inequalities.

An elementary example in the exact saved witness is axis 2, state 4. The
fixed residues are (3,1), (5,3), (6,2), (7,6), leaving resources

    10, 14, 15, 21, 30, 35, 42, 70, 105, 210.

For core points

    8, 13, 59, 68, 89, 97, 139, 173, 179, 209,

the witness assigns probabilities

    .2, .3, .2, .2, .2, .4, .4, .2, .2, .3

to the respective events K_x=10. Their sum is exactly 2.6. At most one
such event can occur: hitting all ten resources includes hitting the single
chosen residue modulo 210, and these ten points are distinct modulo 210.
Every pair of the displayed probabilities sums to at most .8, so checking
pairs separately misses the contradiction.

## Exact separation and independent audit

`code/clique_case_constraints.py` searches maximum-weight cliques in the
event incompatibility graph. The weights are exact integer probability
numerators when processing the rational barrier. It reports one violated
clique per offending fractional state. Only thresholds at positive count
masses are needed: within a constant tail-probability interval, the largest
threshold makes incompatibility easiest to satisfy.

`code/audit_case_cliques.py` does not use the clique optimizer. It enumerates
ordinary residues of each remaining resource to verify every pair in each
reported clique, and recomputes the probability sum in integers. Results:

* 329 of 329 fractional states excluded;
* 46,630 ordinary-residue pair checks;
* largest probability excess: 358333327/100000000;
* corresponding probability sum: 4.58333327, where the valid bound is 1.

Artifacts: `data/nine-prime-three-curve-clique-cuts.json` and
`data/audit-nine-prime-three-curve-clique-cuts.json`.

The old three-direction barrier theorem remains valid for group bounds
admitting its stated fractional states. The new bounds exclude those states,
so that theorem does not limit the strengthened method. Its local mixtures
still do not construct an actual covering.

## Stronger group upper bounds

`fractional_group_scenarios.maximize` now optionally separates violated
cliques during its LP solve. Its numerical dual supplies only suggestions.
`code/clique_event_majorant.py` rounds the nonnegative multipliers to
rational numbers and rebuilds pointwise intercepts exactly. If lambda_C is
the multiplier of a clique, its total contribution to the majorant is
lambda_C because the corresponding event sum is at most one. Nonnegative
resource slopes supply the other terms of the affine majorant.

`code/audit_clique_event_majorant.py` independently checks all clique
incompatibilities by ordinary residues, all pointwise inequalities in
rational arithmetic, and every resource capacity. The existing complete
case-tree auditor now accepts these certificates through the common
capacity dispatcher.

On twelve high-weight cases from the old barrier, using the independently
audited third four-direction reward candidate, all twelve rational clique
bounds strictly improve the previous pair bounds. The independent audit
checks 274 positive clique multipliers. The largest decrease is
3620876565/1000000 = 3620.876565 in the candidate's score units. These twelve
cases do not form a complete partition. See
`data/nine-prime-clique-case-certificates.json` and its audit.

A complete case tree for group axis 6 (high-prime subset {11,17}) is also
verified. With depth limit 3, both trees have 16 nodes. The pair-tree bound
is 6716023195/1000000; the clique-tree bound is 6175864507/1000000, a decrease
of 540.158688. This covers all assignments in that group, including residues
outside the displayed core by monotone domination. It is not a whole-core
exclusion. See `data/audit-nine-prime-clique-axis6-tree.json`.

## Search continuation and tests

All 16,240 numerical scenarios in the old pair-repaired pool violated a
clique inequality. Reducing each nonminimum count mass by a common factor
per case restores the clique constraints while preserving the old resource
mean inequalities and pair constraints. The new pool keeps the same case
order and resource distributions, so selected case indices remain valid.
These are numerical search states, not proof certificates. The largest
previous clique probability sum was approximately 5.13849. The repaired
pool is `data/nine-prime-clique-repaired-scenarios.json`.

With the same 22,783 local states, 4,092 selected group rows and four reward
directions, the numerical mixture optimum changes from -0.009212663465
using the old scenarios to -0.038559545042 using the clique-repaired pool.
This is stronger finite-pool separation, not an unrestricted certificate;
complete local pricing and unrestricted rational group bounds are required.

The new dual was completely priced with rational clique group bounds and
independently audited. Its exact catalogue excess is
-115831466475/1000000; its unrestricted clique-affine excess is
-189251274020/1000000. Both are negative. The audit checks all 9,453,375 plane
tuples, 15,350 exceptional rational cases and 122 integer DPs. The normalized
catalogue gap is -0.115831466480, worse than the saved best weights; the
new candidate is retained for its additional valid states, not as a new
best obstruction. See `data/audit-nine-prime-clique-four-curve-priced.json`.

For exactly these same new weights, separately audited ordinary-affine
bounds are weaker in all 31 groups. The total improvement from clique
bounds is 73109487552/1000000 = 73109.487552. This comparison is in
`data/nine-prime-clique-group-bound-comparison.json`. Changing the weights
between different candidates would not establish that improvement.

The new seed retains 24,653 valid local states (1,870 newly added) and 4,299
group rows. A direct integer check verifies all local states have F<=0,
all preceding states are retained, and every selected group row remains
valid. The numerical pool repair preserves case order and resource
distributions. These checks are recorded in
`data/audit-nine-prime-clique-continuation.json`.

For the same new weights, the gap between prior finite-pool local minima
and exact local minima is 64362.915598. Sharing the old seed states within
each parity reduces this to 34678.157523; adding the larger older known
point pool gives no further reduction. The loss also includes omitted
group rows. This motivates a new probe with all 40,197 actual assignments
and all 16,240 repaired fractional scenarios, rather than just selected
rows. Local states remain a finite sample, so even that probe requires
complete pricing before a covering exclusion could be claimed.

That all-group run reached FEASIBLE, not OPTIMAL, after its 180-second
solver allowance (190 seconds including construction). Its feasible margin
was -0.0419904269421. No optimal dual was exported. This is inconclusive,
not infeasibility or a four-direction barrier. The terminal result is
`data/nine-prime-clique-all-group-probe.json`. Next, close this larger master
using a better solver configuration or removal of duplicate projected
columns before another expensive exact dual-pricing pass.

All 180 regression tests pass. New tests compare weighted clique search to
exhaustive subsets, exercise tiny probability tails, independently audit a
strict improvement over pair bounds, compare shifted bounds to all actual
assignments on small examples, reject false clique certificates, and check
that numerical repair preserves means and pair compatibility. The full
log is `data/tests-clique-full.log`.
