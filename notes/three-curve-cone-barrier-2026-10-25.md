# A verified barrier for the current three reward directions

The global squarefree modulus-2 conjecture remains unresolved. This result
concerns the fixed 22-point core and high-prime pattern 11111 through 23,
the three specified reward directions at each group/point, and group bounds
that admit the explicit fractional case states below. It is not a theorem
that the upper bound 10 is optimal for the distortion method.

Continuation: all 329 fractional states in this witness are now excluded
by independently verified clique inequalities. The theorem below still
applies to bounds admitting these states, but not to the strengthened
bounds. See `notes/clique-case-constraints-2026-10-25.md`.

## The certificate and its implication

Let Q_{h,x,b}(k) denote the cumulative nonnegative threshold reward for
group h, core point x, and direction b. The three directions are the original
baseline and the two basis files `nine-prime-free-threshold-support.json`
and `nine-prime-free-threshold-known-support.json`. Coordinates whose baseline
mass cap is zero are excluded, as in the search. Identically zero directions
need no inequality.

For each core point x, the certificate gives a probability distribution
mu_x on integer local count states satisfying F_x<=0. For each group h,
it gives a probability distribution nu_h on actual group assignments and
fractional resource/count states with specified fixed residues. Every
fractional state satisfies all resource-mean constraints and all conditional
two-point threshold incompatibilities exactly.

The independently checked inequalities are, for every allowed h,x,b,

    E_nu_h Q_{h,x,b}(count_x) - E_mu_x Q_{h,x,b}(count_h)
        >= epsilon Q_{h,x,b}(N_h),

where

    epsilon = 5799742088867890369190237 / 29436385894500000000000000000 > 0.

Therefore any nonnegative combination Q_{h,x}=sum_b a_{h,x,b}Q_{h,x,b}
satisfies the same dominance after weighting and summing. The exact local
minimum is at most its expectation under mu_x. Any group upper bound that
admits the displayed states is at least its expectation under nu_h. Hence

    sum_x (exact local minimum) - sum_h (such group upper bound)
        <= -epsilon sum_{h,x,b} a_{h,x,b} Q_{h,x,b}(N_h) <= 0.

This rules out a positive separation for all nonnegative combinations of
these three directions, even if local pricing is exhaustive. It is more
than a failed coefficient search. It does not apply to group bounds that
exclude some of the displayed fractional states through stronger joint
compatibility or deeper residue case refinement. It also does not rule out
separation against actual assignments alone, since the witness includes
fractional states. Separate local/group mixtures are not a covering system.

## Exact construction and independent verification

The numerical probe used 17,211 local states and 3,476 selected group states.
It found positive normalized surplus 0.000197030765873. This numerical value
was used only to choose a rational witness.

Resource and count probabilities were reconstructed as integer numerators
over 100000000. Mixture weights have denominator 10000000000. Count mass was
moved to the minimum count if needed to restore resource means. If pair
probabilities exceeded their capacity, all nonminimum masses were reduced
and rounded downward. Such reductions preserve the mean constraints. The
independent auditor ignores those repair claims and checks the resulting
integer distributions directly.

Files:

* `data/nine-prime-three-curve-barrier-probe.json`: numerical search only;
* `data/nine-prime-three-curve-barrier.json`: exact witness;
* `data/audit-nine-prime-three-curve-barrier.json`: independent standard-library audit.

The audit verifies 2,030 coordinate inequalities using 691 local integer
states, 200 actual group states, and 329 fractional group states. It checks
3,110 resource distributions and 1,044,562 conditional pair inequalities.
No solver optimality claim is used by this audit.

## A fourth direction escapes this particular witness

There are 4,709 individual threshold coordinates where the group mixture's
tail probability is below the local mixture's. At each group/point having
such a deficit, choose a threshold with the largest deficit and assign its
event weight 1. This gives 679 nonzero curves in
`data/nine-prime-barrier-escape-basis.json`.

Independent comparison with the exact witness gives local minus group score

    156290306497746925237 / 1000000000000000000 > 0.

Thus this specific mixture cannot certify a barrier after that direction is
added. The comparison audit is `data/audit-nine-prime-three-curve-barrier-escape.json`.
Escaping one witness does not by itself produce a covering obstruction or
show that no other dominating mixture exists.

With the fourth direction, the numerical mixture problem on the same finite
state pools has optimum -0.0204220309199743. Its numerical dual supplies a
candidate separating reward. This still uses an incomplete local pool and
only selected group states. It must be checked against the full integer
local domain, the full actual assignment catalogue, and ultimately complete
unrestricted group bounds before any covering exclusion can be asserted.

`code/price_curve_basis_dual.py` rounds the numerical dual to integer event
weights, performs complete local integer pricing, and computes maxima over
all 40,197 actual catalogue assignments. It retains newly found local states
and high-scoring omitted group rows for another mixture probe. Its fractional
scenario scores remain numerical search diagnostics; its integer certificate
requires the separate `audit_refined_support_slopes.py` audit.

## Complete pricing and continuation

The first four-direction candidate passed that independent audit, but its
catalogue excess is negative: -151027982713 / 1000000. Normalizing by total
event weight gives -0.151027982703. The audit checks all 9,453,375 plane
tuples, 1,108 exceptional rational cases and 126 integer DPs. The gap between
the previous finite local pool and the exact local minima alone is
124207.079232. Thus its finite-pool separation was not a covering obstruction.
Files: `data/nine-prime-four-curve-dual-priced.json` and
`data/audit-nine-prime-four-curve-dual-priced.json`.

The next primal-simplex probe used the enlarged 19,120-state seed and 3,716
group rows. It reached numerical optimum -0.01319446454051354. Complete
pricing of its new dual again failed to separate: the exact catalogue excess
is -44932618101 / 1000000 and the ordinary unrestricted affine excess is
-177652262984 / 1000000. Its normalized catalogue gap is -0.044932618107,
still worse than the previously saved best -0.009566797184. Independent
audit checks all 9,453,375 plane tuples, 3,957 exceptional rational cases and
82 integer DPs. It retains 20,877 local states and 3,922 group rows for the
next probe. These new states are useful even though these are not the best
scoring weights. Files: `data/nine-prime-four-curve-corrected-priced.json` and
`data/audit-nine-prime-four-curve-corrected-priced.json`.

A third numerical probe reached optimum -0.010949742475865554 on those
20,877 local states and 3,922 group rows. Its dual was again completely
priced and independently audited. The catalogue excess is
-27766730378 / 1000000; the ordinary unrestricted affine excess is
-161640666298 / 1000000. The normalized catalogue deficit is
-0.0277667303775. Thus the successive four-direction normalized deficits
are -0.151028, -0.044933, -0.027767: improving within this run but still
worse than the previously retained best weights. None is an obstruction.
The third audit checks 9,453,375 plane tuples, 3,889 rational exceptional
cases and 150 integer DPs. Its continuation retains 22,783 local states
and 4,092 group rows; a separate direct check verifies every retained
state has F<=0 and all previous states and row indices are preserved.
Files: `data/nine-prime-four-curve-third-priced.json`,
`data/audit-nine-prime-four-curve-third-priced.json`, and
`data/audit-nine-prime-four-curve-third-continuation.json`.
The next coefficient probe should start from this newest state seed,
while retaining the original baseline and the same three added basis files.

The larger 58,659-state probe with dual simplex terminated ABNORMAL (status
4); this is not infeasibility. A fresh primal-simplex run on those same pools
terminated FEASIBLE at its time limit, with margin -0.01361816026985685 and
no optimality certificate. Neither result establishes a four-direction
barrier or a separation. The latter result is saved in
`data/nine-prime-four-curve-expanded-primal-probe.json`.

## Conditional three-point consistency search

The 329 fractional states in the verified three-direction witness were also
screened for conditional triple-event inequalities. For each three points
and thresholds, if no assignment of the remaining core resources can attain
all three thresholds, the sum of the three event probabilities is at most
two. Fixed residues are included in these tests.

The first search used just the maximum sum of three counts. A stronger
search uses the exact partition of the three points into residue classes
for each remaining resource. Common resources hit all three points;
two-class resources choose a pair or its singleton; three-class resources
choose one point. Enumerating the pair/singleton choices gives the maximum
attainable third count for every pair of required first and second counts.
Only the first infeasible third threshold need be checked, since its tail
probability decreases with the threshold.

Neither search found a violated triple-event inequality in these displayed
states. This is a screening result, not an independent certificate of joint
realizability: weighted marginal constraints, four or more points, and
stronger residue case refinements remain possible. Files:
`data/nine-prime-three-curve-barrier-triple-cuts.json` and
`data/nine-prime-three-curve-barrier-exact-triple-cuts.json`.

The separator is `code/separate_barrier_triples.py`. Any positive cut can be
checked by `code/audit_barrier_triples.py`, which independently enumerates
ordinary residues and computes threshold-capped reachable count triples.
Tests compare the frontier and separation to ordinary residue enumeration, detect a
three-event conflict that pair tests miss, reject a tampered claimed gap,
and preserve marginals of an actual assignment mixture. All four new
tests pass; the preceding complete suite of 170 tests and the six
primal-option regression tests also passed.
