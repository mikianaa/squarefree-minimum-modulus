# Changing individual threshold curves after the rescaling barrier

The general conjecture remains unresolved, as does the fixed 22-point core
and higher-prime pattern 11111. This continuation changes the reward curves
themselves; it does not repeat the coefficient-only search ruled out in
`notes/case-relaxation-barrier-2026-10-25.md`.

## Search model and proof boundary

For each group h and point x, let C_(h,x) be the original curve's total
reward. The new variables u_(h,x,t) are nonnegative and sum to at most one
over thresholds t. The reward at count k is

    Q_(h,x)(k) = C_(h,x) sum_(t <= k) u_(h,x,t).

Thus each curve may change shape while keeping its total mass at most the
old total. Every nonnegative rescaling of the old curve remains available.
The numerical master has point lower-bound variables and group upper-bound
variables. It separates violated constraints from the actual assignment
catalogue and from the saved fractional case states. For the latter it
reconstructs all threshold probabilities from the count distributions;
the old scalar scenario payoffs are not reused with the new curves.

The rounded candidate is evaluated independently of whether all sampled
master constraints have been added. Exact integer local pricing covers the
entire F_x <= 0 domain. Actual assignment scores are evaluated with integer
arithmetic; unrestricted bounds use reconstructed rational affine
majorants. Fractional search rows never enter the actual-assignment proof
pool. Truncating the numerical separation loop can produce a poor candidate,
but cannot turn its independently audited certificate into a false proof.

`code/optimize_threshold_support.py` implements this search. Small tests
show that free curves can improve a fully enumerated two-point obstruction,
and check fixed/unit-resource count shifts and separation of search rows
from actual assignment proof entries.

## Results

The initial unrestricted separation run exhausted its numerical search
budget before producing a candidate. The implementation now retains the
last solved weights and checks them even if sampled constraints remain
unclosed; this is recorded explicitly in each iteration's metadata.

The first complete candidate, in `data/nine-prime-free-threshold-support.json`,
has exact finite-catalogue excess -4003218.127193 and unrestricted affine
excess -8736334.475629. Its audit verifies 170 exceptional plane tuples and
119 integer DPs across all 9,453,375 plane tuples. This is not a successful
covering obstruction.

The previous fixed-curve barrier witness *does* fail for these new curves.
Keeping its two specific mixtures fixed, local expected reward minus group
expected reward becomes approximately +1912765.596681. Direct cumulative
reward evaluation and an independent threshold-probability calculation
agree as exact fractions. This only defeats that particular witness; it
does not prove that no other mixture survives or that the new curves exclude
a covering.

## Reusing known local states

The prior coherent search stored 210,833 local count vectors. Domain and
integer-cost filtering retains 41,375 distinct vectors with F_x <= 0. These
are now available for cheap separation during the master search, rather than
being added all at once or repeatedly rediscovered by integer pricing.

The run using this pool saves `data/nine-prime-free-threshold-known-support.json`:

* 2,270 local states retained in the active support seed, including 146 new
  states from the final exact pricing pass.
* Exact catalogue excess -3783060.183220 on 40,197 actual assignments.
* Unrestricted affine excess -8577714.691432.
* 238 exceptional plane tuples and 146 integer DPs, independently audited.
* Sampled master constraints are still unclosed; no numerical optimality
  claim is made.

This candidate also defeats the old specific barrier mixture, by approximately
+1502904.372409 in its normalized reward units. However, exact local minima
are much smaller than the restricted numerical master predicted. Both new
candidates already fail against the actual finite catalogue. Since the true
group maxima are at least those catalogue maxima, improving unrestricted
group upper bounds alone cannot make either fixed candidate an obstruction.

Artifacts:

* `data/audit-nine-prime-free-threshold-support.json`.
* `data/audit-nine-prime-free-threshold-known-support.json`.
* `data/nine-prime-known-nonpositive-state-counts.json`.
* `data/audit-nine-prime-free-threshold-old-barrier-comparison.json`.
* `data/audit-nine-prime-free-threshold-known-old-barrier-comparison.json`.

The full suite passes 145 tests in
`data/tests-free-threshold-support-final.log`. After adding the independent
new-curve comparison, its four focused tests pass in
`data/tests-free-threshold-barrier-comparison.log`.

## Next search

The fully free threshold master has many more variables than the former
682-coefficient problem and currently spends substantial time adding group
and local constraints before exact pricing. A useful next step is a smaller
curve basis: allow nonnegative combinations of the old curve and the newly
found curves at each group/point, or activate selected independent threshold
directions. This includes changes excluded from the old rescaling family
while reducing the master dimension. Retain the 2,270-state support seed,
the 41,375-state separation pool, the 40,197 actual assignments, and the
16,240 fractional scenarios. Any resulting candidate still needs full-domain
local pricing and unrestricted group-bound verification.

The full objective remains a proof for arbitrary squarefree support or an
independently verified covering counterexample. None of the restricted
searches or mixture comparisons above settles that objective.

Continuation: the three-curve basis now has an independently audited
candidate and 2,895 saved local support states. A true pair incompatibility
excludes the old averaged fractional barrier witness. See
`notes/three-curves-and-case-pairs-2026-10-25.md` for the exact inequality
and the remaining work on rational pair-aware upper bounds.
