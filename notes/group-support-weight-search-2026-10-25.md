# Separate group weights in the nonpositive-support model

The global conjecture and the unrestricted nine-prime case remain open.
Everything here uses the same fixed 22-point core, primes 11,13,17,19,23,
and higher parity pattern 11111. Configuration counts below count group
assignments, not full covering systems.

## Moving beyond the verified reweighting limitation

The preceding certificate in `notes/positive-state-necessity-2026-10-25.md`
rules out a positive excess from any nonnegative combination of 22 fixed
point inequalities, even when their local integer minima are exact.

The new optimization changes the inequalities themselves. At each point x,
it gives each of the 31 group penalties Q_hx its own coefficient a_hx.
For any nonnegative coefficients define

    beta_x(a) = min_{c: F_x(c)<=0} sum_h a_hx Q_hx(c_h).

For any actual covering, or a coherent mixture supported on F_x<=0,
stochastic domination of counts implies

    sum_x beta_x(a)
      <= sum_h max_{actual assignment A} sum_x a_hx Q_hx(hit_hx(A)).

Thus a positive difference is an obstruction when every assignment maximum
is valid over the intended domain. A maximum over a listed catalogue proves
an obstruction only for that catalogue.

The numerical master initially has a finite set of admissible integer count
states. It adds actual minimizing states from exact full-domain point pricing
and separates violated catalogue-assignment constraints. The reported master
upper estimate is exploratory. Every candidate proof instead uses rounded
integer weights, exact integer support minima with attaining states, and
recomputed integer assignment scores. The first search used 1/2<=a_hx<=1
and rounding denominator 1,000. These bounds define the searched family;
they are not an optimality claim over all nonnegative event weights.

## Verified finite-catalogue obstruction

Six master iterations produced a positive lower bound

    91835540333 / 1000000 = 91835.540333

for the shared catalogue of 39,968 configurations. The independent audit
checks all 9,453,375 plane tuples using envelope coverage, 218 exceptional
planes, 36 integer tail DPs, attaining integer states at all 22 points,
and every listed assignment. Files:

* `data/nine-prime-group-support-optimized.json`
* `data/audit-nine-prime-group-support-optimized.json`

This is a rigorous finite-catalogue result. The elementary unrestricted group
upper bounds yield a negative excess, so it is not a proof for unrestricted
assignments, other cores, other patterns, or arbitrary prime supports.

## Omitted assignments matter

Searching all group residue choices found new explicit assignments. After
the first pass, 39,990 configurations are listed and the same point
inequalities still give 12290.738294>0. This expanded result is independently
checked in `data/audit-nine-prime-support-expanded-affine-dual.json`.

A second search found further actual assignments. On the resulting pool
of 40,097 configurations the same proof expression is -76937.043444.
This reevaluation is independently checked in
`data/audit-nine-prime-support-second-pool-dual.json`.
This does not invalidate either earlier certificate: it changes their
explicitly restricted domain. A negative expression also does not provide
a feasible covering or a full coherent-model witness.

The search records CP-SAT upper bounds only as exploratory values. The
independent arithmetic audits do not trust them. Actual candidate residues
are rechecked directly when evaluating the expanded pool. Relevant files:

* `data/nine-prime-support-omitted-assignments.json`
* `data/nine-prime-support-omitted-second.json`
* `data/nine-prime-support-expanded-affine-dual.json`
* `data/nine-prime-support-second-pool-dual.json`

## Improved unrestricted affine majorants

The previous group bound forced the intercept at a point to equal its
base reward. The new bound allows Q_x(k)<=b_x+a_x(k-base), with a_x>=0,
and minimizes the sum of offsets and individual resource capacities.
Numerical slopes are rounded, then each intercept is reconstructed exactly
as max_k[Q_x(k)-a_x(k-base)]. Every resource capacity and final integer floor
is also recomputed with rational arithmetic. The independent auditor checks
all inequalities and capacities without using the LP's objective.

For the current weights this improves the unrestricted excess only slightly,
to -3185341.257374. It still does not prove an unrestricted obstruction.

## Continuation

Reoptimization on the 40,097-entry pool restores a positive excess with
new group weights:

    69871764093 / 1000000 = 69871.764093 > 0.

The independent audit checks the full plane domain, 159 exceptional planes,
47 integer DPs, and the complete expanded pool. This is a new proof with new
weights, not a claim that the preceding negative reevaluation was wrong.
The elementary unrestricted excess remains -3007630.536058.

Latest certificate and audit:

* `data/nine-prime-group-support-reoptimized.json`
* `data/audit-nine-prime-group-support-reoptimized.json`

All runs reported here have finished. The unrestricted group maxima for
these newest weights still need further search and independent upper-bound
certification. No unrestricted nine-prime obstruction has been established.

`data/nine-prime-group-support-resume-input.json` contains the original event
penalty shapes with the 40,097-entry pool. Keeping the original shapes avoids
repeatedly multiplying the certificate's scaling denominators. Use saved
`support_seed` states from the latest optimization as additional initial
constraints. The optimizer now also adds every attaining tail state from
its exceptional integer DPs and preserves generated states even when the
best certificate came from an earlier iteration.
The latest certificate contains 383 saved admissible point states under
`support_seed`; pass it with `--seed` when reoptimizing after the next pool
expansion. First search omitted assignments against the latest certificate,
then retarget that same certificate to the resulting pool to measure any
loss of margin. For another coefficient optimization, retarget the original
`data/nine-prime-integer-support-pricing.json` to the expanded pool so that
the baseline event-weight scaling stays fixed.

Continue alternating exact point pricing, coefficient optimization, and
searches for omitted actual assignments. Positive finite-pool margins must
not be promoted to unrestricted bounds. The full regression suite passes
134 tests in `data/tests-group-support-affine.log`; the optimizer's updated
seed handling was also checked separately afterward.
