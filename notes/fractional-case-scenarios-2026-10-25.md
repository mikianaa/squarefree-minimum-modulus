# Fractional case scenarios for event-coefficient search

The global conjecture remains unresolved. This continuation concerns only
the saved 22-point core and higher-prime parity pattern 11111 through 23.

## The relaxation

Fix one group's resource shapes and a partial assignment of residues. Write
h_x for the number of fixed resources hitting point x, r for the number of
remaining resources, and b=1 if the unit resource is still remaining, else 0.
Let Q_x(k) be the original nonnegative, nondecreasing threshold reward at
total count k. The coefficient being optimized is alpha_x >= 0.

For each remaining resource d != 1, introduce probabilities z_(d,a) over
the residues occupied by the residual points, summing to one. At each point
introduce probabilities u_(x,k), k=b,...,r, also summing to one, satisfying

    sum_k (k-b) u_(x,k) <= sum_(d != 1) z_(d, x mod d).

Maximize sum_(x,k) alpha_x Q_x(h_x+k) u_(x,k). Every ordinary assignment
using occupied residues is feasible by taking point masses. Unoccupied
residues can be replaced by occupied ones without decreasing any count or
score, so the relaxation bounds the maximum over all ordinary assignments.
It need not preserve the full distributional compatibility of the counts.

The LP dual chooses nonnegative slopes a_x and intercepts t_x with

    t_x + (k-b) a_x >= alpha_x Q_x(h_x+k),

and a capacity c_d >= sum_(x: x mod d = a) a_x for each occupied residue a.
Its objective is sum_x t_x + sum_d c_d. This is the free-intercept affine
majorant used by the existing rational certificate, including fixed hits.
The intercepts can be restricted to nonnegative values because k=b and
the rewards are nonnegative. The LP is feasible and bounded, so ordinary
finite-dimensional LP duality gives equality of these relaxation optima.

## Optimization use and proof boundary

A feasible fractional state has point payoff s_x = sum_k Q_x(h_x+k)u_(x,k).
It supplies the linear constraint B >= sum_x alpha_x s_x when searching
coefficients. We generate one maximizing state at each leaf of the saved
case trees under the previous weights, giving 8,176 scenarios. These are
only some states of the case relaxations; changing coefficients can expose
new optimizing states. Closing this finite scenario list does not close
the full relaxation.

`code/fractional_group_scenarios.py` saves the numerical resource and count
distributions as well as their point payoffs. Its output is exploratory;
floating-point scenarios are never used as a proof of an upper bound.
`code/optimize_group_support_weights.py` adds these rows only to its search
matrix. They are not inserted into the actual-assignment catalogue or its
exact maxima. In scenario mode it does not stop merely because the actual
finite pool has positive excess, and ranks candidates by their numerical
scenario margin. Every reported proof still requires exact point pricing
and independently audited rational bounds for all residue cases.

Tests compare the primal relaxation against the affine dual and complete
ordinary-residue enumeration on small instances, including a remaining
unit resource and fully fixed cases. An integration test verifies that
fractional scenarios cannot enter the actual-assignment proof pool.

## Reproduction

The baseline is `data/nine-prime-group-support-third-resume-input.json`;
the old weights are `data/nine-prime-group-support-reoptimized.json`;
the case trees are `data/nine-prime-branched-event-depth6.json`.
The generated scenarios are `data/nine-prime-fractional-case-scenarios.json`.
The original event curves must be retained when generating or consuming
these payoffs; repeatedly scaling already optimized multipliers changes
their units and risks integer overflow.

## Result of the first scenario-aware search

With 0.5 <= alpha_(h,x) <= 1, the search completed eleven full pricing
iterations in 252 seconds. It retained 1,104 valid nonpositive-cost local
support states. The selected candidate's exact catalogue excess is
-137958.820400 on the 40,197 actual configurations. Its numerical margin
against both the catalogue and the original scenarios is -468548.118556;
this numerical quantity is not an unrestricted bound.

New depth-six residue trees were generated for the resulting coefficients.
The full independent audit verifies exact integer support minima across
9,453,375 plane tuples, including 539 exceptional cases and 92 integer DPs,
and all rational tree upper bounds. The resulting unrestricted excess is

    -545194768584 / 1000000 = -545194.768584.

This is still negative. Neither the actual catalogue nor unrestricted
assignments are excluded by these new weights. Changing coefficients
also changes the scale and the point inequalities, so a smaller negative
number by itself is not a scale-independent measure of improvement.

Files:

* `data/nine-prime-scenario-support-weights.json`: optimizer checkpoint,
  including 1,104 support seed states and all eleven iteration summaries.
* `data/nine-prime-scenario-branched-depth6.json`: new complete case trees.
* `data/nine-prime-scenario-branched-dual.json`: combined exact certificate.
* `data/audit-nine-prime-scenario-branched-dual.json`: independent audit.
* `data/tests-fractional-case-scenarios.log`: all 138 regression tests pass.

To diagnose incomplete scenario separation, the original 8,176 leaf LPs
were solved again under the new coefficients, keeping the old tree fixed.
Using only the old fractional scenarios gives margin -468548.031359,
whereas reoptimizing their leaf relaxations gives -607201.324127. Thus the
old scenario list misses approximately 138653.292768 of relaxed score at
these new weights. This comparison is numerical, on the same tree; it does
not confuse newly found fractional states with actual residue assignments.
The small difference from -468548.118556 above comes from also including
actual catalogue assignments in the optimizer's separation matrix.

The second pass is `data/nine-prime-fractional-case-scenarios-second.json`.
The union, deduplicated by the point payoff vectors, has 16,240 scenarios
and is saved as `data/nine-prime-fractional-case-scenarios-union.json`.
The numerical comparison is in
`data/nine-prime-fractional-case-scenario-comparison.json`.

Next, use that union and the 1,104 support states while testing coefficients
down to zero, to distinguish the imposed lower coefficient bound from a
limitation of this relaxation. A zero numerical optimum would not by itself
be an optimality proof. If needed, retain a verifiable dual mixture or
strengthen the residue/count relaxation through further cases or additional
compatibility. The global problem, all nine-prime cores, and even this one
fixed core/pattern remain unresolved.

Update: the zero-floor search and an independent exact mixture certificate
now establish a limitation for all nonnegative rescalings of these fixed
curves with these case relaxations. This closes that coefficient-only route,
but 4,235 individual threshold deficits leave arbitrary changes to the
reward curves open. See `notes/case-relaxation-barrier-2026-10-25.md`.
