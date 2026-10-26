# Individual-state necessity and an integer slope-envelope bound

The global conjecture remains unresolved. The calculations below concern
only the saved 22-point core, the higher primes 11,13,17,19,23 and pattern
11111. They do not exclude this unrestricted core/pattern.

## Stronger necessary condition

Write F_x(c) for the conditional-tail lower bound at core point x. An actual
cover requires F_x(c)<=0 at every point for its actual count state. Therefore
any mixture of covers satisfies both E[F_x]<=0 and

    E[max(0,F_x)] = 0.

The latter is stronger: it does not permit positive and negative states at
the same point to cancel. Thus the previous objective
sum_x max(0,E[F_x]) can vanish without satisfying this new necessary
condition. Neither necessary condition suffices to construct a cover.

The audited latest coherent witness has

    sum_x max(0,E[F_x]) = 656403.066794462,
    sum_x E[max(0,F_x)] = 1022780.941868544.

There are 907 positive-cost and 160 negative-cost states in its support.
This is a diagnostic of this particular witness, not the optimum of either
model. See `data/audit-catalogue-point-cancellation.json` and its separate
resource-feasibility audit `data/audit-nine-prime-pointwise-catalogue-primal.json`.

## A valid separable lower bound

Fix nonnegative cumulative event penalties Q_h(k), taken from the saved
weighted dual. For a plane tuple a write

    F(a,k) = A(a) - sum_h v_h(a) k_h,
    Q(a,k) = Q_plane(a) + sum_h Q_h(k_h).

For a nonnegative slope s define

    L(a,s) = s A(a) + Q_plane(a)
             + sum_h min_{0<=k<=16} (Q_h(k) - s v_h(a) k).

For every tail tuple k, L(a,s)<=s F(a,k)+Q(a,k). Hence:

* If 0<=s<=1, then L(a,s)<=max(0,F(a,k))+Q(a,k).
* If F(a,k)<=0, the stronger conclusion L(a,s)<=Q(a,k) holds
  for every s>=0, without an upper bound on s.

For any finite list T of admissible slopes,

    beta = min_a max_{s in T} L(a,s)

is consequently a valid lower bound in the respective case. Selecting a
slope separately for each plane tuple is allowed: each displayed inequality
holds for every tail tuple. This is a lower bound; interchanging maximum
and minimum does not assert equality with the original integer problem.

At each core point, take expectations and use stochastic domination by the
actual group-assignment hit counts. Because every event penalty is
nonnegative, the expected penalty is at most the corresponding supplied
group penalty. Subtracting upper bounds on those group scores gives either
a lower bound on sum_x E[max(0,F_x)], or a necessary event-score excess for
distributions supported on F_x<=0. A positive excess in the second mode
excludes such distributions. It is not a cost bound for arbitrary states.

## Exact results and scope

The input is the older dual reevaluated on the shared catalogue of 39,968
group configurations. With slopes 0,1/8,...,1 and the original point slope,
the finite-pool lower bound improves from -330462.053616 to
-255478.165546, an increase of 74983.888070. The elementary unrestricted
capacity bound improves to -4158485.577546. Both remain negative.

Allowing slopes 2,3,...,16 in the nonpositive-support mode gives exactly
the same minima for these particular penalties. This is not an optimality
claim over all slopes or event weights.

The generator `code/positive_state_dual.py` uses count grids; the independent
auditor `code/audit_positive_state_dual.py` enumerates free-coordinate grids
and minimizes each of the 17 tail choices explicitly with integer arithmetic
and overflow checks. It also rechecks the original dual and all listed group
assignments. The certificates and audits are:

* `data/nine-prime-positive-state-catalogue-dual.json`
* `data/audit-nine-prime-positive-state-catalogue-dual.json`
* `data/nine-prime-nonpositive-support-catalogue-dual.json`
* `data/audit-nine-prime-nonpositive-support-catalogue-dual.json`

Each audit checks 9,453,375 plane tuples and 39,968 group configurations.
Small exhaustive tests compare the envelope with the actual integer minimum
and reject treating unbounded slopes as positive-cost bounds.
The full suite passes 121 tests in `data/tests-positive-state-coherent.log`.

## Further refinements and their limitations

### A verified limitation of reweighting the current point inequalities

Keeping the current event penalties and the coarse-grid constants beta_x
fixed, optimize a nonnegative coefficient w_x for each point. This only
changes the combination of the 22 already established inequalities. The
finite-catalogue optimization found no positive excess.

There is now an exact certificate for that limited conclusion. For each of
the 31 groups, choose a probability distribution on listed actual residue
assignments. In total the certificate uses 56 assignments. Their expected
aggregate event score at each point is at least beta_x, with a minimum
surplus of 420.006214172363 in the common normalized score units.

Indeed, write s_hx(A) for group h's event score at point x. The certificate
checks sum_h E[s_hx]>=beta_x for all x. For every w_x>=0 it follows that

    sum_x w_x beta_x
      <= sum_h E[sum_x w_x s_hx]
      <= sum_h max_A sum_x w_x s_hx(A).

Therefore no nonnegative reweighting of these fixed 22 inequalities can
give a positive finite-pool excess. The conclusion also holds when more
assignments are allowed. It says nothing about stronger beta_x values,
different event penalties, or additional compatibility constraints.

This mixture is not a cover and is not a feasible witness for all the
coherent count-marginal constraints. It satisfies only these 22 aggregate
score requirements. The certificate and independent integer audit are
`data/nine-prime-reweighted-support-catalogue-dual.json` and
`data/audit-nine-prime-reweighted-support-catalogue-dual.json`.
The auditor also verifies the parent support inequalities; it does not
trust numerical LP optimality as a proof of the limitation.

### Rational refinement of the slope grid

The 1/64 slope grid improves the unweighted catalogue bound to
-252601.352806. An independently audited 57-assignment score mixture still
blocks nonnegative reweighting of those particular point inequalities.
The audit is `data/audit-nine-prime-fine-reweighted-support-catalogue-dual.json`.

The finite slope grid can be refined without enumerating every tuple of
tail counts. For a fixed plane tuple, take the lower convex hull of each
tail's points (k,Q_h(k)), k=0,...,16. Each hull segment supplies an amount
v_h times its length at the segment's penalty-per-unit-volume ratio.
Sorting these ratios and filling the required volume A solves the
fractional tail problem. Its last ratio is an exact rational slope s.
The resulting objective equals L(a,s), giving matching fractional primal
and dual values for that fixed plane tuple. This is a lower bound for
integer tail counts, not an integer feasible state.

Start with the existing integer slope envelope on all plane tuples.
Discard a plane tuple if A>16 sum_h v_h, since no allowed tail counts can
satisfy F<=0. Process the remaining tuples in increasing envelope order.
The best rational value seen supplies an upper bound on the convexified
minimum. Stop once every remaining integer envelope value is at least
the floor of that upper bound. This proves the integer floor of the
global convexified minimum, while needing only a few rational refinements.

For the current 22 points, only 53 plane tuples required refinement.
The resulting unweighted catalogue bound is

    -251558257024698 / 1000000000 = -251558.257024698.

This improves the coarse bound by 3919.908521302. The unrestricted
elementary-capacity bound is still negative, -4154565.669024698.
These values do not exclude a covering.

The independent verifier reconstructs the coarse envelope using free
coordinates, checks every exceptional rational slope against all 17 choices
for each tail, and confirms that every other feasible plane tuple already
lies above the claimed bound. It does not import the fractional knapsack
generator. The stored fractional minimum is not treated as a covering
witness. Files:

* `data/nine-prime-exact-support-slopes.json`
* `data/audit-nine-prime-exact-support-slopes.json`
* `data/nine-prime-exact-reweighted-support-dual.json`
* `data/audit-nine-prime-exact-reweighted-support-dual.json`

The reweighting search again finds a 56-assignment mixture satisfying the
22 refined score requirements. This is a barrier to combinations of these
specific refined inequalities only. Different mixtures used in different
checks need not have comparable score margins.
Both independent audits pass. The full suite now passes 127 tests in
`data/tests-refined-support.log`. All runs in this continuation have finished.

### Exact integer tail minima

The fractional tail relaxation has now been removed for the current event
penalties. For a fixed plane tuple with required volume A, the generator
uses the exact recurrence

    D_0(0)=0,  D_0(r)=infinity for r>0,
    D_h(r)=min_{0<=k<=16} [Q_h(k)+D_{h-1}(max(0,r-v_h k))].

Thus D_h(r) is the minimum integer penalty for covering at least r units
using the first h tail groups. Every minimum carries an attaining integer
count vector. The independent auditor uses a different forward recurrence:
it records exact achieved credit, capped at A, and minimizes penalties after
adding each group. Both use integers with explicit overflow bounds.

The coarse slope envelope still screens the full plane domain. Fractional
pricing supplies an additional valid bound before an integer DP is needed.
For integer penalties, the ceiling of a rational lower bound is valid.
Across the 22 points, 58 exceptional plane tuples needed attention and
only 29 required integer DP. The auditor checks all 9,453,375 plane tuples
via the envelope, exceptional slopes and DPs, and checks an attaining
integer state at every point. Consequently these are exact minima for the
specified event penalties, rather than merely better lower estimates.

The unweighted finite-catalogue excess is now -251084.298, an improvement
of 473.959024698 over fractional-tail pricing. The elementary unrestricted
excess remains negative, -4154091.710. A separately checked mixture of 56
actual group assignments still dominates the 22 exact integer point minima,
with minimum normalized score surplus 226.490030454663.

The exact minima and the score-mixture certificate together prove that,
for these particular point event penalties, neither more accurate local
pricing nor any nonnegative combination of the 22 point inequalities can
produce a contradiction. They do not establish optimality over other
event penalties or other compatibility conditions. The attaining states
at different points are independent and are not a covering construction.

Certificates and independent audits:

* `data/nine-prime-integer-support-pricing.json`
* `data/audit-nine-prime-integer-support-pricing.json`
* `data/nine-prime-integer-reweighted-support-dual.json`
* `data/audit-nine-prime-integer-reweighted-support-dual.json`

All 130 tests pass in `data/tests-integer-support.log`. All these runs have
finished; the global problem and unrestricted nine-prime case remain open.

### Next optimization

Optimize the event weights for this stronger support condition, rather than
only reusing weights from the expectation model. The slope envelope supplies
valid separation bounds but not necessarily an attaining integer state.
Any integration into column generation must distinguish these lower bounds
from the actual reduced costs of candidate columns. Continue independent
full group-assignment pricing; a positive bound for a finite catalogue alone
would still not establish the unrestricted nine-prime obstruction.
The integer minima above close the fractional-tail issue for the present
penalties. A concrete next family lets the 31 group penalties have separate
nonnegative coefficients at each point, instead of a single coefficient per
point. The saved attaining integer states provide initial constraints for
that optimization. Continue pricing omitted actual assignments; a positive
result restricted to the saved catalogue would still not exclude all covers.
