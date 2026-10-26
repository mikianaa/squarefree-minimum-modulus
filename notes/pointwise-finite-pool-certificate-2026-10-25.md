# A certified obstruction for a specified finite configuration pool

The global conjecture and the unrestricted nine-prime case remain unresolved.
This result concerns one fixed 210-core residual, the higher primes
11,13,17,19,23, and higher parity pattern 11111. Its scope is explicitly
restricted to the group configuration lists in the saved certificate.

## What has been proved

There are 31 nonempty higher-prime groups. For each group, restrict its
distribution of actual core residue assignments to the listed finite pool.
Allow arbitrary mixtures within each list, and arbitrary joint count
distributions separately at every core point, dominated by those group
mixtures. The point count distributions are not restricted to the states
generated during the search.

For the expanded pool containing 1,543 listed group configurations in total,
every such feasible collection satisfies

    sum_x max(0,E[F_x]) >= 6574225384000 / 1000000000
                        = 6574.225384 > 0.

Thus this pool cannot support a covering under the stated core and pattern.
The configurations are group assignments, not 1,543 candidate covering
systems. The result excludes all mixtures and combinations drawn from the
listed groups, but says nothing about assignments omitted from their lists.

The older pool of 1,267 configurations has the stronger certified bound
142928152165000/1000000000 = 142928.152165.

## Weighted dual argument

For each point choose 0<=theta_x<=1 and nonpositive event multipliers
lambda_hxt. Full integer point pricing proves, for every local count state,

    theta_x F_x - sum_ht lambda_hxt [c_h>=t] >= beta_x.

Taking expectations and using domination by the supplied group counts gives

    sum_x theta_x E[F_x]
       >= sum_x beta_x
          - sum_h E[sum_xt (-lambda_hxt) [supplied_count_h(x)>=t]].

For each listed group configuration, evaluate the last event score directly.
Its maximum B_h over that group's finite pool is an exact upper bound for
every mixture of the listed configurations. Consequently

    sum_x max(0,E[F_x]) >= sum_x beta_x - sum_h B_h.

All scaling in the certificate is integral. An actual covering would require
F_x<=0 at every point, contradicting a positive bound if its group
assignments all belonged to the stated pools.

## Independent verification

`code/audit_weighted_coherent_dual.py` does not import the optimization model
or use an optimization solver. It independently enumerates 9,453,375 plane
tuples, minimizes over all 17 choices of each tail count, recomputes the
resource hits of every listed assignment, and verifies the final integer
numerator. NumPy is used for exact integer arrays with an overflow check.

Expanded-pool certificate and audit:

* `data/nine-prime-pointwise-expanded-dual.json`
* `data/audit-nine-prime-pointwise-expanded-dual.json`

Older-pool certificate and audit:

* `data/nine-prime-pointwise-weighted-dual.json`
* `data/audit-nine-prime-pointwise-weighted-dual.json`

The same auditor also checks an elementary unrestricted group-score bound.
For each point reward g_x(k), choose a rational slope a_x with
g_x(k)<=g_x(b)+a_x(k-b), where b=1 when the group has the modulus-1 core
resource and b=0 otherwise. Maximize each remaining resource's weighted
hit count over its ordinary residue classes and sum these maxima.
For the expanded certificate this yields the negative unrestricted bound
-4233469.465616. It therefore does not prove an unrestricted obstruction.
CP-SAT group bounds are deliberately not trusted by this auditor.

## Search status and next step

The new pointwise model was validated against independent small full-state
LPs. Compact point transport, equal-volume aggregation, and dual extension
are proved in `notes/coherent-point-transport-2026-10-25.md`. The last three
global-pricing iterations solved all 22 point problems by compact transport;
their numerical primal/dual gaps totaled less than 0.00004. This removes
the earlier large point-oracle gap at those mixtures.

Group pricing found 276 new actual assignments, expanding the pool from
1,267 to 1,543. Its unrestricted lower estimates remained negative and are
not independently proof-checked. The latest feasible rational witness has
positive-point-mean objective 692762003991465/1000000000, compared with the
initial value 890319. Its raw sum is a different objective and must not be
compared with the earlier raw-sum best value 670035.

Latest checkpoint and independent primal audit:

* `data/nine-prime-pointwise-global-pricing.json`
* `data/audit-nine-prime-pointwise-global-primal.json`

All 114 regression tests pass in `data/tests-pointwise-coherent.log`.
All long runs in this continuation have finished.

## Subsequent cross-group catalogue check

The seven saved coherent-group search files contain 308 distinct plane hit
profiles and 1,478 distinct tail hit profiles. An actual residue assignment
is reusable in any group with exactly the same core resource list. Sharing
these profiles across the five plane and twenty-six tail groups gives
39,968 group configurations. This does not identify different groups or
force their distributions to agree.

The preceding positive dual, reevaluated on this larger finite pool, gives
-330462053616000/1000000000 = -330462.053616. The independent integer
auditor checks all 39,968 assignments and 9,453,375 plane tuples in
`data/audit-nine-prime-shared-catalogue-dual.json`. The certificate is
`data/nine-prime-shared-catalogue-dual.json`. Thus the positive certificate
for the older 1,543-entry pool does not exclude this larger known pool.
This negative dual value does not establish a feasible zero-cost mixture
or a covering. All four saved duals in the preceding global-pricing run
were evaluated on this shared pool; none gives a positive bound.

`code/coherent_group_catalogue.py` scores explicit assignments using integer
event rewards, caches their actual point-hit profiles, and transfers them
only between identical resource lists. The decomposed search now accepts
`--catalogue data/nine-prime-group-catalogue.json`. It inserts the four
best known assignments for each group and uses the best as a CP-SAT hint.
Candidates survive solver timeouts, but their scores are never used as
unrestricted upper bounds. Four new tests cover exhaustive scoring,
serialization, timeout handling, and master integration. The full suite
passes 118 tests in `data/tests-catalogue-coherent.log`.

The subsequent catalogue-guided run finished after 259.547 seconds. It added
262 configurations and saved 1,805 group configurations and 210,833 point
states in `data/nine-prime-pointwise-catalogue-pricing.json`. Its best feasible
positive-point-mean value is 656403066794462/1000000000, independently
checked in `data/audit-nine-prime-pointwise-catalogue-primal.json`.
The numerical master lower estimate reached zero; that is not a feasible
zero-cost witness. Resume this newer checkpoint as both pool and warm start.

The stronger individual-state necessity is developed and audited in
`notes/positive-state-necessity-2026-10-25.md`. All 121 tests pass in
`data/tests-positive-state-coherent.log`; the runs above have finished.

For continuation, use this checkpoint as both the pool and warm start, with
the original 327-cut input and integer seed. Enable `--pointwise`,
`--compact-states 16384`, `--level-fraction 0.5`, and group pricing. Prioritize
stronger group pricing: score all already known assignments from the older
coherent-group and pooled marginal-fit files, transfer assignments between
identical resource families, and use their best scores to guide the full
assignment search. Do not interpret a positive finite-pool bound as a
positive unrestricted bound or restart superseded small-prime searches.
