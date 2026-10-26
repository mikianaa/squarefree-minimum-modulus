# Retaining the necessary inequality in each core fibre

The conditional-tail expression F_x is a lower bound for the number of
uncovered higher-coordinate tuples over the core point x. For an actual
covering, every such fibre is covered, so necessarily F_x<=0 for every x.
This is stronger than merely requiring sum_x F_x<=0.

In the existing relaxation, each point x has a probability distribution mu_x
of count states. An actual covering gives a degenerate feasible distribution
at every x, with F_x<=0. Thus each of the following is a necessary condition
in successively stronger relaxations:

1. sum_x E_mu_x[F_x]<=0;
2. E_mu_x[F_x]<=0 for every x;
3. mu_x is supported on states with F_x<=0 for every x.

The latter two can be tested by minimizing, respectively,

    sum_x max(0,E_mu_x[F_x])
    sum_x E_mu_x[max(0,F_x)].

Either minimum being strictly positive would exclude a covering for that
particular core and higher-prime pattern. A feasible zero-cost distribution
would only show that the corresponding relaxation is insufficient; it would
not construct an actual covering. All cores, patterns, and arbitrarily many
primes still need to be addressed for the global conjecture.

## Evidence that the old witnesses use cancellation

`code/audit_point_cost_cancellation.py` checks local cost arithmetic, state
ranges, point masses, and agreement with the already audited raw objective.
It does not repeat the resource-compatibility audits.

For the 327-cut negative witness, raw cost is

    -12450242120427 / 1000000000,

but the positive-state-cost sum is

    704370022644994 / 1000000000.

There are 129 positive-cost and 122 negative-cost support states. Consequently
that saved witness is not a zero-cost witness for condition 3. This calculation
does not establish that reoptimization could not achieve zero.

It also fails condition 2: thirteen point expectations are positive, and
their sum is 149289827122240/1000000000. Thus the simpler strengthened
objective already rejects this particular witness.

For the coherent witness with raw cost 670035, the positive-state-cost sum is
890319. Its distribution has one state per point, so conditions 2 and 3 give
the same value on that witness. Six of its 22 point costs are negative.

See `data/audit-327-cut-point-cancellation.json` and
`data/audit-coherent-point-cancellation.json` for exact per-point values.

## Computational consequence

Condition 2 can reuse the existing integer local-state pricing. In a Benders
master, minimize sum u_x with u_x>=0 and u_x>=z_x, where z_x retains all
supporting inequalities for the optimum point expectation. To certify a
global bound, aggregate the dual into nonnegative point weights theta_x
and event penalties. Reprice theta_x F_x with those event penalties, rather
than incorrectly keeping the coefficient of F_x equal to one. Prove or audit
every group-event upper bound before treating a positive result as an
obstruction. Finite configuration-pool bounds alone are insufficient.

Condition 3 has a harder pricing problem because max(0,F) is not affine in
tail counts. A valid lower-pricing relaxation is to use theta F for
0<=theta<=1 and maximize over theta, separately for each fixed plane tuple.
This would need a new implementation and validation before use.

Condition 2 is now implemented with `--pointwise` in
`code/decomposed_coherent_groups.py`. The independent coherent-witness audit
reports every point numerator and explicitly checks any supplied
`positive_point_means_numerator`. A negative raw sum is not accepted as a
witness that all point means are nonpositive.

The numerical master optionally chooses a nearby mixture at an intermediate
objective level (`--level-fraction`), while retaining the primary optimum as
the numerical finite-pool lower bound. The stabilizing objective is never
reported as an obstruction bound. Compact fixed-mixture point transport and
globally valid dual extensions are documented in
`notes/coherent-point-transport-2026-10-25.md`.

The first pointwise run and a run using diverse integer pricing were
inconclusive. The combined transport/stabilized run improved the independently
auditable feasible objective from 890319 to 766302.317259190. Its numerical
finite-pool lower estimate remained zero. This leaves the stronger relaxation
unresolved; it does not show that an exact zero witness exists.

Condition 3 has not yet been optimized.

Later in the same continuation, full integer dual audits give a positive
bound for the explicitly listed pool of 1,543 group configurations. This
excludes that restricted pool only. The current feasible pointwise value is
692762.003991465, and the unrestricted group bounds remain inconclusive.
See `notes/pointwise-finite-pool-certificate-2026-10-25.md` for the exact
certificate, scope, completed tests, and next direction.
