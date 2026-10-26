# Exact penalties for coherent group marginals

The global conjecture is unresolved. This lemma concerns optimization of the
existing conditional-tail relaxation; it does not assert that its minimum is
positive, or that a relaxed negative point distribution is a covering.

## A one-hit inequality

For a fixed core point of parity e, put D_i=q_i-1-[b_i=e], f_i=D_i-c_i.
The plane counts satisfy 0<=c_i<=min(14,D_i), and tail counts h_S satisfy
0<=h_S<=16. The local expression is

    F(c,h) = product_i f_i - sum_{|S|>=2} h_S product_{i not in S} f_i.

Removing one plane hit c_j changes F by

    product_{i!=j} f_i
       - sum_{|S|>=2, j not in S} h_S product_{i not in S, i!=j} f_i
    <= product_{i!=j} D_i = B_j.

Removing one tail hit h_S changes F by product_{i not in S} f_i,
which is at most B_S=product_{i not in S} D_i. All factors and tail
counts are nonnegative. Hence for any coordinate g, decreasing its count
by k increases F by at most k B_g, even after other counts are decreased.

## Integer quantile repair

Let a point distribution have total integer mass m, demand survival counts
d_g(t), and allowed survival counts s_g(t) supplied by an actual group
assignment mixture. Couple the demand count and supply count in decreasing
quantile order, and replace the demand count by their minimum. This only
decreases that coordinate and leaves every other coordinate unchanged.

For threshold t, the quantile intervals of lengths d_g(t) and s_g(t) are
nested. The repaired survival count is therefore min(d_g(t),s_g(t)). The
number of removed hits, counted with their integer masses, is exactly

    Delta_g = sum_t max(d_g(t)-s_g(t),0).

Apply this independently to every coordinate. Changing coordinate g leaves
the marginals of the other coordinates intact, so earlier repairs remain
valid and the same formula applies at every step. The result satisfies all
threshold constraints and its cost increases by at most sum_g B_g Delta_g.
Splitting integer weights at quantile endpoints needs no floating arithmetic.

## Consequence for point optimization

Replace each threshold constraint demand_g(t)<=s_g(t) by

    demand_g(t) <= s_g(t) + slack_g(t),   slack_g(t)>=0,

and add sum_{g,t} B_g slack_g(t) to the objective. Every feasible original
distribution is feasible with zero slacks. Conversely, the preceding repair
turns any penalized distribution into an original feasible distribution
whose objective is no greater than its penalized objective. Thus the two
complete point LPs have exactly the same optimum.

In the dual, the new slack variables impose -B_g<=lambda_g(t)<=0. These
bounds therefore lose no optimal value. They prevent unconstrained large
negative event prices at zero-supply thresholds, the main source of weak
cuts in the first decomposed calculation.

The same proof also justifies converting intermediate, incomplete-column
point solutions with nonzero slacks into feasible primal witnesses. It does
not justify a lower bound from a restricted state pool: full integer pricing
is still required for every supporting inequality.

## Implementation and scope

`code/coherent_transport_repair.py` implements the one-hit bounds and integer
quantile repair. Its witnesses go through the existing separate standard-
library auditor `code/audit_coherent_groups.py`. The generator uses capped
duals in `code/decomposed_coherent_groups.py`; integer pricing still searches
every local plane tuple and optimizes all tail counts exactly.

Tests compare the capped point optimization with an independently constructed
full-state LP, exhaustively check one-hit inequalities for a two-prime model,
and verify quantile removal counts in 100 randomly generated integer examples.
These checks supplement the algebraic proof above; they are not a proof of
the global covering conjecture.

For a deterministic actual group assignment, the initial point problem can
also be solved at plane box vertices: F is affine in each plane count, and
is nonincreasing in every tail count. Use each tail's full available count
and, for each plane, either zero or its clipped available count. Adding these
at most 2^5 vertices supplies the known feasible value 670035 immediately.
