# Finite-support point transport and globally valid dual extensions

This is an exact reduction for the existing coherent-group relaxation at a
fixed group mixture. It does not assert that the global covering problem or
the nine-prime instance has been settled. Numerical solutions remain subject
to rational reconstruction and independent audit.

## Restricting plane counts to supplied endpoints

Fix the group mixtures and a core point x. The supplied count distribution
for group h has a finite set of atoms. Clip a plane supply at its available
higher-coordinate count D_h. Its relevant endpoints are zero and its atom
counts after clipping.

Between two consecutive endpoints a<b, the supplied survival probability is
constant on thresholds a+1,...,b. Conditional on every other coordinate,
replace an intermediate demand count k by a or b with probabilities preserving
its mean. F is affine in that plane coordinate, so the expected cost does not
change. At every threshold in this gap, the new demand is at most the old
demand at threshold a+1, hence at most the common supplied probability.
Demand outside the gap is unchanged. Repeating for all plane coordinates
preserves the marginals of all other coordinates and all domination bounds.

Thus the point optimum occurs on the Cartesian product of these endpoint
sets. For a sparse mixture this can be much smaller than the full plane grid.

## Tail transport

For a fixed plane state j of mass m_j, let f_i=D_i-c_i and
v_jS=product_{i not in S} f_i. The coefficient of a tail count h_S is -v_jS<=0.
Tail supply below its smallest atom need never be used: increasing such a
count to that atom preserves domination and cannot increase F.

Write the tail atom counts as a_0<...<a_k. Introduce
y_jSl=P(plane state j and tail count >=a_l), for l>=1, with

    0<=y_jSk<=...<=y_jS1<=m_j,
    sum_j y_jSl <= supplied_survival_S(a_l).

The tail contribution is

    -v_jS [a_0 m_j + sum_l (a_l-a_(l-1)) y_jSl].

For each plane state, place the nested threshold events on initial intervals
[0,y_jSl] in [0,m_j]. The union of their endpoints, across all tails, constructs
one joint distribution realizing every transport variable. Therefore this
linear model is exact, not just a relaxation of the fixed-mixture point LP.

For each tail S, plane states with the same coefficient v_jS can be combined.
Let M_v be the sum of their masses. Use one transport chain bounded by M_v
for that value v. Split its resulting threshold masses back among the plane
states in proportion to m_j/M_v. This preserves every global tail marginal
and the objective, because the coefficient v is constant within the group.
It also realizes compatible conditional tail distributions at each plane
state, which the nested-interval construction then combines. Thus this
aggregation is exact. Zero-volume groups need no transport variables.

## Extending its dual to all count states

The reduced LP has a mass dual beta and nonpositive multipliers on global
survival bounds. Let P_h(c) be their cumulative sum at each count endpoint.
For planes, zero is an endpoint. For tails, assign potential zero from zero
through the smallest atom, which has no nontrivial threshold constraint.

Interpolate P_h linearly between consecutive endpoints. Then F-sum_h P_h
is affine in each coordinate separately on every endpoint box, and its minimum
on that box occurs at a vertex. All vertices obey the reduced LP dual bound.

Beyond the largest supplied count, extend P_h with slope -B_h, where the
one-hit lemma gives F(c+1)-F(c)>=-B_h. Consequently F-P_h is nondecreasing
beyond that endpoint. Below a tail's smallest atom, P_h is zero and F is
nonincreasing in the tail count, so the minimum is at the smallest atom.
These facts prove the extended dual inequality on the entire count domain,
including count values absent from the current mixture.

For a plane the admissible domain ends at min(14,D_h). Event multipliers above
that limit must be zero, even if the un-clipped group supply contains larger
counts. Equivalently, the potential is constant beyond the clipping limit.
An initial implementation incorrectly extended the slope through those
impossible thresholds, producing valid but needlessly weak inequalities.
The fixed-mixture boundary regression detects this case and verifies dual
attainment after correction.

The extension does not change the expected potential under the supplied
distribution. Hence it attains the same point optimum before numerical
rounding. This gives a supporting inequality tight at the current mixture
and valid for every other mixture, which is what the Benders master needs.

## Implementation and verification

`code/coherent_point_transport.py` solves the reduced numerical LP and
reconstructs joint count states by nested intervals. It also constructs the
extended event potentials. Multipliers are rounded to denominator 10^6 and
then the entire point-state domain is repriced with the existing integer
routine. Thus the saved supporting inequality remains valid independently
of numerical dual accuracy or the endpoint interpolation argument.

The numerical generator skips this reduction when its configured state or
variable budget is exceeded, or its short LP solve does not finish. It then
uses the existing full-state column-generation oracle. Small examples with
three or more supplied count levels agree with independent full-state LPs,
including all intermediate counts omitted by the reduction. Their extended
duals also attain the same values under exhaustive integer point pricing.

The separate `code/audit_coherent_groups.py` recomputes every resource hit,
threshold domination constraint, raw cost, and positive part of each point
expectation from the rational witness. No numerical point value is itself
treated as a covering obstruction.
