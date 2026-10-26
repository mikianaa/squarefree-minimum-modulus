# Quantitative residual bounds retaining core intersections

The global minimum-2 conjecture remains unresolved. A new integer-only chain
gives u(2310)>=119, u(30030)>=501, and u(510510)>=510. Here u(N) is the minimum
number of uncovered residues modulo N among families with one congruence for
every divisor of N greater than 2. Missing congruences can be added, so the
lower bounds apply to incomplete families too.

The final lower bound has an independent exact dual certificate. It gives
another proof of the existing seven-prime obstruction, with a quantitative
margin. It does not exclude eight primes. The original plane-load proof is
preserved in `seven-prime-obstruction-2026-10-25.md`.

An independently checked actual family gives u(510510)<=7436, with residual
parity counts (2779,4657). Thus the current verified bracket is

    510 <= u(510510) <= 7436.

## Quantitative use of the core catalogue

Use the previously audited normalization and 210-core pair masks. For a pair
residual S, record its two layer sizes and the maximal intersections C_d(e)
for odd d in {3,5,7,15,21,35,105}. The remaining classes 105 and 210 remove at
most one point per layer and one additional point in one selected layer. This
gives two lower size pairs l_e; replace each capacity at this lower pair by
min(l_e,C_d(e)). Keep every raw profile. The normalized dominance used in the
earlier seven-prime proof cannot be used for absolute residual counts.

For a fixed allocation of the seven even classes 2qd, the resulting layer
deficit has the form

    max(0, D*r - sum_i min(r,C_i)),                     (1)

where D=q-1-[b_q=e]. Before clipping, this is convex and equals zero at zero.
Its positive part is nondecreasing. The actual core residual has size at least
l_e and its capacities are bounded by the original C_i. Therefore evaluating
(1) at l_e gives a valid quantitative bound even if the actual size is larger.

There is also a shared-label restriction for the odd class 3q. Since b_3=0,
every core point in parity 0 has label 2 modulo 3. Write u,v for the numbers
of parity-1 pair-residual points at labels 1,2. The odd 3q class has capacity
pair bounded by either

    (0, min(r_1,u))  or  (r_0, min(r_1,v)).

The two independent marginal maxima cannot be granted simultaneously.
Label 0 meets no core point and can be discarded. Both choices preserve the
convex monotonicity argument: in the second branch the contribution r_0 simply
reduces D by one. The even 6q class still chooses one parity independently.

Enumerate these two odd-label choices and all 128 even allocations. Pareto
pruning now uses componentwise absolute lower pairs. Intersect their upward
closures with the independently audited compression frontiers from the previous
note. The numbers of enhanced raw profiles for core patterns 000,001,010,011
are 4026,4480,5137,4876 respectively.

At q=11 the resulting total lower bounds are:

| Pattern | 0000 | 0001 | 0010 | 0011 | 0100 | 0101 | 0110 | 0111 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Lower | 186 | 169 | 150 | 133 | 131 | 122 | 121 | 119 |

Applying the exact compressed parity extension through 13 gives global lower
501 at 30030. Including earlier solver optima raises it to 510, but the chain
used below uses only the integer-only rows. A separate ordinary-residue and
histogram implementation verifies all 40 frontiers in 133.797 seconds.

## A scalable extension bound

Let the two current anchor boxes have coordinate lengths D_i(e), and let
r_0,r_1 be a lower pair. For each layer, let F_e(r_e) be the sharp maximum of
the sum of all nonempty cylinder capacities, given by coordinate compression
in `compressed-residual-capacity-2026-10-25.md`. Each individual even shape S
has capacity bound

    c_S(e) = min(r_e, product_{i not in S} D_i(e)).

All odd qd classes together have capacity at most F_e. If the even shapes
assigned to layer e have total individual allowance s_e, their actual capacity
sum is at most min(F_e,s_e). Therefore a lower deficit pair is

    U_e >= max(0, (q-1-[b_q=e])*r_e - F_e - min(F_e,s_e)).       (2)

For a fixed allocation, F_e and s_e are concave in r_e and vanish at zero.
Their minimum is concave too. Thus the clipped expression in (2) is
nondecreasing, justifying propagation of componentwise lower pairs.

A knapsack program allocates each even shape to one layer. It uses only
2^k-1 items for k odd prime coordinates, rather than enumerating all
2^(2^k-1) assignments. Its computed global lower at 510510 is 517. The full
knapsack frontiers have small exhaustive regression checks; the independently
audited theorem below uses the slightly weaker bound 510 instead.

## Independent exact dual certificate

Put D_e=q-1-[b_q=e], A_e=D_e*r_e-2F_e, and

    B_e = F_e + min(A_e,0),     K = sum_e max(A_e,0).

For any t_0,t_1 in [0,1], every allocation satisfies

    U_0+U_1 >= K + t_0*B_0+t_1*B_1
                   - sum_S max(t_0*c_S(0),t_1*c_S(1)).        (3)

For A_e>=0, take a convex combination of U_e>=A_e and
U_e>=D_e*r_e-F_e-s_e. For A_e<0, combine the latter inequality with U_e>=0.
This gives U_e>=max(A_e,0)+t_e*(B_e-s_e). Each even resource contributes at
most the maximum of its two weighted capacities, proving (3).

The independent audit constructs F_e by explicitly sorting point weights in
the anchor boxes. It tries rational t pairs obtained from individual capacity
ratios and the endpoints, then rounds each valid lower bound upward because
the actual uncovered count is an integer. It checks 11332 input pairs across
all 32 anchor patterns at 510510 in 1.532 seconds. Every resulting integer
lower is at least 510.

The weakest case is pattern 011110, starting with counts (322,208) at 30030.
The weights (t_0,t_1)=(208/288,1) give the lower bound

    146656/288 = 4583/9 > 509,

so at least 510 points remain. These dual certificates do not rely on the
knapsack solver's optimality or on the earlier CP-SAT optima.

## What remains unresolved

At the next prime 19, the knapsack relaxation accepts zero deficit from the
prefix lower pair (102,470), for both choices of b_19, in pattern 011111.
This is a feasible capacity relaxation, not an actual residual or covering.
The computation therefore does not exclude eight primes. New restrictions on
realizable higher-dimensional residuals remain necessary for this route.
Explicit even-resource assignments realizing zero deficit in this relaxation
are stored in `data/hybrid-9699690-relaxation-witnesses.json`; their capacity
sums are checked directly, independently of the knapsack objective value.

Two exploratory weighted CP-SAT optimizations at 2310, pattern 0111, ran for
60 seconds each. They returned feasible witnesses but lower bounds 0 and 2;
they are not inputs to the proof. Their independently counted weighted upper
values were 176 for weights (1,2) and 153 for weights (2,1).

The literature check did not supply an unrestricted density inequality for this
purpose. Theorem 1.9 of [Harrington, Klein, Lowrance and Trifonov](https://arxiv.org/html/2605.18644v1)
assumes at most three prime divisors; their Problem 4 asks about broader
conditions. The [BBMST theorem](https://arxiv.org/html/1901.11465) supplies the
even-modulus reduction, not the presence of modulus 2.

## Files and reproduction

The previously audited core and compressed-frontier certificates are inputs.
After reproducing those as described in their notes, run:

    python code/core_capacity_frontier.py
    python code/audit_core_capacity_frontier.py
    python code/hybrid_slice_frontier.py --output data/hybrid-frontier-510510.json
    python code/audit_hybrid_slice_bound.py
    python code/audit_partial_families.py

The two producer scripts require NumPy; the three audit scripts use only the
standard library. Their outputs are `data/core-capacity-frontiers.json`,
`data/audit-core-capacity-frontiers.json`, `data/hybrid-frontier-510510.json`,
`data/audit-hybrid-slice-510510.json`, and `data/audit-partial-families.json`.
All 34 repository regression tests pass in 62.384 seconds. No proof-assistant
verification or novelty claim is made.
