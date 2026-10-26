# The bound 7 remains optimal when distortions above 1/2 are allowed

Date: 26 October 2026.

## The precise conclusion

Keep the independent Bernoulli squarefree catalogue comparison, increasing
prime order, and the total-stage-bound-less-than-1 criterion used in the
upper-bound-7 certificate. Allow every real parameter 0<=delta_p<1, with
the natural probability cap

    r_p = min(1, 1/[p(1-delta_p)]).

At cutoff C0=6, for every such schedule, the sum of the catalogue estimates
through prime 10007 is at least

    1029993468 / 10^9 = 1.029993468 > 1.

This includes taking the minimum of the first-moment, stop-loss, and
quadratic estimates at every stage. A separately implemented audit gives
the stronger lower bound 103219703/10^8 = 1.03219703.

Thus allowing delta_p>1/2 cannot improve the uniform minimum-modulus bound
below 7 within this specified comparison. The existing bound-7 certificate
and this lower certificate establish optimality at 7 for the full parameter
range. The numerical lower bounds above are not claimed to be the exact
minimum of the parameter optimization problem.

This is a limitation of an upper-estimate expression, not a construction of
a covering and not a lower bound on its actual covered mass. It does not
exclude another comparison distribution, different prime ordering,
conditional distortions, family-specific information, or bounds using
actual congruence compatibility and overlap. The modulus-2 conjecture
remains unresolved.

## 1. The extension to all delta<1 is well defined

In a fibre with covered fraction alpha, the distorted density multiplier is
zero on the covered part and 1/(1-alpha) on the uncovered part when
alpha<=delta. When alpha>delta, it is

    (alpha-delta)/[alpha(1-delta)] on the covered part,
    1/(1-delta) on the uncovered part.

For every 0<=delta<1 these are nonnegative, preserve the fibre's total mass,
and are bounded by 1/(1-delta). The degenerate fibres alpha=0 and alpha=1
have their evident values (all uncovered or all covered). Consequently the
conditional atom bound is min(1,1/[p(1-delta)]), and the usual conditional
convex comparison applies with Bernoulli variables of those probabilities.
The exact stage loss is E[(alpha-delta)_+]/(1-delta).

The first-moment and quadratic alternatives remain valid for this whole
range. The quadratic bound is always at least the stop-loss bound. Whenever
a schedule selects the first moment at a stage, changing that stage's
delta to zero turns its stop-loss into the first moment and weakens every
later Bernoulli probability. Replacing every such choice gives a pure
stop-loss schedule with no larger total. It therefore suffices to lower-bound
all pure stop-loss schedules on the finite horizon.

At cutoff 6 no allowed modulus occurs at primes 2 or 3, so setting their
distortions to zero only improves the estimate. The initial count law is
(1/3,1/2,1/6). At prime 5 put A_5(n)=2^n-1, and at all later primes put
A_p(n)=2^n.

## 2. A complete finite control set for the unrestricted range

Let a_p be the largest possible A_p(n) strictly less than p: a_5=3, and
a_p is the largest power of two below p when p>=7. Once delta>=a_p/p,
every term with A_p(n)<p is zero. Every remaining term, with A_p(n)>=p,
is nondecreasing in nu=1/(1-delta). Reducing delta to a_p/p therefore
does not increase the current loss or any later comparison probability.
This reduction also covers the region where the original probability bound
was capped at 1.

We may hence restrict to delta<=a_p/p, where r=1/[p(1-delta)]<=1. The
stage charge is

    g_p(n,r) = max(0, 1+r*(A_p(n)-p)).

With the other controls fixed, the present charge is piecewise affine in r,
and every future expected charge is affine in r. Its coordinate minimum
is at a breakpoint or endpoint. Successively replacing the controls of any
finite schedule yields a schedule no worse with the atoms

    p=5: a in {0,1,3},
    p>=7: a in {0,1,2,4,...}, a<p,
    delta=a/p, r=1/(p-a).

For such an atom the stage charge is simply

    g_p(n,a) = (A_p(n)-a)_+/(p-a).

No real delta in [0,1) has been omitted by these reductions.

## 3. Certified action dominance, without imposing a half-range cutoff

Take horizon P=10007. At a stage p let a be the largest atom and b the
preceding atom. For a fixed earlier count n and fixed future controls, the
difference between the full objective using a and that using b equals

    (1/(p-a)-1/(p-b)) * [d_p(n)+D_{>p}(n)],

where

    d_p(n) = (A_p(n)-p) * 1{A_p(n)>=a},
    D_{>p}(n) = F_{>p}(n+1)-F_{>p}(n).

Here F_{>p} is the expected remaining cost through P. The first expression
is the exact slope of the current piecewise-affine cost between consecutive
breakpoints. The second is the exact change in future cost when the Bernoulli
variable at p changes from zero to one.

For a later stage q define the nonnegative marginal cost

    h_q(n,a) = g_q(n+1,a)-g_q(n,a).

A backwards dynamic program lower-bounds every D_{>p}. It even lets the
future controller choose its action after seeing the count n, a freedom
not available to the original fixed schedule. With terminal value zero,

    W_q(n) = min_a [ h_q(n,a)
                  + (1-1/(q-a))*W_next(n)
                  + (1/(q-a))*W_next(n+1) ].

The count is capped at 24 in the primary calculation. This remains a lower
bound because h_q is increasing in n. At that cap, d_p(n) is also increasing,
so capping the terminal expression used below is safe. The cap is above
every negative part of d_p for this horizon.

Now lower-bound E[d_p(N)+W_{>p}(N)] over all earlier controls. A second
backwards recursion again permits state-dependent controls, using the two
endpoints of each earlier probability interval. It uses signed integer
floor divisions, including when the terminal expression is negative.
If the resulting lower bound is nonnegative, atom a is dominated by b for
all the other currently permitted controls, and is removed.

This calculation certifies 326 removals. Importantly, the comparison is
for the FULL objective through 10007: a replacement can raise a shorter
prefix while reducing its continuation. All conclusions below retain this
distinction.

Many larger distortions remain permitted. For example delta_7=4/7 is
retained, as is delta_4861=4096/4861. The proof does not replace the original
problem by the earlier delta<=1/2 problem.

The first three surviving control sets are

    delta_5  in {0,1/5},
    delta_7  in {0,1/7,2/7,4/7},
    delta_11 in {0,1/11,2/11,4/11}.

They produce 32 cases. No global optimizer's reported optimality or
floating-point dominance test is trusted.

## 4. Lower-bound all remaining controls by mean-budget majorization

For each later prime q, use the interval

    1/q <= r_q <= 1/[q-largest_remaining_atom(q)].

Allowing the entire interval is a relaxation of the surviving discrete
control set. Before a stage p, let l_i be the decreasing list of the lower
bounds and u_i the decreasing list of the upper bounds, sorted separately.
Then u_i>=l_i. Put x=sum(r_i-l_i).

Starting at l, fill the first coordinate up to u_1, then the second, and
so on, using total excess x. The resulting decreasing vector v(x) has
prefix sums

    min(sum_{i<=j}u_i, sum_{i<=j}l_i+x).

These dominate the sums of the j largest entries of every admissible
probability vector: one bound follows from the j largest upper bounds;
the other follows by subtracting the smallest lower bounds of the
complement. Sorting the upper and lower bounds separately can introduce
additional unattainable configurations, but only makes the lower relaxation
more favorable.

The Bernoulli convex-order argument from the half-range proof therefore
applies. At fixed x, the sum of independent Bernoulli variables with
probabilities v(x) gives a lower expected stage cost. The exact initial
count law in each of the 32 cases is convolved with this law.

The excess state changes by x_next=x+r_p-1/p. It is rounded DOWN to a
grid of denominator 4000. Control intervals cover every possible value
of floor(4000*(r_p-1/p)); the current piecewise-affine cost is minimized
over the entire interval, including all hinge breakpoints, and the future
value is evaluated at the lower next grid state. The exact relaxed value
is increasing in x, which justifies this direction of state rounding.

The terminal continuation after prime 5003 is zero. Thus this part proves
a lower bound on the prefix of each PRUNED schedule. Since pruning never
increases the full objective through 10007 and all its remaining tail costs
are nonnegative, the result is a lower bound on the ORIGINAL full objective
through 10007.

## 5. Exact arithmetic and separate verification

Primary action elimination:

* Integer scale 10^12, count cap 24.
* 326 eliminations, all with positive certified lower bounds.
* Smallest certified dominance bound: 0.000008779191.

Separate action audit:

* Integer scale 10^10, count cap 28.
* Different rounding decomposition for current and future marginal costs.
* All 326 eliminations verified; smallest audit lower bound 0.0000087317.

Primary final lower bound:

* Probability and mixing scales 10^9, state grid denominator 4000.
* 20 retained terms for the unknown count; omitted nonnegative contributions
  are discarded.
* Directed rational-sum enclosures at scale 10^18 determine interpolation
  weights, with exact-fraction fallback at ambiguous boundaries.
* Explicit integer-overflow checks for the vectorized arithmetic.
* Minimum over all 32 cases: 1029993468/10^9 > 1.

Separate final audit:

* Probability scale 10^8, state grid denominator 5000, 18 count terms.
* Explicit small-prime subset enumeration.
* Capacity breakpoints padded UP at scale 10^12, producing a stochastically
  smaller probability vector, instead of the primary interpolation enclosure.
* Convolution with the full prefix law before cost evaluation.
* Minimum over all 32 cases: 103219703/10^8 > 1.

The weakest case in both calculations is
(delta_5,delta_7,delta_11)=(1/5,2/7,4/11), with prefix loss 31/105 and
count law (144,324,260,95,16,1)/840. The primary continuation lower bound is
734755373/10^9. Different grids and rounding rules explain the different
final lower bounds; they need not agree exactly.

## Reproduction

Run from the repository root:

```text
python code/certify_fullrange_action_pruning.py data/fullrange-action-pruning.json
python code/audit_fullrange_action_pruning.py data/fullrange-action-pruning.json data/audit-fullrange-action-pruning.json
python code/certify_fullrange_stoploss_barrier.py data/fullrange-action-pruning.json data/fullrange-stoploss-barrier.json
python code/audit_fullrange_stoploss_barrier.py data/fullrange-action-pruning.json data/audit-fullrange-action-pruning.json data/audit-fullrange-stoploss-barrier.json
python -m unittest discover -s code -p test_fullrange_stoploss.py
```

The pruning programs use standard-library integer arithmetic. The final
barrier programs use NumPy only for bounded integer array operations. No
solver is needed. The separate audits do not import the corresponding
primary certificate generators.

The exact tests check large-delta domination, marginal continuation bounds
against exhaustive small rational schedules, separately sorted unequal
probability bounds, both grid constructions against exact rational vectors,
and independent enumeration of the 32 prefix cases.

The compact manifest `data/fullrange-stoploss-barrier-summary.json` records
the exact conclusions, certificate dependencies, and file hashes.
