# An all-schedule barrier below 7 in the original distortion range

Date: 26 October 2026.

Later update: the larger parameter range is now also excluded for the same
fixed catalogue comparison. For all 0<=delta_p<1, the original objective
through prime 10007 is at least 1.029993468. See
`stoploss-full-range-barrier-2026-10-26.md`. This older half-range certificate
remains valid; its statement about unexamined larger distortions records
the scope of that certificate, not the current research status.

## Exact scope

Fix the squarefree catalogue comparison used in the upper-bound-7
certificate: increasing prime order, the independent Bernoulli comparison
with probabilities r_p=1/[p(1-delta_p)], and the minimum of first-moment,
positive-part, and quadratic stage bounds. Allow every prime to have its
own real parameter in the entire interval 0<=delta_p<=1/2.

This is the parameter interval specified in Section 3 of
[Cummings--Filaseta--Trifonov](https://arxiv.org/html/2211.08548v1#S3)
and enforced by both existing covering-certificate verifiers.

For cutoff C0=6, the sum of these catalogue stage bounds through prime
5003 is at least

    1010143243 / 10^9 = 1.010143243 > 1

for **every** such parameter schedule. Thus the existing comparison and
parameter range cannot certify a uniform minimum-modulus upper bound 6,
or any smaller cutoff, using the total-loss-<1 criterion. Combined with
the successful bound-7 certificate, this makes 7 the best integer upper
bound attainable by this specified method.

This is a lower bound on an upper-estimate expression, not a lower bound
on the actual covered mass, and not the construction of a covering.
It does NOT exclude delta_p>1/2, a different prime ordering, different
comparison distributions, actual-family-specific bounds, or bounds using
compatibility and overlap information across stages. The minimum-modulus-2
conjecture is still open. These distinctions are part of the theorem.

## 1. Reduction to pure stop-loss

The quadratic bound is at least the stop-loss bound, pointwise in the
comparison variable. If the first moment is used at a stage with positive
delta, replace that delta by zero. Its new stop-loss is the first moment;
all later Bernoulli probabilities are no larger, and every later stage
functional is increasing in the count. Replacing all stages that selected
the first moment therefore produces a pure stop-loss schedule with no
larger total. Consequently the infimum of the minimum-of-three objective
equals the infimum of pure stop-loss.

No allowed modulus occurs at primes 2 or 3 when C0=6. Their distortions
can be zero without increasing any later cost. Hence the initial count
law is (1/3,1/2,1/6) on counts (0,1,2).

For p=5 put A(n)=2^n-1; for every p>=7 put A(n)=2^n. With
nu=1/(1-delta) and r=nu/p, the stage cost is

    E g_p(N,nu),   g_p(n,nu) = max(0, 1+nu*(A(n)/p-1)).

## 2. Only 60 initial cases are needed

Keep all other controls fixed and vary one control r. The current stage
cost is convex and piecewise affine in r; the expected sum of all later
costs is affine in r, because adding a single independent Bernoulli
variable is affine in its probability. A coordinate minimum occurs at a
hinge breakpoint or an endpoint. The finite-horizon objective is continuous
on a compact product of intervals. Applying this observation successively
to its first three active controls yields a schedule no worse than the
original with

    delta_5  in {0,1/5,1/2},
    delta_7  in {0,1/7,2/7,1/2},
    delta_11 in {0,1/11,2/11,4/11,1/2}.

Thus 3*4*5=60 cases exhaust the possible improvement in those controls.
This does not discretize any of the later controls: they remain continuous.
For each case, the initial cost and the independent count law through
prime 11 are calculated with exact fractions.

## 3. A lower comparison distribution at fixed total probability

Let q_1<...<q_k be the earlier primes at least 13. Their Bernoulli
probabilities obey l_i=1/q_i <= r_i <= u_i=2/q_i. Put

    x = sum_i(r_i-l_i).

At fixed x, form v(x) by starting at all l_i and adding the excess x to
the first coordinate until it reaches u_1, then to the second, and so on.
At most one coordinate is partially filled.

Both l_i and u_i decrease with i, so v(x) is decreasing. Its first j
entries sum to

    min(sum_{i<=j}u_i, sum_{i<=j}l_i + x).

For ANY admissible r, the sum of its j largest entries is at most this:
it is at most the sum of the j largest upper bounds, and at most the
total sum minus the sum of the k-j smallest lower bounds. Therefore
v(x) majorizes r.

For independent Bernoulli variables, making a fixed-sum pair of
probabilities more unequal lowers the expectation of any convex function
of their sum. Indeed, conditional on the other count W, the part that
depends on the product ab is

    ab * [phi(W+2)-2phi(W+1)+phi(W)],

whose coefficient is nonnegative. A mass-spreading pair transfer reduces
ab. Successive such transfers give the majorization comparison. Adding
the same independent prefix count preserves the inequality. The functions
n -> g_p(n,nu) are convex, so replacing the unknown earlier probabilities
by v(x) gives a valid LOWER stage cost.

This is the elementary Bernoulli convex-order principle underlying
Hoeffding's work on independent trials. The proof above supplies the
particular direction and the unequal coordinate bounds needed here;
the numerical certificate does not rely on an external theorem lookup.

## 4. A one-dimensional relaxed control problem

After the next prime p, the total excess changes by

    x_next = x + (nu-1)/p,   1<=nu<=2.

Replace the stage distribution at every step by the extremal distribution
from Section 3 with that x. This can only lower every schedule's cost, even
though these extremal distributions need not be simultaneously realizable
by one history. It is deliberately a more favorable problem for a
hypothetical successful schedule.

Its exact optimal continuation value is nondecreasing in x: the extremal
probabilities v_i(x) are nondecreasing in x, the stage functions are
increasing in the count, and the same future controls are available.

The primary certificate rounds x DOWN to j/H, H=4000. Split the entire
continuous control interval into intervals for which

    floor(H*(nu-1)/p) = h.

For such an interval, any next state is at least (j+h)/H. A lower bound is
therefore the independently minimized present cost on that whole interval
plus the saved continuation lower bound at j+h. The present expectation
is piecewise affine in nu; its minimum is attained at an interval endpoint
or a hinge breakpoint nu=p/(p-2^n). All are enumerated with integer
comparisons. This proves that the dynamic program covers continuous
controls, not just a sampled parameter grid.

The terminal continuation after prime 5003 is zero. All omitted later
stage costs are nonnegative, so no infinite-tail approximation is involved
in this impossibility certificate.

## 5. Arithmetic and result

`code/certify_stoploss_half_barrier.py` uses probability scale 10^9,
mixing scale 10^9, grid denominator 4000, and the first 20 terms of the
unknown count distribution. Every arithmetic and state approximation
rounds DOWN. The omitted count probabilities contribute nonnegative cost
and are simply discarded. Exact fractions determine grid locations,
prefix cases, and interpolation coefficients.

Each extremal law is a convolution of an upper-probability prefix and a
lower-probability suffix. Their probability polynomials are calculated by
integer recurrences. Interpolation uses lower enclosures of both mixing
weights. Combining these subdistributions with the exact prefix law and
nonnegative stage costs preserves the lower-bound direction. Explicit
moment and scale bounds prevent int64 overflow.

All 60 cases exceed 1. The weakest case is

    (delta_5,delta_7,delta_11) = (1/5,2/7,4/11).
    Prefix count law = (144,324,260,95,16,1)/840.
    Exact prefix loss = 31/105.
    Prefix lower at scale 10^9 = 295238095.
    Continuation lower at scale 10^9 = 714905148.
    Total lower = 1010143243 > 10^9.

The approximately 1.107 numerical candidate found at the start of this
investigation is not used in the proof. Neither solver optimality nor
floating-point tolerances are trusted.

## Reproduction and checks

```text
python code/certify_stoploss_half_barrier.py data/stoploss-half-barrier.json
python code/audit_stoploss_half_barrier.py data/audit-stoploss-half-barrier.json
python -m unittest discover -s code -p test_stoploss_barrier.py
```

The separate audit imports neither the primary generator nor the existing
covering verifiers. It enumerates the initial subsets and controls afresh,
uses grid denominator 5000 and probability scale 10^8, retains 18 unknown
count terms, and convolves the full prefix law before evaluating costs.
It independently gives the stronger lower bound

    101255014 / 10^8 = 1.01255014 > 1

over all 60 cases. The two numbers need not coincide because the grids,
truncations, and scales differ; both are rigorous lower bounds for the same
continuous optimization problem. `data/stoploss-half-barrier-summary.json`
records both exact values, the scope, and hashes of the programs and results.

Additional exact tests compare the small-prime case enumerations, verify
the convex-order comparison on rational three-coordinate examples, and
check all 60 lower bounds through prime 17 against exact exhaustive
coordinate minima. The generalized upper-bound search and its existing
tests also remain available.
