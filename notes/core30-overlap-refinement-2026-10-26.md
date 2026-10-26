# Retaining the uncovered core: modulus 10 is necessary at minimum 7

Date: 26 October 2026.

## Result and scope

Every finite covering system with distinct squarefree moduli, all greater
than 6, must contain modulus 10. Together with the existing upper bound 7,
any remaining system of minimum modulus at least 7 must contain **both
moduli 7 and 10**. This does not prove the upper bound 6.

The new comparison has certified total loss below 0.989675744340. The
primary integer calculation gives 0.98967574433903727787657..., and a
separate implementation gives 0.98967574433903962056... after converting
to the same normalization. Both include the infinite prime tail. Their
difference is less than 2.4e-15, from different count truncations and
outward rounding.

This changes the comparison distribution using the actual uncovered
set. It is outside the previously certified all-parameter barrier.

## Complete reduction of the initial configurations

Suppose a covering has no modulus <=6 and no modulus 10. Among moduli
whose prime factors belong to {2,3,5}, only 15 and 30 are allowed.
Add either if absent. If the class modulo 30 is contained in the class
modulo 15, it is redundant and may be moved outside that class without
destroying the covering. Thus these classes may be assumed to cover
exactly three points of the 30-point CRT core.

Write the core as Z_2 x Z_3 x Z_5. Normalize the modulus-15 class to
x_3=x_5=0. Coordinate permutations reduce the extra point modulo 30 to
(0,0,1), (0,1,0), or (0,1,1). Arbitrary coordinate permutations preserve
all coordinate cylinders; they need not be affine. Let U be the 27
uncovered points and start with the uniform probability on U. The
initial progressions have zero mass.

## Sharper comparison, with an elementary upper proof

For every divisor d of 30 choose a residue a_d modulo d and set

    f(x) = sum_{d|30} 1[x = a_d (mod d)].

The d=1 term is 1, so 1<=f<=8. Missing cylinders may be added in an
increasing upper bound. For H(j)=max_(a_d) sum_(x in U) (f(x)-j)_+,
the following table is valid:

| j | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H(j) | 68 | 41 | 20 | 12 | 4 | 3 | 2 | 1 | 0 |

An elementary proof of the upper bounds avoids enumeration:

* A cylinder for d=1 contains 27 points of U. A parity cylinder contains
  at most 14, since the removed points have parity counts (2,1). For
  divisors 3,5,6,10,15,30 use the full-core sizes 10,6,5,3,2,1. Hence
  sum_U f <=27+14+10+6+5+3+2+1=68.
* Since f>=1, sum_U(f-1)=sum_U f-27<=41.
* For j>=2 enlarge U to the full core. The original uniform-coordinate
  product comparison gives 2^K with unnormalized K-masses (8,14,7,1).
  Its stop-loss sums for j=2,...,8 are (20,12,4,3,2,1,0).

The polygon through these values is the stop-loss transform of the new
normalized comparison variable F:

| F | 1 | 2 | 4 | 8 |
| --- | --- | --- | --- | --- |
| Probability | 6/27 | 13/27 | 7/27 | 1/27 |

Thus E_U phi(f)<=E phi(F) for every convex nondecreasing phi. Indeed,
expand phi on {1,...,8} into its value at 1, its nonnegative first
difference times (f-1)_+, and its nonnegative second differences times
(f-j)_+ for j=2,...,7. Both measures have total mass 1. Linear
interpolation handles noninteger stop-loss thresholds.

The independent integer enumeration `audit_core30_no10.py` checks all
810000 cylinder catalogues in each of the three configurations and
rechecks maximizing witnesses point by point. Every entry of the table
is sharp. Only the elementary upper proof is needed for the theorem.

## Transfer to later primes

At a prime p>=7, apply the usual conditional product comparison only to
earlier coordinates q>=7, with caps r_q=1/[q(1-delta_q)]. Keep the core
intact. For N successful comparison coordinates the numerator becomes

    sum_{J subset of successful coordinates} f_J(x).

There are t=2^N terms, each bounded by a core catalogue. Jensen gives

    E_U phi(sum_J f_J(x)/p)
      <= (1/t) sum_J E_U phi(t f_J(x)/p)
      <= E phi(t F/p).

The valid comparison is therefore F*2^N/p, with F independent of the
tail Bernoulli count. Initialize the count with masses (6,13,7,1)/27
and start adding primes at 7. The initial moments are

    E 2^K=68/27, E 3^K=5, E 4^K=234/27.

The first-moment, stop-loss, quadratic, and CFT infinite-tail bounds
then apply. Independence is asserted only for the comparison variables,
not for the actual coordinates of U.

## Arithmetic certificates

The schedule has denominator 10^6, delta<=1/2, exact distribution range
through 3000000, and 1000000 finite primes, ending at 15485863.

| Part | Primary upper bound (approximate) |
| --- | --- |
| Distribution range | 0.98800906798623694 |
| Finite quadratic continuation | 0.00088249744327022 |
| Infinite tail | 0.00078417890953012 |
| Total | 0.98967574433903728 |

The separate audit uses unnormalized core masses (6,13,7,1)/30, initial
covered mass 1/10, and a tail-only count. It sums four core contributions
separately, uses scale 10^50 and count cap 100, and uses the earlier
independent logarithm implementation. Its unnormalized total is below
0.990708169906. The primary verifier folds the core into its count law,
using scale 10^60 and cap 120. Both use directed integer rounding. No
numerical optimization status is used in either proof.

```
python code/audit_core30_no10.py data/audit-core30-no10-envelope.json
python code/verify_core30_continuation.py data/schedule-core30-no10.json data/cert-core30-no10.json
python code/audit_core30_continuation.py data/schedule-core30-no10.json data/audit-core30-no10-continuation.json
python -m unittest discover -s code -p test_core30_comparison.py
```

This schedule requires the absence of modulus 10 and the new core lemma.
It is not an unconditional cutoff-6 certificate. Hashes and exact totals
are in `data/core30-overlap-summary.json`.

## Numerical investigation and remaining work

These observations are not impossibility proofs:

* Deleting the modulus-7 plane and trying all 24 orders of 2,3,5,7 did
  not improve the tested schedule. In the usual order, deleting 7 alone
  reproduces the old bound after shifting the hinge by one.
* Keeping the actual 30-point core and removing stage-5 covered points
  improves a difficult tested case from about 1.107 to 1.096 with
  distribution cutoff 100000. Allowing delta_5=0 gives about 1.094.
* At cutoff 3000000, reoptimizing the difficult case with delta_5=1/5
  gives 1.082355496564312. This is an upper-estimate candidate for that
  case, not a lower bound on all schedules or initial measures.
* A floating LP searched arbitrary probability measures on that core,
  with exhaustive cylinder-envelope separation. For one difficult case
  and the fixed continuation it reached about 1.094121 at cutoff
  100000. It does not address all measures and all continuations.

Upper bound 6 remains unproved. Remaining systems must contain both 7
and 10. Their joint geometry and later-stage overlaps remain to be
controlled. The old all-parameter barrier cannot be applied to these
new comparison methods.
