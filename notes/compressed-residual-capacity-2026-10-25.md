# Sharp residual capacity bounds by coordinate compression

Subsequent core-intersection refinements give stronger quantitative bounds;
see `quantitative-core-extension-2026-10-25.md`. The results below are retained
as the original compression calculation and its verified inputs.

The squarefree minimum-2 conjecture remains unresolved. The seven-prime
obstruction remains the strongest restriction on prime support in this
repository. The bounds here strengthen quantitative residual estimates:

* Integer-only propagation gives u(2310) >= 106 and u(30030) >= 411.
* Including the previously obtained CP-SAT optima gives u(30030) >= 488.
  This second statement retains its dependence on solver optimality.
* An actual family with all 254 distinct divisors >=3 of 9699690 leaves
  79979 residues uncovered. It is a partial covering, not a counterexample.

Here u(N) is the minimum number of uncovered residues among all families using
one class for each divisor >=3 of squarefree N. Adding missing classes only
decreases the residual, so the lower bounds also hold for incomplete families.

## Sharp bound using only cardinality

Let X be a product of finite sets of sizes p_1,...,p_k, and let R be any subset
of cardinality r. For nonempty S contained in {1,...,k}, write

    c_S(R) = max_a |{x in R : x_i = a_i for every i in S}|.

Give x the weight 2^z(x), where z(x) is the number of its zero coordinates.
Let W(r) be the sum of the r largest weights in X. Then

    sum_{nonempty S} c_S(R) <= W(r) - r.                 (1)

This is sharp over arbitrary sets of size r. It does not say that an extremal
set is realizable as the residual of a distinct congruence family.

To prove (1), compress every coordinate line to an initial interval. A cylinder
that does not fix that coordinate keeps its intersection size. For a cylinder
that fixes it, the new zero slice contains at least as many points as any old
slice: every formerly contributing line remains nonempty. Thus no c_S decreases.
Compression preserves downset properties already established in other
coordinates, since nested line fibers have nested initial intervals.

After all coordinates are compressed, each maximum occurs at the all-zero
fixed labels. Including S empty, the sum of intersections is exactly
sum_{x in R} 2^z(x), which is at most W(r). To attain it, take the first r
points ordered by decreasing number of zero coordinates and then by increasing
coordinate lexicographic order. This is a downset and has the r largest weights.

The number of points having z zero coordinates is the coefficient of t^z in
the polynomial product_i ((p_i-1)+t). Thus W(r) uses only k+1 integer buckets.

## Extension through one new prime

Let N be squarefree, q a new prime, and R the residual of the q-free classes
in a family modulo Nq. The q class removes one of the q slices. Each remaining
class qd, d>1 dividing N, covers at most c_d(R) points in one other slice.
Consequently the total remaining uncovered count is at least

    max(0, q*r - W(r)).                                (2)

With only a lower bound r>=L, minimize q*r-W(r) over that interval before
clipping at zero. Its slopes are q minus the decreasing point weights, so it
is convex. A minimizing size is max(L, number of points of weight >q).

Starting from the independently checked u(210)=22, (2) gives 90 at 2310
and then 242 at 30030. Starting instead from the solver-certified 121 at
2310 gives 397 at 30030. The first lower counts that make (2) positive are:

| Base period | New prime | Required residual lower count |
| --- | --- | --- |
| 30030 | 17 | 511 |
| 510510 | 19 | 6305 |
| 9699690 | 23 | 79978 |

## Retaining both parity layers

Normalize each odd prime class to label 0 and its double to label 1 in parity
b_p, as in the existing finite obstructions. In parity e, the residual lies
in a box with side lengths D_p(e)=p-1-[b_p=e]. The two alphabets can be ordered
with common labels first and the one exclusive label last. Each coordinate
has a common first label, since p>=3.

Compress inside each layer. Every odd cylinder capacity is bounded by the sum
of its two zero-slice counts; an even cylinder chooses one layer. Let A assign
each nonempty odd-coordinate subset to the parity of its even double. For a
point x in layer e, with zero-coordinate set Z(x), define

    w_e,A(x) = sum_{nonempty S subset Z(x)} (1 + [A(S)=e]).

Let T_e,A(r) be the sum of the r largest such weights in that layer's box.
For fixed sizes r_0,r_1, the exact maximum of the sum of odd and even cylinder
capacities over arbitrary residual subsets is

    max_A (T_0,A(r_0) + T_1,A(r_1)).

The upper bound follows by assigning each even cylinder to a maximizing layer
and compressing. For sharpness, for each fixed A choose the largest-weight
points in each box, breaking ties to make downsets. Their zero slices realize
the odd-cylinder maxima simultaneously, and the selected even assignments
give the stated weight sum. Maximizing over A proves equality.

For extension through q with its double in parity b, a valid pair of lower
bounds, for each A, is

    ell_e = max(0, (q-1-[b=e])*r_e - T_e,A(r_e)).         (3)

For a fixed A, the expression before clipping is convex in r_e and equals
zero at r_e=0. Its positive part is therefore nondecreasing. This justifies
substituting componentwise lower pairs for the unknown actual sizes. Enumerate
A and retain the Pareto-minimal output pairs. This is quantitative dominance,
unlike the normalized cover-feasibility dominance of the seven-prime proof.

The starting pairs are the independently audited 210 certificate. At 2310,
the total lower bounds by anchor pattern are:

| Pattern | 0000 | 0001 | 0010 | 0011 | 0100 | 0101 | 0110 | 0111 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Integer lower | 174 | 158 | 141 | 127 | 121 | 112 | 109 | 106 |

Propagating the full frontiers through 13 gives global lower 411. Intersecting
the 2310 frontiers with the previously solver-certified total bounds raises
this to 488. An arithmetic audit cannot remove dependence on those solver inputs.

Exploratory changes of the initial core to (3,5,11) or (3,7,11), followed by
the missing prime, give frontiers dominated by the (3,5,7) core in all eight
2310 patterns. Their intersection leaves the current frontiers unchanged.
The reversed 13-then-11 extension from 210 is weaker in its global total
bound. These exploratory outputs are not inputs to the audited certificates.

## A precise limitation of cardinality-only iteration

The independently checked 9699690 family leaves 79979 points, with parity
counts (30484,49495). Thus the exact u(9699690) is at most 79979. Since the
propagated bound (2) is nondecreasing in its input lower bound, even knowing
that exact minimum cannot make this particular rule certify more than 13
uncovered points after adjoining 23. Applying the same rule through 29 then
gives zero. This is a limitation of sequential use of (2), starting at this
period, without a later independent lower bound. It does not limit methods
that retain parity, geometry, or compatibility between the congruence classes.

## Reproduction and checks

* `code/extremal_slice_capacity.py`: exact size-only formulas and examples.
  Examples using witness upper counts are explicitly conditional, not lower bounds.
* `code/compressed_parity_capacity.py --frontiers`: NumPy integer implementation
  of (3); outputs `data/compressed-parity-frontiers.json`.
* `code/audit_compressed_parity.py`: separate standard-library implementation
  using ordinary residues and integer weight histograms; all 40 frontiers
  passed in 101.094 seconds (`data/audit-compressed-parity.json`).
* `code/audit_extension_witness.py`: verifies all 254 congruences over the
  full period and records the witness hash and cardinality-method ceiling.
* `code/test_extremal_slice_capacity.py` and
  `code/test_compressed_parity_capacity.py`: exhaustive small-set sharpness,
  extension, and attaining-construction checks.
* `code/test_extension_search.py`: checks reconstruction of actual congruences
  from the search coordinates against ordinary full-period counting.

The extremal formulas have mathematical proofs above and finite regression
checks, not proof-assistant verification. All 29 repository regression tests
passed in 55.610 seconds. No novelty claim is made.
