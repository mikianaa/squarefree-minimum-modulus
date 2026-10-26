# A solver-free obstruction with at most six prime factors

## Result and scope

Every finite covering system with distinct squarefree moduli greater than 2 must
have at least seven distinct prime factors in the lcm of its moduli.

The finite statement proved below is that the divisors d >= 3 of 30030 cannot
form a covering system, for any residue choices. BBMST and prime replacement
extend this to every squarefree lcm with at most six prime factors. This does
not settle Selfridge's global minimum-modulus problem. No novelty claim is made
relative to the full literature.

The proof is a finite integer enumeration with explicit counting inequalities,
not a CP-SAT infeasibility assertion. The small base enumeration was checked
both with and without symmetry reduction. A second implementation checks every
subsequent residual-count pair and obtains the same numerical bounds.

## Normalization and parity boxes

Write L = 2*product(P), where P is a set of odd primes. Add unused moduli d|L,
d>=3, if necessary. Normalize the class of each p in P to CRT label 0. A class
of modulus 2p at label 0 is redundant and may be moved to a nonzero label;
coordinate permutations then put its label at 1. Let b_p in {0,1} denote its
parity. Flip parity to put b_3=0.

The points not covered by these anchor classes form two boxes, with dimensions

    D_p(e) = p - 1 - 1[b_p=e],    e in {0,1}.

All remaining uncovered points lie in these boxes. Every odd divisor d>1 of
product(P) fixes the coordinates for primes dividing d. Its class can meet at
most

    C_d(e) = product_{p in P, p not dividing d} D_p(e)

points of box e. If the residual in that box has size r_e, the sharper bound is

    c_d(e,r_e) = min(r_e, C_d(e)).

These bounds apply to arbitrary positions of the residual points.

## Base: exact three-coordinate enumeration at 210

After anchoring p and 2p for p=3,5,7, the moduli 15,21,35 give three axis lines
with shared labels across the parity layers. Moduli 30,42,70 give one additional
line each in a chosen layer. Moduli 105 and 210 add at most three points total.

The three shared lines have only 100 configurations up to permutations fixing
the anchor labels. Enumerate those, then all placements of the even lines. The
following bounds result; see `three_axis_bound.py` and the accompanying audit
note for the precise enumeration. An enumeration of all 2304 configurations
without quotienting label permutations gives the same table.

| b_3,b_5,b_7 | Residual after anchors | Maximum covered by pair moduli and doubles | Minimum final residual |
| --- | ---: | ---: | ---: |
| 000 | 63 | 30 | 30 |
| 001 | 58 | 30 | 25 |
| 010 | 56 | 30 | 23 |
| 011 | 54 | 29 | 22 |

## Propagation through a new prime q

Consider a normalized family at Lq, and let (r_0,r_1) be the residual counts of
its q-free part at L. The class q clears one complete q-slice. The class 2q
clears at most r_{b_q} more points in parity b_q. For each odd divisor d>1 of
product(P), the modulus qd clears at most c_d(e,r_e) in each parity. The modulus
2qd can clear at most c_d(e,r_e) in one chosen parity e.

Allocate these even moduli to parities, writing A_e for the odd divisors whose
2qd class is assigned to e. The final uncovered count in parity e is at least

    G_e(r_e) = max(0,
        (q-1-1[b_q=e])*r_e
        - sum_d min(r_e,C_d(e))
        - sum_{d in A_e} min(r_e,C_d(e))).

The bound generously ignores overlaps, q-slice conflicts, and consistency of
the locations in the two parities. Therefore it cannot exclude a real family
by being too restrictive.

For L=210, q=11, there are only seven odd divisors d>1 and 128 allocations of
their doubles. For each anchor pattern, minimize G_0(r_0)+G_1(r_1) over all
integer counts satisfying the base total bound and the anchor-box size limits.
This yields the U values in the table below.

The primary verifier performs the minimization by piecewise-linear convexity:
every pair with total at least u dominates some (a,u-a), and each G_e is
minimized over all counts above that coordinate lower bound. The independent
auditor enumerates the bounded integer ranges directly. Both obtain identical
U values.

## Exclusion at q=13

For L=2310, a cover at Lq would require

    (q-1)*(r_0+r_1) <= r_{b_q}
      + sum_d (c_d(0,r_0)+c_d(1,r_1)+max(c_d(0,r_0),c_d(1,r_1))).

The maximum term permits every even modulus to choose whichever layer gives
the larger capacity. The table gives the smallest deficit (left minus right)
for q=13, for each b_13. The auditor checks every pair of integer counts in the
anchor boxes with r_0+r_1 >= U. All deficits are strictly positive.

| b_3,b_5,b_7,b_11 | U: uncovered modulo 2310 at least | Minimum deficit, b_13=0 | Minimum deficit, b_13=1 |
| --- | ---: | ---: | ---: |
| 0000 | 156 | 498 | 414 |
| 0001 | 138 | 387 | 309 |
| 0010 | 114 | 203 | 143 |
| 0011 | 106 | 170 | 124 |
| 0100 | 103 | 137 | 106 |
| 0101 | 90 | 63 | 45 |
| 0110 | 96 | 107 | 83 |
| 0111 | 90 | 53 | 69 |

The primary verifier only needs pairs with total U. To justify this reduction,
suppose an allocation could cover a larger pair (R_0,R_1). For each positive
coordinate, divide its necessary capacity inequality by R_e. Every ratio
min(R_e,C)/R_e is nonincreasing in R_e. The same allocation would consequently
pass the necessary inequalities at any componentwise smaller pair; a zero
coordinate imposes no requirement. Every larger pair dominates a pair of sum U,
all of which already fail the summed inequality. The full-grid audit avoids
this monotonicity reduction entirely and reaches the same deficits.

This exhausts all sixteen anchor patterns, so no cover exists at 30030.

## From this lcm to arbitrary prime supports

BBMST ensures that 2 divides the lcm. Simpson--Zeilberger's prime-replacement
theorem permits replacing an odd prime by a smaller absent odd prime. These
replacements preserve squarefreeness, distinctness, and the restriction d>=3.
Any counterexample with at most six prime factors would therefore give a cover
whose moduli divide 2*3*5*7*11*13, contradicting the finite obstruction.

References and the threshold condition are recorded in
`minimum-three-audit-2026-10-25.md`.

## Reproduction and limitations

Both commands below use only the Python standard library:

    python code/verify_six_prime_obstruction.py
    python code/audit_parity_capacity.py

The first reproduces `data/six-prime-obstruction.json`. The second independently
reproduces the integer-range minimization and the full-grid deficits in
`data/audit-six-prime-obstruction.json`. It also audits a stronger intermediate
two-coordinate frontier, recorded separately.

The test suite checks the unquotiented base enumeration, an actual 22-point
residual witness, and the auxiliary solver changes. Its solver-model checks
require OR-Tools; the obstruction verifiers do not. No proof assistant has
checked the argument.

An initial propagation to seven primes does not settle that case: for anchor
pattern 01111 at 30030, the current frontier only gives a total lower bound 74,
and all 62 lower pairs pass the capacity tests for adjoining 17, for either
parity of the modulus 34. Passing this relaxation is not evidence of a cover.
