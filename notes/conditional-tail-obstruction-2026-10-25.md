# Conditional tail capacities for one eight-prime core

## Scope and result

The global minimum-modulus-2 question remains open in this project. The
previous independently audited obstruction still excludes **all supports of
at most seven primes**, not all eight-prime supports.

This note proves a more limited eight-prime result. Fix the following classes
for the 14 divisors of 210 greater than 2 (each pair is modulus, residue):

```
(3,0), (5,0), (6,4), (7,0), (10,1), (14,1), (15,2),
(21,2), (30,26), (35,3), (42,37), (70,53), (105,83), (210,158).
```

Their residual is

```
8,13,19,59,67,68,74,89,97,103,104,109,134,139,157,164,
173,179,187,194,199,209.
```

It has seven even and fifteen odd points. No extension of this core using
distinct divisors greater than 2 of 9699690 can cover. After normalizing the
classes q and 2q for q=11,13,17,19, all 16 possible anchor parity patterns
are excluded by independently checked integer dual inequalities. The
uniform residual lower bound supplied by these certificates is 130.
This is not asserted optimal. For pattern 1111, retaining individual core
point identities gives the stronger audited lower bound 15257.

This is the core of the existing 79979-uncovered-point partial family. The
result therefore closes that fixed-core construction route through prime
19, but says nothing conclusive about other cores or larger prime supports.

## 1. Augmenting the plane exclusions

Let R denote the specified core residual, with parity layers R_e. Normalize
the classes q_i and 2q_i so the remaining alphabet in the i-th prime
coordinate at a point of parity e has size

    D_i(e) = q_i - 1 - [b_i=e].

This normalization changes no core coordinate. If the 2q_i class was
redundant with q_i, it can be moved to another label, only enlarging the
covered set. Missing classes can also be added. Thus lower bounds after
normalization apply to the original family.

For each i, the other classes involving exactly q_i outside the core are
q_i d for the 14 core divisors d>=3. Let n_i(x) count their core projections
containing x. They delete at most n_i(x) distinct labels in coordinate i.
At each x independently, enlarge those deleted-label sets until they have
size min(D_i(e),n_i(x)). This artificial enlargement can only decrease the
uncovered set. It need not itself be a congruence family.

Write a_i(x)=min(D_i(e),n_i(x)) and f_i(x)=D_i(e)-a_i(x). The residual fibre
after this enlargement is a product with exactly product_i f_i(x) points.

For a subset S of at least two higher coordinates, each remaining modulus
has the form d product_{i in S}q_i with d any divisor of 210, including 1
and 2. If its core projection contains x, its intersection with the enlarged
residual fibre has at most product_{i not in S} f_i(x) points. This upper
bound deliberately ignores whether its fixed higher labels survive.

Let h_S(x) count how many of the 16 core projections in group S contain x.
The union bound now gives a lower bound for the original uncovered count:

    F = sum_{x in R} [ product_i f_i(x)
                      - sum_{|S|>=2} h_S(x) product_{i not in S} f_i(x) ].

This is the key improvement: tail capacities are evaluated **after** the
plane deletions. A positive lower bound for F excludes a cover.

Classes whose core projections miss R can first be moved to projections
that meet R; this only adds coverage on the relevant core fibres. Hence
the projection models may restrict to nonempty core intersections.

## 2. Exact concentration inequalities

For the 14 core shapes define

    H_14(m0,m1) = max_{T subset R, |T intersect R_e|=m_e}
                 sum_{d|210, d>=3} max_a |T intersect (a mod d)|.

All 4194304 subsets of R were enumerated. A second implementation, using
ordinary residue counts and Gray-code updates with Python integers,
independently reproduced the complete 8-by-16 table.

Adding shapes 1 and 2 gives the exact identity

    H_16(m0,m1) = H_14(m0,m1) + m0 + m1 + max(m0,m1).

For any choice of one projection per shape, with point hit counts z(x),
and nonnegative thresholds t0,t1,

    sum_x (z(x)-t_{parity(x)})_+
        <= K_j(t0,t1)
         = max_{m0,m1} [H_j(m0,m1)-t0*m0-t1*m1].

To prove this, choose T to be the points with z(x)>t_{parity(x)}, then bound
the sum of projection intersections on T by H_j. The same inequalities
hold for the truncated plane counts a_i, since truncation only decreases
each summand.

Use j=14 for the four plane groups and j=16 for the eleven tail groups.
These inequalities retain restrictions on concentration which ordinary
maximum-slice capacities alone lose.

## 3. Local states and exact dual proof

A local state in parity e consists of four counts
0<=a_i<=min(14,D_i(e)) and eleven counts 0<=h_S<=16. Its cost is the
bracketed expression defining F. The aggregate relaxation allows a
distribution over these states in each parity layer, with masses 7 and 15,
and imposes all mixed-threshold concentration inequalities above.

For each group and threshold pair choose a nonpositive integer multiplier
y, and a positive scale L. Define

    p_e = min_local_states [ L*cost(state)
                            - sum y*(group_count-t_e)_+ ].

Then every actual extension satisfies

    L*F >= 7*p_0 + 15*p_1 + sum y*K.

The proof is direct: sum the pointwise inequalities, then multiply each
upper concentration bound by its nonpositive multiplier. No solver
optimality assertion is needed. Round the rational bound upward because
the actual uncovered count is an integer.

Although a full local state has fifteen coordinates, pricing is small.
For each of the four-plane count tuples, the eleven tail counts minimize
independently: each appears only linearly in the cost and in its own
threshold penalties. Enumerate all 17 choices for each such count. The
independent audit uses Python integers throughout, including these minima.

The stronger point-identity model additionally links each point's hit-count
marginals to one projection choice per resource. Its nonpositive link
duals are priced over all ordinary-residue projections. This gives
15256349133/1000000, hence 15257, for pattern 1111. The all-pattern result
uses the simpler aggregate model and does not need those links.

## 4. Audited pattern bounds

Bits are ordered by primes 11,13,17,19. Most optimizations stopped at the
first positive certificate; these numbers should not be compared as optima.

| Pattern | Residual lower bound | Pattern | Residual lower bound |
| --- | ---: | --- | ---: |
| 0000 | 644 | 1000 | 1564 |
| 0001 | 1981 | 1001 | 1808 |
| 0010 | 669 | 1010 | 130 |
| 0011 | 3246 | 1011 | 705 |
| 0100 | 543 | 1100 | 2852 |
| 0101 | 487 | 1101 | 952 |
| 0110 | 2954 | 1110 | 615 |
| 0111 | 729 | 1111 | 11214 |

The independent all-pattern audit passed in about 19 seconds, including
the complete subset enumeration. All 40 regression tests passed. Tests
include actual small grids with colliding plane labels, direct checks of
the augmentation/conditional-capacity inequality, and exhaustive local
state comparisons against the separable pricing calculation.

## 5. Files and reproducibility

* `code/core_hit_concentration.py`: NumPy subset-capacity producer.
* `code/audit_fibre_concentration.py`: independent subset enumeration and
  audit of the earlier, insufficient, fixed-tail bound 104242 versus 132663.
* `code/conditional_tail_lp.py`: conditional-tail LP and exact dual producer;
  OR-Tools is used to find multipliers, not to validate the bound.
* `code/audit_conditional_tail.py`: independent standard-library integer
  audit of one conditional-tail certificate.
* `code/conditional_tail_patterns.py`: all-pattern producer and auditor.
* `data/core-hit-concentration-7-15.json`: exact subset-capacity table.
* `data/conditional-tail-eight-lp.json` and
  `data/audit-conditional-tail-eight.json`: point-identity certificate 15257.
* `data/aggregate-conditional-tail-eight-lp.json`: aggregate 1111 certificate.
* `data/conditional-tail-core-pattern-certificates.json` and
  `data/audit-conditional-tail-core-patterns.json`: all 16 certificates and
  their independent audit.

Recheck the all-pattern theorem without NumPy or an optimization solver:

```
python code/conditional_tail_patterns.py --audit
```

Recheck the stronger single-pattern bound:

```
python code/audit_conditional_tail.py
```

The proof is not proof-assistant verified, and no novelty claim is made.

## 6. Remaining work

Other core residuals need their own valid concentration bounds. The method
is general, but this exact H table belongs only to the specified R. A
relaxation using only the anchor-box compression bound and the per-parity
maximum slice capacities is weaker and has not given a positive bound in
the preliminary test, even for this core. Thus replacing the exact subset
geometry by the older capacity profiles is not yet justified as an
effective way to close all remaining cases.

An eight-prime obstruction would still leave arbitrary larger prime
supports unresolved. The global theorem requires an argument controlling
those supports, or a counterexample.
