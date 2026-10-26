# Paired projection capacities exclude the entire 010 branch

Subsequent update: branch 011 has also been independently excluded. The
complete result is in `eight-prime-obstruction-2026-10-25.md`. Statements
below about 011 remaining open describe the intermediate stage of this note.

## Result and scope

For N=9699690=2*3*5*7*11*13*17*19, every normalized family of distinct
divisors greater than 2 with (b3,b5,b7)=010 leaves at least **15** residues
uncovered. This includes every 210-core and all sixteen higher-prime anchor
patterns. An independent arithmetic audit checks the complete enumeration,
capacity inequalities, and selected local dual inequalities. No solver
optimality assertion is used. This is not a proof-assistant verification or
a claim of novelty.

Together with the earlier 000 and 001 certificates, only 011 remains at this
period. Thus any hypothetical cover by divisors of N omitting 2 must contain
6, 10, and 14, with the class of modulus 6 in the opposite parity from both
10 and 14. Missing classes could otherwise be added in an excluded parity
pattern. For arbitrary eight-prime supports, this statement applies after
prime compression, without asserting that the original moduli literally
include 6, 10, and 14. Larger prime supports and the global conjecture remain
open.

## Why the earlier marginal relaxation was insufficient

The complete joint-slice catalogue for 010 contains 4340 histogram profiles,
7945 final-count capacity tables, and 491 pointwise dominating tables at
fixed counts. An ordinary-residue audit checks 351244 pair residuals and
independently reproduces these tables. The pair residual is taken before
placing the classes modulo 105 and 210.

One input with final counts (7,16) and higher pattern 1110 has an exactly
feasible rational aggregate-relaxation distribution with cost

    -6521997797 / 1000000000.

The independent checker verifies its 51 positive-mass states, both parity
masses, and all 4079 concentration constraints. This rules out obtaining a
positive lower bound from that particular relaxation and capacity table,
regardless of further optimization. It is not a congruence covering system.
The capacity table also occurs as catalogue case 584, so merely separating
the count envelopes into their individual profiles does not remove this
obstruction.

## Coupling two projection resources

Let R be an actual pair residual and R_e its parity layers. For two core
moduli d and d', choose their residue classes A and B. If T has m_e points
in each parity layer, then

    |T intersect A| + |T intersect B|
      <= sum_e [ min(m_e, |R_e intersect (A union B)|)
                 + min(m_e, |R_e intersect A intersect B|) ].

This is exact for fixed A and B: take the largest weights of 1_A+1_B in each
layer. Maximizing over both residue choices gives a table P_{d,d'}(m0,m1).
For every perfect matching M of the fourteen core moduli d>=3 dividing 210,

    H_14(m0,m1) <= sum_{(d,d') in M} P_{d,d'}(m0,m1).

The generator minimizes this sum by dynamic programming, recording a
matching for each table entry. The auditor needs only the recorded matching:
it independently enumerates ordinary residue choices, sorts the point
weights, and checks that every modulus occurs once. Matching optimality is
not a proof requirement. Taking the minimum with the earlier joint-slice
and compression upper bound remains valid.

If the pair residual has counts (r0,r1), the final two core classes delete
at most one point in each parity through 105, and one additional point in
one parity through 210. Every final residual therefore contains a subset
with one of the two count pairs

    (max(0,r0-2), max(0,r1-1)),
    (max(0,r0-1), max(0,r1-2)).

The pair-residual capacity table bounds every such subset; exact placement
of these last classes need not be enumerated. The concentration transform
and conditional-tail integer dual are as in
`conditional-tail-obstruction-2026-10-25.md`.

## Local validation and complete branch proof

An initial 23-point final residual was checked by enumerating all 8388608
subsets, with a separate Gray-code audit. Its exact-capacity dual gives
3212 uncovered points for higher pattern 1110. A paired-capacity dual gives
1324 for the same input. Both have independent integer dual audits.

The difficult histogram profile 302 has three pair-residual orbit types.
An ordinary CRT-permutation audit establishes completeness for that profile
and a uniform bound 68 for all higher patterns. The coupled bounds thereby
overcome the negative-cost witness of the weaker marginal relaxation.

For the whole 010 branch, the original joint bounds leave 26 profiles to
refine. These have 82 pair-residual orbit representatives. Every ordinary
residual in those profiles is explicitly checked to be an image of a listed
representative under permutations fixing labels 0 and 1 in coordinates 5
and 7. No inclusion reduction is used. The three new local duals found in
the profile-302 investigation, together with the existing pool, already
supply positive bounds for all these refined inputs.

The complete audit verifies:

* all 351244 ordinary pair residuals and all 4340 histogram profiles;
* all 7945 final-count tables and their 491 recorded pointwise dominators;
* all 82 paired capacity certificates;
* explicit label-permutation coverage of the 314 ordinary residuals in the
  26 refined profiles;
* 10505 selected capacity/pattern inequalities, using 28 distinct local
  dual checks and 1653300 plane-count tuples;
* a minimum positive integer bound of 15.

The audit took about 124 seconds on this machine. The numerical optimizer
is not imported by it. NumPy accelerates integer table operations; final
selected dual numerators are also recomputed with Python integers.

## Reproduction

Using a Python environment with NumPy, run from the repository root:

```
python code/audit_paired_tail_branch.py --pattern 010
```

Primary artifacts:

* `data/paired-tail-refinement-010.json`
* `data/joint-tail-profiles-010.json`
* `data/joint-tail-envelopes-010.json`
* `data/audit-paired-tail-branch-010.json`
* `data/joint-tail-010-primal-probe.json`
* `data/audit-joint-tail-010-primal.json`

The 58-test regression suite passed after the new paired-capacity,
rational-witness, and independent concentration-transform tests were added. The full branch audit is an additional
mathematical certificate check, not a substitute for those tests.
