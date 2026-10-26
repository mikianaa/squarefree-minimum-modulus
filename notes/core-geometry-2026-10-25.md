# Core geometry and the seven-prime search

Subsequent result: the seven-prime case is now excluded by the audited plane-load
argument in [the seven-prime obstruction note](seven-prime-obstruction-2026-10-25.md).
The unsuccessful relaxations below are retained as research history.

The global Selfridge problem remains unresolved. These computations refine
necessary conditions at 510510; they do not exclude that period or construct
a covering system avoiding modulus 2.

## A checked finite catalogue

Normalize the prime and twice-prime classes as in the six-prime proof. For the
210 core, first place the three odd pair classes 15, 21, 35 and their doubles.
The core geometry script enumerates their residual sets S in the two parity
layers. It fixes labels 0 and 1 and introduces unused labels in first-use order.

For each S, record its two sizes and, for every d in
{3,5,7,15,21,35,105}, the largest intersection of S with a class modulo d in
each parity layer. Denote these intersection maxima by C_d(e).

The final core classes 105 and 210 remove at most one point in each parity,
plus one further point in the parity chosen for 210. For each of those two
choices, use the optimistic counts

    r_e = max(0, |S_e| - 1 - 1[e = parity of 210])

and capacities min(r_e,C_d(e)). These are necessary-capacity profiles, not
claims that a residual set with those statistics exists.

The normalization is sound for COVER FEASIBILITY. If T is the actual residual,
then t_e=|T_e| >= r_e, and its actual class capacity c_d(e) satisfies

    c_d(e)/t_e <= min(t_e,C_d(e))/t_e
                <= min(r_e,C_d(e))/r_e.

The zero-count case imposes no constraint. A profile with larger capacity/count
ratios in every coordinate is easier to cover. Hence a dominating catalogue
suffices for testing feasibility. This argument must not be used to claim
absolute uncovered-count bounds or impose absolute prefix bounds on the
dominating profiles.

| Core pattern | Canonical pair residual sets | Without even-line quotient | Distinct optimistic profiles | Dominating catalogue |
| --- | ---: | ---: | ---: | ---: |
| 000 | 61604 | 481937 | 2298 | 60 |
| 001 | 45415 | 377945 | 2734 | 66 |
| 010 | 44836 | 354214 | 3255 | 115 |
| 011 | 39174 | 327848 | 3332 | 126 |

The complete profile SETS agree in the two enumerations, not only their counts.
Every profile is dominated by a retained one, checked with exact integer cross
products. The odd-line quotient is common to these two enumerations; its
separate unquotiented audit is in the existing three-axis tests.

Reproduce the catalogue and its audit with:

    python code/audit_core_geometry.py

Results are in `data/core-geometry-profile-audit.json`. The geometric profile
and triple-removal inequalities also have separate ordinary-residue tests in
`code/test_core_geometry.py`.

## Weighted total capacities still leave all anchor patterns open

Adjoin the coordinates 11,13,17. For each core profile and each high-prime parity
pattern, count the total possible contribution of every remaining modulus.
Odd moduli contribute in both parity layers; their even doubles choose one.
For nonnegative weights w0,w1, a necessary condition is

    w0 * demand0 + w1 * demand1
      <= weighted odd capacity + sum_i max(w0*c_i0, w1*c_i1).

All relevant weights can be tested exactly: normalize w0+w1=1 and check the
endpoints and the breakpoints where w0*c_i0=w1*c_i1. Between breakpoints the
capacity deficit is linear. This tests all nonnegative TWO-LAYER weights;
it does not test weights that vary within a parity layer.

Of the 2936 combinations of retained profile and high parity pattern, 298
fail this necessary inequality. Every one of the 32 complete anchor patterns
still has surviving profiles. For core pattern 011, all 126 profiles survive
all eight high-prime patterns. Surviving profiles need not be realizable, and
passing the inequality is not a feasible placement certificate.

Reproduce with:

    python code/block_profile_capacity.py

Detailed tests and their integer weights are saved in
`data/block-profile-capacity-510510.json`.

## Constructive checks of the uncovered-count threshold

`extend_residual_witness.py` extends an explicit partial family by a new prime;
`improve_cover_witness.py` then changes the non-anchor residue classes by local
search. Each output is verified independently on its full ordinary-integer
period. These routines supply upper witnesses only, never nonexistence proofs.

Starting from the optimal 2310 family at pattern 0111, the extension to pattern
01111 modulo 30030 first left 824 residues uncovered. Allowing changes in the
smaller classes improved this to 799 and then 792. The known lower bound for
this pattern remains 282. In particular, the experiment has not determined
whether the sufficient threshold 786 can be proved: 792 is above it.

The 62-class witness is saved in `data/min-uncovered-30030-witness.json`.
An additional direct check over all 30030 integers confirms the count, and a
regression test preserves it. All seven new geometry/witness tests pass.

The next obstruction must retain more placement information, or establish a
substantially stronger residual lower bound. The finite catalogues above do
not by themselves provide an argument for arbitrarily many prime factors.
