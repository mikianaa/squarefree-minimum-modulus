# An elementary proof of the sharp event bound 6

The unrestricted minimum-modulus conjecture remains unresolved. This note
replaces one 16,583-step SAT/RUP proof in the fixed-core investigation by a
short counting argument with three main cases. It also identifies a reusable
weighted-capacity principle and the scope of transferring existing bounds.

## Counting principle

Put

    D = {3,5,6,7,10,14,15,21,30,35,42,70,105,210}.

Choose one residue class modulo each d in D, and let z(x) count how many of
these fourteen classes contain x. For any finitely supported weights w(x)>=0,

    sum_x w(x) z(x) <= C(w)
        := sum_(d in D) max_(a mod d) sum_(x == a mod d) w(x).

This follows by summing the contribution of each chosen residue class.
Write C(T) for weights 1 on a finite set T and 0 elsewhere. In particular,
thresholds t(x) cannot all hold if sum_x w(x)t(x)>C(w). For k such thresholds,
their event indicators then sum to at most k-1. This reasoning applies to
arbitrary moduli and nonnegative weights; only the numerical examples below
use D and the particular core points.

## Statement

Let [P] be 1 when P holds and 0 otherwise. Define

    A = [z(8)>=9],   A' = [z(8)>=10],   B = [z(68)>=10],
    C = [z(74)>=9],  V  = [z(89)>=8],   E = [z(104)>=9],
    L = sum_(x in {19,67,97,139,164,209}) [z(x)>=4].

Then, for every choice of the fourteen residue classes,

    2A + A' + 3B + C + V + E + L <= 6.                 (1)

The notation C for an indicator is distinguished from the capacity C(T) by
its argument. This is a statement about one finite resource family, not a
covering-system obstruction by itself.

## Capacity calculations

The following small tables are obtained directly from the counting formula.

| T | C(T) |
|---|---:|
| {8,74}, {8,89}, {8,104} | 16, 15, 16 respectively |
| {68,74}, {68,89}, {68,104} | 16, 17, 16 respectively |
| {8,68,x}, x=19,67,97,139,164,209 | 21,20,20,20,22,21 respectively |

For three disjoint pairs partitioning the six points counted by L:

| Pair P | C({8} union P) | C({68} union P) |
|---|---:|---:|
| {19,97} | 16 | 17 |
| {67,209} | 15 | 15 |
| {139,164} | 17 | 17 |

The remaining unweighted calculations are:

| T | C(T) |
|---|---:|
| {74,19,67,97,164} | 24 |
| {104,19,67,97,164} | 24 |
| {74,89,67}, {74,89,97} | 18,17 respectively |
| {74,104,19}, {74,104,67} | 21,21 respectively |
| {89,104,67}, {89,104,97} | 17,18 respectively |
| {74,89,104} | 23 |

Finally, take weights (2,2,1,1,2) on (67,89,97,139,164).
The fourteen per-modulus maxima, in the displayed order of D, are

    4,5,4,2,3,2,4,2,3,2,2,2,2,2,

whose sum is 39. Consequently

    2z(67)+2z(89)+z(97)+z(139)+2z(164) <= 39.          (2)

If V=1 and all four other points in (2) contribute to L, the left side is
at least 40. This is the one weighted calculation needed in the proof.

## Three cases

**Case 1: A=1.** The first capacity table gives C=V=E=0: for example,
z(8)>=9 and z(89)>=8 would require 17 hits on a set of capacity 15.

If B=1 as well, every x counted by L has
z(8)+z(68)+z(x)>=23 if its event occurs, whereas C({8,68,x})<=22.
Thus L=0, and the score in (1) is 5+A'<=6.

If B=0 and A'=0, each of {19,97} and {67,209} contributes at most one
event to L, since 9+4+4=17 exceeds its capacity with point 8.
The third pair contributes at most two. Hence L<=4 and the score is <=6.
If A'=1, the same argument also excludes two events in {139,164},
because 10+4+4=18>17. Now L<=3, again giving a score <=6.

**Case 2: A=0 and B=1.** Again C=V=E=0 by the first table, and A'=0.
Each of the three disjoint pairs contributes at most one event to L,
because 10+4+4=18 exceeds every capacity in the last column of the
pair table. Thus L<=3 and the score is <=6.

**Case 3: A=B=0.** Also A'=0. Put k=C+V+E. We have k<=2, because
all three events would require 9+8+9=26 hits on a set of capacity 23.

If k=0, use L<=6. If k=1, use L<=5: for C or E, all four low events
at {19,67,97,164} would require 9+16=25>24; for V, use (2).

If k=2, at least two of the low events are absent. When C=V=1,
the events at 67 and 97 are impossible; when C=E=1, those at 19 and 67
are impossible; when V=E=1, those at 67 and 97 are impossible.
Each assertion follows from the corresponding three-point capacity above.
Thus L<=4, and in every subcase k+L<=6. This proves (1).

## Sharpness and independent checks

Equality is attained by choosing, in the displayed order of D, residues

    2,4,1,5,7,13,14,5,7,24,25,19,89,89.

Direct substitution gives V=1 and five of the six low events, while
A=A'=B=C=E=0. Hence the score is exactly 6.

The human proof uses 27 unweighted capacities and (2). Independently,
`code/audit_event_capacity_bound.py` recomputes every capacity with integer
arithmetic and checks all 3,072 event-threshold profiles. All 2,079 profiles
of score greater than 6 violate one of these 28 inequalities. The audit
also checks the attaining residue assignment. This exhaustive check is
additional validation; the three cases above are a standalone proof.

Files:

* `data/small-event-bound-human-certificate.json`
* `data/audit-small-event-bound-human.json`
* `data/small-event-bound-capacity-certificate.json` (a separate 17-inequality
  certificate, verified over the same 3,072 profiles)
* `data/audit-small-event-bound-capacity.json`

## What this changes for the larger problem

The bound applies to every single higher-prime group using these same fourteen
core resources, not only the group where it was discovered. Decreasing hit
counts by clipping also preserves the bound. Similarly, inequalities for the
sixteen-resource family including moduli 1 and 2 transfer between all groups
using that family. They must not be transferred between the two families.

Applying this observation to the fifteen previously proved event bounds gives
327 distinct group inequalities (3 times 5 plane groups plus 12 times 26
remaining groups). Forty of them exclude the previously audited round-three
fractional witness. This is progress on the fixed-core relaxation; it does not
yet provide an exclusion for this core or for arbitrary nine-prime supports.

After adding all 327 constraints, a new feasible fractional witness has
251 support states and exact objective

    -12450242120427 / 1000000000 < 0.

The independent audit checks all 486 resource distributions, 682 point links,
8,639 concentration bounds, 38,180 rank bounds, 86,461 pair bounds,
9,772,351 triple bounds, eight stored subset bounds, and all 327 event bounds.
It replays all 3,417,086 proof additions from the fifteen distinct SAT proofs;
identical proof inputs are reused within the audit, avoiding 327 replays.
The cache key includes the core points, actual resource moduli, events, bound,
and proof identity. Tests reject changed thresholds and transfers between
the fourteen-resource and sixteen-resource families. All 90 tests pass.

Files are `data/nine-prime-transferred-group-model.json` and
`data/audit-nine-prime-transferred-group-primal.json`. This witness proves
that these particular relaxation constraints still do not force a positive
objective. It is neither a covering system nor evidence that the unrestricted
conjecture is false. The reported value is not asserted to be the LP optimum.

A separate search allowing point-capacity weights 1,2,3 on all 31 groups found
no additional violated threshold incompatibility for that old witness. All
31 CP-SAT searches returned INFEASIBLE, but those negative search outcomes
have not been independently proof-checked and are not used as a theorem.

The unrestricted conjecture still requires control of arbitrary core geometry
and arbitrarily many prime factors. The counting principle is general; the
specific capacities and case split proved here are not uniform in those data.
