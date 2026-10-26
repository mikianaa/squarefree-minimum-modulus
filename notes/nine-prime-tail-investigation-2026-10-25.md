# Conditional-tail models through prime 23

## Scope

The eight-prime obstruction remains complete. This investigation concerns
the next period, 223092870, and does not yet exclude it. The global
minimum-2 conjecture remains unresolved.

The conditional-tail inequality extends formally to five higher primes
11,13,17,19,23 over the same 210-core. There are five plane groups and
twenty-six groups involving at least two higher primes. The implementation
now includes all 31 nonempty groups. Depending on anchor parity, its exact
local pricing enumerates about 400000 to 480000 plane-count tuples per
parity. Tail counts still minimize independently over 0,...,16.

`multitail_lp.py` generates aggregate certificates. `multitail_identity.py`
additionally retains each core point and one distribution over ordinary
core projections per resource. Separate Python-integer auditors verify
dual bounds and rational primal witnesses. The latter are distributions in
relaxations, **not covering systems**.

## Verified limitations

The following exact negative-cost feasible distributions have passed the
independent checkers. These are witnesses, not claimed optima; their
numerical magnitudes should not be compared as optimal values.

| Model | Core counts | Higher pattern | Exact cost numerator, denominator 10^9 |
| --- | --- | --- | ---: |
| Paired-capacity aggregate, critical 011 input | (10,13) | 01111 | -5844334159857 |
| Exact-subset aggregate, actual fixed core | (7,15) | 11111 | -8855820644859 |
| Point identities plus exact subset capacities | (7,15) | 11111 | -9789010947435 |
| Exact-subset aggregate plus weighted threshold counts | (7,15) | 11111 | -7042922505060 |
| Point identities plus weighted threshold counts | (7,15) | 11111 | -14773587089740 |
| Point identities, weighted threshold counts, and all pair-event cuts | (7,15) | 11111 | -21568836819427 |

The actual fixed core is the 22-point core in
`conditional-tail-obstruction-2026-10-25.md`. Its exact subset capacity table
was already independently checked over all 4194304 subsets. The paired
input is exactly case 635 of the independently audited 011 refinement.

Each five-prime aggregate audit checks 8639 concentration inequalities.
The point-identity audits also check all 486 projection resources and all
682 point/group linking inequalities. The threshold-enhanced audits check
38180 additional inequalities for the five parity-weight pairs specified
below. Thus further optimization of any one of these same relaxations
cannot yield a positive obstruction for its listed input.

## Integer threshold-count inequalities

Let z(x) be a group's core hit count, clipped to the available alphabet for
a plane group. For positive thresholds t0,t1 let

    m_e = #{x in R_e : z(x)>=t_e}.

Every actual projection family satisfies

    t0*m0 + t1*m1 <= H_j(m0,m1),

where j is 14 for a plane and 16 for a tail group. This follows by summing
the actual hit counts on the superlevel set; clipping only decreases them.
Consequently, for any nonnegative integer weights w0,w1,

    w0*m0+w1*m1 <= max {w0*a+w1*b : t0*a+t1*b <= H_j(a,b)}.

The maximum is over integer pairs within the core parity sizes. This yields
linear constraints on probabilities of threshold events. The initial
implementation uses weights (1,0),(0,1),(1,1),(2,1),(1,2); it does not claim
to include every possible weighted facet.

The first point-identity witness violates 120 of these constraints. Adding
them excludes that witness, but further independently verified witnesses
show that this enhancement alone is insufficient. Exhaustive small-resource
tests check these inequalities against every actual projection assignment.

## Pairwise compatibility of threshold events

For distinct core points x,y and a group with k shapes, let c be the number
of shapes whose moduli divide x-y. Each of those c resources can hit both
points; every other resource can hit at most one. Thus

    z(x)+z(y) <= k+c.

Clipping the plane counts preserves the inequality. If t+u>k+c, the two
events z(x)>=t and z(y)>=u cannot both occur. Hence

    1[z(x)>=t] + 1[z(y)>=u] <= 1.

It suffices to keep the boundary pairs t+u=k+c+1. Larger thresholds yield
weaker inequalities. These cuts retain the identities of two specific
points, information lost by the whole-core threshold-count constraints.

The verified point-identity-plus-threshold witness violates 5095 pair-event
inequalities when all incompatible threshold pairs are counted. The new
separation code retains only the nonredundant boundary pairs. The combined
model adds 3210 boundary pair cuts and 161 rank cuts before producing a new
negative-cost witness. Its independent audit checks **all 86461** boundary
pair inequalities, not merely those added to the optimizer, as well as the
486 resources, 682 links, 8639 concentration bounds, and 38180 rank bounds.
This particular combined relaxation is therefore also insufficient.

## Three-point compatibility: the next refinement

For any specified core point set T, define

    C(T) = sum_d max_a |T intersect (a mod d)|.

Every actual projection family satisfies sum_{x in T} z(x)<=C(T). If
positive thresholds t_x sum to more than C(T), the simultaneous threshold
events are impossible, giving

    sum_{x in T} 1[z(x)>=t_x] <= |T|-1.

Checking all 47740 triple/group combinations finds five violated triple
cuts in the witness that passed every pair cut. Each returned violation is
also directly recomputed with Python integers. The largest has tail group
axis 16 (higher primes 11,13,19), core points **68,97,134**, and thresholds
**12,4,6**. Their total threshold is 22, while C(T)=21; hence the three
event probabilities must sum to at most 2. The witness exceeds 2 by
139297013/1000000000. All three pairwise threshold sums are within their
pair capacities, explaining why pair cuts alone miss this incompatibility.

These triple cuts have now been added to the optimizer. A new rational
negative-cost witness survives every triple capacity inequality; see the
continuation below. They do not constitute a nine-prime obstruction.

## Integer projection search

An independent heuristic search chooses actual integer core projections
for all 486 resources. For fixed plane projections, each tail projection is
optimized directly by maximizing its weighted intersection with the core.
A sixty-second coordinate search found expression **678850**, independently
recomputed from the stored choices. This is an **upper** bound on the minimum
of that expression, not a lower bound or an obstruction. Higher-prime labels
remain unspecified, so it is also not a covering system.

The large difference from the fractional witnesses motivates stronger
compatibility constraints; it does not prove an integrality gap or positivity
of the true integer minimum.

## Reproduction and files

All six rational witnesses were re-audited after the code changes; the
actual-core inputs were checked against the saved 14-class core and its
previously audited exact capacity table. The consolidated report and source
hashes are in `data/nine-prime-relaxation-audits-summary.json`. All **67**
regression tests pass, and every eight-prime certificate hash is unchanged.

Verified rational witnesses and reports:

* `data/nine-prime-critical-paired-lp.json`
* `data/audit-nine-prime-critical-paired-primal.json`
* `data/nine-prime-fixed-core-exact-lp.json`
* `data/audit-nine-prime-fixed-core-exact-primal.json`
* `data/nine-prime-fixed-core-identity-lp.json`
* `data/audit-nine-prime-fixed-core-identity-primal.json`
* `data/nine-prime-fixed-core-rank-lp.json`
* `data/audit-nine-prime-fixed-core-rank-primal.json`
* `data/nine-prime-fixed-core-identity-rank-lp.json`
* `data/audit-nine-prime-fixed-core-identity-rank-primal.json`
* `data/nine-prime-fixed-core-pair-rank-lp.json`
* `data/audit-nine-prime-fixed-core-pair-rank-primal.json`
* `data/nine-prime-core-triple-threshold-violations.json`
* `data/nine-prime-fixed-core-integer-projections.json`
* `data/audit-nine-prime-integer-projections.json`

Use `audit_multitail.py` for aggregate witnesses and
`audit_multitail_identity.py` for point-identity witnesses, with explicit
`--source` and `--output` paths. Both use standard-library integer arithmetic.
Capacity validity remains a separately audited premise.

The point-identity generator saves its local state columns and added cuts.
Passing its output as the next `--source`, with the same cut flags, restores
these objects and the best prior valid bound. No timeout is interpreted as
infeasibility or as evidence that the relaxation can never work.


## Triple and arbitrary-subset continuation

The strengthened point-identity model still has a rational negative-cost
witness, with numerator **-18686968416170** and denominator 10^9. A separate
standard-library audit checks all **9,772,351** boundary triple inequalities,
in addition to the previously checked resource, point-link, concentration,
rank, and pair constraints. This is a rigorous limitation of that relaxation,
not a covering construction.

Exhaustive four-point separation checks 226,765 subset/group combinations
and finds three new violations. All use core points **8,97,109,194** and
thresholds **8,4,4,8**. Their sum is 24, while the exact projection-sum capacity
is 23, so the four event probabilities must sum to at most 3. The largest
violation is 61578118/10^9 (group axis 8, higher primes 11 and 19).

Adding those cuts produces another independently audited negative-cost witness:
**-18549671763562/10^9**, with 274 local support states. The audit checks every
triple capacity inequality and all eight stored subset constraints. It does
not claim to check every arbitrary-subset inequality.

A new CP-SAT separator searches subsets of any size. It minimizes the sum
of event deficits subject to total thresholds exceeding the exact subset
projection capacity. On the final witness it reports INFEASIBLE for every
one of the 31 groups: no further violated cut of this form. These solver
infeasibility reports have not been independently certified. The three
violations it found on the preceding witness exactly match exhaustive
four-point separation; each proposed violation is recomputed in integers.

A stronger, exact simultaneous three-point feasibility check also finds no
violation among all 47,740 group/triple combinations (148 partition profiles).
For three points each modulus partitions them as one common class, three
singletons, or one singleton plus a pair. Enumerating how many resources of
each pair type choose their singleton determines feasibility exactly; the
all-singleton resources fill remaining demands. This check goes beyond the
sum-capacity criterion. Its formula is tested against direct enumeration of
resource assignments on small examples. This is still only a check of one
fractional witness, not a nine-prime exclusion.

New files:

* `data/nine-prime-fixed-core-triple-rank-lp.json`
* `data/audit-nine-prime-fixed-core-triple-rank-primal.json`
* `data/nine-prime-core-four-point-threshold-violations.json`
* `data/nine-prime-core-arbitrary-subset-violations.json`
* `data/nine-prime-fixed-core-subset-rank-lp.json`
* `data/audit-nine-prime-fixed-core-subset-rank-primal.json`
* `data/nine-prime-core-subset-final-separation.json`
* `data/nine-prime-core-joint-triple-violations.json`

The verified global support bound remains nine. Further work must retain
stronger joint projection compatibility, strengthen other inequalities, or
change the conditional-tail bound. Do not rerun the unchanged pair/triple
relaxations expecting a positive certificate.

Validation after this continuation: all **72 regression tests pass**.


## Whole-group continuation

The surviving witness passes exact four-point compatibility and the expanded
weighted rank check. Whole-group event separation nevertheless finds five
new bounds; all five have checked propositional proofs. Adding them leaves
a new independently audited negative-cost witness. See
[the whole-group investigation](whole-group-event-bounds-2026-10-25.md) for
the inequalities, 602,394 checked proof steps, and the next candidate cuts.
The global support bound remains nine, and all 76 regression tests pass.
