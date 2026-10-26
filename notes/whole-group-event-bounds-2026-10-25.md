# Whole-group compatibility in the nine-prime relaxation

The global squarefree minimum-modulus question remains unresolved. The verified
support bound remains nine distinct primes. This continuation finds and proves
five new inequalities for one actual 22-point core; it does not exclude that
core's nine-prime extension yet.

Update: the bound 6 below now has an elementary proof using three main cases,
given in `notes/elementary-event-bound-2026-10-25.md`. The same note explains
how fifteen proved inequalities transfer to 327 group constraints. The SAT
proofs here remain valid historical certificates.

## Why small subsets did not suffice

The previously audited fractional witness in
`data/nine-prime-fixed-core-subset-rank-lp.json` passes exact simultaneous
four-point feasibility: 226,765 subset/group combinations, including 2,310,450
candidate threshold-event tuples whose probabilities could violate a cut.
The check enumerates achievable resource hit vectors by Boolean convolution,
then takes their downward closure. A separate tuple-set implementation checks
proposed infeasibilities, and small examples are compared with direct residue
assignment enumeration. No four-point violation was found.

The same witness passes 557,428 weighted rank inequalities using all 73 coprime
nonnegative weight pairs in the rectangle 0<=w0<=15, 0<=w1<=7. This is a finite
check, not a claim that every possible nonlinear rank constraint was tested.

## New whole-group inequality

For a group of higher-prime factors, let z(x) be the number of its core
resources hitting core point x. Each resource chooses one residue class.
For nonnegative integer event weights w(x,t), every consistent choice obeys

    sum_(x,t) w(x,t) 1[z(x)>=t] <= B,

where B is the largest possible score over all shared residue choices.
Consequently every mixture of consistent choices obeys the same inequality.
The bound also applies if some plane hit counts have been clipped downward.

A separating LP selects the event weights; a CP-SAT model maximizes their
score over actual shared resource choices. Its reported bounds alone are
exploratory. The first five bounds were simplified to small integer weights
and then proved using a different SAT encoding and checked proof traces.

For the plane group corresponding to prime 11, one particularly small example
uses the following events on the actual core residues:

| Core residue x | Threshold t | Weight |
| --- | --- | --- |
| 8 | 9 | 2 |
| 8 | 10 | 1 |
| 19 | 4 | 1 |
| 67 | 4 | 1 |
| 68 | 10 | 3 |
| 74 | 9 | 1 |
| 89 | 8 | 1 |
| 97 | 4 | 1 |
| 104 | 9 | 1 |
| 139 | 4 | 1 |
| 164 | 4 | 1 |
| 209 | 4 | 1 |

Their weighted score is at most **6** for any choice of the 14 core residues
(one class for each divisor d>=3 of 210). The previous fractional witness has
score **6.733028127**, so it cannot be a mixture of consistent group assignments.
This incompatibility was not detected by the preceding subset-event checks.

## Proof checks and scope

For each proposed inequality, the independent encoding asks whether its score
can exceed B. It uses exactly one residue choice per resource, threshold flags
that imply enough actual hits, and fresh equivalent copies of those flags for
integer weights. Restricting to residues hitting at least one core point is
valid: replacing an empty choice by a nonempty one can only increase these
nonnegative, monotone event scores.

Glucose produces an UNSAT proof. Every added clause is checked for reverse
unit propagation using a separate Minisat instance before being accepted.
This is a checked propositional proof, not a proof-assistant formalization of
the mathematical encoding. The encoding is also tested against complete
resource assignment enumeration on small weighted examples, including multiple
thresholds at the same point.

| Higher-group axis | Bound | RUP additions checked |
| --- | --- | --- |
| 8 | 19 | 61,842 |
| 10 | 31 | 345,402 |
| 17 | 61 | 105,038 |
| 0 | 6 | 16,583 |
| 16 | 52 | 73,529 |

Total: **602,394** checked additions. CNF hashes, proof hashes and compressed
proof paths are stored with the cuts. The new point-identity auditor rebuilds
the encoding, checks both hashes, and replays these proofs; it does not rely
only on a stored `verified` flag. The group-proof branch requires PySAT, while
the old arithmetic-only branches remain standard-library checks.

## Result after adding the bounds

All five inequalities were added to the local model, including their weighted
threshold contributions to pricing, rational primal repair and integer dual
certificates. A new rational negative-cost witness survives:

    -9963858788173 / 1000000000.

Its independent audit checks 486 resource distributions, 682 point links,
8,639 concentration bounds, 38,180 rank bounds, all 86,461 boundary pair cuts,
all 9,772,351 boundary triple capacity cuts, eight stored subset cuts and the
five newly proved weighted event bounds. It has 260 local support states.
This is still not an actual covering system, nor an optimal relaxation value.

The separator finds another candidate inequality in each of the same five
groups on this new witness. These second-round bounds have **not yet been
simplified, proof-checked or added**. Their data is retained as a concrete next
step. A systematic whole-group separation loop is more promising than repeating
unchanged small-subset searches, but no success of that loop is assumed.

## Files and validation

* `code/joint_subset_events.py`: exact small-subset compatibility check.
* `code/group_event_coherence.py`: whole-group event separation.
* `code/audit_group_event_cut.py`: SAT encoding and checked proof production.
* `code/multitail_identity.py`: weighted event inequalities in the main model.
* `code/audit_multitail_identity.py`: integer witness/dual checks and proof replay.
* `data/nine-prime-core-joint-four-point-violations.json`
* `data/nine-prime-full-rank-violations.json`
* `data/nine-prime-group-coherence-separation.json`
* `data/nine-prime-group-event-certified-cuts.json`
* `data/audit-nine-prime-group-event-cuts.json` and its five `.rup.gz` proofs
* `data/nine-prime-fixed-core-group-events-lp.json`
* `data/audit-nine-prime-fixed-core-group-events-primal.json`
* `data/nine-prime-group-coherence-second-separation.json`

All **76 regression tests pass**. A new edge case found during integration was
fixed: arbitrary-subset separation may return a two-point cut, which must be
routed into the pair-cut machinery rather than the three-or-more-point handler.


Validation also exposed a Windows Glucose proof-logging shutdown failure on
some satisfiable toy formulas. The test checks SAT/UNSAT outcomes with Minisat
and requests Glucose proof logs only for known UNSAT toy instances. This keeps
both sides of the encoding test and avoids that backend failure; the five
research bounds were UNSAT and their proof checks exited normally.


## Further continuation

The second and third rounds are now complete: fifteen event bounds have
proofs, and a new fractional witness has passed a separate audit replaying
all 3,417,086 proof additions. A model retaining explicit whole-group
assignment mixtures has also been implemented and remains inconclusive.
See [the coherent-group continuation](coherent-group-mixtures-2026-10-25.md).
All 82 regression tests pass with a clean exit.
