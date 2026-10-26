# Unrestricted event bounds by residue case splitting

The global conjecture remains open. These bounds apply to every assignment
inside each of the 31 groups for the fixed 22-point core and higher pattern
11111. Unlike a catalogue maximum, their validity is not restricted to
assignments already found by the search. The combined covering obstruction
is still insufficient for this fixed core/pattern.

## Case reduction

For a group, fix the residues of several core resources. Let h_x be the
number of hits already forced at point x. An event with threshold t is
already satisfied if t<=h_x; otherwise its threshold on the remaining
resources becomes t-h_x. Events requiring more hits than there are remaining
resources can be omitted from this case's score.

Bound the remaining score with the exact affine majorants described in
`notes/group-support-weight-search-2026-10-25.md`, and add the score of the
already satisfied events. This gives an unrestricted upper bound for the
particular partial assignment.

Branch on a remaining resource d, using every residue occupied by a point
of the core residual. A residue missing all core points is dominated by any
occupied residue: replacing it can only increase hits, and every event
weight is nonnegative. Thus these cases suffice to bound the unrestricted
maximum, even though unoccupied residues are not explicitly branched on.

Every node keeps its exact local bound and the minimum of that bound with
its inherited ancestor bound. Splitting replaces a node by all occupied
residue cases. Unsplit leaves retain valid bounds. The maximum leaf bound
is consequently valid even if the search stops early or reaches its node
budget. The generator refines the leaf with greatest current bound first.

## Independent check

`code/audit_branched_event_bound.py` reconstructs every shifted event list
from the actual fixed residues. It checks all affine inequalities and
resource capacities with rational arithmetic, every inherited bound, and
every complete child-residue list. It rejects missing cases, repeated nodes,
unreachable nodes and false bounds. No LP objective or CP-SAT upper bound is
trusted. Small tests compare fully split bounds with exhaustive ordinary
residue enumeration, including assignments using unoccupied residues.

## Current results

At maximum depth six, the 31 case trees contain 10,341 nodes and 8,176
leaves. For group 0, its unrestricted score upper bound decreases from
1501501511806 to 1038235694140. Applying all 31 improved bounds to the
current point inequalities changes the unrestricted excess from
-3007630.536058 to -795482.326767. Both are negative, so neither is a
covering obstruction.

Meanwhile a further search found 100 explicit omitted group assignments.
The expanded 40,197-entry catalogue has exact excess -31153.213220 under
these same weights. This is distinct from the unrestricted bound. It also
does not produce a covering or a full coherent-model feasible witness.

The combined certificate and audit are:

* `data/nine-prime-support-third-pool-branched-dual.json`
* `data/audit-nine-prime-support-third-pool-branched-dual.json`

The audit passes, including all local integer support minima and every
residue case. Tree-generation data are in
`data/nine-prime-branched-event-depth6.json`; new actual assignments are in
`data/nine-prime-support-omitted-third.json`. The suite passes 136 tests in
`data/tests-branched-event-bounds.log`. All runs here have finished.

## Next use

The next weight optimization should account for unrestricted group bounds,
not only maximize the margin against the known catalogue. A promising
separation relaxation for a fixed residue case uses fractional choices of
remaining resource residues and count distributions at each point, with
mean-count demand bounded by the fractional resource supply. Its dual is
the affine majorant LP. A maximizing fractional assignment provides a
linear score constraint for the group coefficients; final claims must still
be checked by the independent rational tree and integer point auditors.
Such fractional assignments are optimization constraints, not coverings.

For catalogue-based continuation, the original event curves and expanded
pool are in `data/nine-prime-group-support-third-resume-input.json`. The
latest optimized certificate retains 383 admissible support seed states.
Do not re-scale the already multiplied event weights repeatedly. Deeper
case splits and better weight choices remain possible; the current bound
is not a method-wide optimality result.
