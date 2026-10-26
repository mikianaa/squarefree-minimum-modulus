# Seven-prime obstruction for the squarefree minimum-2 problem

A later, separate quantitative proof leaves at least 510 uncovered residues
modulo 510510; see `quantitative-core-extension-2026-10-25.md`. The original
plane-load proof and its audits below remain valid.

## Result and scope

The integer certificates exclude a covering using distinct divisors >=3 of
510510 = 2*3*5*7*11*13*17. With the BBMST theorem that a distinct squarefree
covering has an even modulus, and prime compression, this says that any
counterexample to Selfridge's minimum-2 question requires **at least eight
distinct prime factors**. The global problem is not settled by this finite
obstruction. No novelty claim or proof-assistant verification is made.

Verification status: producer certificates complete; both independent arithmetic
and ordinary-residue completeness audits passed. The arithmetic audit checked
2936 base combinations and 116 additional cases. The final audit
reports are `data/audit-plane-capacity-510510.json` and
`data/audit-core-completeness-510510.json`.

This argument does not depend on the CP-SAT optimum u(2310)=121, solver
infeasibility reports, or any timed-out search.

## Normalization and core catalogue

Augment a hypothetical covering by adding every missing divisor >=3 of the
period; this preserves distinctness and coverage. Independently permute the
CRT labels so the class modulo each odd prime p is label 0 and its double 2p
is label 1 in a selected parity b_p. If the 2p class originally used the same
label as p, it was redundant and can first be moved. Flip parity to set b_3=0.

There are 32 parity patterns for the six odd primes. In layer e, the remaining
labels in coordinate p have size D_p(e)=p-1-1[b_p=e].

The 210 core consists of coordinates 2,3,5,7. Enumerate the three odd pair
classes 15,21,35 and their doubles. Permutations fixing labels 0,1 reduce their
possible residual sets S to a finite catalogue. The preceding full-even-line
audit checks that this reduction preserves the entire capacity-profile set.

For d in {3,5,7,15,21,35,105}, let C_d(e) be the maximum intersection of S_e
with one class modulo d. The classes 105 and 210 can remove at most one point
in each parity plus one further point in the parity selected for 210. This
gives two optimistic size/capacity profiles per S, with sizes r_e and caps
min(r_e,C_d(e)). If the actual residual has size t_e>=r_e, its normalized class
capacities are bounded above by those of the optimistic profile:

    actual_capacity/t_e <= min(t_e,C_d(e))/t_e
                           <= min(r_e,C_d(e))/r_e.

A zero-sized layer imposes no constraint. Hence profiles can be discarded when
another has at least as large a capacity/size ratio in every coordinate. The
four dominating catalogues have 60,66,115,126 profiles. This reduction is for
cover feasibility; it is not an absolute residual-count bound.

## Plane loads and the elementary hinge inequality

Fix a core profile, a high-prime pattern for 11,13,17, and the parities chosen
for the even plane classes 2dq, where d is a nontrivial odd divisor of 105 and
q is a high prime. For each q there are seven such choices, thus 128 allocations.
The three high coordinates give 128^3 allocation triples.

The odd classes dq are granted their maximum core capacities independently in
both parity layers. Their even doubles contribute only in their chosen layer.
In a fixed layer, each class puts one indivisible weight into one of the D_q(e)
coordinate labels. Granting independent placements in the two layers only
makes covering easier.

For nonnegative integer weights w_i packed into D bins with loads x_j, and any
nonnegative integer t,

    sum_j (t-x_j)_+ >= max(0, D*t - sum_i min(t,w_i)).               (1)

Indeed, sum_j min(t,x_j) <= sum_i min(t,w_i), as well as <=D*t.

A single integer vector v realizes the right side of (1) simultaneously for
all integer t. Construct it by keeping the largest weights that exceed the
remaining average, then distribute the remaining sum as evenly as possible
over the remaining bins. If there are fewer weights than bins, append zeros.
Its hinge sums satisfy

    sum_j (t-v_j)_+ = max(0, D*t - sum_i min(t,w_i)).               (2)

The identity follows by considering t below or above the final balancing
level; it is also checked independently in the audit via second differences
of the right side, without constructing v.

After the three sets of plane classes, the capacity deficit at a high cell is
at least (r_e-x-y-z)_+. Apply (1)-(2) successively in each coordinate. The total
layer deficit is therefore at least

    F_e = sum_{x in v_11,e; y in v_13,e; z in v_17,e} (r_e-x-y-z)_+.

The producer computes these values for all allocation triples with integer
histogram tensors. The independent audit instead builds the hinge functions
directly, obtains their integer distributions by second differences, and uses
polynomial convolution. All arrays use exact integer arithmetic, with values
far below the int64 limit.

## Allocating the remaining lines and points

Remaining classes contain two or three high primes. If their odd core part is
d, grant capacity c_d(e) times the product of high-coordinate sizes not fixed
by the class; use c_1(e)=r_e. Each odd class gets both layer capacities, and its
even double selects one layer.

First compare min(F_0+F_1) with the sum of the two odd capacities plus the larger
even capacity for every shape. This excludes 2920 of 2936 profile/pattern
combinations.

For the remaining 16, retain the full componentwise minimal frontier of
(F_0,F_1). Allocate the 32 even line/point shapes exactly by a two-layer integer
knapsack. This excludes 14 further combinations. The two survivors are:

| Core pattern | High pattern | Catalogue index |
| --- | --- | ---: |
| 010 | 111 | 1 |
| 011 | 111 | 2 |

The arrays are necessary lower frontiers; passing the test is not a covering.

## Explicit triple classes close the penultimate case

For each of the two surviving core cases, discard every pair residual S whose
optimistic profiles are dominated by already-excluded catalogue entries. This
leaves 53 pair residuals for 010 and 45 for 011.

Now explicitly place 105 and 210. A redundant class can be moved to gain
coverage, so it suffices to remove a nonempty class modulo 105 and then one
remaining point modulo 210. This enumerates 28940 and 19819 distinct residual
sets, respectively. An independent enumeration uses ordinary residues modulo
105 and 210 rather than CRT-coordinate groups.

After removing profiles already dominated by excluded ones, the 010 case has
three dominating actual profiles. All three fail the plane bound. The 011
case has fourteen actual profiles needing further treatment, all dominated
by one marginal-capacity profile. The distinction lost in that domination is
exactly which residue modulo 3 attains the capacity in each parity layer.

## Shared modulo-3 choices close the last case

There are 1108 critical actual residual sets in the 011 case, grouped into
fourteen enhanced profiles. In parity 0, every remaining point has label 2
modulo 3. In parity 1, let u and v be the counts with labels 1 and 2.

An odd plane class 3q has one common residue modulo 3. Its capacity pair is
therefore either (0,u) or (r_0,v). Label 0 misses the entire residual and is
dominated by either of these choices. The earlier estimate allowed
the larger capacity from each layer simultaneously as a relaxation; that
relaxation is now tightened to the two valid alternatives.

For q=11,13,17 there are eight independent choices of these pairs. Keep all
other class capacities at their optimistic values and repeat the elementary
plane bound. Every one of the fourteen enhanced profiles fails in all eight
choices. The smallest capacity deficit is 312. Thus every normalized anchor
pattern at 510510 is excluded.

## From one primorial to arbitrary supports of seven primes

BBMST ensures that a hypothetical distinct squarefree covering has 2 among its
prime factors. Replace any larger odd prime p by a smaller absent odd prime q:
restrict its CRT coordinate to q labels and relabel them. Classes fixing an
excluded label disappear; classes that remain have p replaced by q in their
moduli. Coverage of the restricted product is preserved. Since q was absent,
the remaining moduli are still distinct and squarefree. Since q>=3, no modulus
2 is introduced. Repeating compresses any support of at most seven primes into
{2,3,5,7,11,13,17}. Add unused coordinates if needed. The finite obstruction
therefore excludes all such supports.

BBMST reference: [The Erdős–Selfridge problem with square-free moduli](https://arxiv.org/abs/1901.11465).
The prime-compression argument is also the principle behind the Simpson–Zeilberger
replacement theorem discussed in the six-prime note.

## Reproduction and audit

Run from the repository root with NumPy available. The obstruction and audits
use no SAT, CP-SAT or linear-programming solver.

To recheck the saved certificates in one command:

    python code/verify_seven_prime_obstruction.py

To regenerate the producer certificates first:

    python code/audit_core_geometry.py
    python code/plane_capacity_bound.py
    python code/plane_capacity_bound.py --refine-survivors
    python code/refine_core_triples.py --pattern 010 --output data/core-triple-refinement-010.json
    python code/refine_core_triples.py --pattern 011 --output data/core-triple-refinement-011.json
    python code/finish_coupled_core3.py

The saved refined-frontier and explicit-triple certificates are then checked
by the following independent audits:

    python code/audit_plane_capacity.py
    python code/audit_core_completeness.py

`refine_core_triples.py` and `finish_coupled_core3.py` regenerate the two final
refinement stages. `test_plane_capacity.py` checks the load inequality against
all small bin packings and compares tensor entries with direct cell sums.

The earlier orbit linear programs, integer orbit searches, local projection
searches, and one-coordinate CP-SAT probes are exploratory only and are not
inputs to this obstruction.
