# Expanded catalogue and simultaneous core projections

The general squarefree minimum-modulus conjecture remains unresolved. This
work concerns one fixed 22-point core, higher primes 11,13,17,19,23, and
pattern 11111. Neither the new finite-catalogue obstruction nor the search
below excludes arbitrary nine-prime cores or arbitrarily many primes.

## Revalidated previous computation

The preceding expanded-pool coefficient process has terminated. Its former
process handle is missing, and its final output records FEASIBLE, not
OPTIMAL. The restricted-master margin is -0.026959970316238136; repricing
every supplied state gives -0.022494400220835348, leaving numerical duality
gap 0.004465570095402788. It ran 247.766 seconds including preprocessing.
The positive finite-pool closure flag is **false**. No inference of
infeasibility or whole-method optimality follows from this run.

The previous research turn made progress: it enlarged the explicit catalogue
and saved a new coefficient suggestion. The intervening billing lookup did
not change the mathematical state. This continuation uses the completed
output rather than restarting a process on the basis of a timeout.

## Exact pricing gives a larger finite-catalogue obstruction

Even a nonoptimal coefficient suggestion can be rounded to nonnegative
integer weights and independently evaluated. For the saved suggestion,
complete local pricing gives finite-catalogue excess

    4431904473 / 1000000 = 4431.904473 > 0

on the 42,532 actual group entries. The independent audit checks 9,453,375
plane tuples, 32,537 rational exceptional planes, and 182 integer DPs.
Every local minimum is exact and has an attaining state. The unrestricted
clique-affine excess remains -81644513580/1000000, so the result only
excludes the specified finite group catalogue.

Certificate: `data/nine-prime-clique-expanded-priced.json`.
Audit: `data/audit-nine-prime-clique-expanded-priced.json`.

Two subsequent searches retain the same weights and local proof:

| Actual entries | Newly found distinct profiles in this search | Exact catalogue excess |
| --- | --- | --- |
| 42,532 | initial input | +4431.904473 |
| 43,378 | 39, with one second per group | +1132.903354 |
| 45,581 | 92, with five seconds per group | -5493.839687 |

Entry increases include sharing profiles among groups with identical core
divisor resources. They are not counts of distinct coverings. Solver
upper bounds are saved as diagnostics only and are not proof inputs.

For the last pool, an independent integer-only failure witness checks
22 valid local states and 31 actual group residue assignments. The local
states' total reward is 143904845471 and the actual assignments' reward
is 149398685158. Their difference is -5493839687. Exact local minima can
only be smaller, and unrestricted group maxima can only be larger. Thus
these particular weights cannot prove exclusion even with exact
unrestricted group maxima. This says nothing about all possible weights.

Witness and audit:
`data/nine-prime-clique-third-expanded-failure-witness.json` and
`data/audit-nine-prime-clique-third-expanded-failure.json`.

On the **same** 45,581-entry pool, normalized gaps are:

| Saved weights | Excess divided by total event mass |
| --- | --- |
| Earlier cached-slope best | -0.013348374299 |
| Closed old-master weights | -0.015497388176 |
| Current expanded-master suggestion | -0.005493839687 |

The latest weights are better among these three, but are still insufficient.
The comparison is `data/nine-prime-clique-third-expanded-weight-comparison.json`.

## A stronger simultaneous-assignment test

The existing marginal relaxations select local states independently at
different core points. A different search now chooses an actual residue
for every core divisor in every one of the 31 groups simultaneously.
For its resulting counts k_h(x), put

    f_i(x) = max(0, q_i-1-[b_i=x mod 2]-k_i(x)),
    F_x = product_i f_i - sum_{|S|>=2} k_S(x) product_{i not in S} f_i.

The target is F_x<=0 for all 22 points. This still omits the high-prime
residue labels, so success would be a relaxation witness, not a covering.
It would nevertheless give a stronger method barrier than a single failed
weight vector: cap the plane counts at D_i to obtain admissible local
states. For every nonnegative threshold reward, the sum of local minima
is at most the rewards of these states, and those rewards are at most the
sum of the actual group maxima. Consequently no nonnegative separable
threshold-reward inequality could exclude this core using only these
local F conditions and independent group maxima.

No such witness was found. A 120-second coordinate search performed 2,057
passes, reaching sum_x max(F_x,0)=792460, with 20 violated points.
An additional pass optimizing each entire group found no improvement.
All 31 group subproblems returned numerical OPTIMAL, but these solver
statuses are not used to assert a rigorous local or global optimum.
Independent ordinary-residue arithmetic verifies the displayed assignment
and its objective. The positive value is an incumbent objective, **not a
lower bound** on other assignments.

The group optimization has a simple exact reduction. Fix every group except
h. For a tail group, F_x is affine and nonincreasing in k_h. For a plane
group, write F_x=A max(D_h-k_h,0)-B with B>=0. If A<=0, its positive part
is zero; if A>0, its positive part is nonincreasing in k_h. Therefore
max(F_x,0) has nonnegative successive decreases. These are precisely the
threshold rewards accepted by the existing actual-group optimizer. The
search uses explicit returned assignments and retains the incumbent when
a proposed assignment does not improve it.

Code: `code/search_joint_core_projections.py` and
`code/audit_joint_core_projections.py`.
Search: `data/nine-prime-joint-projections-search.json`.
Polished assignment and audit:
`data/nine-prime-joint-projections-polished.json` and
`data/audit-nine-prime-joint-projections-polished.json`.

## Continuation and verification

The latest search source is
`data/nine-prime-clique-third-expanded-source.json`; the corresponding seed
is `data/nine-prime-clique-third-expanded-continuation.json`. They retain
30,339 local states, 45,581 actual group entries, and 4,675 cached group
rows. Continue using the 16,240 clique-repaired scenarios and the original
baseline plus three added reward directions. All old actual columns keep
their indices; scenario cache indices have been shifted by the enlarged
actual prefixes. The search source preserves the original baseline curves.

`data/audit-nine-prime-clique-third-expanded-continuation.json` records an
independent recomputation of every enlarged-pool maximum, unchanged local
proof payloads, preserved prefixes, and remapped scenario rows. It explicitly
depends on the already completed parent exact-minimum audit; that expensive
audit is not rerun. SHA-256 hashes bind the parent certificate, parent audit,
and continued certificate. Reproduce the transformation check with
`tmp/audit_expanded_support_continuation.py` after running the parent audit.

Three new tests compare incremental F updates with direct arithmetic,
including capped planes, check monotone whole-group costs, verify a
feasible toy projection example, and reject corrupted audit data. The three
existing fixed-weight/retarget tests also pass again. See
`data/tests-joint-core-projections.log` and
`data/tests-support-failure-continuation.log`. The earlier 184-test full
suite was not rerun for these isolated additions.

All processes launched in this continuation terminated normally. The next
coefficient search can use the larger seed and saved active states from
`data/nine-prime-clique-expanded-cut-probe.json`. The simultaneous search
can instead be strengthened by coupled group moves or a complete nonlinear
model; its present failure is not evidence that a feasible projection
assignment is impossible. Neither route resolves the missing argument for
arbitrarily many prime factors.
