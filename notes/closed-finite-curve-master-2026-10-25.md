# Closing the finite four-direction search after clique repair

The global squarefree modulus-2 conjecture remains unresolved. This work
concerns the fixed 22-point core, high-prime pattern 11111 through 23, and
four specified reward directions. A finite coefficient optimum is only a
search result; a proposed covering obstruction still needs complete local
integer pricing and unrestricted group bounds with independent audits.
Even an exclusion of this fixed core would not settle all nine-prime cores
or the case of arbitrarily many prime factors.

## The complete finite problem

The input contains 24,653 valid local states and 56,437 group states:
40,197 actual core assignments and 16,240 clique-repaired fractional cases.
There are 2,709 allowed normalized reward coordinates. For each coordinate,
write its cumulative reward as Q_{h,x,b}, normalized to have total mass one.

The mixture formulation maximizes a common margin epsilon such that the
group expectation minus the local expectation is at least epsilon in every
coordinate. Its coefficient dual chooses nonnegative w_{h,x,b}, with sum w=1,
to minimize

    sum_h max_{g in G_h} sum_{x,b} w_{h,x,b} Q_{h,x,b}(g_x)
      - sum_x min_{k in L_x} sum_{h,b} w_{h,x,b} Q_{h,x,b}(k_h).

This dual has one group upper-bound variable and one local lower-bound
variable for each group/point. It can begin with one state inequality per
group/point. After each solve, every supplied group and local state is
priced; violated inequalities are added, and the LP is solved again. A
numerical optimum for the complete finite pools is reported only when the
final scan finds no violation above the stated tolerance.

`code/probe_curve_basis_cuts.py` implements this procedure. Positive
constraint duals recover the group mixtures, and negative local constraint
duals recover the local mixtures. Tests check this sign convention and
independently audit both a positive rational mixture and a separating
integer-weight certificate on small examples.

## Terminal statuses and numerical tolerances

The original all-state GLOP run stopped FEASIBLE at margin
-0.0419904269421. Switching the mixture model to CLP stopped FEASIBLE at
-0.0316451141781 after 180 solver seconds. Neither run established optimality.

Identical numerical projected columns were removed without rounding.
This eliminated only 13 local and 8 group columns, leaving 24,640 and 56,429
respectively. The benefit was small. This preserves the numerical LP, not
a new rational equivalence theorem; all later arithmetic certificates are
checked against the original integer rewards and actual states.

A first GLOP constraint-generation run stopped FEASIBLE at
-0.0299887839089. A CLP continuation then found a useful coefficient
suggestion, but failed the retained-row check because the default solver
tolerance was looser than the 1e-8 separation test. Its last saved candidate
is explicitly marked NUMERICAL_CONSTRAINT_VIOLATION; it is not an optimal
solution. The program now asks for primal and dual tolerances of 1e-10 and
records any retained-row violation as a terminal numerical failure rather
than mistaking it for missing constraints.

With those tolerances and the saved active rows, CLP closes the finite
problem in eight iterations, taking 80.531 seconds including preparation:

* numerical optimum: -0.02973483488210175;
* independent full-pool repricing within the same numerical model:
  -0.029734834882095662;
* numerical primal/dual gap: 6.09e-15;
* largest local constraint violation: 7.56e-16;
* largest group constraint violation: 1.18e-15;
* active state inequalities: 1,287 local and 1,684 group, plus normalization.

The saved result is `data/nine-prime-clique-all-group-cut-tight-probe.json`.
It has `finite_pool_closed: true` and `verified: false`: numerical closure
does not replace an arithmetic covering certificate.

## Keeping suggestions distinct from proofs

Optimality is not necessary for a useful weight suggestion. Any nonnegative
weights can be rounded to integers, priced over the complete local domain,
and compared with rigorous actual-group upper bounds.
`code/price_curve_basis_dual.py` now accepts saved feasible coefficient
suggestions as well as optimal ones, while requiring explicit coefficient
weights. A feasible mixture status alone supplies no such weights. The
certificate records the probe status and does not assume solver optimality.

The constraint-generation search saves active states and current coefficient
suggestions for continuation. Its progress checkpoint is not evidence that
a process is live; process handles and terminal results remain authoritative.
Warm starts reuse actual count states and group identifiers, recomputing
every coefficient for the current problem rather than reusing an old matrix.

The finite optimum uses only the known local pool. Complete pricing can
still find lower-cost local states, and the finite actual assignment
catalogue is not the unrestricted group domain. Those are the two remaining
checks before even this fixed core could be excluded.

## Complete arithmetic pricing and the next actual assignments

Both suggested weight sets were completely priced and independently
audited. The early feasible suggestion has catalogue excess
-6567499078/1000000, normalized by total event mass to -0.006567402638.
That improves the preceding best on the same old 40,197-entry catalogue.
Its unrestricted clique-affine excess remains -99561349948/1000000.
The audit checks 9,453,375 plane tuples, 9,354 rational exceptional cases
and 181 integer DPs. See `data/audit-nine-prime-clique-all-group-priced.json`.

For the closed finite optimum's weights, complete local pricing still
leaves a positive old-catalogue excess:

    148548042 / 1000000 = 148.548042 > 0.

The independent audit verifies 9,453,375 plane tuples, 16,037 rational
exceptional cases and 167 integer DPs. Thus all local states are covered
by this certificate; the remaining restriction is the named group catalogue.
The unrestricted clique-affine excess is still -94223840099/1000000.
See `data/nine-prime-clique-closed-priced.json` and
`data/audit-nine-prime-clique-closed-priced.json`. This does not exclude
arbitrary group assignments.

The small positive margin was immediately tested against new actual group
assignments. One second of CP-SAT search per group, with a known catalogue
assignment as a hint, found 43 new distinct profiles. Sharing these and
previously unshared profiles across identical core resource families expands
the catalogue from 40,197 to 42,532 entries. The same weights then give

    -11819536974 / 1000000 = -11819.536974.

No solver upper bound is used for this conclusion. A separate integer-only
failure witness supplies one valid F<=0 local state per core point and one
actual residue assignment per group. Their respective scores are
135643669186 and 147463206160. Exact local minima cannot exceed the former;
any unrestricted group upper bounds must admit the latter. Therefore these
particular weights cannot exclude the unrestricted core, regardless of how
tightly their group maxima are computed. This is a fixed-weight failure,
not an all-weight barrier or an actual covering.

The independent check is `code/audit_support_failure.py`, with report
`data/audit-nine-prime-clique-expanded-failure.json`. It checks only 22 local
states and 31 actual assignments and does not rely on either an optimizer
or the expensive exact-minimum audit. Candidate generation and catalogue
sharing are in `code/expand_support_catalogue.py`.

## Current continuation

Use `data/nine-prime-clique-expanded-source.json` as the search source and
`data/nine-prime-clique-expanded-continuation.json` as the seed, retaining
`data/nine-prime-clique-repaired-scenarios.json` and the same three added
basis files. The source preserves the original baseline reward directions;
it is not a stale old proof with a changed pool.

There are 28,430 retained local states and 42,532 actual group entries,
plus 16,240 repaired fractional cases. All old actual columns retain their
indices. The 4,537 cached group row indices are explicitly shifted where
they refer to scenarios, because new actual columns precede the scenarios.
The local proof payload and event weights are unchanged by retargeting.
See `data/audit-nine-prime-clique-expanded-continuation.json`.

On this same expanded catalogue, normalized gaps of the previously saved
weights are -0.012032462597 (old cached-slope best), -0.017057213993 (early
all-group suggestion), and -0.011819536974 (closed-optimum weights). The last
is best among these three, but remains negative. The comparison is
`data/nine-prime-clique-expanded-weight-comparison.json`; scores from
different catalogues must not be compared as if they had the same scope.

All 184 existing and new master regressions pass, including backend
agreement, exact numerical-column deduplication, full-pool closure,
warm-start equivalence, and exact audit of a feasible coefficient suggestion.
Three additional tests cover the independent failure witness, false witness
rejection, preservation of old columns and remapping of scenario indices.
Logs: `data/tests-closed-curve-master-full.log` and
`data/tests-support-failure.log`.
