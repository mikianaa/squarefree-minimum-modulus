# Full group mixtures and continued weighted event refinement

The unrestricted squarefree minimum-modulus conjecture remains unresolved.
The verified global support threshold is still nine primes. This continuation
adds ten proved inequalities and starts a stronger model that retains complete
resource assignments within each higher-prime group.

## Repeated proved inequalities

`code/refine_group_event_bounds.py` repeats separation, integer simplification,
proof production, proof replay and model strengthening. It saves each stage
and can resume completed proof checks after interruption. Only bounds with
successfully replayed UNSAT proofs enter the model.

Starting from the five previously proved event bounds:

| Round | New bounds | New RUP additions checked | Total bounds |
| --- | ---: | ---: | ---: |
| 2 | 5 | 1,463,081 | 10 |
| 3 | 5 | 1,351,611 | 15 |

Together with the first round, this is **3,417,086** checked RUP additions.
The second-round model generated a negative rational witness with numerator
-8608648524844 and denominator 10^9. The third-round model saves numerator
-22688857858436 with denominator 10^9 and 255 point states. These are particular
feasible distributions, not relaxation optima; their objectives should not
be compared as an optimal-bound sequence. The separate final witness audit passed: it replayed all fifteen proofs and
checked every arithmetic constraint. Its report is
`data/audit-nine-prime-group-refinement-round3-primal.json`. The intermediate
second-round witness was not separately re-audited.

Resumption exposed valid pair cuts with thresholds strictly beyond the minimal
boundary. They are now accepted using the defining inequality t+u>capacity,
rather than requiring exact membership in the boundary-only separator list.
A regression test checks that these resume correctly and invalid cuts fail.

## Stronger model: mixtures of actual assignments

For each higher group h, let A_h be the finite set of complete choices of one
core residue per resource. For a in A_h, let z_a(x) count its hits at core point
x. Introduce a probability distribution lambda_h on A_h.

Each core point x also has a probability distribution mu_x over its entire
local count vector v. Impose, for every h,x and threshold t>=1,

    Pr_(v~mu_x)[v_h >= t] <= Pr_(a~lambda_h)[z_a(x) >= t].

These are 10,692 threshold inequalities for the present 22-point core and
31 higher groups. They imply every nonnegative weighted event inequality for
each group, including all fifteen individually proved cuts. They also imply
the previously used concentration bounds and mean resource links.

The objective is the expected conditional-tail expression summed over core
points. Every actual set of resource choices gives feasible point masses and
group masses: take the actual counts, clipped at the available coordinate
sizes for plane groups. Thus a strictly positive lower bound would exclude
this fixed core and parity pattern.

This remains a relaxation. Different groups' mixtures need not share one
joint distribution, and the model still omits consistency of actual labels
in the higher-prime coordinates. A negative feasible mixture, if found, would
be a limitation of this stronger relaxation, not a covering system.

## Computation and present limits

`code/coherent_group_model.py` uses column generation in both directions:
exact point-state pricing, and pricing of consistent group assignments.
All intermediate assignments found by the group solver are available as new
columns. Saved point states and group assignments can be resumed without
discarding previous work.

The previously independently checked integer assignment with expression
**678850** supplies a feasible starting point. Its conversion to the stronger
model is independently tested, including all 10,692 threshold inequalities.
This is an upper bound on the minimum expression, not an obstruction.

Four finite runs remain **inconclusive**. The latest file retains 3,523
point columns and 2,068 group columns. Its restricted master remains near
678850, while its best reported pricing lower bound is -13379460. The latter
uses CP-SAT upper bounds without separate proof replay and is exploratory.
No negative full-group mixture or positive audited certificate has been found.

`code/audit_coherent_groups.py` is ready to check a rational mixture directly
using only Python integers: it recomputes every hit count from explicit residue
assignments, verifies masses and all threshold domination inequalities, and
recomputes the objective. Unlike the event-cut model, a negative mixture would
need no optimization proof. Positive dual bounds still need separately checked
group-pricing proofs.

The fourth run seeds event prices from an earlier dual and stabilizes them.
It completes one pricing iteration in 37 seconds, but the next restricted LP
uses the remainder of its 180-second allowance without a second iteration.
Its best lower bound is unchanged. This is an inconclusive terminal run,
not a process waiting in the background.

The pricing signs and dual calculation are tested on a small one-group model
whose exact optimum is 5: fourteen resources each hit at most one of two
points with total demand 19. The model recovers that value. Stabilized pricing
also recovers 5, and every old penalty value is checked after conversion to
event prices. All **90** regression tests pass with a clean Python exit after
adding arithmetic event proofs and same-resource proof transfer.

## Current files and next work

* `data/nine-prime-group-refinement-summary.json`
* `data/nine-prime-group-refinement-round2-cuts.json`
* `data/nine-prime-group-refinement-round2-model.json`
* `data/nine-prime-group-refinement-round3-cuts.json`
* `data/nine-prime-group-refinement-round3-model.json`
* `data/nine-prime-coherent-groups-model.json`
* `data/nine-prime-coherent-groups-resumed.json`
* `data/nine-prime-coherent-groups-expanded.json`
* `data/nine-prime-coherent-groups-stabilized.json` (latest full-group state)

Resume the latest full-group state rather than rebuilding from scratch. Its
slow dual progress suggests seeding or stabilizing prices using the earlier
verified event-bound duals, or enriching the configuration pool. Any resulting
certificate must still be audited at its actual scope. Neither completing
this one core nor completing the nine-prime case would by itself settle the
unrestricted conjecture.

The latest finite-event model is separately strengthened by transferring each
proved inequality to every group with identical core resources. See
`notes/elementary-event-bound-2026-10-25.md`. Its 327 inequalities are still
weaker than the full coherent-group model described here.

One possible smaller search within the full-group model is to fix two actual
assignments for each group and a rational mixing probability. At each core
point there are at most 32 choices of the five plane counts. Tail high-count
mass can be transported among those plane states using linear variables:
the objective is linear in each tail count and contains no product of tail
counts. After solving, conditional tail marginals can be realized by nested
intervals on a common uniform variable, with at most 27 joint tail states
per plane state. This could avoid the slow large restricted master. It is a
proposed search reduction, not yet implemented or validated.
