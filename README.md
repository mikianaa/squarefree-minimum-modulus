# Minimum modulus of distinct squarefree covering systems

How large can the smallest modulus of a finite covering system be when the moduli are
distinct and squarefree? Write M_sf for the supremum. Known before this project:
M_sf >= 2 (Krukenberg 1971: explicit systems with least modulus 2 and lcm 210) and
M_sf <= 118 (Cummings, Filaseta, Trifonov 2025, building on the method of Balister,
Bollobas, Morris, Sahasrabudhe and Tiba, who proved that some modulus must be even).
Dalton and Trifonov (2022, Problem 4) ask whether M_sf = 2, i.e. whether every distinct
squarefree covering system contains the modulus 2. This repository aims to determine M_sf.
Progress and status live in [PROGRESS.md](PROGRESS.md); this file only says what the
project is and how it is laid out.

## What is being done

* **Upper side.** The convex stop-loss comparison developed for the general minimum-modulus
  problem (repository `minimum-modulus`) is specialised to squarefree moduli: coordinates
  Z_p, one modulus per subset of primes, Euler factors 1 + nu/p. An integer-arithmetic
  verifier certifies cutoffs C0 ("all moduli > C0 cannot cover"). Certified on 26 October:
  C0 = 8, hence **M_sf <= 7** (8 is not squarefree). The candidates left are 2, 3, 5, 6, 7.
  Per-prime breakpoint optimization and distribution evaluation through prime 3000000
  give total loss <0.999087538091, checked by two integer-arithmetic implementations.
  This overturns the earlier numerical stop at 10. A separate all-schedule lower
  certificate proves that **7 is optimal for this fixed comparison even with every
  0<=delta<1 allowed**: at cutoff 6, the original loss estimate through prime 10007
  is at least 1.029993468. All 326 action eliminations and 32 remaining initial cases
  have separate integer verification. Different comparison bounds and geometry-aware
  arguments are not excluded. See [the full-range barrier](notes/stoploss-full-range-barrier-2026-10-26.md).
  See [the bound-7 proof and reproduction](notes/upper-bound-seven-2026-10-26.md).
  A subsequent refinement retains the uncovered 30-point core and proves
  that any system with minimum modulus at least 7 must also contain modulus
  10. The no-10 branch has total <0.989675744340 in two integer checks.
  Upper bound 6 remains open; see [the core refinement](notes/core30-overlap-refinement-2026-10-26.md).
  The subsequent [joint 7-and-10 investigation](notes/core210-joint-overlap-2026-10-26.md)
  proves the sharp 104-point residual bound modulo 210, audits all 288846
  equality-case orbits, and gives a scoped fixed-schedule obstruction.
  A later certified refinement incorporates stage-11 geometry and excludes
  one branch even with six stage-11 residues unrestricted: total loss
  <0.999530236097, with independent integer core and infinite-tail checks.
  Further transport and case-split arguments leave only four small
  stage-11 projection conditions for this core, with all residues modulo
  11 unrestricted. Other initial configurations remain, so this does not
  establish upper bound 6.
  The latest adaptive proof removes all four of those conditions for the
  difficult core: 2038 audited leaves cover every later assignment, with
  maximum <0.999996892707. A checked profile transfer covers 7 equality-case
  orbits. Other initial cores remain; the global upper bound is still 7.
* **Lower side.** Exact SAT / CP-SAT searches for distinct squarefree covering systems with
  least modulus >= 3 (moduli dividing a squarefree lcm). An exact integer plane-load
  proof now excludes all four normalized branches at lcm 9699690, leaving at
  least 2 residues uncovered. By BBMST and prime compression, any counterexample
  needs at least nine distinct prime factors. The first remaining primorial is
  223092870, through prime 23. See [the eight-prime proof](notes/eight-prime-obstruction-2026-10-25.md).
* **Determination.** M_sf = 2 if no system avoids the modulus 2 (Dalton-Trifonov Problem 4);
  otherwise an explicit system with least modulus in {3, 5, 6, 7} must be found.

## Layout

| Path | Contents |
| --- | --- |
| `PROGRESS.md` | Current state, timeline, closed routes, open directions. |
| `code/verify_squarefree.py` | The certificate verifier (standard library only, about 20 s; 120 s with `--p-exact 1000000`). Usage: `python code/verify_squarefree.py data/schedule10.json out.json`. |
| `code/audit_verify.py` | Independent re-implementation of the verifier, written from the report alone; reproduces the C0 = 10 total exactly. |
| `code/schedule_opt.py` | Floating-point optimiser of the distortion schedule (input to the verifier, not a proof step). |
| `code/exact_two_stage.py`, `code/stage11.py` | Exact worst case of the stage-5 and stage-7 losses at C0 = 8 over all families, and the stage-11 adversary (toward "<= 7"). |
| `code/cpsat_cover.py` | CP-SAT search for distinct covering systems with moduli dividing N and minimum modulus >= M (squarefree N gives squarefree moduli). |
| `code/verify_seven_prime_obstruction.py` | Re-run both independent audits excluding all supports of at most seven primes; requires NumPy, no optimization solver. |
| `code/verify_eight_prime_obstruction.py` | Re-run all four branch audits excluding every support of at most eight primes; requires NumPy, no optimization solver. |
| `code/verify_six_prime_obstruction.py` | Standard-library-only verifier excluding minimum modulus >= 3 with at most six prime factors. |
| `code/audit_parity_capacity.py` | Independent integer-range and full-grid audit of the six-prime obstruction. |
| `code/three_axis_bound.py`, `code/parity_capacity.py` | Exact base enumeration at 210 and propagation of residual counts in the two parity layers. |
| `code/extremal_slice_capacity.py`, `code/compressed_parity_capacity.py`, `code/audit_compressed_parity.py` | Sharp coordinate-compression bounds and independent integer audit of quantitative parity frontiers; integer-only u(30030)>=411. |
| `code/core_capacity_frontier.py`, `code/audit_core_capacity_frontier.py` | Quantitative core intersections and shared odd-3q labels strengthen the integer-only bound to u(30030)>=501. |
| `code/hybrid_slice_frontier.py`, `code/audit_hybrid_slice_bound.py` | Scalable even-resource allocation and independent exact dual certificates proving u(510510)>=510. |
| `code/conditional_tail_lp.py`, `code/conditional_tail_patterns.py`, `code/audit_conditional_tail.py` | Conditional tail capacities and independent integer audits excluding all 16 extensions through 19 of one specified 210-core; not a global eight-prime obstruction. |
| `code/eight_prime_anchor_branch.py`, `code/audit_eight_prime_anchor_branch.py` | Independent integer audit excluding every core in the eight-prime low-anchor branch 000, with at least 6696 uncovered residues. |
| `code/compressed_tail_branch.py`, `code/audit_compressed_tail_branch.py` | Count-only conditional-tail certificates and independent audit excluding the entire eight-prime branch 001, with at least 305 uncovered points. |
| `code/core_orbit_reduction.py`, `code/portable_tail_certificates.py` | Core orbit/inclusion reduction and portable-dual screening; the final eight-prime proof independently checks all orbit coverage it uses and does not use inclusion reduction. |
| `code/paired_slice_capacity.py`, `code/audit_paired_tail_branch.py` | Pairwise projection capacity bounds and complete independent audits excluding eight-prime branches 010 and 011. |
| `code/audit_tail_primal.py` | Exact rational check of a negative-cost feasible distribution in a weaker aggregate relaxation; not a covering system. |
| `code/multitail_lp.py`, `code/multitail_identity.py`, `code/audit_multitail.py`, `code/audit_multitail_identity.py` | Five-higher-prime conditional-tail investigation and independent rational audits; no nine-prime obstruction yet. |
| `code/rank_capacity.py`, `code/point_pair_capacity.py`, `code/core_subset_thresholds.py`, `code/separate_subset_events.py`, `code/joint_triple_events.py` | Threshold-count, subset-event, and exact joint-triple compatibility checks for the nine-prime relaxations; negative-cost witnesses still survive. |
| `code/group_event_coherence.py`, `code/audit_group_event_cut.py` | Whole-group weighted event bounds, with SAT proof production and separate RUP checks; fifteen bounds verified, global problem still open. |
| `code/audit_event_capacity_bound.py`, `code/weighted_event_capacity.py` | Integer-only capacity certificates and weighted-threshold separation. The sharp event bound 6 also has a three-case elementary proof in `notes/elementary-event-bound-2026-10-25.md`. |
| `code/transfer_group_event_bounds.py` | Transfer the fifteen proved bounds to all groups with identical resources, producing 327 valid inequalities. |
| `code/coherent_group_model.py`, `code/audit_coherent_groups.py` | Stronger explicit per-group assignment mixtures and integer audit; current large-core run is inconclusive. |
| `code/decomposed_coherent_groups.py`, `code/coherent_transport_repair.py` | Decomposed group-mixture optimization with exact one-hit penalties, integer quantile repair, and new assignment pricing; no nine-prime obstruction yet. |
| `code/coherent_point_transport.py`, `code/audit_weighted_coherent_dual.py` | Exact compact point transport and independent weighted-dual verification; certifies exclusion of the explicitly listed finite configuration pool, not the unrestricted nine-prime case. |
| `code/coherent_group_catalogue.py` | Reuses actual assignments across identical resource families and guides group pricing. The prior positive dual becomes negative on the larger shared catalogue; this does not establish a covering. |
| `code/positive_state_dual.py`, `code/audit_positive_state_dual.py` | Stronger individual-state necessity and independent integer slope-envelope audit; improves the shared-catalogue bound but does not yet exclude it. |
| `code/reweight_nonpositive_support.py`, `code/audit_support_reweighting.py` | Optimize combinations of fixed point inequalities; an explicit audited score mixture certifies that reweighting the current coarse inequalities alone cannot exclude the catalogue. |
| `code/support_tail_knapsack.py`, `code/refine_support_slopes.py`, `code/audit_refined_support_slopes.py` | Exact fractional-tail refinement at critical plane tuples and independent rational-slope coverage checks; integer tail counts remain a stronger condition. |
| `code/integer_tail_knapsack.py`, `code/refine_integer_support.py` | Exact integer support pricing for fixed event penalties; independently audited by a forward-credit recurrence in `audit_refined_support_slopes.py`, with an attaining integer state at each point. |
| `code/optimize_group_support_weights.py` | Separate group coefficients at every point, with exact integer support pricing and auditable finite-catalogue certificates. |
| `code/fractional_group_scenarios.py` | Fractional residue/count relaxations of complete residue cases, used as separate coefficient-search constraints; these are not actual covering assignments or proof bounds. |
| `code/case_relaxation_barrier.py`, `code/audit_case_relaxation_barrier.py` | Exact mixture barrier for all group/point rescalings of the specified event curves and case relaxations; independent rational verification also identifies individual threshold deficits that this barrier does not cover. |
| `code/optimize_threshold_support.py` | Changes threshold reward curves, shares valid local states within each parity, reuses checked group-constraint indices, and saves candidates with exact integer local pricing and independently auditable actual-pool bounds. |
| `code/diagnose_support_gap.py` | Locates the discrepancy between a prior local state pool and audited exact local minima for fixed weights; compares point-specific and parity-shared pools. See [the diverse-state continuation](notes/diverse-support-search-2026-10-25.md). |
| `code/cached_support_slopes.py`, `code/benchmark_cached_slopes.py` | Reuse nonnegative slope suggestions, recompute every bound for current weights, and benchmark unchanged exact minima. The optimizer can retain all known states while activating only needed constraints; see [the continuation](notes/cached-support-slopes-2026-10-25.md). |
| `code/probe_curve_basis_barrier.py`, `code/certify_curve_basis_barrier.py`, `code/audit_curve_basis_barrier.py` | Search, rationally reconstruct, and independently audit scoped barriers for nonnegative combinations of supplied reward curves. The verified three-direction barrier admits explicit fractional group states and does not prove whole-method optimality. |
| `code/probe_curve_basis_cuts.py` | Close the supplied finite coefficient problem by repricing all known local and group states and adding only violated inequalities. Separates numerical closure from exact certificate verification and saves active states for continuation. |
| `code/expand_support_catalogue.py`, `code/audit_support_failure.py` | Search new actual group profiles, share them without losing old indices, and independently check a fixed-weight failure witness using only valid local states and explicit residue assignments. |
| `code/search_joint_core_projections.py`, `code/audit_joint_core_projections.py` | Search simultaneous actual core assignments for every group, then audit their local F values. Whole-group moves use exact monotone threshold costs. A feasible result would be a projection relaxation witness, not a covering; heuristic failure proves no obstruction. |
| `code/escape_curve_basis_barrier.py`, `code/price_curve_basis_dual.py` | Find threshold directions escaping a displayed witness, then price numerical candidate weights over the complete local integer domain. Three four-direction candidates remain negative after independent audit. |
| `code/separate_barrier_triples.py`, `code/audit_barrier_triples.py` | Conditional triple-event inequalities with bound two, and independent ordinary-residue DP checks for reported cuts. That screen found no violation; the stronger bound-one clique inequalities below exclude the fractional states. |
| `code/clique_case_constraints.py`, `code/audit_case_cliques.py`, `code/repair_clique_scenarios.py` | Mutually exclusive event-family constraints exclude all 329 displayed fractional barrier states; the independent audit checks 46,630 pairs. Repair the numerical scenario pool while preserving resource means. |
| `code/clique_event_majorant.py`, `code/audit_clique_event_majorant.py` | Rational affine bounds using clique multipliers, independently checked by ordinary residue enumeration; strengthen selected cases and complete group case trees. |
| `code/separate_barrier_cases.py` | Whole-assignment numerical cut search against selected fractional cases. Positive solver bounds require independent certification; the current simpler clique cuts have such certificates. |
| `code/audit_case_pair_check.py` | Independently verifies that specific fractional case states and their averaged witness violate true pair-threshold compatibility; motivates stronger case relaxations without claiming a covering obstruction. |
| `code/pair_event_majorant.py`, `code/audit_pair_event_majorant.py` | Rational affine upper bounds with pair-event multipliers and independent verification; also usable in complete residue case trees. |
| `code/repair_pair_scenarios.py` | Replaces incompatible fractional search states with numerically pair-compatible distributions by reducing count means; not a proof generator. |
| `code/price_support_group_assignments.py`, `code/retarget_support_certificate.py` | Search omitted actual assignments, reevaluate certificates on expanded pools, and improve unrestricted affine event majorants with exact rational checks. |
| `code/branched_event_bound.py`, `code/audit_branched_event_bound.py` | Unrestricted residue-case upper bounds, with independent tree coverage and rational arithmetic checks; improves the global group bounds but does not yet exclude the fixed nine-prime core. |
| `code/parity_cover.py` | Exact CP-SAT search on the residual CRT boxes, with arbitrary label-permutation symmetries removed. |
| `code/core_projection.py`, `code/layered_core_projection.py` | Necessary capacity models retaining the placements of classes in several additional prime coordinates. Feasible outputs are not covering constructions. |
| `code/sat_uncovered_audit.py` | Independent residue-based SAT encoding and optional RUP proof audit for fixed-pattern residual bounds. |
| `code/core_residual_geometry.py`, `code/audit_core_geometry.py` | Enumerate and audit dominating core-intersection capacity profiles; these are necessary conditions, not realizable residual sets. |
| `code/block_profile_capacity.py` | Exact weighted total-capacity tests of those profiles at 510510. |
| `code/extend_residual_witness.py`, `code/improve_cover_witness.py` | Construct and verify partial covering families; their uncovered counts are upper bounds only. |
| `code/search_core_extensions.py`, `code/audit_extension_witness.py` | Actual eight-prime partial-family search and independent full-period audit; current witness leaves 79979 residues. |
| `code/audit_partial_families.py` | Full-period ordinary-integer checks of the saved six-, seven-, and eight-prime partial families. |
| `code/mc_stage.py` | Monte Carlo check of single stage bounds. |
| `data/cert*.json`, `data/schedule*.json` | Certificates (proof objects) and their schedules for C0 = 14, 12, 10; `cert10-full.json` uses exact stages through p = 15,485,863. |
| `data/cft118-*.json`, `data/audit-*.json` | Validation against Cummings-Filaseta-Trifonov's published numbers at C0 = 118, and the independent audit outputs. |
| `notes/` | Dated reports. |
| `meta/` | Meta-analysis log of the agents' progress: `agent-progress.jsonl` (one record per milestone or attempt: bounds, method, candidates, failure cause), `meta_stats.py` (token usage and cost from the local Claude Code and Codex transcripts), rendered `agent-progress.md` and `token-usage.json`. Schema in `meta/README.md`. |
| `literature/SOURCES.md` | References. |

## Provenance

The upper-bound method, code and certificates were produced in October 2026 by Claude (Fable 5.1)
agents working with the repository owner; the verifier was audited by an independent
re-implementation but not by a proof assistant. The mathematical inputs are the
Balister-Bollobas-Morris-Sahasrabudhe-Tiba / Cummings-Filaseta-Trifonov distorted measure,
two comparison lemmas proved in the `minimum-modulus` repository
(`src/conditional-stoploss-improvement.tex`), and two inequalities of Rosser and Schoenfeld (1962).

The solver-free six-prime obstruction and its independent integer audit were developed
by Codex on 25 October 2026. The argument and reproducible certificates are documented
in `notes/six-prime-obstruction-2026-10-25.md`; it has not been checked by a proof assistant.

The seven-prime extension, including an independent plane-arithmetic audit and an
ordinary-residue completeness audit, is documented in
`notes/seven-prime-obstruction-2026-10-25.md`. It is independent of the solver-certified
2310 optimum and has not been checked by a proof assistant.

The complete eight-prime obstruction was developed by Codex on 25 October
2026. Paired projection capacities, ordinary-residue completeness checks,
and independent integer dual audits exclude all four normalized branches.
See `notes/eight-prime-obstruction-2026-10-25.md`; the global conjecture
remains unresolved.
