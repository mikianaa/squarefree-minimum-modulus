# Progress

Status of the determination of M_sf, the largest possible minimum modulus of a finite
covering system with distinct squarefree moduli. Published state of the art: 2 <= M_sf <= 118.
Last updated 26 October 2026.

## Current state

Latest upper result: **M_sf <= 7**, certified by two integer-arithmetic
implementations with total loss <0.999087538091. The same convex stop-loss
comparison succeeds after per-prime breakpoint optimization and distribution
evaluation through prime 3000000. Thus the earlier numerical stop at 10 was
not a method barrier. See `notes/upper-bound-seven-2026-10-26.md`.

Latest method limitation: **7 is optimal for the current catalogue comparison
even with all 0<=delta_p<1 allowed**. For cutoff 6, the original loss-estimate
sum through prime 10007 is at least 1.029993468 for every real schedule. A
separate integer audit gives at least 1.03219703. The proof verifies 326
action eliminations and all 32 remaining initial cases, then bounds every
remaining continuous control. Different comparisons and geometry-aware
arguments remain outside its scope. See `notes/stoploss-full-range-barrier-2026-10-26.md`.

New geometry refinement: **a system with minimum modulus at least 7 must
contain modulus 10**, hence both 7 and 10. Retaining the actual uncovered
30-point core yields a new comparison law with masses (6,13,7,1)/27 on
{1,2,4,8}. Two integer continuation implementations give normalized total
<0.989675744340 for the no-10 branch. The core inequality has an elementary
capacity proof and a complete separate enumeration of all three patterns.
The general upper bound remains 7; a difficult tested case containing 10
still gives a numerical estimate about 1.082. This is no impossibility
claim for the improved method. See `notes/core30-overlap-refinement-2026-10-26.md`.

Joint 7-and-10 investigation: retaining the full core modulo 210 gives a
sharp **104-point residual bound**. All **288846 equality-case orbits**
(35490420 assignments after base normalization) have an independent
completeness audit. Joint measure/tail searches improve a difficult tested
case to about **1.00294**, still only a heuristic, incompletely priced
score. A separate rational dual proves a value >=**1.001192376** for every
initial measure on that core with one specified stop-loss schedule; this
does not exclude other schedules. Upper bound 6 remains unresolved.
See `notes/core210-joint-overlap-2026-10-26.md`.

Later-stage refinement: the previously difficult fixed core now has a
**certified branch exclusion**, after retaining stage-11 geometry and
optimizing both the residual measure and later parameters. The fully fixed
core has loss <0.970803496584. A broader certificate fixes only ten of its
stage-11 residues and allows arbitrary choices or omissions for moduli
330,385,462,770,1155,2310; its total is **<0.999530236097**. The core comparison
has an independent exhaustive integer audit of 810000 catalogues, and the
infinite continuation has two integer implementations. The sharp residual
787 after stage 11 is also audited for this one core. Other initial
configurations remain, so **the global upper bound is still 7**, with no
impossibility result for this stronger method. See the later-stage section
of the same note and `data/audit-core2310-six-free-branch.json`.

The branch has since been enlarged: **only four stage-11 base projections
remain fixed**, at moduli 22,33,55,66. Every component modulo 11, all twelve
other stage-11 moduli, and every later prime are unrestricted; omissions
are allowed. A rational transport audit reduces all 70 pairs of 77/110
projections to ten integer core certificates, with total
**<0.998351837531**. This uses a complete 7-plane in the specified core and
measures with density depending only on x modulo 30. See
`data/audit-core2310-density30-transport.json` and the final section of the
joint-overlap note. The separate complete audit of all **70 projection
pairs** also passes, with the sharper maximum **<0.997109604548**; see
`data/audit-core2310-split77-110-branch.json`. These are conditional branch
results, not M_sf<=6.

Latest continuation: **all four remaining stage-11 conditions have now
been removed for the difficult fixed core-210 family**. A complete adaptive
tree covers 1800 basic projection cases with 2038 integer-certified leaves,
including 17 cases split further by the projection of modulus 165. Every
stage-11 residue and omission, and every later prime, is unrestricted.
The maximum rational total is **<0.999996892707**. Every selected core has
an independent integer enumeration audit; the infinite tail uses the two
independent common coefficient audits. See
`data/audit-core2310-complete-projection-tree.json` and the final sections
of `notes/core210-joint-overlap-2026-10-26.md`.
A separately checked residual-profile transfer excludes **7 equality-case
orbits** in the existing catalogue. Other core-210 profiles and larger
residuals remain, so **the global upper bound is still 7**.

On the lower side, **all supports of at most eight primes are excluded**. Every
counterexample omitting modulus 2 needs at least **nine distinct prime
factors** in its lcm. The four branch bounds at 9699690 are 6696, 305, 15,
and 2. The first remaining primorial is **223092870**, through prime 23.
See `notes/eight-prime-obstruction-2026-10-25.md`. The general conjecture
remains unresolved.

| Side | Result | Status | Where |
| --- | --- | --- | --- |
| Expanded catalogue and simultaneous projections | Independent complete local audit excludes the 42,532-entry catalogue with excess +4431.904473 | Two actual-assignment searches add 39 and 92 profiles, growing the shared pool to 45,581 entries; a separate integer witness proves the same weights fail with excess -5493.839687. Retains 30,339 local states and 4,675 cached group rows. Simultaneous actual core-projection search remains infeasible at its incumbent, with no impossibility claim. | `notes/expanded-catalogue-and-joint-projections-2026-10-25.md`, `data/audit-nine-prime-clique-expanded-priced.json`, `data/audit-nine-prime-clique-third-expanded-failure.json` |
| Closed finite master and new actual profiles | The four-direction finite master closes numerically; independent exact pricing excludes the old 40,197-entry catalogue with excess +148.548042 | Searching actual assignments finds 43 new profiles. Shared catalogue grows to 42,532 entries and the same weights give -11819.536974. A separate integer-only witness proves these weights cannot exclude unrestricted group assignments. Retains 28,430 local states and 4,537 group rows; all 184 master regressions plus three new witness tests pass. | `notes/closed-finite-curve-master-2026-10-25.md`, `data/audit-nine-prime-clique-closed-priced.json`, `data/audit-nine-prime-clique-expanded-failure.json` |
| Conditional event cliques | All 329 fractional states in the preceding three-direction barrier violate a verified mutually exclusive event-family inequality | Independent integer audit checks 46,630 incompatible pairs. Rational clique majorants improve all 12 selected pair-aware case bounds and a complete group {11,17} tree by 540.158688. All 16,240 old scenarios repaired numerically; 180 tests pass. The covering problem remains open. | `notes/clique-case-constraints-2026-10-25.md`, `data/audit-nine-prime-three-curve-clique-cuts.json`, `data/audit-nine-prime-clique-axis6-tree.json` |
| Three-direction cone barrier | Exact mixtures prove that every nonnegative combination of the three specified reward directions fails with group bounds admitting the displayed fractional states | Independent audit verifies 2,030 coordinates and 1,044,562 pair inequalities. A fourth direction escapes this witness, but three completely priced candidates still have negative catalogue excess. Latest state pool: 22,783 local states and 4,092 group rows. Triple-event screening finds no violation in the 329 fractional states. This is not whole-method optimality or a covering. | `notes/three-curve-cone-barrier-2026-10-25.md`, `data/audit-nine-prime-three-curve-barrier.json`, `data/audit-nine-prime-four-curve-third-priced.json` |
| Cached slopes and compact search | Better audited fixed-core candidate has catalogue excess -47976.248035 and normalized excess -0.009567; no obstruction | Same-weight repricing reduces two difficult points from 35,477 / 27,597 exceptional planes to one each without changing their exact minima. Retains 17,211 local states and 3,476 group constraints. A newer compact candidate has a smaller raw deficit but worse normalized gap, so the better weights are preserved. All 164 tests pass. | `notes/cached-support-slopes-2026-10-25.md`, `data/nine-prime-compact-master-resume.json` |
| Diverse local states and parity sharing | Audited fixed-core catalogue excess improves to -80607.833091; no obstruction | 13,185 local states and 2,235 cached group constraints retained. Local states now shared within parity. Normalizing by total event weight improves the gap from -0.037614 to -0.010636, still negative. Independent audit checks 9,453,375 plane tuples, 63,589 exceptional cases and 226 integer DPs. All 160 tests pass. | `notes/diverse-support-search-2026-10-25.md`, `data/audit-nine-prime-parity-shared-support.json` |
| Certified pair-aware upper bounds | Rational certificates strictly improve 38 of 40 selected cases; a complete group {11,17} case tree improves its bound by 2786.819499 in normalized units | Independent checks verify every pair inequality and rational bound. Repaired 2,793 of 16,240 fractional scenarios. New search retains 3,381 local states, but its audited catalogue excess -546600.730794 is still negative. All 156 tests pass. | `notes/pair-affine-certificates-2026-10-25.md` |
| Three-curve and pair-compatible continuation | Audited three-curve candidate gives catalogue excess -1109526.907486; no obstruction | A specific averaged fractional barrier witness violates the true pair bound at core points 67 and 164 by exactly 4961018949/50000000000. Independent rational check passes. Pair-aware LP pricing improves 38 of 40 selected cases numerically; rational pair-dual upper bounds are the next step. All 150 tests pass. | `notes/three-curves-and-case-pairs-2026-10-25.md` |
| Free threshold curves | New curves defeat the preceding specific fixed-curve barrier mixture, but still fail against the actual catalogue | Independent exact point audits verify catalogue excess -3783060.183220 for the latest candidate. Reused 41,375 valid nonpositive local states and retained a 2,270-state seed. The numerical master is not closed; no new obstruction or optimality result. | `notes/free-threshold-support-2026-10-25.md`, `data/audit-nine-prime-free-threshold-known-support.json` |
| Fixed-curve case-relaxation barrier | All 682 nonnegative group/point rescalings are rigorously insufficient with the specified old reward curves and 8,176 case relaxations | Exact mixtures of 243 local integer states and 214 fractional case states dominate every curve coordinate. Independent standard-library audit passes. The witness has 4,235 individual threshold deficits, so arbitrary changes to event weights remain open. This is not a covering or a general method limitation. | `notes/case-relaxation-barrier-2026-10-25.md`, `data/audit-nine-prime-case-relaxation-barrier.json` |
| Fractional case scenarios | Optimized group coefficients against 8,176 fractional case states; new unrestricted excess -545194.768584 | Independent audit checks exact integer point minima and all new rational case bounds. Still no fixed-core obstruction. Repricing the same old cases finds omitted fractional states; 16,240 scenarios and 1,104 local support states are saved for continuation. All 138 tests pass. | `notes/fractional-case-scenarios-2026-10-25.md`, `data/audit-nine-prime-scenario-branched-dual.json` |
| Unrestricted group-score bounds | Residue case splitting improves the fixed-core unrestricted excess from -3007630.536058 to -795482.326767 | Independent rational audit checks 31 complete case trees with 10,341 nodes. The bound is still negative, so this does not exclude the fixed core/pattern. All 136 tests pass. | `notes/branched-event-majorants-2026-10-25.md`, `data/audit-nine-prime-support-third-pool-branched-dual.json` |
| Nonpositive-support finite-pool obstruction | The specified 40,097-entry pool is excluded for the fixed core and pattern 11111 | Certified excess 69871.764093. A later search adds 100 assignments; the same weights give -31153.213220 on the new 40,197-entry pool. The earlier scoped certificate remains valid; the expanded pool and unrestricted assignments remain unresolved. | `notes/group-support-weight-search-2026-10-25.md`, `data/audit-nine-prime-group-support-reoptimized.json` |
| Coherent-mixture search | Catalogue sharing and solver hints found 262 further group assignments | Latest checkpoint has 210,833 point states and 1,805 group configurations. Independently audited feasible positive-point-mean objective 656403.066794462. General problem unresolved. | `data/nine-prime-pointwise-catalogue-pricing.json`, `data/audit-nine-prime-pointwise-catalogue-primal.json` |
| Individual-state necessity | Exact integer event-penalty minima over every F_x<=0 local count state are now verified | Independent forward-DP audit checks 29 integer subproblems, 58 exceptional plane tuples, and the full 9,453,375-tuple plane domain. Catalogue excess -251084.298; still no exclusion. | `data/audit-nine-prime-integer-support-pricing.json` |
| Scoped baseline reweighting limitation | Nonnegative combinations of the specified baseline 22 exact integer point inequalities cannot exclude the catalogue | An independently audited mixture of 56 actual assignments dominates those 22 event-score requirements. Separate group coefficients now overcome this limited barrier on the stated finite pool. This was never a general method optimality theorem. | `data/audit-nine-prime-integer-reweighted-support-dual.json` |
| Pointwise necessity | A covering requires F_x<=0 separately at every core point | The new objective disallows cancellation between fibres. The old negative 327-cut witness has thirteen positive point expectations, so it does not settle feasibility of this strengthened relaxation. | `notes/pointwise-cover-necessity-2026-10-25.md` |
| Nine-prime investigation | Fifteen proved bounds transferred to 327 group constraints; the sharp bound 6 now has an elementary three-case proof | A new negative-cost witness passes independent arithmetic audit of all 327 constraints and replay of 3,417,086 proof steps, with identical proofs reused across identical resource families. Stronger full-group model has 3,523 point columns and 2,068 group configurations but remains inconclusive. No nine-prime exclusion or cover. | `notes/elementary-event-bound-2026-10-25.md`, `data/audit-nine-prime-transferred-group-primal.json`, `notes/coherent-group-mixtures-2026-10-25.md` |
| Lower (obstructions) | Any counterexample omitting 2 needs at least nine distinct prime factors in its lcm | All four branches at 9699690 independently excluded. BBMST plus prime compression extends this to every support of at most eight primes. The finite uniform bound is 2; no solver optimality or proof-assistant claim. | `notes/eight-prime-obstruction-2026-10-25.md`, `data/eight-prime-obstruction.json` |
| Aligned construction test (direction 3 of the direction-analysis note) | Pure counting reproves k<=6 obstructed; k=7,8 (the open case) are undecided by counting alone | For the "G-residues zero" slice extension, R_allowed := L minus the size of {x: omega(gcd(x,L))<=ceil(log2 q)-1} is below the certified u(L) lower bound at k=4,5,6 (210,2310,30030), giving a one-line reproof, but exceeds it at k=7 (1074 vs 510) and k=8 (30682 vs 2), so counting cannot settle the open nine-prime case. The naive residue-0 realization misses the whole coprime-to-L class (1,658,880 points at k=8, 54x the allowed residual), showing the restriction is not nearly tight; resolving k=8 needs the static two-family model (direction 2), not this cheap test. No change to primes_needed (still 9) or the first remaining primorial (still 223092870). | `notes/aligned-construction-test-2026-10-26.md`, `code/aligned_construction_test.py`, `data/aligned-construction-test.json` |
| Static two-family model, exact test (direction 2) at k=5,6,7 | INFEASIBLE at k=5 (true optimum) and k=6 (heuristic witness); UNKNOWN at k=7 even after a 16.7x variable-count compression (104,367 vars, ~14.5 min CP-SAT time) -- the scaling wall at k=7 is about constraint structure, not variable count | CP-SAT needing q-1-fold coverage of each witness's real residual by the divisors >1 of L. The full non-aligned model buys nothing over the aligned special case on the k=5 data. A lossless compression (one boolean per occurring (divisor,residue) pair instead of per all residues) cuts k=7 16.7x and k=8 an estimated 25.3x, but does not make k=7 solvable, so k=8 (~1.38M compressed vars) was not attempted. k=6, k=7 use heuristic (non-optimal) witnesses, so their verdicts are calibration data, not new obstructions. | `notes/static-two-family-test-2026-10-26.md`, `code/static_two_family_test.py`, `data/static-two-family-2310-13.json`, `data/static-two-family-30030-17.json`, `data/static-two-family-510510-19-compressed.json` |
| Coherent-residual G-side, exact reduction at k=8->9 (q=23, the correct open case) | The G-side (slice-coverage) of the single most achievable coherent target (m=4199=13*17*19, \|R\|=2310, the smallest omega(m) whose capacity ratio clears 22) is exactly SAT, found via a new CRT reduction to a tiny equivalent problem on Z_2310, with the real solution independently re-verified by direct integer arithmetic on all 2310 real points of Z_9699690. Not every omega=3 target that clears the 22-fold threshold is confirmed SAT: testing all four to completion gives a clean monotonic trend (margin 0.94->SAT 22s, 0.63->SAT 238s, 0.23->UNKNOWN, 0.10->UNKNOWN after ~15min), with margin above the threshold, not just clearing it, driving solve difficulty -- and an inverse relationship with F-friendliness (easier G-side cases have smaller, less F-achievable |R|). F-side (can an actual family leave residual near this target) remains open. Widening the coherent-target search past omega=3 to omega=4 found m=770 (|R|=12597, SAT in 58s, independently re-verified): the best EXISTING general-purpose witness has 0 overlap with multiples of 770 too, giving the best confirmed F-side comparison found tonight for the primary case: **6.35x short**. Two methods that deliberately steered toward a (smaller, omega=3) target did WORSE (53x and 2056x), an informative finding on its own. The same pattern repeats one level up (9->10, q=29): G-side SAT at multiple targets, with the F-side comparison improved from 525x (m=96577, best margin, smallest |R|=2310) down to 27.2x by choosing larger, still-confirmed-SAT coherent targets (m=5005, |R|=44574, SAT in 272s, every result independently re-verified). The practical frontier for this exact method sits between |R|=44574 (SAT, m=5005) and 74290 (UNKNOWN after 311s, m=3003). Either way, F-side achievability remains open and gets harder as primes are added; trying to steer a construction toward any of these targets still does not help on the F-side. Does not change primes_needed (9) or any certified bound; settles one half of one specific route, not the global question. | `notes/coherent-virtual-reduction-2026-10-26.md`, `code/coherent_virtual_reduction.py`, `data/coherent-virtual-reduction-4199-23.json` |
| Spread-vs-coherent capacity tradeoff, quantified at nine primes | A single full-class residual (coherent) clears the 22-fold requirement already at omega(m)=4 (m=390, ratio 22.2), vs the spread greedy witness's 7.88 at comparable \|R\|; confirms numerically why the "second induction on duplicates" is load-bearing, not optional, at this scale. A naive actual construction leaning toward that coherence (safe-residue-restricted greedy) does NOT close the gap: it comes out worse (7.02 at m=390) or barely better (8.15 at m=30) than the plain spread greedy, at 6-11x the residual size. | Closed form cap(R)=Sum_d L/lcm(d,m) for R={m\|x}, validated by brute force. Coherent residuals are NOT directly realizable by a distinct-moduli F except via the already-known duplicated-modulus mechanism (24 Oct slice-extension entry, mod-105 case), itself excluded one level down there; whether an analogous obstruction holds at nine primes is open and not attempted here. An actual CRT-safe-residue construction (bug in the first version caught by its own uncovered->contains-multiples-of-m invariant, then fixed) shows naively leaning toward coherence is not a shortcut, and local search on top of it only improves 5% (still 11.7x the coherent target), reasonably strong evidence the smaller-scale duplicate-modulus trick (2310/105) has no easy nine-prime analogue, though not a proof. Supplies a cheap building block for the existing maximiser search, does not run it. | `notes/coherence-capacity-tradeoff-2026-10-26.md`, `code/coherent_residual_capacity.py`, `code/safe_residue_construction.py`, `data/coherent-residual-capacity-table.json` |
| Direct greedy construction at the first open primorial (direction 4) | New heuristic upper witness u(223092870)<=1213060 (0.544% density) after local-search tuning; cap(R)/\|R\|=7.88 stays at 28% of the 28-fold coverage a tenth-prime extension would need, essentially unmoved by tuning | A one-pass numpy greedy (small divisors first), then a proper Gauss-Seidel coordinate-descent local search (coverage-count array, O(L) reshape per modulus, 5 sweeps to convergence) over all 510 moduli>=3 of the nine-prime primorial 223092870. Tuning shrank the residual 8.5% (1325141->1213060) but barely moved the capacity ratio (7.78->7.88), the strongest evidence yet that the 28% figure is not a crude-method artifact. Self-verified (independent coverage rebuild; the one-pass witness also passed a blind auditor re-derivation). No actual 9-prime cover found or expected; evidence, not a proof. | `notes/nine-prime-greedy-construction-2026-10-26.md`, `code/greedy_cover_nine_primes.py`, `code/local_search_nine_primes.py`, `data/local-search-nine-prime.json` |
| Eight-prime branch | Entire 011 branch excluded at 9699690 | Independent audit checks 317869 ordinary pair residuals, 4323 profiles, 332 paired-capacity certificates, 19422 inequalities, and 1181475 plane tuples. Uniform lower bound 2. | `data/audit-paired-tail-branch-011.json`, `notes/eight-prime-obstruction-2026-10-25.md` |
| Eight-prime branch | Entire 010 branch excluded at 9699690 | Complete joint-profile enumeration and 82 paired-capacity certificates leave at least 15 points uncovered. Independent audit checks 351244 ordinary pair residuals, 4340 profiles, 10505 inequalities, and 1653300 local plane tuples. | `notes/eight-prime-anchor010-2026-10-25.md`, `data/audit-paired-tail-branch-010.json` |
| Relaxation limitation | One joint-marginal 010 input has an exact negative-cost feasible distribution | Counts (7,16), higher pattern 1110; rational cost -6521997797/1000000000, 51 states, all 4079 inequalities checked. Not a covering system. Paired projection constraints overcome this particular limitation. | `data/audit-joint-tail-010-primal.json`, `notes/eight-prime-anchor010-2026-10-25.md` |
| Upper | **M_sf <= 7** (C0 = 8) | Certified 26 Oct by two integer-arithmetic implementations: total <0.999087538091. Per-prime breakpoint optimization, distribution through prime 3000000, explicit infinite tail. All exported delta<=1/2. No new comparison lemma required. | `notes/upper-bound-seven-2026-10-26.md`, `data/upper-bound-seven-summary.json`, `data/cert7-breakpoints.json`, `data/audit-cert7-breakpoints.json` |
| Scoped upper-method optimality | **No bound below 7 from the same catalogue comparison with 0<=delta<=1/2** | At cutoff 6, the finite loss-estimate sum through prime 5003 is at least 1.010143243 for every real schedule. A separate computation gives 1.01255014. Majorization, 60 exhaustive prefix cases, and a downward-rounded scalar-state dynamic program; no solver optimality assumptions. Does not exclude delta>1/2 or stronger bounds. | `notes/stoploss-half-range-barrier-2026-10-26.md`, `data/stoploss-half-barrier-summary.json` |
| Full-parameter upper-method optimality | **No bound below 7 from the same catalogue comparison even for all 0<=delta<1** | For cutoff 6, the original loss-estimate sum through prime 10007 is at least 1.029993468; separate audit lower bound 1.03219703. All 326 action eliminations and 32 initial cases verified, with no floating-point or solver-optimality assumptions. Stronger comparisons and actual-family information remain outside the theorem. | `notes/stoploss-full-range-barrier-2026-10-26.md`, `data/fullrange-stoploss-barrier-summary.json` |
| Earlier upper bound | M_sf <= 10 | Superseded by 7. The previous certificates remain valid: total 0.92457, or 0.92265 with a longer distribution range. | `data/cert10.json`, `data/cert10-full.json`, `data/audit-cert10.json` |
| Earlier numerical limitation | Old C0 = 8 schedule gave 1.017877 | Superseded by the successful per-prime schedule. The fixed-parameter small-stage sharpness experiments remain valid but did not prove a global limitation. | `notes/scope-of-upper-bound-ten-2026-10-25.md`, `data/schedule8-float-only.json` |
| Lower | M_sf >= 2 | Krukenberg 1971 (lcm 210); 14- and 15-class systems mod 210 verified | literature, `data/cpsat-lcm210-min2-sat.json` |
| Earlier obstruction | All supports of at most seven primes excluded | Two integer-only routes exclude 510510: the original plane-load proof and a quantitative core/compression proof leaving at least 510 residues uncovered. Both have independent arithmetic audits. Superseded as a support threshold by the eight-prime obstruction above. | `notes/seven-prime-obstruction-2026-10-25.md`, `notes/quantitative-core-extension-2026-10-25.md`, `data/quantitative-seven-prime-obstruction.json` |
| Residual minimum | u(2310)=121 | CP-SAT OPTIMAL for five anchor patterns; audited counting bounds already exceed 121 for the other three. All five witnesses checked on the complete period. Independent SAT proof-log validation of the 121 lower bound remains unfinished. The six-prime obstruction does not depend on this new solver result. | `data/min-uncovered-2310-five-patterns-optimal.json`, `notes/seven-prime-investigation-2026-10-25.md` |
| Residual lower bounds | Integer-only u(2310)>=119, u(30030)>=501, u(510510)>=510 | Retain all raw core intersection profiles and the shared modulo-3 label of 3q, then use compressed parity propagation. Independent audit verifies 40 frontiers; a separate rational dual calculation checks 11332 input pairs in all 32 final patterns. Including solver optima gives u(30030)>=510, but those inputs are not used for the seven-prime bound. | `notes/quantitative-core-extension-2026-10-25.md`, `data/audit-core-capacity-frontiers.json`, `data/audit-hybrid-slice-510510.json` |
| Residual bracket | 510<=u(510510)<=7436 | Upper bound from an actual 126-class family, independently checked over all 510510 integers. Its parity counts are (2779,4657). Neither endpoint is asserted optimal. | `data/min-uncovered-510510-witness.json`, `data/audit-partial-families.json` |
| Fixed eight-prime core | All 16 high-prime anchor patterns excluded for the 22-point core of the current partial family | Conditional tail capacities, exact subset concentration, and independently audited integer duals give at least 130 uncovered points for every pattern. Pattern 1111 has the stronger audited bound 15257 when point identities are retained. This earlier fixed-core result is subsumed qualitatively by the complete eight-prime obstruction. | `notes/conditional-tail-obstruction-2026-10-25.md`, `data/audit-conditional-tail-core-patterns.json` |
| Eight-prime branch | All cores with b3=b5=b7 are excluded at 9699690 | Portable conditional-tail duals and arbitrary-subset compression give at least 6696 uncovered residues for all 16 higher patterns. Independent standard-library audit checks 15819570 ordinary core combinations and all 144 frontier/pattern bounds. | `notes/eight-prime-anchor-branch-2026-10-25.md`, `data/audit-eight-prime-anchor000.json` |
| Eight-prime branch | Entire 001 branch excluded at 9699690 | The count-only certificate leaves at least 305 points uncovered for every core and all 16 higher patterns. Independent audit checks 12495350 core combinations and 1004175 local plane tuples. This is one of the four branches in the complete obstruction. | `notes/eight-prime-anchor001-2026-10-25.md`, `data/audit-compressed-tail-branch-001.json` |
| Residual bracket | 2<=u(9699690)<=79979 | Lower bound from the four branch audits. Upper bound from an actual 254-class family checked independently over the complete period; parity residual counts (30484,49495). Neither is asserted optimal. The witness also limits the unconstrained cardinality-only compression recurrence to at most 13 through 23, then zero through 29. | `data/eight-prime-obstruction.json`, `data/min-uncovered-9699690-witness.json`, `data/audit-extension-witness.json` |
| Lower target | A system with least modulus >= 3, or a global impossibility proof | Unresolved. All supports of at most eight primes are excluded. The first remaining primorial is 223092870, through prime 23. The finite obstruction does not control arbitrarily many prime factors. Older searches and weaker relaxations are retained as historical evidence. | `notes/eight-prime-obstruction-2026-10-25.md`, `PROGRESS.md` |

## Timeline

| Date (2026) | Side | Result | Files |
| --- | --- | --- | --- |
| 26 Oct | lower | Direction 3 of the direction-analysis note (aligned construction test, "G-residues zero") pinned down exactly: pure integer counting reproves the already-known k<=6 obstructions in one line, but the allowed residual exceeds the certified u(L) lower bound at k=7 and at the open k=8 case, so it is undecided there. The naive residue-0 realization of the construction is off by two orders of magnitude (misses the whole coprime-to-L class). Cross-checked by an independent DP re-derivation and the `auditor` subagent from the note's prose alone. | `notes/aligned-construction-test-2026-10-26.md`, `code/aligned_construction_test.py` |
| 26 Oct | lower | Direction 4 (take M_sf>=3 seriously, search constructively): a one-pass numpy greedy at the first open primorial 223092870 (nine primes) gives a new heuristic upper witness u(223092870)<=1325141 and a slice-capacity ratio cap(R)/\|R\|=7.78 against the 28-fold coverage a tenth-prime extension would need -- 28% of the requirement, the lowest fraction yet in the widening-shortfall table (compare 43%, 38%, 40%, 36% at five through eight primes). No construction found; evidence continues to favor the conjecture, not a proof. | `notes/nine-prime-greedy-construction-2026-10-26.md`, `code/greedy_cover_nine_primes.py` |
| 25 Oct | upper | Stochastic-orders check of the comparison lemma: it is the Bernoulli case of the comonotonic convex-order theorem (Meilijson-Nadas 1979, Kaas et al. 2002) and is attained within its information class, so no result of that literature sharpens the per-stage bound (variance constraints give no convex-order improvement, Cheung-Vanduffel 2013). The only discarded information is the exact union identity P_n(U B_k) = sum_k (charge on fibres not yet covered). Exhaustively over stages 5 and 7 at C0 = 8 it is worth 0.0015 against the deficit 0.0179; adapted stage-11 adversaries lose 2-5 % of their charge. Lead for a multi-stage lemma, not a route to 7 with present tools. | `notes/stochastic-orders-and-union-accounting-2026-10-25.md`, `code/union_accounting.py`, `data/union-accounting-C0-8.json` |
| 25 Oct | theory / lower | Replaced the 16,583-step SAT proof of the sharp event bound 6 by three main cases using 28 elementary capacities. Independent integer enumeration checks all 3,072 event profiles and an equality witness. Transferred fifteen bounds to 327 group constraints; forty reject the old fractional witness, but a new 251-state negative witness passes the full independent audit. All 90 regression tests pass. | `notes/elementary-event-bound-2026-10-25.md`, `data/audit-small-event-bound-human.json`, `data/audit-nine-prime-transferred-group-primal.json` |
| 25 Oct | lower | Completed all four eight-prime branches. New paired-capacity audits exclude 010 and 011 with bounds 15 and 2, respectively; every counterexample now needs at least nine prime factors. Complete ordinary-residue and integer dual audits pass, along with all 58 regression tests. | `notes/eight-prime-obstruction-2026-10-25.md`, `data/eight-prime-obstruction.json` |
| 25 Oct | lower | Excluded the entire 001 low-anchor branch with uniform uncovered lower bound 305, using only the core count frontier and compression capacities. The independent integer audit covers every core and all higher patterns without orbit reduction. All 50 regression tests pass. The old 167-input 001 screening is superseded for this theorem. | `notes/eight-prime-anchor001-2026-10-25.md` |
| 25 Oct | lower | Excluded the entire normalized 000 low-anchor branch at 9699690, with uniform uncovered lower bound 6696 and an independent ordinary-residue/integer audit. Core orbit reduction and portable screening prepare the remaining branches. Eight additional generated duals reduce 167 unresolved 001 capacity inputs to 41; these new duals still await independent audit. | `notes/eight-prime-anchor-branch-2026-10-25.md` |
| 25 Oct | lower | New conditional-tail inequality evaluates higher-codimension capacities after plane deletions. Independent enumeration of all 4194304 core subsets and integer dual audits exclude all 16 extensions through 19 of one specified 210-core. Uniform lower bound 130; stronger 1111 bound 15257. All 40 regression tests pass. The global eight-prime case remains open. | `notes/conditional-tail-obstruction-2026-10-25.md` |
| 25 Oct | lower | Quantitative core intersections and shared odd-3q choices raise integer-only residual bounds to 119 at 2310 and 501 at 30030. A scalable parity-allocation bound plus independent rational duals proves u(510510)>=510. An actual seven-prime family leaves 7436 points. All 34 regression tests pass. The next-prime relaxation still permits zero from the unverified residual pair (102,470). | `notes/quantitative-core-extension-2026-10-25.md` |
| 25 Oct | lower | Proved a sharp arbitrary-set cylinder-capacity formula by coordinate compression. Parity propagation gives integer-only u(2310)>=106 and u(30030)>=411; solver inputs raise the latter to 488. Independent ordinary-residue/histogram audit verifies all 40 frontiers. Actual eight-prime partial family leaves 79979 uncovered points, independently counted. All 29 regression tests pass. | `notes/compressed-residual-capacity-2026-10-25.md` |
| 25 Oct | lower | Solver-free obstruction at 510510: all 32 normalized anchor patterns excluded, so every counterexample needs at least eight prime factors. Independent arithmetic audit checks 2936 base and 116 additional cases; independent ordinary-residue audit checks completeness. All 21 regression tests pass. | `notes/seven-prime-obstruction-2026-10-25.md`, `data/seven-prime-obstruction.json` |
| 25 Oct | lower | Audited exact core-geometry catalogues reduce four patterns to 60, 66, 115, 126 optimistic capacity profiles. Weighted capacity excludes 298 of 2936 profile/pattern combinations, but leaves all 32 seven-prime anchor patterns open. A full-period-verified partial family at 30030 leaves 792 residues uncovered; no optimality claim. | `notes/core-geometry-2026-10-25.md` |
| 25 Oct | lower | Solver-certified exact u(2310)=121; the propagated bound for one 30030 pattern rises from 74 to 282. Linked projection models remain inconclusive at 510510. A separate SAT/RUP path verifies the 210 bound; its 2310 audits time out. | `notes/seven-prime-investigation-2026-10-25.md` |
| 25 Oct | lower | Solver-free obstruction at 30030, hence every counterexample needs at least seven prime factors. Two independent integer implementations agree. Exact u(210)=22 verified without a solver; u(2310)>=90. Corrected an unsound timeout-as-infeasible exclusion in the auxiliary residual loop. | `notes/six-prime-obstruction-2026-10-25.md`, `data/audit-six-prime-obstruction.json` |
| 23 Oct | upper | The convex stop-loss comparison of the general problem, applied to squarefree moduli, gives minimum modulus <= 14 numerically (published 118). Found by the "transfer" agent of the multi-agent run in the `minimum-modulus` repository. | `notes/sqfree-agent-report-2026-10-23.md` |
| 24 Oct | upper | Certificates C0 = 14 (total 0.8234), C0 = 12 ("<= 11", 0.8856), C0 = 10 (0.9246). Verifier validated on CFT's own delta table at C0 = 118 (stages 1-7 exactly 194/1365, tail 0.000402960697 vs CFT 0.000402960685). C0 = 8 fails numerically at 1.018-1.020. | `data/cert14.json`, `data/cert12.json`, `data/cert10.json`, `data/cft118-check.json` |
| 24 Oct | upper | Independent audit: a second verifier written from the report alone reproduces total 0.924567866690 (12 digits) and the part sums; outward roundings differ beyond the 20th digit. | `code/audit_verify.py`, `data/audit-cert10.json` |
| 24 Oct | upper | Exact adversary over all stage-5 families (orbits under Aff(Z_2) x Aff(Z_3) x Aff(Z_5)) and all 810,000 stage-7 residue configurations at C0 = 8: max loss_5 + loss_7 = 409327/1932805 = 0.2118 at (D5, D7) = (0.195, 0.1425). Run interrupted (machine out of memory) before the stage-11 continuation. | `code/exact_two_stage.py`, `code/stage11.py`, `data/exact-two-stage-*.json` |
| 24 Oct | lower | Brute CP-SAT over squarefree lcm with minimum modulus 3: 210 UNSAT (0.6 s; min 2 SAT with all 15 moduli, Krukenberg type), 2310 and 30030 UNKNOWN after 900 s with 8 workers. | `data/cpsat-lcm210-*.json`, `code/cpsat_cover.py` |
| 24 Oct | repo | Project split out of `minimum-modulus` into this repository. | |
| 24 Oct | lower | Slice extension: a cover with lcm Lq (q a new prime) consists of a q-free family over L with residual R plus the shapes qd (d | L) distributed over the q slices; one slice is covered by q itself, each other slice must cover R. Minimal residuals at 210 (22 points) and 2310 (121 points) are spread, so Sum_d max_a |R cap (a + dZ)| is about 4.5 |R|, far below the 12 |R| needed at q = 13. The most concentrated residual with >= 22 points at L = 2310 is a full class mod 105 (Codex's bound 22 = 2310/105 is exactly this size), and it is the ONLY single-class residual that passes the count Sum_d gcd(m,d)/d >= 12 (12.09); the slice assignment for it exists (CP-SAT, 4 s). Hence a cover with least modulus 3 and lcm 30030 exists if Z_2310 can be covered by all divisors >= 3 with the modulus 105 used twice. That test is UNKNOWN after 300 s, but it is impossible by the same slice count one level down: both copies of 105 are 11-free, so the 11-free part is 15 progressions over Z_210 missing >= 22 - 2 = 20 points, and f(210, 20) = 183 < 10 * 20. |
| 24 Oct | lower | Trivial slice count (`code/slice_count.py`): a cover at L q needs f(L, u(L)) = Sum_{d>1, d | L} min(u, L/d) >= (q-1) u(L). With u(210) = 22 (CP-SAT OPTIMAL): f = 194 < 220, a one-line reproof that lcm 2310 has no cover with least modulus 3. For 2310 -> 30030 the count would need u(2310) >= 212, but u(2310) <= 121, so 30030 needs the direct search or a sharper residual bound. All single-class residuals at 30030 passing the count for q = 17 (m = 210, 330, ..., 1001) admit slice assignments. | `code/slice_count.py` | `code/slice_extend.py`, `code/duplicate_test.py`, `data/duplicate-2310-105-300s-unknown.json` |
| 24 Oct | upper | Exact worst case at p = 5, 7 (C0 = 8) equals the comparison value: stage 5 exactly 2/23 both, stage 7 0.1248222 exact vs 0.1248289 comparison. Route closed. Truncation ceiling at C0 = 8 computed: 0.905 total, savings per stage 0.005 (11), 0.008 (13), 0.004 (17), 0.005 (19). | `code/truncation_ceiling.py` |
| 24 Oct | upper | Stage-11 adversary (CP-SAT 540 s) for the aligned stage-5/7 measure at C0 = 8: a spread family realises 0.055092 against comparison 0.055291 and the attempted truncated model 0.050363. This defeats that proposed improvement under this prefix measure; it is not an absolute optimality theorem for cutoff 10. | `data/stage11-relaxation-tower-C0-8.txt`, `code/stage11.py` |
| 24 Oct | lower | Minimum uncovered count for least modulus 3: lcm 210 exactly 22 of 210 (10.5 %, OPTIMAL 0.3 s); lcm 2310 at most 121 of 2310 (5.2 %, FEASIBLE after 540 s, solver bound 19; Codex's rigorous bound 22). | `code/cpsat_min_uncovered.py`, `data/min-uncovered-lcm2310-feasible121.json` |
| 24 Oct | lower | Pigeonhole: some non-free slice holds at most floor((tau(L)-1)/(q-1)) shapes, so the q-free residual lies in a union of that many progressions. As a WLOG constraint in the direct model (`--small-slice k`) it proves 2310 in 0.6 s (k = 1). Phase-1 test at 2310 -> 30030 (k = 2): families with residual inside two progressions exist (e.g. 0 mod 2 union 2 mod 3, 156 points), so the pigeonhole alone does not settle 30030. | `code/few_progression_residual.py`, `data/few-progression-2310-13-k2-phase1-feasible.json` |
| 24 Oct | lower | Lossless normal form added to the CP-SAT search; lcm 2310 with least modulus 3 proved infeasible in 27 s. | `code/cpsat_cover.py`, `data/cpsat-lcm2310-min3-normalized-unsat.json` |

## Closed routes (do not retry)

* Exact adversary at the stages p = 5 and 7 for C0 = 8: the comparison is already tight there (the aligned tower realises it), gain 7e-6.
* Truncation alpha <= 1 at p >= 11: spread families realise 99.6 % of the untruncated comparison at stage 11 (same phenomenon as in the general problem). No successful route below 10 has been demonstrated with these tested variants.

* Parameter tuning of the convex comparison below C0 = 10: the floating-point optimum at C0 = 8 is 1.018 with per-stage bands; full exactness (-0.002) and a larger N before the tail (-0.0003) do not close the gap.
* The tested low-cutoff moment estimates remain insufficient. This is not a
  proof that all moment methods fail. A separate exact barrier applies only to
  the factored quadratic continuation specified in `notes/quadratic-tail-barrier-2026-10-25.md`.

* Iterating only the sharp arbitrary-set cardinality bound from period 9699690:
  even the exact input minimum cannot yield more than 13 through 23, then the
  next extension through 29 gives zero. This limitation is specific to that
  recurrence; it does not apply to bounds retaining residual geometry or parity.

## Open directions

* Five concrete next steps for the lower side, not yet exhausted, are listed
  in `notes/direction-analysis-2026-10-25.md`: (1) a uniform inductive
  invariant via the slice-capacity relaxation's maximiser structure
  (Codex is actively iterating on this one, see the uncommitted
  `code/extremal_slice_capacity.py`/`paired_slice_capacity.py`; do not
  duplicate), (2) the static two-family model (F covers Z_L\R, a separate G
  using every divisor of L covers R at least q-1 times), (3) the aligned
  (G-residues-zero) special case of (2), now resolved as far as pure counting
  goes -- see `notes/aligned-construction-test-2026-10-26.md`: obstructs
  k<=6, undecided at the open k=8 case, and needs the full model (2) to go
  further, (4) constructive heuristic search for an actual nine/ten-prime
  counterexample (take M_sf>=3 seriously) -- a first one-pass greedy at
  223092870 found no construction and continues the widening-shortfall
  trend (cap(R)/|R| = 28% of the ten-prime requirement, the lowest fraction
  yet in the table), see `notes/nine-prime-greedy-construction-2026-10-26.md`;
  a tuned local search at this scale (reusing the O(L) reshape trick instead
  of `improve_cover_witness.py`'s per-residue dictionary, which is not
  memory-tractable here) is the natural next step, (5) a polyhedral characterization
  of the nine-prime relaxation's valid inequalities instead of one-witness-at-
  a-time cuts.
* Upper: no tested stage-wise comparison certifies below 10. The examined
  adversaries nearly attain those estimates; this does not exclude stronger
  comparisons that retain additional residue compatibility or geometry.
* Lower: extend the conditional-tail argument beyond eight prime factors or
  develop an argument controlling arbitrarily many primes. The 9699690
  obstruction is complete; the first remaining primorial is 223092870.
  Do not restart the superseded 30030, 510510, or eight-prime branch searches.
* Residual geometry: the compression formulas are sharp for arbitrary sets of fixed
  size, so improving them requires restrictions on sets that can actually arise
  as congruence-family residuals. The parity frontiers are quantitatively valid;
  the old normalized capacity-profile dominance only preserves cover feasibility.
  Retaining core capacities and shared modulo-3 labels now gives a positive
  quantitative seven-prime bound. The old eight-prime relaxation admitting
  counts (102,470) is superseded by the complete conditional-tail proof.
* Conditional-tail bounds retaining paired core projections now cover every
  eight-prime core. For nine primes, additional high-prime groups must be
  included, and neither the existing four-prime local inequalities nor their
  positive margins transfer without a new proof and audit.
* All four eight-prime branches are independently excluded. The 010 and 011
  proofs retain pairwise projection geometry and have complete ordinary-residue
  orbit audits; they do not depend on the older inclusion reduction. A possible
  next finite step is a five-higher-prime conditional-tail model over the same
  210-core, with 32 higher parity patterns. New multi-prime generators and
  independent auditors include all 31 groups through 23. Exact negative-cost
  witnesses rule out several unchanged relaxations even for one actual core.
  Even the model with all pairwise threshold events has an audited negative-cost
  witness. Triple cuts and three four-point cuts have now been added, but a
  new independently audited negative-cost witness survives. The arbitrary-
  subset separator reports no remaining cuts of this sum-capacity form; an
  exact joint-triple check also finds none. Whole-group event separation now
  gives fifteen proof-checked inequalities. A new negative-cost witness still
  survives them, with all proofs and arithmetic independently audited. The
  stronger explicit group-mixture model remains inconclusive. Its latest
  decomposed checkpoint is `data/nine-prime-pointwise-catalogue-pricing.json`;
  pass it as both the pool and warm start to preserve assignments and cuts.
  Use `--catalogue data/nine-prime-group-catalogue.json` for shared assignment
  pricing. Only the older 1,543-entry pool is rigorously excluded; the expanded
  shared catalogue is not. Fixed-point reweighting now has an explicit audited
  limitation, even with exact integer local pricing. Separate group coefficients
  give a positive obstruction for the specified 40,097 configurations, but
  the expanded 40,197-entry pool remains unresolved. Certified residue case
  trees substantially improve unrestricted maxima. Fractional case scenarios
  are now incorporated into coefficient optimization, but the audited new
  unrestricted excess remains -545194.768584. Repricing exposes omitted
  fractional states, not covering assignments. Testing coefficients down to
  zero now leads to an exact fixed-curve case-relaxation barrier, independently
  verified by rational mixtures. Do not repeat coefficient-only optimization
  for these curves/cases. Free threshold curves now defeat that specific
  witness, but exact local pricing exposes new weak states and both new
  candidates fail against the actual catalogue. A three-curve basis now
  closes sampled constraints and saves 2,895 active local states, but exact
  catalogue excess remains negative. The old fractional witness violates
  an independently checked pair incompatibility. Pair-aware affine upper
  bounds are now independently certified, and 2,793 incompatible scenarios
  have been repaired. The repaired-pool run retains 3,381 local states but
  still has negative exact catalogue excess. Continue with diverse valid
  local states, the 41,375-state separation pool, and pair-aware case pricing.
  See `notes/pair-affine-certificates-2026-10-25.md`.
  Catalogue continuation uses `data/nine-prime-group-support-third-resume-input.json`
  and the support seed in `data/nine-prime-three-curve-pair-repaired-support.json`.
  Use `data/nine-prime-pair-repaired-scenarios.json` in place of the old
  pair-incompatible scenario union.
* Write-up: a short note (118 -> 10) with the argument, the inequality list and the certificates.
