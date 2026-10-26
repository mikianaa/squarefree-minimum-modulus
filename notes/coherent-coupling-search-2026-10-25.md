# Decomposed coherent-mixture searches

The global conjecture and the nine-prime obstruction remain unresolved.
The previous continuation proved an elementary event bound and audited a
negative witness satisfying 327 event inequalities. This continuation tests
the stronger requirement of explicit coherent resource assignments per group.

## Two assignments per group

Fix two actual residue assignments A_h,B_h for each higher-prime group h,
with a fixed mixing probability. At a core point x these supply two possible
hit counts. Local counts need only be stochastically dominated by that supply.

For each of the five plane counts it suffices to consider 0 and the two
supplied counts, after clipping at the number of available higher residues.
Indeed, the cost is affine in each plane count when all other counts are
fixed. Counts between consecutive supplied endpoints can be split between
the endpoints without changing the conditional mean or cost, and without
violating domination. This leaves at most 3^5=243 plane states per point.

Each tail count has a nonpositive coefficient, so its full two-point supply
distribution can be used. For every plane state j introduce its probability
m_j, and, for each tail h, the mass y_jh using the higher tail count. Impose
0<=y_jh<=m_j and sum_j y_jh equal to the supplied high-count probability.
The cost is linear in these variables: it contains no products of tail counts.
For a fixed plane state, the tail marginals can be realized simultaneously
by nested intervals [0,y_jh] within [0,m_j]. Thus at most 27 joint tail
states suffice per plane state.

The point problems are independent once the group mixtures are fixed.
`code/two_configuration_mixtures.py` implements this reduction. Its small
three-plane test agrees with a separate LP enumerating every admissible count
vector, including intermediate counts omitted by the reduced model. Exact
rational witnesses are repaired and checked by `code/audit_coherent_groups.py`.

A 120-second run examined 482 choices of second assignments and common
mixing probabilities. None improved the baseline. Allowing valid downward
changes to plane counts already lowers the earlier integer-projection
expression 678850 to **670035**. The attaining rational witness is audited
using 31 groups and all 10,692 threshold links. It has 22 point states.
Its two assignments coincide, so this improvement does not require a
nontrivial mixture. This is an upper bound on the relaxed minimum, not an
obstruction or a claimed optimum.

* `data/nine-prime-two-configuration-mixtures.json`
* `data/audit-nine-prime-two-configuration-mixtures.json`

## Fitting the negative point witness

`code/fit_coherent_marginals.py` pools configurations across groups with the
same resources. The saved pool contains 237 distinct plane configurations
and 1,185 distinct tail configurations. For each group, a small LP maximizes
gamma<=1 such that its configuration mixture dominates gamma times the event
marginals of the negative 327-cut witness.

The restricted-pool gamma values range from about .453 to .627. After exact
per-point repair, an explicit coherent witness has positive cost

    7224066529508926 / 1000000000.

Its 1,100 group assignments and all 10,692 threshold links pass the separate
integer audit. Failure of this finite pool to realize the target is not a
proof that no other assignments could realize it.

* `data/nine-prime-coherent-marginal-fit.json`
* `data/audit-nine-prime-coherent-marginal-fit.json`

One-second CP-SAT pricing of the marginal-fitting duals supplies additional
actual configurations and twenty candidate separating event inequalities.
These candidates are saved in `data/nine-prime-coherent-marginal-pricing.json`.
Their upper bounds have not been independently proof-checked, so they are
exploratory and have not been inserted into the certified 327-cut model.
All **94** regression tests pass with a clean exit.

## Independent mixing probabilities

The initial two-configuration search used a common fixed mixing probability.
The implemented stronger approach optimizes a separate theta_h in [0,1]
for each group, with the actual pair A_h,B_h fixed. Each point LP has a
right-hand side affine in theta, so its optimum is a convex piecewise-linear
function of theta. Its dual gives a supporting affine function.

A Benders master needs only 31 mixture variables and 22 point-cost
variables, adding supporting inequalities after solving the independent
point LPs. Plane and tail marginal duals give the coefficient of theta_h,
with sign determined by which assignment supplies the larger count. This
can improve mixtures that uniform probabilities made unattractive, while
avoiding the slow full configuration master. Final probabilities are rounded,
the primal witness reconstructed, and the existing integer audit applied.

`code/optimize_two_configuration_mixtures.py` implements this step. For the
selected pair in each group, ten iterations converge numerically to 670035,
choosing the original assignment in every group. The rational witness passes
the independent audit. This is a conclusion about those fixed pairs only.
See `data/nine-prime-optimized-two-configurations.json`.

## Larger configuration pools

`code/decomposed_coherent_groups.py` permits arbitrary finite assignment
pools and generates local point columns by full integer pricing. The first
120-second calculation had a persistent large oracle gap: its finite-pool
lower estimate stayed at -11803760 while its best feasible value was 678850.
This exposed poor numerical dual choices, not a covering obstruction.

A new exact-penalty lemma bounds every event dual by the cost of deleting
one hit. The associated integer quantile repair constructs feasible point
distributions without the earlier global shrink-to-zero step. The proof and
implementation are in `notes/coherent-transport-penalty-2026-10-25.md`.
The first bounded run improved the numerical finite-pool lower estimate to
about -1.013 million. It still had no negative coherent witness.

The generator now also combines the Benders duals into weighted group-event
objectives and searches for missing actual assignments. It rounds event
weights upward, prices every local point state with integer arithmetic, and
records CP-SAT upper bounds for the group objectives. These produce exploratory
global lower estimates only: independent group-bound proofs would still be
required for an obstruction. Any discovered actual assignments can safely be
added without relying on the CP-SAT upper bounds.

For continuation, pass the latest complete checkpoint as both `--pool` and
`--warm-start`; this preserves the group assignments, point states, and cuts.
The original source and integer assignment input remain unchanged.

The latest complete checkpoint is `data/nine-prime-decomposed-priced-model.json`:
66,247 point states and 1,267 group configurations, including 71 new actual
assignments found by group pricing. Its numerical finite-pool lower estimate
is -602445.9820693743 and its audited feasible value remains 670035. No
complete global lower estimate was obtained because some short group-pricing
runs did not return an upper bound. The separate standard-library audit is
`data/audit-nine-prime-decomposed-priced-model.json`.

An intermediate unscaled continuation stopped on a numerical solver assertion
and produced no checkpoint. Objective rescaling removed the observed failure;
the final 189-second run completed with zero numerical retries. All 101
regression tests pass in `data/tests-coherent-transport.log`.

The next structural refinement is to require the necessary inequality in
each core fibre, rather than permitting negative expectations at some points
to cancel positive expectations at others. Exact diagnostics show that the
old negative 327-cut witness has thirteen positive point expectations.
See `notes/pointwise-cover-necessity-2026-10-25.md` before extending the raw
sum-cost search again.

The subsequent pointwise continuation is complete and supersedes this next
step: compact transport and integer dual audits now exclude an explicitly
listed pool of 1,543 group configurations for the fixed core/pattern. The
unrestricted case is still open. Resume the newer checkpoint
`data/nine-prime-pointwise-global-pricing.json` and follow
`notes/pointwise-finite-pool-certificate-2026-10-25.md`.
