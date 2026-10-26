# Exact limitation of the fixed-curve residue-case relaxation

The global covering conjecture remains unresolved. The result below is a
limitation of one specified relaxation for one 22-point core and higher-prime
pattern 11111. It is not a limitation of the distortion method, a proof that
the minimum-modulus upper bound 10 is optimal, or a covering construction.

## What is now proved

Keep the original threshold reward curves Q_(h,x)(k), for 31 groups and 22
points, and the 8,176 residue cases of
`data/nine-prime-branched-event-depth6.json`. Within each case use the
fractional resource/count relaxation defined in
`notes/fractional-case-scenarios-2026-10-25.md`.

Allow a separate arbitrary nonnegative coefficient alpha_(h,x) on each of
the 682 reward curves. These coefficients cannot give a positive covering
obstruction with these group relaxations, even if the local integer
nonpositive-cost state minima and the case LP optima are computed exactly.
The lower restriction alpha >= 0.5 from the preceding experiment is not
responsible for this failure.

The exact certificate consists of distributions on 243 valid integer local
states and 214 valid fractional group-case states. These distributions give
local expected scores D_(h,x) and group expected scores S_(h,x) satisfying

    S_(h,x) - D_(h,x) >=
        (128598143 / 20000000000) Q_(h,x)(maximum count) > 0

at all 682 coordinates. The normalized surplus is approximately 0.0064299.
All probabilities, inequalities, and local costs are checked as rational
numbers. No solver optimum is trusted.

## Why the certificate blocks every such coefficient choice

Let A_x(alpha) be the exact minimum weighted reward over integer local
states with F_x <= 0. Every local state in the certificate satisfies this
condition, so its mixture bounds that minimum from above:

    A_x(alpha) <= sum_h alpha_(h,x) D_(h,x).

Let B_h(alpha) be the maximum of the group-case LP optima. Each fractional
state in the certificate belongs to one listed leaf case, hence

    B_h(alpha) >= sum_x alpha_(h,x) S_(h,x).

Summing gives

    sum_x A_x(alpha) - sum_h B_h(alpha)
       <= sum_(h,x) alpha_(h,x) [D_(h,x)-S_(h,x)] <= 0.

The inequality is strict for any nonzero coefficient vector in this
instance. Zero coefficients achieve margin zero. This proves the optimum
margin is zero for this family; it is stronger than observing a numerical
optimizer return zero. The proof does not need to enumerate every integer
local state, because exhibiting valid states already supplies the needed
upper bound on its minimum.

The current implementation rounds each group's real affine upper bound
down to an integer. This does not rescue it. The smallest absolute surplus
at any active coordinate is

    14843440860113154396329 / 750000000150000000000 > 19.

For integer coefficients, each active group contributes at least this much
surplus, whereas unit-floor rounding saves less than one per active group.
Common-denominator scaling covers rational coefficients with that rounding
convention. Each child-case relaxed state also embeds in its ancestor case
by making newly fixed resources deterministic, so taking minima of inherited
affine bounds does not evade this argument. Different stronger arithmetic
bounds or deeper cases are outside this limitation.

## Independent verification and files

`code/case_relaxation_barrier.py` constructs a mixture numerically, then
rationalizes its probabilities. If rounding would make a count distribution
consume too much resource supply, it moves mass toward its minimum count.
The resulting exact score dominance is checked before saving.

`code/audit_case_relaxation_barrier.py` uses only the Python standard library.
It independently rebuilds all rewards and local costs, checks each probability
sum and resource-demand inequality, checks the complete topology of all
10,341 tree nodes and 8,176 leaves, and verifies the 682 exact score balances.
It imports neither the generator nor an optimization solver.

* Certificate: `data/nine-prime-case-relaxation-barrier.json`.
* Audit: `data/audit-nine-prime-case-relaxation-barrier.json`.
* Zero-floor search: `data/nine-prime-zero-floor-support-weights.json`.
  It completed two iterations, returned margin zero, and retained 1,126
  valid local support states.
* Regression suite: 141 tests pass in
  `data/tests-case-relaxation-barrier.log`. The three barrier tests were
  also rerun after adding the integer-rounding and threshold diagnostics.
  They reject missing cases, invalid local states, invalid probability mass,
  and false coordinate surpluses.

## The barrier identifies a route that remains open

This witness does **not** dominate each individual threshold probability.
The audit finds 4,235 coordinates (group, point, threshold) where its group
mixture has smaller probability than its local-state mixture. For example,
at residual point 13 in the group for higher primes {13,17,19}, the score
balance for the event "count at least 7" is

    supply - demand = -12185291157 / 12500000000.

The original curve assigns that event weight zero. Thus the mixtures that
defeat all rescalings of the old curves do not defeat arbitrary new curves.
This is a concrete deficit in this witness, not a proof that changing the
curve will succeed: a different feasible mixture may survive the new event.
All exact threshold deficits are saved in the audit file.

The next search should allow individual threshold-event weights to change,
starting from these deficits, and retain exact local-state pricing plus
unrestricted group bounds. Alternatively it can strengthen count/residue
compatibility within the cases. Repeating coefficient-only optimization
with these same curves and cases cannot settle even the fixed-core case.

Continuation: free threshold curves now defeat this specific mixture in an
independent exact comparison, but fail after full local-state pricing.
This does not contradict the fixed-curve theorem above. See
`notes/free-threshold-support-2026-10-25.md` for the new certificates and
the smaller-basis search proposed next.
