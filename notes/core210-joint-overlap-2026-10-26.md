# Joint geometry after moduli 7 and 10

Date: 26 October 2026.

## Status

The unconditional upper bound remains **7**. Upper bound 6 has not been
proved, and no impossibility theorem for all improved comparisons or
all parameter schedules has been proved.

This investigation retains the entire uncovered CRT core modulo 210,
including the actual residues of moduli 7 and 10. It produces:

1. A sharp elementary theorem: at least **104 of the 210 core points**
   remain after all allowed moduli dividing 210 are used.
2. A complete, independently audited catalogue of the equality cases:
   **288846 coordinate-permutation orbits**, representing **35490420**
   assignments after the five-class normalization below.
3. A numerical joint-core comparison that is substantially better on
   tested configurations, but still has a difficult case near 1.003.
4. An exact obstruction for **one specified parameter schedule**: every
   initial probability measure on one actual 104-point core gives the
   defined stop-loss comparison functional at least **1.001192376**.
   This is not an all-schedule obstruction.

None of the numerical scores below 1 in this investigation is claimed
as a covering exclusion: the catalogue-pricing searches were heuristic,
and the independent CP-SAT pricing attempts did not close their bounds.

## Why the residual bound is exactly 104

The allowed divisors of 210 are

    7,10,14,15,21,30,35,42,70,105,210.

Missing classes may be added. Remove the class modulo 7, leaving 180
points. The remaining classes have the following maximum sizes there:

| Modulus | 10 | 14 | 15 | 21 | 30 | 35 | 42 | 70 | 105 | 210 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Outside the 7-class | 18 | 15 | 12 | 10 | 6 | 6 | 5 | 3 | 2 | 1 |

Their sum is 78. If either the 14-class or the 21-class lies inside the
7-class, losing its whole capacity already leaves more than 104 points.
Otherwise the classes modulo 10 and 21 meet in exactly one point
outside the 7-class, and so do those modulo 15 and 14. Group these as
two disjoint pairs of sets before applying the union bound. This gives

    uncovered >= 180 - (78 - 1 - 1) = 104.

The argument subtracts the two intersections safely even if they occur
at the same point, since the pairs use disjoint sets. An explicit family
in `data/core210-dense-prefix.json` attains 104. The integer audit checks
the full witness, capacities, and all 360 relevant normalized pair
intersections. The initial CP-SAT search had only FEASIBLE status; the
sharpness theorem comes from this elementary bound and the witness.

## Equality cases and their complete enumeration

Equality requires the two pair-unions and all six remaining parts to
be disjoint outside the 7-class. In CRT coordinates this forces:

* 10 and 14 to use different parity values;
* 15 and 21 to use different modulo-3 values;
* 10 and 15 to use different modulo-5 values;
* 7,14,21 to use three different modulo-7 values.

Coordinate permutations therefore normalize the five classes to

    7:0, 10:0, 14:1, 15:6, 21:16.

Their union has 83 points. The remaining moduli are
30,35,42,70,105,210, whose outside parts must have sizes 6,6,5,3,2,1
and be pairwise disjoint and disjoint from the normalized base.
The residual symmetry group has order 3!*4!=144.

| Through remaining modulus | none | 30 | 35 | 42 | 70 | 105 | 210 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Orbits | 1 | 3 | 3 | 9 | 114 | 3132 | 288846 |
| Assignments without quotienting | 1 | 7 | 60 | 420 | 9060 | 338004 | 35490420 |

The independent audit counts normalized assignments without symmetry
reduction, checks every representative for validity and lexicographic
canonicity, checks every stabilizer, rejects duplicate orbits, and
verifies that orbit sizes sum to the unquotiented count.

**This is a catalogue of equality cases only.** Prefixes leaving 105 or
more points are not covered by this enumeration and cannot be discarded
merely because their residual set is larger.

## The joint-core comparison

For a chosen actual residual set U, allow any probability Q supported
on U. For every divisor d of 210 choose one residue a_d and define

    f_a(x) = sum_{d|210} 1[x = a_d (mod d)],
    H_Q(t) = max_a E_Q (f_a-t)_+.

Here 1<=f_a<=16. At later primes p>=11, the earlier coordinates >=11
are compared with independent Bernoulli variables of probabilities
1/[p(1-delta_p)]. Keep the core itself intact. As in the core-30 proof,
Jensen over the 2^N core catalogues gives the valid stage upper estimate

    sum_n P(N=n) 2^n/[p(1-delta_p)] H_Q(p*delta_p/2^n).

This accounts for all initial overlaps by restricting Q to U. It still
allows the maximizing future catalogue to change with the stage and
the comparison count. It does not track overlaps between two later
actual stages.

Two parts of H_Q are easy:

    H_Q(1) = sum_{d|210,d>1} max_r Q(x=r mod d),
    H_Q(t) = (16-t) max_x Q(x)   for 8<=t<=16.

For the second identity, the uniform full-core comparison has masses
(48,92,56,13,1) at loads (1,2,4,8,16), divided by 210. It implies
sum_x(f_a-t)_+ <=16-t for t>=8. Multiply by max Q for the upper bound;
align every cylinder at a point maximizing Q for equality.

The intermediate thresholds require catalogue optimization. A linear
program chooses Q, and coordinate-ascent pricing adds actual catalogue
constraints. These LP optima and stopping tests are numerical only.

## Numerical experiments and their limits

The search generated 120 actual prefixes and selected difficult examples
for measure optimization. These samples are not the full equality
catalogue. A first 104-point example gave a priced score near 0.981534
at cutoff 100000, but other 104-point configurations were harder.

The difficult configuration retained in `data/core210-hard-prefix.json`
belongs to catalogue orbit 168749 (zero based), with normalized residues
(12,3,13,18,93,193) for (30,35,42,70,105,210).

| Calculation for this difficult case | Numerical result |
| --- | --- |
| Initial fixed continuation, optimized core measure, cutoff 100000 | 1.029632267475 |
| Alternating the core measure and controls through 10000, distribution range 3000000 | about 1.004151 |
| Use actual dyadic thresholds instead of interpolating integer envelope values | 1.003106638983 |
| Finite-pool continuous stop-loss hull and noninteger control breakpoints, fixed final measure | 1.002937636754 |

The later controls after 10000 were held at the stored core-30 schedule
in the fast alternating search. The final hull is exact for the finite
pool used, not for the unrestricted set of catalogues. All scores in
this table remain numerical, incompletely priced candidates.

Capping the first subsequent fibre load at 11 was also tried. At the
tested difficult optimum its maximizing catalogue already avoids the
affected high loads, so this did not reduce that objective. It is not
a proof that all collision-aware improvements fail.

Independent 50-second CP-SAT pricing runs for the first example found
matching witnesses but retained substantial gaps (threshold 2: about
1.303 versus upper 1.488; threshold 4: about 0.499 versus upper 0.705).
Timeout and FEASIBLE status are not treated as completed proofs.

## Exact fixed-schedule obstruction

Fix the schedule encoded by the integer atoms a_p in
`data/core210-fixed-schedule-lower-coefficients.json`: delta_p=a_p/p,
for the 216812 primes from 11 through 3000000. Define

    B(Q) = sum_p sum_n P(N_p=n) 2^n/(p-a_p) H_Q(a_p/2^n),

where N_p counts earlier primes starting at 11. This is the pure
stop-loss catalogue functional above, with the specified schedule. It
is not the actual covered probability, nor a minimization over all
possible comparison methods or parameter schedules.

Grouping its nonnegative terms gives a constant, a coefficient of
H_Q(1), coefficients of intermediate dyadic H_Q(t), and a coefficient
of max Q. The generator computes lower coefficients by downward-rounded
count probabilities at scale 10^30, count cap 60, discarding the
nonnegative remaining count/prime tails. A separate audit recomputes
lower coefficients at scale 10^24, count cap 50, rounding the two incoming
probability flows separately.

A rational dual has 123 positive terms, each a point indicator, a
cylinder indicator, or the stop-loss value of an explicitly listed
actual catalogue. Its coefficient budgets fit inside the separately
recomputed lower coefficients. At every one of the 104 core points,
the constant plus the weighted features is at least

    25029809419968797730091 / 25000000000000000000000
      = 1.00119237679875190920364 > 1.

Averaging this pointwise inequality over any Q proves B(Q)>1 for every
initial measure on this core. No solver optimum or successful pricing
termination is trusted in this audit; the listed catalogues give lower
constraints regardless of other omitted catalogues.

This shows that changing initial weights alone cannot rescue **this
fixed schedule and functional**. It does not exclude different schedules,
including the subsequently tested noninteger schedule, or joint
accounting of later stages. The unconditional upper bound remains 7.

## Reproduction

```
python code/core210_maximum_prefixes.py audit data/core210-dense-prefix.json data/audit-core210-residual104.json
python code/audit_core210_maximum_prefixes.py data/core210-maximum-prefixes.json data/audit-core210-maximum-prefixes.json
python code/audit_core210_fixed_schedule.py data/core210-fixed-schedule-dual.json data/audit-core210-fixed-schedule.json
tmp/research-venv/Scripts/python.exe -m unittest discover -s code -p test_core210_joint.py
```

The equality catalogue and both exact obstruction audits are retained
alongside numerical searches. Hashes and scopes are summarized in
`data/core210-joint-overlap-summary.json`.

## Later-stage overlap: a certified branch modulo 2310

Further work on 26 October retains the actual stage-11 classes, then optimizes
the probability measure on the uncovered points modulo 2310. **The global
upper bound remains 7.** The results here concern the specific core-210 family
in `data/core210-hard-prefix.json` and a specified stage-11 subfamily.

For this fixed 104-point core, the maximum sizes of cylinders modulo
2,3,5,6,7,10,14,15,21,30,35,42,70,105,210 are respectively

    55,46,30,24,24,17,14,11,9,6,6,5,3,2,1.

Their sum is 253. The modulus-11 class removes 104 of the 1144 lifted points;
each other modulus 11d removes at most the corresponding cylinder capacity.
Thus every extension at stage 11 leaves at least 10*104-253=787 points.
The saved extension attains every capacity with disjoint classes outside the
11-class. `audit_stage11_capacity.py` independently verifies this sharp
conditional minimum. It is not a minimum theorem over all core-210 families.

### A comparison with exhaustive finite verification

Let Q be any probability measure on the 787-point residual support. Consider
the catalogue F of all 32 divisor types of 2310, one cylinder per type. To
eliminate a prime coordinate p, write F=f0+S, where f0 is the catalogue of
types not containing p and S is the sum of active types containing p. After
forgetting their p-coordinates, these latter types give a second catalogue f1.
At each point in the smaller core let T be the total weight over p lifts and
C the largest weight of one lift. For increasing convex phi,

    sum_y w_y phi(k*(f0+S_y))
      <= (T-C) phi(k*f0) + C phi(k*(f0+f1))
      <= (T-C) phi(k*f0) + C/2 phi(2*k*f0) + C/2 phi(2*k*f1).

Indeed 0<=S_y<=f1 and sum_y w_y S_y<=C*f1; the secant bound proves the first
inequality, and midpoint convexity the second. The f1=0 case is immediate.
Applying this first to 11 and then to 7 produces four groups of weighted
catalogues on the core modulo 30, with scales 1,2,4. In each group, scale
terms share the SAME catalogue; only different groups are maximized
separately. All 2*3*5*6*10*15*30=810000 catalogues on this core are enumerated.

For phi(t)=(t-j)+, the resulting maxima define a convex stop-loss upper
envelope H(j). Chords between thresholds 0,1,2,4,8,16,32 give a probability
law supported on {1,2,4,8,16,32} dominating every possible core catalogue.
The second differences of the finer integer-threshold envelope give an
alternative law on {1,...,32}, used for numerical optimization.

The core proof uses rational point weights and integer enumeration. A
separate audit rebuilds the groups with explicit modular indexing, generates
catalogues by mixed-radix decoding, and directly evaluates all seven chord
knots. It does not import the generator's reduction or catalogue code.

The numerical search initially gave 0.903 with heuristic, incomplete pricing;
that number is **not a certificate**. The safe conditional comparison on
those same weights was 1.044. Optimizing this safe comparison gives about
0.971879 with the original tail parameters. A cutting-plane LP uses exact
exhaustive pricing but floating-point arithmetic; its optimum is only a
candidate. Rationalizing its point weights and reoptimizing the per-prime
parameters yields two integer continuation certificates, each below

    0.970803496584 < 1.

They process primes 13 through 3000000 by stop-loss comparison, continue
quadratic estimates to the millionth prime 15485863, then bound the infinite
tail. The primary verifier convolves the power-of-two core seed into a count
law at scale 10^60; the separate audit leaves the six core loads outside a
tail-only count law at scale 10^50. Both round upwards, bound truncation
remainders, and use separate prime/logarithm helpers.

This excludes all distinct squarefree covering extensions with moduli >6 of
the fixed 27-class reference family. It is compatible with the older
fixed-schedule obstruction: the starting core, comparison, and parameters
have changed.

### Six stage-11 residues can be left unrestricted

The branch was broadened by another measure optimization. Keep the core-210
family fixed and fix only these ten stage-11 classes:

| Modulus | Residue |
| --- | --- |
| 11 | 0 |
| 22 | 3 |
| 33 | 10 |
| 55 | 17 |
| 66 | 7 |
| 77 | 41 |
| 110 | 12 |
| 154 | 97 |
| 165 | 13 |
| 231 | 62 |

The six moduli 330,385,462,770,1155,2310 may each have ANY residue, or be
absent. For each free modulus d charge max_a Q(x=a mod d). Future distortion
steps preserve the core marginal Q, so these union-bound charges can be
added to the later-prime loss. Q is allowed to vanish on extra, artificially
chosen points; it only needs to be supported outside the fixed subfamily.
The reference family's other six residues therefore impose no assumption
on the actual covering system in this enlarged branch.

The rational six-free candidate has free-class charge about 0.028821198
and later-prime charge below 0.970709038377. The combined exact audit is
`data/audit-core2310-six-free-branch.json`; it checks the hashes of both
independent continuation certificates and the independent core audit.
The resulting exact total is

    624707097231637813275842198349656998003569289449038045633
    / 625000700000000000000000000000000000000000000000000000000
      < 0.999530236097 < 1.

This is a branch exclusion, not a proof for all stage-11 configurations,
all 288846 dense core-210 orbits, or the less dense core-210 configurations.
There is no new impossibility theorem for further overlap refinements.

### Reproduce the later-stage certificates

```
python code/audit_stage11_capacity.py
python code/audit_core2310_power_envelope.py data/core2310-six-free-envelope.json data/audit-core2310-six-free-envelope.json
python code/verify_core2310_continuation.py data/schedule-core2310-six-free.json data/cert-core2310-six-free-tail.json
python code/audit_core2310_continuation.py data/schedule-core2310-six-free.json data/audit-core2310-six-free-tail.json
python code/audit_core2310_free_classes.py
```

Candidate generators are `later_stage_overlap_probe.py`,
`optimize_conditional_lp.py`, `core2310_conditional_envelope.py`, and
`optimize_core2310_schedule.py`. Scalar survivor-mass tracking was also
tested (`data/core210-survivor-mass-probe.json`): it saved only about
1.7e-7 before that trial lost its surviving mass. This numerical failure
does not establish a method barrier. The successful branch calculation
uses actual residual geometry, not only the amount of surviving mass.

## Removing residue conditions by fresh-colour domination

A subsequent continuation enlarges the allowed initial support to the
complement of the actually fixed subfamily. The old 787-point support
unnecessarily excluded points belonging only to the artificially chosen
residues of free moduli. Releasing those zero-weight restrictions improves
the seven-free numerical bound from about 1.010273 to 0.984243. Exhaustive
LP pricing gives approximately 0.989861 with eight free classes, 0.991911
with nine, and 1.003985 with ten. These numbers alone are not certificates
or impossibility statements.

The ten-free case can instead be split by the residue of the modulus-110
class modulo 10. All ten cases have integer certificates; the largest is

    1994219089441354292935240698399114291973196526185554256949
    / 1999999880000000000000000000000000000000000000000000000000
      < 0.997109604548 < 1.

The full audit is `data/audit-core2310-split110-branch.json`. Its broader
projection statement requires only the following conditions at stage 11,
in addition to the fixed core-210 family:

| Modulus, if present | Required projection of its residue |
| --- | --- |
| 22 | 1 modulo 2 |
| 33 | 1 modulo 3 |
| 55 | 2 modulo 5 |
| 66 | 1 modulo 6 |
| 77 | 6 modulo 7 |

All components modulo 11 are unrestricted. Modulus 11 itself, all ten
remaining stage-11 moduli, and every later prime are unrestricted; omissions
are allowed. Thus the result is stronger than excluding just the ten
reference families used by the certificates.

### Why only fresh, distinct colours need to be checked

In the CRT representation, call the residue modulo 11 a colour. Fix base
projections of k<11 stage-11 classes. Assigning them distinct colours
removes exactly t(x) points from the 11 lifts of each base point x, where
t(x) is the number of active projections. Every other colour assignment,
or omission of classes, leaves at least s(x)=11-t(x) points.

The certified reference measure is uniform on the s(x) remaining lifts of
each x, with per-point mass w(x) and base mass q(x)=s(x)w(x). For another
colour assignment, spread the same base mass q(x) uniformly on its s'(x)
surviving lifts. Its largest fibre weight is q(x)/s'(x)<=w(x).

The conditional comparison above depends on the base mass T and a largest
fibre weight C. For each nonnegative convex phi with phi(0)=0, its
coefficient of C is

    [phi(2*f0)+phi(2*f1)]/2 - phi(f0) >= 0.

Consequently the reference envelope also bounds the new measure. For a
free class of modulus 11m, the new measure of any residue class is at most
the maximum, over a modulo m, of sum_{x=a mod m} w(x). A wholly unused
reference colour exists because k<11, so this projection bound is exactly
the reference maximum class mass. The direct audit checks fibre uniformity,
s(x)=11-t(x), distinct reference colours, and equality of the direct and
projection capacities. This is an analytic domination proof for every
colour assignment, not a sampling argument.

For modulus 110 there remain only ten base projections. A missing class
may be added arbitrarily, or treated using the same larger-fibre argument.

### Shared, independently checked infinite-tail coefficients

The same rounded parameter sequence is used for all cases. Instead of
repeating the million-prime computation for each initial law, two integer
implementations compute upper coefficients G(f) for f=1,2,4,8,16,32.
They use pure stop-loss through prime 3000000 and pure quadratic bounds
after it, including the infinite tail. For any certified seed probabilities
b_f, sum_f b_f G(f) is a valid continuation bound. Using a separate minimum
of different estimates for each f would not preserve convexity and is
deliberately avoided.

`certify_core2310_tail_basis.py` uses tail-count suffix sums at scale 10^60;
`audit_core2310_tail_basis.py` uses direct summation, a different count cap,
and independent arithmetic helpers at scale 10^50. The case audit takes
the larger of the two rational coefficients in each coordinate.

```
python code/certify_core2310_tail_basis.py data/schedule-core2310-conditional-final.json data/core2310-common-tail-basis.json
python code/audit_core2310_tail_basis.py data/schedule-core2310-conditional-final.json data/audit-core2310-common-tail-basis.json
python code/verify_core2310_split110.py build
python code/verify_core2310_split110.py audit
```

The next investigation also varies the projection of modulus 77, giving
70 pairs of base residues. All 70 now have independent integer core audits,
and `verify_core2310_split77_110.py audit` passes. The certified maximum is
again <0.997109604548, attained in the certificate for projections 6 modulo
7 and 5 modulo 10. The complete record is
`data/audit-core2310-split77-110-branch.json`. This removes the condition
on modulus 77, leaving only the four projections at 22,33,55,66. Other
core-210 families and other small base projections remain outside this
branch analysis. The numerical checkpoint is retained separately in
`data/core2310-split77-110-summary.json`.

```
tmp/research-venv/Scripts/python.exe code/split_core2310_mod77_110.py
python code/verify_core2310_split77_110.py build
python code/verify_core2310_split77_110.py audit
```

### A compact transport certificate also removes the modulus-77 projection

`data/audit-core2310-density30-transport.json` now verifies the same branch
with only FOUR projection conditions, those for 22,33,55,66 in the table
above. Moduli 11,77,110,154,165,231,330,385,462,770,1155,2310 are completely
free, and every component modulo 11 is free. Omissions are allowed. The
core-210 family itself is still the specified one.

The transport proof uses ten reference measures, one for each projection
modulo 10 of the 110-class, with point density depending only on x modulo
30. Its certified maximum is

    6655678650640544783918707135582152405762848091403651209293
    / 6666666400000000000000000000000000000000000000000000000000
      < 0.998351837531 < 1.

The reference 77-class projects to 6 modulo 7. The core-210 family leaves
this 7-plane complete over its core-30 residual R. Write w(x) for the point
density at x in R and q(y)=s(y)w(y mod 30) for the base-210 mass. To handle
any other projection a modulo 7, within each 30-fibre swap base masses
between the 7-planes a and 6 if both are present. If a is absent, leave the
base masses unchanged. The transported masses fit in the actual surviving
11-fibres with per-point density at most w(x); collisions among 11-colours
can only increase their available size.

In the first conditional reduction use the upper cap C=w(x), whether or
not it is attained. Base masses are at least this cap. For each 30-fibre,
the first reduction's arrays are permutations of the reference arrays,
since w(x) is constant there. The sums and maxima used by the subsequent
7-coordinate reduction are consequently identical. Thus the same core-30
envelope applies without enumerating seven separate 77 projections.

For a free class 11m with 7 dividing m, its mass is bounded by the maximum
sum of w(x) over a class modulo m/7 in R. If 7 does not divide m, use the
same sum with the number of surviving 7-lifts as multiplicity. These
bounds are attained in the reference by an unused 11-colour, using the
complete 7-plane when needed. The audit checks these equalities, every
transported base mass, all seven permutations for each of the ten cases,
and all rational continuation totals. It does not infer the theorem from
numerical optimization.

```
tmp/research-venv/Scripts/python.exe code/core2310_density30_transport.py search
python code/core2310_density30_transport.py build
python code/core2310_density30_transport.py audit
```

This compact proof and the separate 70-pair calculation have the same
scope. Neither removes the four conditions at 22,33,55,66 or deals with
the other core-210 families. In particular, the global upper bound stays 7.

## Complete removal of the four stage-11 conditions for the difficult core

The subsequent adaptive calculation **removes all four remaining stage-11
conditions** for the family in `data/core210-hard-prefix.json`. No distinct
squarefree covering system with every modulus >6 can extend this core-210
family, regardless of any later classes. All sixteen stage-11 moduli,
their residues and omissions, and all later primes are unrestricted.
This is not an unconditional upper bound 6: other initial core-210
families remain.

The final certificate is `data/audit-core2310-complete-projection-tree.json`.
It verifies 1800 basic projection cases and 2038 leaves of a complete
adaptive tree. Its largest rational bound is

    9999963527085141756663694041855110948303124054395703104521
    / 9999994600000000000000000000000000000000000000000000000000
      < 0.999996892707 < 1.

The small margin is certified with upward-rounded infinite-tail bounds;
it is not a floating-point estimate. Every selected core certificate has
an independent overflow-checked integer enumeration audit.

### Case coverage and why the additional split works

The density-30 transport handles all seven projections of modulus 77.
The other initial projections range independently over

    (a22,a33,a55,a66,a110) in Z2 x Z3 x Z5 x Z6 x Z10.

There are exactly 1800 choices. The earlier colour-domination argument
handles every component modulo 11, so these five projections suffice.
The fixed representatives use distinct colours modulo 11.

Initially 1660 candidates passed the exact primary core bound. Shared
exhaustive cutting planes raised this to 1783. For the remaining 17 basic
cases, condition additionally on the projection of modulus 165 modulo
15, using a fresh eighth colour. All 15 projections must be included.
Thus the complete proof has 1783+17*15=2038 leaves. One candidate required
a final optimization before its integer certificate passed.

The final audit checks that each basic case appears either once with 165
free, or exactly fifteen times with its distinct projections fixed. It
also checks the reference colour sets, rational densities, all seven
transports, every free-class projection capacity, the independent core
audit hashes, the common tail coefficients, and strict total <1 at every
leaf. No candidate-generation solver status is used to establish coverage.
Missing moduli can be added arbitrarily before applying the argument, so
omissions do not introduce additional uncovered branches of the proof.

Intermediate complete results are retained:

* `data/audit-core2310-release66-transport.json`: 60 cases remove the
  condition at 66, leaving conditions at 22,33,55.
* `data/audit-core2310-release33-66-transport.json`: 180 cases also remove
  the condition at 33, leaving only 22 and 55, with maximum <0.999013601374.

Simply subtracting an overlap credit for the free stage-11 classes was
insufficient at one difficult measure. An independent integer witness in
`data/audit-core2310-free-union-sharp.json` exhibits nine actual classes
that are disjoint on that measure's positive support and attain every
individual class capacity. Their union mass is exactly 6834860/100000010.
This is a limitation for that specified measure, not a covering or a
general method barrier. The successful refinement instead conditions on
165 and reoptimizes the residual measure used for all subsequent primes.
Reoptimizing the parameter schedule alone did not resolve that test.

### Candidate generation versus exact verification

`density30_branch_lp.py` reduces the measure optimization to 30 density
variables and five power-knot envelopes. Because the density is constant
on each 30-fibre, the maxima inside both conditional reductions scale
linearly with the density. The reduced optimization can generate
finite-pool candidates rapidly, but their LP objectives are NOT upper
bounds. Every candidate selected for the proof is then exhaustively
evaluated over all 810000 core-30 catalogues.

`batch_core2310_power.py` accelerates this enumeration using the integer
features C, (C-2)+, (C-4)+. Its matrix operands are nonnegative integers,
and an explicit bound ensures every product and partial sum is below
2^52. Binary64 matrix arithmetic therefore represents these integer
operations exactly; no numerical-error tolerance is used in the proof.
The subtraction of the integer constant is performed in integer arithmetic.

`audit_batch_core2310_power.py` independently generates catalogues by
mixed-radix decoding, reconstructs the conditional reductions with explicit
modular indexing, and sums direct lookup-table values using unsigned
32-bit or signed 64-bit integer arithmetic. It proves an overflow bound
before choosing the narrower representation. It reproduces every power
envelope and checks the seed law using rational arithmetic. Regression
checks reproduce three earlier independently verified certificates exactly.

The 1800 original candidates and 378 additional selected candidates all
have independent core audits; the final proof uses the specified 2038-leaf
subset. The two independently computed common infinite-tail coefficient
sets remain the ones described above. The parameter schedule is unchanged.

### Reproduce the final proof from saved certificates

```
python code/certify_core2310_tail_basis.py data/schedule-core2310-conditional-final.json data/core2310-common-tail-basis.json
python code/audit_core2310_tail_basis.py data/schedule-core2310-conditional-final.json data/audit-core2310-common-tail-basis.json
python code/batch_all_core2310_projections.py audit
python code/finalize_core2310_projection_tree.py extra
python code/finalize_core2310_projection_tree.py audit
python code/audit_core210_profile_transfer.py
```

The audit commands accept saved core certificates; rerunning floating
optimization is unnecessary. Audit batches may be split using `--shard`
and `--shards`. The final tree manifest is
`data/core2310-complete-projection-tree.json`.

### Transfer to other equality-case orbits

The full density-30 proof depends on the residual set U modulo 210 only
through the counts |{y in U:y=z mod 30}|, provided one 7-plane is complete
over the support of those counts. For a reference stage-11 choice, each
30-fibre has one distinguished 7-lift with one fewer surviving colour;
all other surviving 7-lifts have equal sizes. The sums and maxima in the
second reduction therefore depend only on the count, not on which other
7-lifts are present. Free-class capacities also depend only on these
counts, using the complete plane when the modulus contains 7.

CRT-coordinate permutations of Z2, Z3 and Z5 preserve all cylinder types
and permute the complete basic case tree. An independent audit reconstructs
the residual sets of the previously audited equality catalogue, checks the
count profiles up to these permutations, checks a complete 7-plane, and
compares the relevant seven-coordinate multisets at every possible fibre
size. It finds **7 equality-case orbits** to which the full exclusion
transfers; see `data/audit-core210-profile-transfer.json`.

This does not settle the other equality profiles or prefixes leaving more
than 104 residues modulo 210. The general upper bound remains 7. No claim
of impossibility for further improvements is made.
