# Diverse integer support states

The global modulus-2 conjecture remains unresolved. This continuation concerns
the same fixed 22-point core and high-prime pattern 11111 through prime 23.
It does not exclude all nine-prime systems or address arbitrary supports.

## Why the extra states are valid

For fixed plane counts c_i, put f_i=D_i-c_i, A=product_i f_i, and
v_S=product_{i outside S} f_i. A necessary local covering condition is

    sum_S v_S h_S >= A,       0 <= h_S <= 16.

All v_S are nonnegative integers. The convexified tail minimization returns
rational counts satisfying this inequality. Replacing every count h_S by
ceil(h_S) preserves the inequality and the upper bound 16. This constructs
an actual integer count state in the necessary local domain. Its reward is
evaluated at the rounded counts; the interpolated convexified objective is
not its reward. Such a state need not arise from a covering.

`refine_integer_support.diverse_support_states` selects the 32 cheapest
feasible plane tuples under the coarse lower envelope, plus cheap tuples
on each fixed-coordinate slice. It constructs states by the rounding above.
The exact local-minimum calculation and its independent auditor are unchanged.
Extra states only constrain the numerical master used to choose event weights.

The optimizer also retains group-constraint row indices. It reuses them only
when SHA-256 hashes of every normalized search matrix agree. This includes
the actual assignment pool, fractional scenarios, reward directions, and
normalization as represented by those matrices. A different matrix discards
the cached indices. Point states are separately checked against the integer
domain and the nonpositive local cost before insertion.

## Validation and reproduction

All 159 tests pass (`data/tests-diverse-support.log`, 45.281 seconds), including
integer feasibility and reward evaluation after rounding, unchanged exact
minima, slices with no feasible states, a one-plane problem with no tails,
and audited resumed optimization with matching and changed search matrices.

The first continuation uses the original baseline, the three saved curve
directions, the repaired fractional scenarios, and the known point pool:

```text
tmp/research-venv/Scripts/python.exe code/optimize_threshold_support.py
  --source data/nine-prime-group-support-third-resume-input.json
  --seed data/nine-prime-three-curve-pair-repaired-support.json
  --scenarios data/nine-prime-pair-repaired-scenarios.json
  --point-pool data/nine-prime-pointwise-catalogue-pricing.json
  --basis data/nine-prime-free-threshold-support.json
  --basis data/nine-prime-free-threshold-known-support.json
  --output data/nine-prime-diverse-support.json
  --time 180 --iterations 5 --diverse 1
```

Seed enrichment adds 1,811 valid local states to the previous 3,381, for
5,192 states before further separation. The preceding seed has no cached
group indices, so this first run reuses zero group constraints.

The time parameter limits numerical optimization; initialization and the
complete exact pricing pass can take additional time. A positive sampled
margin would still need complete independently audited group bounds for the
same weights. Even that would only resolve the specified core and pattern.

## First completed continuation

Two exact pricing passes retain 8,992 local states. The second candidate is
selected; independent verification in `data/audit-nine-prime-diverse-support.json`
gives

* finite actual-catalogue excess: -232040632834 / 1000000;
* ordinary unrestricted affine excess: -1763916023814 / 1000000;
* 9,453,375 plane tuples, 20,501 exceptional rational-slope cases, and
  83 exact integer subproblems checked.

The finite-catalogue excess improves from -546600.730794 to -232040.632834.
It remains negative. The numerical sampled master had positive objective
37243.0416316804 and was not closed at timeout; that value is not a proof.
The full integer local oracle still finds much cheaper states than those
represented by the master's current local constraints.

The first pricing pass gave -334515.060462 on the actual catalogue and
-344665.3436429634 against the numerical scenario sample. The second gives
-256881.61008213917 against that sample. Fractional scenario margins are
search diagnostics and are not independently certified group maxima.

The resume command replaces `--seed` with `data/nine-prime-diverse-support.json`
and `--output` with `data/nine-prime-diverse-support-resumed.json`.
It reuses exactly 1,675 group constraints; re-enrichment adds no repeated states.

## Verified resumed candidate and the missing-state gap

The resumed run completes one exact pricing pass and retains 10,946 states.
`data/audit-nine-prime-diverse-support-resumed.json` independently verifies
the catalogue excess -144572.600201 and the ordinary unrestricted affine
excess -1549921.013518. It checks 9,453,375 plane tuples, 42,090 exceptional
rational-slope cases, and 242 exact integer subproblems. Both excesses remain
negative, so there is no obstruction.

The sampled master was closed and had numerical objective 23847.294791213;
the integer-pricing candidate's sampled-scenario margin is -165601.053604845.
For these same rounded weights, `code/diagnose_support_gap.py` compares the
preceding saved local states plus the known point catalogue against the
audited exact local minima. Their total difference is exactly
189448348123/1000000. Thus almost all the observed discrepancy is explained
by omitted local count states, rather than by failure to close the sampled
group constraints. This remains a diagnostic for these weights, not a claim
about the optimal value of the complete relaxation.

The five largest local differences are at core points 164, 74, 59, 134, and
109, respectively 18530.474317, 18323.416085, 15154.123408, 12528.509844, and
11570.892854. Every one of the 22 points has a positive difference. The
gap is distributed across the core; it is not caused by a single point.
See `data/nine-prime-diverse-support-gap.json`.

## Sharing states within a parity

The local count domain and F depend on x only through x mod 2. Therefore a
state found for one core point is admissible for every other point of the
same parity. The count penalties differ by point, but feasibility does not.
This sharing must not cross parities: the D_i values and plane limits can
change. A shared state is still only a necessary-domain state, not a common
assignment of residues at two points.

Sharing the preceding known and saved pools gives 14,906 distinct even-point
states and 33,418 distinct odd-point states. Without changing any weight,
this reduces the total missing-state gap from 189448.348123 to 113931.874915.
The diagnostic is `data/nine-prime-parity-sharing-probe.json`, reproduced by
adding `--share-parity` to the gap diagnostic. This is not a newly optimized
or positive certificate.

The optimizer's `--share-parity` option shares the filtered known pool and
saved seed, then refreshes it from newly found states after each pricing
pass. It separates violated rows from that pool instead of inserting every
shared state as a master constraint. Group matrices are unchanged, so cached
group constraints remain usable. All 160 tests pass in
`data/tests-parity-support-full.log` (45.695 seconds), including domain
separation between parities and independent audit of resumed optimization
with sharing enabled.

The subsequent run uses seed `data/nine-prime-diverse-support-resumed.json`,
output `data/nine-prime-parity-shared-support.json`, and `--share-parity` in
addition to the earlier command. It reuses 1,896 group constraints.

## Verified parity-shared continuation

The last run retains 13,185 local states and 2,235 group constraints. It
separates against 630,015 point/state pairs after sharing (this is not a
count of distinct count vectors). The independently verified certificate
and audit are `data/nine-prime-parity-shared-support.json` and
`data/audit-nine-prime-parity-shared-support.json`:

* actual finite-catalogue excess: -80607833091 / 1000000;
* ordinary unrestricted affine excess: -1061225868783 / 1000000;
* 9,453,375 plane tuples, 63,589 exceptional rational-slope cases, and
  226 exact integer subproblems checked.

The sampled master was not closed when its 180-second optimization budget
ended. Its last numerical objective 6882.713693375 is not a certificate.
Completing exact pricing took the run to 304.312 seconds, excluding
initialization. The saved exact candidate's numerical scenario margin is
-94596.712542982. A further diagnostic using the preceding parity-shared
local pool gives a total missing-state gap 99024.255486; its largest terms
are at core points 74, 164, 104, 157, and 68.

Raw negative gaps can shrink merely by scaling all weights toward zero.
To avoid interpreting that as progress, `data/nine-prime-diverse-margin-comparison.json`
also divides each catalogue gap by its total event weight. The successive
values are -0.037613916894, -0.019919796181, -0.013292836473, and
-0.010635672729. This comparison is invariant under uniform rescaling.
The directions improve by this diagnostic too, but a negative value is not
an obstruction and its magnitude is not a distance to proving the conjecture.

## Continuation and scope

Resume with `data/nine-prime-parity-shared-support.json` as seed and keep
`--share-parity --diverse 1`, the original baseline, both basis files, the
repaired scenarios, and the known point catalogue. The matching group matrices
allow all 2,235 saved group constraints to be reused. Retain the full support
seed even when the best candidate comes from an earlier iteration.

The local gap remains spread across many core points. Further local-state
separation is needed before spending more time tightening group upper bounds
for these weights: an actual finite-pool negative gap already prevents these
fixed weights from proving an unrestricted obstruction. Large numbers of
exceptional planes at core points 97 and 179 suggest reusing additional
nonnegative slope values from previous rational certificates to strengthen
the coarse envelope; any such values must be recorded and independently
audited. This is a prospective optimization, not implemented here.

If the three-direction master converges to zero, investigate its precise
restricted barrier or expand the directions rather than declaring the entire
method optimal. A positive finite-pool margin would still require complete
group bounds for its exact weights. Resolving this one core/pattern would
still leave the other nine-prime cases and arbitrary prime supports. Neither
a counterexample nor a proof of the global modulus-2 conjecture has been found.
