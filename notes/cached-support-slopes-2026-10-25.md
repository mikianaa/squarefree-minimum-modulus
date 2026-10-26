# Cached rational slopes in exact support pricing

The general squarefree modulus-2 conjecture is still unresolved. This work
continues the fixed 22-point core and high-prime pattern 11111 through 23.

## Validity of the additional coefficients

For fixed plane counts c, write A=product_i(D_i-c_i) and
v_S=product_{i outside S}(D_i-c_i). If Q is the integer event penalty and
the local covering necessity is F=A-sum_S v_S h_S <= 0, then every s>=0 gives

    Q(c,h) >= Q(c,h)+s F(c,h)
            >= s A + Q_plane(c) + sum_S min_k (Q_S(k)-s v_S k).

Thus any nonnegative slope gives a lower bound. A slope from a previous
candidate is only a suggestion: the expression is recalculated using the
new penalties, and the new coefficient is saved for independent verification.
No old minimum or old bound is asserted for new weights.

`code/cached_support_slopes.py` selects up to 12 additional integer
coefficients from prior rational slopes and prior coarse coefficients. It
rescales by certificate denominators, tries the neighboring integers, and
retains frequent values and values spread across the range. A conservative
bound matching the vectorized evaluator filters coefficients that could
overflow signed 64-bit arithmetic. Tests check that the maximal admitted
coefficient passes the evaluator's guard and the next integer fails it.

The optimizer records the actual coefficient list separately at each core
point. Its independent auditor recomputes every corresponding lower bound,
every exceptional case, and each attaining integer state as before. The
`--cached-slopes 12` option enables the new selection; the default zero keeps
the prior selection. `--dual-simplex` separately permits the numerical
master to use dual simplex when adding constraints; it changes no proof rule.

## Same-weight benchmark

Using the already audited `data/nine-prime-parity-shared-support.json`,
`code/benchmark_cached_slopes.py --points 8 17` gives:

| Core point | Old exceptional planes | New exceptional planes | Old / new integer DPs | New time |
| --- | ---: | ---: | ---: | ---: |
| 97 | 35,477 | 1 | 96 / 1 | 0.781 s |
| 179 | 27,597 | 1 | 77 / 1 | 0.797 s |

The exact local minima remain 726848680 and 3753181952, respectively, in
denominator 1000000. This is an equality check against the previously
audited minima, not a new obstruction. The benchmark and coefficient lists
are saved in `data/nine-prime-cached-slopes-benchmark.json`.

These two points benefited from small coefficients absent from the old
coarse grid: its smallest positive value was 125000, whereas useful saved
coefficients at point 97 range from 3523 to 18367 and at point 179 from
21912 to 89247. The coarse grid was valid but inefficient in these cases.

## Search continuation

The first new run uses `data/nine-prime-parity-shared-support.json` as seed,
the original baseline `data/nine-prime-group-support-third-resume-input.json`,
both existing nonbaseline curve directions, the repaired scenario pool,
the known point catalogue, and options

    --time 240 --iterations 5 --diverse 1 --share-parity --cached-slopes 12

It writes `data/nine-prime-cached-slope-support.json` and reuses 2,235 group
constraints. These runs still select only combinations of the specified
three directions at each group/point. They are not an optimization over
every possible inequality or a test of all squarefree covering systems.

## Verified first candidate

The first new candidate and its independent audit are
`data/nine-prime-cached-slope-support.json` and
`data/audit-nine-prime-cached-slope-support.json`. The audit verifies:

* finite actual-catalogue excess -47976248035 / 1000000;
* ordinary unrestricted affine excess -714839589096 / 1000000;
* all 9,453,375 plane tuples, 1,763 exceptional rational cases, and
  114 exact integer subproblems.

The sampled numerical master was not closed; its last objective was
1351.585454223, whereas the exact local candidate's numerical scenario
margin is -57729.020898063. The saved seed has 15,037 local states and
2,441 group constraints. The negative exact excess does not exclude the
actual catalogue, and complete group bounds alone cannot make these fixed
weights an obstruction.

All 163 tests pass in `data/tests-cached-slopes-final.log` (113.370 seconds
while a research process was also running). An integration test deliberately
adds an arbitrary nonnegative rational slope suggestion to the old cache,
then checks that the new coefficient is saved and that independent audit
of the new certificate succeeds. Correctness depends on recomputing and
verifying the new bounds, not on trusting old cache entries as proofs.

The continuation uses this first candidate as seed and adds `--dual-simplex`,
writing `data/nine-prime-dual-simplex-support.json`. It retains the original
baseline and the same three directions, scenarios, and known point pool.
Timings from runs with different seeds are not a controlled comparison of
the two numerical simplex algorithms.

That unscaled-cache dual-simplex comparison was deliberately interrupted
after five point results, once the scale-adjusted implementation had produced
a complete candidate. Its session ended with exit code 1 after Ctrl-C;
`data/nine-prime-dual-simplex-support.json` was not written. Its partial log
is retained, but it is neither an infeasibility result nor a certificate.

## Adapting to changed weight scales and keeping the master small

The first dual-simplex continuation chose much smaller weights. Reusing
denominator-scaled slopes alone then became inefficient: its first three
points required 78,863, 49,163, and 53,276 exceptional-plane checks. This
does not invalidate those lower bounds, but it motivates adapting their
scale. The updated selector uses the ratio of the new and old total integer
penalty masses at the core point. Under uniform rescaling this is exactly
the factor by which an optimal slope rescales. For nonuniform changes it is
only a heuristic; every selected coefficient is still nonnegative and every
new bound is recomputed and audited. A zero old mass falls back to the
denominator ratio.

The `--compact-master` option retains every saved valid support state, but
initially activates only the saved attaining states, maximal fallback states,
and the selected diverse states. All retained states remain in the known
pool, including through parity sharing. Violated constraints are inserted
as needed. Thus closing separation checks all known states even though they
were not all installed at initialization; closing separation still says
nothing about local states outside the pool until exact pricing runs.

The complete retained pool is written back to `support_seed`. An integration
test confirms both independent audit of a compact-run candidate and retention
of every inherited seed state while fewer state constraints are active.
In compact runs the historical `new_states` diagnostic counts newly
activated constraints, not necessarily newly discovered count vectors.

All 164 tests pass in `data/tests-compact-master-full.log` (75.950 seconds
while research computations were also running). The compact continuation
uses the same seed as the dual-simplex continuation, adds `--compact-master`,
uses the revised mass-scaled cache, and writes
`data/nine-prime-compact-master-support.json`. It initially retains all 15,037
states while activating only 1,818 local constraints and reusing 2,441 group
constraints. Its first solve took 1.234 seconds, but had a weaker restricted
master; this is not a comparison of equally closed optimization problems.

## Verified compact candidate and retained continuation

`data/audit-nine-prime-compact-master-support.json` verifies the completed
compact candidate: catalogue excess -8554427561/1000000, ordinary
unrestricted affine excess -98140321528/1000000, all 9,453,375 plane tuples,
11,773 exceptional rational cases, and 126 exact integer subproblems.
The sampled master was not closed. The run retains 17,211 local states,
of which 4,993 were active constraints, and 3,476 group constraints.

The smaller raw deficit is misleading if compared without normalization.
The first cached-slope candidate has total event weight 5014870401318 and
normalized catalogue excess -0.009566797184. The compact candidate has total
event weight 677279443230 and normalized excess -0.012630573165, which is
worse. Neither excludes the actual catalogue. The earlier parity-shared
candidate's normalized value was -0.010635672729, so the first cached-slope
candidate remains the best of these three by this diagnostic.

`data/nine-prime-compact-master-resume.json` keeps the first candidate's
audited proof payload and weights, and imports the compact run's full
17,211-state pool and 3,476 group-constraint rows. The inputs, actual group
catalogue, and group-matrix hashes were checked equal; every original seed
state is retained; all merged states were checked in their integer local
domains. The input, group catalogue, pricing history, positive-state bounds,
and exact integer refinements are unchanged from the audited weight source.
This is a continuation checkpoint, not a stronger proof.

Resume from that checkpoint with the original baseline and both basis files,
the repaired scenario pool, the known point catalogue, and

    --diverse 1 --share-parity --cached-slopes 12 --dual-simplex --compact-master

Next work should close separation against the saved local pool and then
apply exact local pricing. If the specified three-direction family has no
positive separation, establish a precisely scoped barrier or expand the
directions; zero numerical objectives or shrinking raw negative deficits do
not prove that the method is optimal. All general-modulus and arbitrary-prime
requirements of the original goal remain in force. No full covering or
global obstruction has been obtained.
