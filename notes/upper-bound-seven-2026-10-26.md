# The same stop-loss comparison certifies minimum modulus at most 7

Date: 26 October 2026.

Later update on the same date: separate certificates now prove that 7 is
optimal for this fixed catalogue comparison, first for 0<=delta_p<=1/2
and then for the full range 0<=delta_p<1. See
`stoploss-full-range-barrier-2026-10-26.md`. Different comparison bounds and
geometry-aware extensions remain outside that limitation theorem.

## Result and scope

Every finite covering system with distinct squarefree moduli has a modulus
at most 7. The existing distortion / convex stop-loss comparison gives a
uniform total loss strictly less than 1 under the hypothesis that all
moduli exceed 8. Since 8 is not squarefree, this proves the stated bound.

This disproves the proposed barrier at 10 for this specific comparison.
No stronger comparison lemma or compatibility inequality is required.
It does not settle the conjecture that every such system contains modulus 2,
and by itself does not prove that 7 is optimal for the method. This is an internally
checked computer-assisted result, not a claim of external peer review.

The primary verifier certifies the rational upper bound

    999087538090134249032793704260561781812927932344492427708399 / 10^60
      < 0.999087538091 < 1.

The independent implementation also certifies a bound below 1:

    999087538090134249030821201485083204038153377114994384897372 / 10^60.

The difference is less than 2e-21 and is fully accounted for by the
primary verifier's conservative truncation remainders at primes 2,3,5,7.
The independent implementation evaluates those small stages exactly.
The later-stage and tail contributions agree as integers.

## Why the previous search missed this

The saved earlier search used banded distortions and coordinate steps of
prescribed sizes. Its 1.0178771888087672 was an attained numerical upper
estimate, not a lower bound over all schedules. Moreover, it replaced the
distribution calculation by quadratic moments beyond prime 100000.

The new search makes each early prime's distortion an independent control,
uses exact one-coordinate breakpoint choices, optimizes the quadratic
continuation, and retains the count distribution through prime 3000000.

For p=5,7, put A(n)=2^n-1; for p>=11 put A(n)=2^n. Set

    r_p = 1/[p(1-delta_p)].

The stage cost in state n becomes

    g_p(n,r_p) = max(0, 1 + r_p(A(n)-p)).

The distribution update is

    w_next(n) = (1-r_p)w(n) + r_p w(n-1).

For a fixed prefix distribution and fixed future controls, all future
expected costs are affine in r_p. The present cost is piecewise affine,
with breakpoints delta_p=A(n)/p. Hence an optimal choice for this single
coordinate occurs at a breakpoint or an endpoint. Backward continuation
vectors allow a forward coordinate sweep to evaluate these choices without
recomputing the whole schedule at every trial. This is coordinate descent;
no global-minimality conclusion is asserted or needed.

The initial distribution, after primes 2 and 3 with zero distortions, is

    P(N=0)=1/3, P(N=1)=1/2, P(N=2)=1/6.

For the later quadratic continuation, if K is the future coefficient of
E[4^N], its backward update is

    F_p(K) = K + (1+sqrt(1+12pK))^2/(4p^2),
    delta_p = 1/(1+sqrt(1+12pK)).

It ends in the same CFT infinite-tail estimate already used by the original
verifier. This optimization affects only candidate generation; the final
certificate verifies each explicit rounded parameter directly.

Numerical search milestones, not impossibility bounds:

| Distribution retained through | Numerical total after optimization |
| --- | ---: |
| 100000 | 1.011793571033 |
| 1000000 | 1.000653162504 |
| 3000000 | 0.999010207148 |

The last number precedes parameter clipping and rounding. The exploratory
search allowed delta up to 1-1/p. Its 1402 parameters above 1/2 were clipped
to 1/2, and all parameters were rounded to multiples of 10^-6. The resulting
schedule satisfies the original verifiers' unchanged delta<=1/2 assumption
and still certifies a total of 0.999087538090... .

## Actual certificate

The proof uses exactly the comparison inequalities recorded in Sections
3.2-3.3 of `sqfree-agent-report-2026-10-23.md`. Only the parameter schedule
and the division between distribution and quadratic evaluation changed.

| Contribution | Upper bound, displayed approximately |
| --- | ---: |
| p<=8 | 0.208333437500 |
| 8<p<=3000000, distribution calculation | 0.789177808757 |
| Remaining primes through stage 1000000 | 0.000834639257 |
| Infinite prime tail | 0.000741652576 |
| Total | 0.999087538090 |

All proof arithmetic uses integer or rational upper enclosures. The
Poisson-binomial distribution is retained through N=120 at scale 10^60.
The omitted count tail is bounded explicitly using

    E[2^N 1{N>120}] <= (2/3)^120 E[3^N].

There are 1000000 explicitly evaluated prime stages, ending at 15485863.
Beyond that stage delta=1/2, and the existing CFT / Rosser-Schoenfeld
infinite-tail bound applies. Thus the result has no restriction on the
number of prime factors of the covering system's lcm.

The artifact `data/schedule7-breakpoints.json` has SHA-256

    74ebcb30fc6eb0b4cca72e2e4edecf9f8729a550547bee1ac1bf8c593d5520af

and 254572 consecutive constant-parameter bands after rounding.
`data/upper-bound-seven-summary.json` records the certificate references,
the exact integers, and the scope of the conclusion.

## Reproduction

From the repository root, using ordinary Python (standard library only):

```text
python code/verify_squarefree.py data/schedule7-breakpoints.json data/cert7-breakpoints.json --p-exact 3000000
python -I code/audit_verify.py data/schedule7-breakpoints.json data/audit-cert7-breakpoints.json --p-exact 3000000
```

The primary verifier was unchanged. The only change to the independent
verifier was replacing its linear band lookup by binary search and checking
strictly increasing band endpoints, to support a large per-prime schedule.
No probability, loss, rounding, or tail calculation in either verifier was
changed for this result. Boundary tests check the lookup change.

The numerical generator, requiring NumPy, is
`code/optimize_stoploss_breakpoints.py`; it starts from the old saved
schedule. The certified schedule was generated by its exploratory precursor
`tmp/probe_stoploss_control.py`. The optimizer's floating-point state
truncation is not trusted by the proof; both verifiers independently rebuild
the probability upper bounds with the larger N=120 limit and an explicit
remainder. Candidate-generation tests run with:

```text
python -m unittest discover -s code -p test_stoploss_breakpoints.py
```

## Correction to earlier status

The earlier stop at 10 reflected the search and evaluation choices, not an
established method barrier. The numerical local-sharpness experiments at
primes 5 and 7 remain valid at their stated fixed parameters, but they did
not exclude improvements from optimizing the rest of the schedule and
evaluating more of its distribution. The current upper bound in this
repository is 7; the independent lower bound remains 2.
