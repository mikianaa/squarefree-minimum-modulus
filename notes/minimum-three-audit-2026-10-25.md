# Minimum modulus 3: audit and exact finite bounds

The global question remains open: does every finite covering system with distinct
squarefree moduli greater than 1 contain modulus 2? No counterexample or global
proof is obtained here. Later in this investigation, the 30030 case was excluded
by the solver-free argument in `six-prime-obstruction-2026-10-25.md`. Statements
below describing its earlier unresolved status are historical.

## What a finite obstruction implies

BBMST proves that a squarefree distinct covering has an even modulus, so its lcm
contains the prime 2. The prime replacement theorem of Simpson and Zeilberger,
Theorem 1, permits replacing a prime by a smaller absent prime. Klein's Lemma 6
states the threshold-preserving version: the new prime must be at least the
minimum-modulus threshold. For threshold 3, keep 2 and compress all odd primes
to the first available odd primes. Squarefreeness and exclusion of modulus 2
are preserved. Adding unused moduli does not destroy coverage.

Therefore nonexistence with moduli dividing 2310 rules out every counterexample
whose lcm has at most five prime factors, whatever their sizes. Nonexistence at
30030 would rule out every counterexample with at most six prime factors. The
latter is established in the subsequent six-prime obstruction note.

Sources:
- BBMST, https://arxiv.org/abs/1901.11465
- Simpson--Zeilberger, Theorem 1,
  https://sites.math.rutgers.edu/~zeilberg/mamarimY/Zeilberger_y1991_p59.pdf
- Klein, Lemma 6, https://math.colgate.edu/~integers/aa38/aa38.pdf
- Selfridge's problem is restated as Problem 4 of Dalton--Trifonov:
  https://cs.uwaterloo.ca/journals/JIS/VOL25/Trifonov/trif3.pdf

## An exact, solver-free certificate for u(210) = 22

Let u(L) be the minimum number of uncovered residues when choosing at most one
class for each divisor d >= 3 of the squarefree integer L.

Normalize the classes modulo 3, 5, and 7 to label 0 in their respective CRT
coordinates. Move a class modulo 2p off label 0 if necessary (it was redundant),
and normalize its odd coordinate to label 1. Write b_p for its parity. Flipping
parity lets us take b_3 = 0. These operations cannot increase the uncovered set's
cardinality, so lower bounds for normalized systems apply to all systems.

After these six anchor classes the two parity boxes have side lengths
p - 1 - 1[b_p = parity]. The remaining pair moduli 15, 21, 35 are axis lines,
each with the same fixed labels in both parity boxes. Their doubles 30, 42, 70
give one additional line each in a chosen parity. The final moduli 105 and 210
can cover at most two and one additional points, respectively.

Each odd coordinate occurs in exactly two odd pair classes. Permuting labels
2,...,p-1 fixes the anchor normalization. Up to these permutations, the two labels
are one of (1,1), (1,2), (2,1), (2,2), (2,3), with the last unavailable at p=3.
Thus only 4*5*5 = 100 odd-line configurations need be inspected per parity
pattern. For each, enumerate every possible placement of the three even lines,
using integer bit masks to count the union exactly.

| Anchor parities b_3,b_5,b_7 | Anchor residual | Maximum pair union | Further capacity | Uncovered at least |
| --- | ---: | ---: | ---: | ---: |
| 000 | 63 | 30 | 3 | 30 |
| 001 | 58 | 30 | 3 | 25 |
| 010 | 56 | 30 | 3 | 23 |
| 011 | 54 | 29 | 3 | 22 |

`code/three_axis_bound.py` implements this enumeration using only the Python
standard library. `data/three-axis-210-certificate.json` records its output.
The test suite also enumerates all 2304 odd-line configurations without the
symmetry quotient, obtaining the same table. A separate full-period verification
of `data/min-uncovered-lcm210-witness22.json` attains 22, proving equality. That
file's solver status is FEASIBLE: the lower bound comes from this enumeration,
not from its solver's optimization bound.

## Why u(2310) >= 26

Split modulo 11. The 11-free classes leave r >= 22 points modulo 210 in every
slice. The class of modulus 11 covers one slice. The remaining shapes 11d,
d|210, d>1, together cover at most

    f(r) = sum_{d|210, d>1} min(r, 210/d)

points of the residual across the other ten slices. Hence at least 10r-f(r)
points remain. At r=22, f(r)=194 and 10r-f(r)=26. For r>=22, at most five
summands of f can increase by one per unit increase in r, so 10r-f(r) is
increasing. This proves the bound 26, without a solver. Combined with prime
replacement, it excludes counterexamples on at most five prime factors.

## Search correctness and new model

The auxiliary loop in `few_progression_residual.py` previously excluded a residual
whenever its slice-assignment solver failed to produce a solution, including a
timeout. This was unsound: UNKNOWN is not an infeasibility proof. It now stops
with an undecided result and saves the candidate on any non-INFEASIBLE failure.
Only a proven INFEASIBLE residual may be excluded. A mocked-timeout regression
test verifies that it is not skipped. This issue does not affect the independent
direct CP-SAT searches or the enumeration certificate above.

`parity_cover.py` uses the anchor boxes directly and imposes first-use ordering
on labels 2,...,p-1 in a fixed resource order. This quotients arbitrary CRT
coordinate permutations fixing 0 and 1, rather than just affine symmetries.
Coverage is encoded exactly; all produced systems are checked independently on
the full integer period. A fixed parity pattern is one case, not a global proof.
All eight cases at 2310 were refuted in under one second each in the initial run.
Small optimization problems with and without the symmetry reduction agree.

The optional pair relaxation retains only pair moduli and treats the omitted
higher-order classes as a budget equal to their individual maximum capacities.
Infeasibility is sufficient to exclude a case. Feasibility of this relaxation
does not construct a covering system.

## Previously unfinished searches

The old six-hour small-slice run finished UNKNOWN; its result is now preserved
in `data/cpsat-lcm30030-min3-small2-6h-unknown.json`.
The old 30030 case split recorded 300 distinct cases: 4 UNSAT, 12 SKIP by
counting, 284 UNKNOWN. Another 196 cases have no result. The partial log is
preserved in `data/case-split-30030-slice1-partial.jsonl`. This is not an exhaustive
nonexistence proof.
