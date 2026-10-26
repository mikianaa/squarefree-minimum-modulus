# A counterexample omitting modulus 2 needs at least nine primes

## Result

Every finite covering system with distinct squarefree moduli greater than 1,
whose least common multiple has at most eight distinct prime factors,
contains modulus 2.

The finite computational statement is stronger quantitatively: every family
of distinct divisors greater than 2 of

    N = 9699690 = 2*3*5*7*11*13*17*19

leaves at least **2** residues modulo N uncovered. All four branches of its
normalized core have independently audited positive integer certificates.
The audits use exact arithmetic and do not require an optimization solver.
They are not proof-assistant verifications; no novelty claim is made.

This does **not** settle the general minimum-2 conjecture. Counterexamples
with nine or more prime factors remain open. The first remaining primorial
is **223092870**, through prime 23. The global minimum-modulus upper bound
certified elsewhere in this repository remains 10.

## The complete case split

Add missing admissible moduli if necessary, normalize each prime class to
coordinate label 0, and each twice-prime class to label 1 in its chosen
parity. A redundant twice-prime class can be moved to label 1, only enlarging
coverage. Reverse the parity coordinate so that the class modulo 6 is even.
Let b_p denote the parity of the class modulo 2p. The three low-prime bits
are therefore (b3,b5,b7)=(0,b5,b7).

| Core pattern | Parities of classes modulo 6, 10, 14 | Audited uncovered lower bound |
| --- | --- | ---: |
| 000 | even, even, even | 6696 |
| 001 | even, even, odd | 305 |
| 010 | even, odd, even | 15 |
| 011 | even, odd, odd | 2 |

Each row includes every placement of the other core classes modulo divisors
of 210 and all sixteen choices of (b11,b13,b17,b19). Thus these are four
exhaustive branches, not four particular covering families.

The 000 and 001 proofs use a core count frontier and arbitrary-subset
compression capacities. The 010 and 011 proofs also retain pairwise
projection geometry; see `eight-prime-anchor010-2026-10-25.md` for the
inequality and its independent weighted-point audit. Every final dual bound
is checked by enumerating all local plane-count tuples and minimizing each
tail count independently. Solver optimality and floating-point signs play
no role in accepting a certificate.

## Completing branch 011

The independent ordinary-residue enumeration gives 317869 pair residuals,
4323 joint histogram profiles, and 8027 final-count capacity tables. Its
533 recorded dominating tables are checked pointwise at equal counts.

The initial joint-capacity screen requires refinement of 113 profiles.
These reduce to 332 actual pair-residual orbit representatives. The auditor
checks every paired-capacity table independently and explicitly verifies
that the 1614 ordinary residuals in the refined profiles belong to the
listed CRT-label orbits. There is no dependency on the earlier inclusion
reduction.

After reusing the 010 dual pool and adding six further duals, the remaining
unresolved inequalities all belong to one pair residual, with parity sizes
(11,15), histogram profile 612. Its
optimistic final count pairs are (9,14) and (10,13). Twenty further rounded
dual certificates complete the screen. The complete independent audit
checks 19422 final inequalities, using 20 selected local dual/pattern
checks and 1181475 plane-count tuples. The minimum selected integer bound
is 2. The audit takes about 252 seconds on this machine.

The 010 audit separately checks 351244 ordinary pair residuals, 4340
profiles, 82 paired-capacity certificates, 10505 final inequalities, and
1653300 local plane-count tuples. Its minimum bound is 15.

## From the primorial to arbitrary prime supports

[BBMST](https://arxiv.org/abs/1901.11465) proves that a distinct squarefree
covering system has an even modulus. Consequently 2 belongs to its prime
support. This theorem does not itself require the modulus 2 to occur.

Suppose a cover omitting 2 has at most eight prime factors in its least
common multiple. If an odd prime p in its support is larger than an absent
odd prime q, restrict the p-coordinate of the CRT product to any q labels
and relabel them by the residues modulo q. Classes fixing an excluded label
disappear. In every remaining modulus containing p, replace p by q.

The restricted product remains covered. The resulting moduli stay
squarefree and distinct because q was absent. No modulus 2 is created,
since the replacement prime is odd. Repetition embeds the support into
{2,3,5,7,11,13,17,19}; add unused coordinates if needed. This contradicts
the finite obstruction above. Hence a counterexample needs at least nine
distinct prime factors.

## Reproduction and artifacts

From the repository root, using a Python environment with NumPy:

```
python code/verify_eight_prime_obstruction.py
```

This reruns all four independent audits before writing the combined
summary. The optional `--record-existing` flag merely summarizes previously
saved audit reports; it does not reverify them. The current combined summary
was recorded after the four branch audits had passed in this research run.

* `data/eight-prime-obstruction.json`: combined statement, branch bounds,
  source-file hashes, and explicit unresolved global scope.
* `data/audit-eight-prime-anchor000.json`
* `data/audit-compressed-tail-branch-001.json`
* `data/audit-paired-tail-branch-010.json`
* `data/audit-paired-tail-branch-011.json`

The 58-test regression suite passes. An actual partial family still gives
u(9699690)<=79979, so the current independently checked bracket is
2<=u(9699690)<=79979. Neither endpoint is asserted optimal.

The attribution of the general question is to John Selfridge, going back
to at least 1981, as recorded by
[Cummings, Filaseta, and Trifonov](https://arxiv.org/html/2211.08548v1).
