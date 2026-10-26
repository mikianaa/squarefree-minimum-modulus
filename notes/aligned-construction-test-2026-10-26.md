# The aligned construction test (direction 3) decides k <= 6 only, by pure counting

Written 26 October 2026 (Claude Code, Sonnet 5), following up on direction 3 of
`notes/direction-analysis-2026-10-25.md`. This adds no new impossibility result;
it pins down exactly how far the "cheap, decisive" test from that note actually
goes, and why it does not reach the open k = 8 case without the heavier
machinery of direction 2 (the static two-family model).

## The test

For primes p_1 < ... < p_k, L = prod p_i, and next prime q, the slice-extension
model for a cover with lcm Lq needs a q-free family F (distinct moduli >= 3,
dividing L) over Z_L with residual R, plus q-1 further copies ("slices") each
covering R using shapes qd (d | L, d > 1). Direction 3 asks about the special
case where every slice uses the trivial residue 0 mod d for every d it holds
("G-residues zero"): then a point x of R gets covered by exactly
2^omega(gcd(x,L)) - 1 of the available shapes (every nonempty subset of its own
divisors), so sufficiency against q-1 slices needs 2^omega(gcd(x,L)) - 1 >= q-1,
i.e. omega(gcd(x,L)) >= s := ceil(log2 q). So this construction requires F to
cover T := {x : omega(gcd(x,L)) <= s-1} exactly, leaving R inside the
complementary set {x : omega(gcd(x,L)) >= s}.

`code/aligned_construction_test.py` computes |T| and R_allowed := L - |T|
exactly (integer arithmetic, no solver) for k = 4..8, and compares R_allowed
against the best already-certified rigorous lower bound on u(L) (the true
minimum uncovered count for ANY distinct family with moduli >= 3 of L, not
just this restricted construction). If R_allowed is below that lower bound,
no family at all -- let alone this restricted one -- can confine its residual
to R_allowed, so the construction is impossible by a one-line counting
argument. Output: `data/aligned-construction-test.json`.

| k | L | q | s | R_allowed | density | known u(L) >= | decided infeasible by counting |
| - | - | - | - | - | - | - | - |
| 4 | 210 | 11 | 4 | 1 | 0.48% | 22 (exact) | yes |
| 5 | 2310 | 13 | 4 | 24 | 1.04% | 121 (exact) | yes |
| 6 | 30030 | 17 | 5 | 36 | 0.12% | 501 (integer-only) | yes |
| 7 | 510510 | 19 | 5 | 1074 | 0.21% | 510 (integer-only) | **no** (1074 > 510) |
| 8 | 9699690 | 23 | 5 | 30682 | 0.32% | 2 (weakest branch) | **no** (30682 > 2) |

So the counting argument alone reproves (via a new, much simpler route) that
k = 4, 5, 6 already fail -- consistent with, and no stronger than, the
previously established six/seven-prime obstructions, which used far more
machinery. It does **not** reach k = 7 or k = 8: at both, R_allowed exceeds
every already-certified u(L) lower bound, so pure counting cannot exclude a
family whose residual happens to land entirely inside the omega >= s region.

This matches the diagnosis in section 1(a) of the direction note: the gap
between "needed density" and "achievable density" closes somewhere around
five to six primes, and after that point only the specific shape of the
target, not its size, can decide the question.

## Why the naive realization of this construction fails badly anyway

A natural way to realize "F covers T" is: for every nonempty subset S of the
k primes with |S| <= s-1 (other than S = {the single banned prime}), use the
modulus d_S = prod(S) with residue 0. This covers every x divisible by at
least the primes in S -- in particular every x with exact divisor set S'
(|S'| <= s-1) that is itself a candidate modulus, since x is then ≡ 0 modulo
its own exact divisor product.

This completely misses one specific, large piece of T: the points coprime to
L (omega = 0), since no residue-0 modulus (d > 1) ever contains a point
coprime to its own modulus. Their count is exactly prod_p (p-1):

| k | L | coprime-class size | R_allowed |
| - | - | - | - |
| 7 | 510510 | 92,160 | 1,074 |
| 8 | 9,699,690 | 1,658,880 | 30,682 |

At k = 8 the single coprime class is 54 times larger than the entire allowed
residual. So the "residue 0 everywhere" realization of the aligned
construction is not a near miss at k = 7, 8; it is off by orders of magnitude,
and the coprime class has to be covered by *other* composite moduli using
non-trivial residues.

That is in fact possible in principle: a modulus d with prime set P_d and
residue coprime to d intersects the coprime class in exactly
prod_{q not in P_d} (q-1) points, and -- crucially -- every modulus that
shares some of P_d's primes with another (e.g. 3, 15, 21, 105 all involve the
prime 3) gets an *independent* residue choice modulo the primes they share,
since they are different integers. This is exactly the same "plane capacity"
flexibility that the already-certified six/seven/eight-prime obstructions
exploit (composite moduli choosing independent residues per shared prime).
Whether that flexibility is enough to actually tile the coprime class (and
the rest of T) without exceeding the 254 (k=8) available distinct moduli is
not decided by this note; it is precisely the question that the static
two-family model (direction 2) is set up to answer, and is not a "cheap"
counting check.

## Conclusion and scope

* The aligned/G-residues-zero construction is cleanly excluded by counting at
  k <= 6 (reproof, not new).
* It is genuinely undecided by counting at k = 7, 8; the k = 8 case is the
  first open primorial (223092870, nine primes) and remains the only
  interesting target.
* The naive residue-0 realization fails by two orders of magnitude at k = 8,
  but more flexible composite-residue realizations are not excluded by this
  note. Resolving k = 8 needs the static two-family model (item 2 of
  `notes/direction-analysis-2026-10-25.md`), reusing the plane/core capacity
  machinery (`code/core_capacity_frontier.py`, `code/plane_capacity_bound.py`)
  already built for the eight-prime obstruction, adapted to the omega-
  threshold target T instead of a fixed normalized core. That adaptation is
  not done here and is the natural next step.
* This does not change the global upper or lower bound. `primes_needed`
  remains 9; the first remaining primorial remains 223092870.

## Reproduction

```
python code/aligned_construction_test.py
```

Pure integer arithmetic, no solver, no floating point; runs in well under a
second. `data/aligned-construction-test.json` records the full table and the
coprime-class sizes.
