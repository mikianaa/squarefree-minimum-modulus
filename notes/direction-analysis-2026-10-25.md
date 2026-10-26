# Why the lower-side search is stalling, and directions that are not "nine primes again"

Written 25 October 2026 (Claude Code, Fable 5.1) at the owner's request, after the
nine-prime relaxations produced a sixth audited negative-cost witness. This note
contains no new theorem; its numbers are reproducible from the saved witnesses.

## 1. Diagnosis

**(a) Counting died at five primes and the gap widens by a factor ~2 per prime.**
The slice count excludes a cover at Lq only if the q-free residual has density at
least about 5-9 % (largest u with f(L,u) >= (q-1)u):

| k | L | q | counting needs u(L)/L > | known residual density |
|---|---|---|---|---|
| 4 | 210 | 11 | 7.1 % | 10.5 % (exact) |
| 5 | 2310 | 13 | 9.1 % | 5.2 % (exact) |
| 6 | 30030 | 17 | 7.6 % | 1.7-2.6 % |
| 7 | 510510 | 19 | 8.0 % | 0.1-1.5 % |
| 8 | 9699690 | 23 | 7.0 % | <= 0.82 % |
| 9 | 223092870 | 29 | 5.5 % | unknown |
| 10 | 6469693230 | 31 | 5.6 % | unknown |

So from k = 5 on, every exclusion has to extract geometry from the residual, and the
amount of geometry needed roughly doubles per prime. That is why the proofs went
counts -> parity -> planes -> paired projections -> group events -> pointwise
mixtures, and why at nine primes each strengthened relaxation is defeated by a new
witness that is not a covering system (relaxation gap). Column generation against
one witness at a time is the slowest possible way to find the right inequality
class, and each round needs its own independent audit.

**(b) Even a nine-prime exclusion carries no induction hypothesis.** The next step
has 1023 shapes and 31 slices; nothing proved at k primes is reused at k+1 except
the residual bound, and (a) shows the residual bound is far too weak to reuse.

**(c) The two sides cannot meet in the middle.** Any distortion / density argument
treats stages p >= 7 by a loss comparable to their density. But the six moduli
3, 6, 5, 10, 15, 30 alone can cover 27 of the 30 classes mod 30, so after the
exact stages 3 and 5 at most 10 % may remain, while stage 7 alone has density
(1/7)(3/2)(4/3)(6/5) = 0.34. Hence no stop-loss certificate can go below the
small-prime structure, and M_sf = 2 (if true) needs a structural proof. The upper
side saturating at 10 is consistent with this.

## 2. The decisive measurement: actual residuals are spread, not concentrated

For a cover with lcm Lq, the shapes qd across the q-1 non-free slices give one
progression per divisor d > 1 of L, and every point of the q-free residual R lies
in at least q-1 of them. So cap(R) := Sum_{d | L, d > 1} max_a |R cap (a mod d)|
must be >= (q-1)|R|. For the saved near-optimal witnesses:

| L | q | |R| | cap(R)/|R| | required | shortfall |
|---|---|---|---|---|---|
| 2310 | 13 | 121 | 5.17 | 12 | 6.8 |
| 30030 | 17 | 792 | 6.04 | 16 | 10.0 |
| 510510 | 19 | 7436 | 7.15 | 18 | 10.9 |
| 9699690 | 23 | 79979 | 7.87 | 22 | 14.1 |

The ratio grows like the divisor sum (about Prod(1+1/p)), i.e. doubly
logarithmically, while the requirement grows like the next prime. So the residuals
that small families actually leave are nowhere near "coverable by the next prime",
and the gap is widening. A cover can only arise from a family whose residual is
far more concentrated than the near-optimal ones; the repository's own 2310
experiments show the most concentrated residuals are single classes a mod m,
and those force a duplicated-modulus cover one level down.

## 3. Directions

1. **Uniform inductive invariant instead of another primorial.** The target
   statement is of the form: for every distinct family F over L (moduli >= 3)
   either cap(R_F) < (q-1)|R_F| for all primes q > max(L), or R_F contains a
   coset that forces a cover with a duplicated modulus over a proper divisor
   of L; and no such duplicated cover exists (second induction on duplicates,
   t <= floor((tau(L)-1)/(q-1)) extra distinct shapes, modulus 2 allowed).
   Test cheaply first: maximise cap(R_F) - (q-1)|R_F| over F at 2310 and 30030
   (the existing slice_capacity.py relaxation) and record the structure of the
   maximisers. If the maximisers are always coset-like, the two-level
   induction is the proof shape; if a spread maximiser reaches the threshold,
   the conjecture itself is in doubt.

2. **Static two-family model.** Replace the slice bookkeeping by: F (distinct,
   >= 3) covers Z_L \ R and a single family G (one residue per divisor d > 1 of
   L) covers every point of R at least q-1 times. This is a necessary condition
   with 2^k residue variables and no partition into slices; it is the right
   object for both SAT-style construction and LP-style exclusion. A point with
   multiplicity >= q-1 under G has >= q-1 pairwise coherent divisors through it,
   so R splits into coherent blocks (cosets a mod m with omega(m) >= log2 q,
   all residues a_d = a mod d) and sporadic points costing q-1 incidences each
   out of a budget of L (Prod(1+1/p) - 1) ~ 2.6 L. This block/sporadic
   dichotomy is a case split that is uniform in k.

3. **Aligned construction test (cheap, decisive either way).** Take all
   G-residues zero. Then R must lie in {x : omega(gcd(x, L)) >= s}, s = ceil(log2 q)
   (s = 5 for q = 23, 29, 31; s = 6 for 37..61). Question: can a distinct family
   with moduli >= 3 dividing 9699690 cover every x with at most four of the
   eight primes dividing it? The residual allowed has density about 6e-4, far
   below the best witness (0.82 %), so a negative answer is expected and would
   be a clean lemma; a positive answer would be the first serious step toward a
   counterexample, since the slice side (22 disjoint hitting families of
   divisors) looks borderline feasible by Turan-type counts.

4. **Take the possibility M_sf >= 3 seriously and search constructively at 9-11
   primes with heuristics, not exact solvers.** Residual-first: choose a target
   residual shape (a few coherent blocks), cover its complement by greedy /
   local search over distinct divisors, then try to partition the remaining
   divisors into q-1 covers of the residual. The quantity to track is
   cap(R)/|R| against q-1 (table in section 2); if it stays at a third of the
   requirement at 9 and 10 primes, that is the strongest available evidence for
   the conjecture and tells a prover exactly what to bound.

5. **Stop generating cuts against single witnesses.** If the nine-prime
   relaxation route is kept at all, characterise the full family of valid
   inequalities (a polyhedral description of the actual projection patterns of
   one group) once, rather than separating one witness at a time; the notes
   show the separators already report no more cuts of the forms tried.

Items 1 and 3 can each be started in an afternoon with existing code
(slice_capacity.py, cpsat_cover.py, the full-period audit scripts).
