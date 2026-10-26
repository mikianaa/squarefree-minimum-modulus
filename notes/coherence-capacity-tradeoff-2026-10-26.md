# The spread-vs-coherent tradeoff in cap(R)/|R|, quantified at nine primes

Written 26 October 2026 (Claude Code, Sonnet 5). A follow-up to
`notes/nine-prime-greedy-construction-2026-10-26.md`, prompted by item 2 of
`notes/direction-analysis-2026-10-25.md` ("R splits into coherent blocks...
and sporadic points") and the slice-extension history already in
`PROGRESS.md` (24 Oct entry: the most concentrated residual at a smaller
scale was a full class mod 105, and realizing it as an actual residual
needed the modulus 105 *twice* -- a duplicate, which a further count then
excluded one level down). This note does not resolve anything; it puts
numbers on the tension that the rest of the program is already built
around, at the scale of the open nine-prime case.

## The two residuals compared

The greedy/local-search witness in `nine-prime-greedy-construction-2026-10-26.md`
is a *spread* residual (every one of the 510 moduli contributes to covering
"the rest", so what's left behind is scattered): |R| = 1,213,060,
cap(R)/|R| = 7.88, badly short of the 28 a tenth-prime extension would need
(28% of the requirement), consistent with the already-established slow
growth of this ratio across k = 5..9 (5.17 -> 6.04 -> 7.15 -> 7.87 -> 7.88).

At the opposite extreme, a *coherent* residual R = {x : m | x} for some
divisor m of L = 223092870 (a single full residue class, size L/m) has an
exact, cheap closed form for its capacity, proved here and checked against
brute-force enumeration at a small scale (L=2310):

    cap(R) = sum_{d | L, d > 1} L / lcm(d, m)

(derivation: for divisor d, the points of R = {m | x} that also lie in a
fixed residue class a mod d exist, by CRT, iff a is divisible by
gcd(d, m); whenever that holds the count of solutions is exactly
L / lcm(d, m), independent of which valid a -- so every modulus's best
choice contributes the same amount regardless of which compatible residue
is picked, including a = 0).

## Numbers

| omega(m) | example m | \|R\| = L/m | cap(R) | cap(R)/\|R\| |
| --- | --- | --- | --- | --- |
| 1 | 2 | 111,546,435 | 445,837,245 | 4.00 |
| 2 | 6 | 37,182,145 | 241,509,695 | 6.50 |
| 3 | 30 | 7,436,429 | 85,460,851 | 11.49 |
| **4** | **210** | **1,062,347** | **22,161,973** | **20.86** |
| 4 | 390 | 572,033 | 12,699,007 | **22.20** |
| 5 | 2310 | 96,577 | 3,774,143 | 39.08 |
| 6 | 30030 | 7,429 | 545,531 | 73.43 |
| 7 | 510510 | 437 | 61,003 | 139.60 |
| 8 | 9699690 | 23 | 6,121 | 266.13 |
| 9 | 223092870 (m=L) | 1 | 511 | 511.00 |

At omega(m) = 4 (m = 210, |R| = 1,062,347 -- almost the same size as the
spread witness's 1,213,060) the ratio is already 20.86, and by m = 390 it
has already crossed **22**, the exact threshold a tenth-prime extension
needs. A fully coherent residual of comparable *size* to the spread witness
clears the requirement comfortably; the spread witness, at the same size,
is at barely a third of it.

## Why this is not a route to a counterexample

A residual R = {m | x} is only useful if some actual distinct-moduli family
F realizes it, i.e. F covers every x with m not dividing x and nothing
else. This is exactly the situation already studied at a smaller scale in
`PROGRESS.md`'s 24 Oct slice-extension entry: the single most concentrated
residual found at L = 2310 for the (then relevant) extension was a full
class mod 105, and *achieving* that residual needed modulus 105 twice --
once as part of the base family F over Z_2310, and again (doubled, i.e.
qd = 13*105 appearing effectively as a repeated resource across slices) in
the continuation, which is not allowed for a distinct-moduli cover. That
duplicate route was itself excluded one level down by a further count
(f(210, 20) = 183 < 10*20), not simply assumed impossible.

This note does not repeat that duplicate-counting argument at nine primes;
it only establishes that if an analogous duplicate-realizability obstruction
holds at this scale too, it is doing real work -- without it, coherent
residuals like m = 210 or m = 390 would already have enough capacity to
pass the one necessary condition this whole research program keeps
measuring. Equivalently: this is a quantitative illustration of exactly why
direction-analysis item 1's "second induction on duplicates" is not an
optional refinement but the load-bearing part of that proof shape at this
scale, and it pins down numerically how coherent a residual would need to
be (roughly m with omega(m) = 4, i.e. concentrated on four of the nine
primes) before the duplicate question even becomes relevant.

## A natural next question (not attempted here)

Between "fully spread" (greedy, ratio 7.88) and "fully coherent" (single
class, ratio >= 22 by omega(m) = 4) there is presumably a continuum of
partially-coherent residuals -- a few large coset blocks plus scattered
points, as item 2 of the direction-analysis note already anticipates. The
real question for the existing maximiser search (Codex's
`code/extremal_slice_capacity.py` / `paired_slice_capacity.py`, already in
progress in this tree) is whether any such *partially* coherent residual is
both (a) actually realizable by a distinct family F without a duplicate,
and (b) has cap(R)/|R| anywhere near 22-28. This note supplies the
closed-form coherent-residual capacity formula above as a cheap building
block for exploring that continuum, but does not explore it.

## Follow-up: a naive "lean toward coherence" construction does not pay off

Tried the obvious next step: build an *actual* family biased toward the
coherent target instead of just asserting its theoretical capacity.
`code/safe_residue_construction.py` restricts every candidate modulus d to
residues a with a mod gcd(d, m) != 0; by the standard CRT solvability
criterion (x ≡ a (mod d) and x ≡ 0 (mod m) have a common solution iff
a ≡ 0 (mod gcd(d, m))), such a progression can never contain a multiple of
m, so the resulting family's uncovered set is guaranteed to contain every
multiple of m, automatically.

**Bug caught by this guarantee itself**: the first version of this script
treated every divisor d coprime to m (gcd = 1) as fully safe (no residues
excluded) and got uncovered = 2,854,002 for m = 390, *below* the coherent
target of 572,033 -- impossible under the stated guarantee, which is
exactly what caught the bug. The CRT criterion with gcd(d,m) = 1 gives
"a ≡ 0 (mod 1)", true for every a, so a divisor coprime to m has *no* safe
residue at all (every one of its residue classes, being equidistributed
modulo m by CRT, contains some multiple of m) -- the opposite of what the
first version assumed. Fixed by skipping every divisor coprime to m
entirely; an assertion that the construction's uncovered set actually does
contain every multiple of m now guards against a repeat.

With the fix, both m tried come out *worse*, not better, than the
unconstrained spread greedy (7.88 at |R| = 1,213,060):

| m | omega(m) | usable moduli | uncovered | cap(R)/\|R\| | fraction of 28 |
| - | - | - | - | - | - |
| 30 | 3 | 447 / 510 | 13,524,536 | 8.15 | 0.291 |
| 390 | 4 | 479 / 510 | 7,071,195 | 7.02 | 0.251 |
| (unconstrained greedy, for comparison) | -- | 510 / 510 | 1,213,060 | 7.88 | 0.281 |

Excluding every divisor coprime to m, and every unsafe residue of the
rest, costs so much raw covering power (residual 6-11x larger) that the
mild coherence gain is a wash at best (m=30) or a net loss (m=390). This is
a negative result for this specific naive approach, not for the broader
"partially coherent" idea: it shows that crudely forbidding "risky"
residues is not how to approach the coherent end of the spectrum, since it
throws away far more capacity than it protects. A real attempt would need
to choose *which* points to sacrifice more surgically (e.g. an actual LP
over the two objectives) rather than blanket-excluding whole divisors --
exactly the kind of search `code/extremal_slice_capacity.py` is already
doing, not duplicated here.

## Reproduction

The closed form and the table above were computed directly (no file
written), validated against brute-force enumeration at L = 2310:

```
python -c "
import math
L = 2310; m = 30
divs = [d for d in range(2, L+1) if L % d == 0 and d > 1]
cap = sum(L // (d*m//math.gcd(d,m)) for d in divs)
print(cap, L//m, cap/(L//m))
"
```

reproduces cap(R) = 691, |R| = 77, ratio 8.974, matching an independent
brute-force check over every divisor and residue.
