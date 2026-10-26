"""Aligned construction test (direction 3 of notes/direction-analysis-2026-10-25.md).

For primes p_1 < ... < p_k with L = prod p_i and next prime q, the "G-residues
zero" slice model needs every point of the q-free residual R to have at least
q-1 incidences from the divisors d | L, d > 1, each contributing 1 if d | x
(choosing the trivial residue 0 mod d for every d). A point x has exactly
2^omega(gcd(x,L)) - 1 such divisors, so sufficiency needs
2^omega(gcd(x,L)) - 1 >= q - 1, i.e. omega(gcd(x,L)) >= s := ceil(log2 q).
Hence this construction needs a distinct family F (moduli >= 3 dividing L)
whose uncovered residual R is contained in {x : omega(gcd(x,L)) >= s}, i.e. F
covers T := {x : omega(gcd(x,L)) <= s-1} exactly.

This script computes, by exact integer arithmetic (no floats, no solver):

* |T| and |R_allowed| = L - |T| for k = 4..8 (the rows of the table in
  section 1 of the direction note, q = 11..23);
* for each k, the best known rigorous lower bound on u(L) (the true minimum
  uncovered count for ANY distinct family with moduli >= 3 of L), taken from
  already-certified results recorded in PROGRESS.md;
* the verdict: if |R_allowed| < (that lower bound), no family can confine its
  uncovered residual to R_allowed, so the aligned construction is impossible
  at that k by a one-line counting argument -- independent of q-1-incidence
  geometry. If |R_allowed| >= the lower bound, counting alone does not decide
  it (ambiguous; needs the full non-aligned "static two-family" treatment,
  direction 2 of the same note).

Also computes, for k = 8 (the first ambiguous case, L = 9699690, q = 23), the
exact size of the single largest pattern class (points coprime to L) to show
why the trivial "residue 0 on every small divisor" construction fails badly:
that one class alone is already far larger than |R_allowed|, and no modulus
using residue 0 can ever intersect it.
"""
import json
import math
from itertools import combinations
from pathlib import Path

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]

# Certified / best-recorded rigorous lower bounds on u(L) for the q-free
# period L, taken from PROGRESS.md (exact values where the true minimum is
# known; otherwise the strongest published integer-only or audited lower
# bound). Not asserted optimal except where noted.
KNOWN_U_LOWER = {
    210: (22, "exact (CP-SAT OPTIMAL)"),
    2310: (121, "exact (CP-SAT OPTIMAL, five anchor patterns)"),
    30030: (501, "integer-only lower bound, quantitative-core-extension-2026-10-25.md"),
    510510: (510, "integer-only lower bound, quantitative-core-extension-2026-10-25.md"),
    9699690: (2, "minimum branch bound (011 pattern), eight-prime-obstruction-2026-10-25.md; "
                 "other branches give 15, 305, 6696"),
}


def agg(remaining, budget, qm1):
    """Sum over S subset of `remaining` with |S| <= budget of prod_{q in remaining\\S} (q-1)."""
    total = 0
    rl = list(remaining)
    for size in range(0, min(budget, len(rl)) + 1):
        for subset in combinations(rl, size):
            sset = set(subset)
            prod = 1
            for q in rl:
                if q not in sset:
                    prod *= qm1[q]
            total += prod
    return total


def table():
    rows = []
    for k in range(4, 9):
        ps = PRIMES[:k]
        q = PRIMES[k]
        qm1 = {p: p - 1 for p in ps}
        L = 1
        for p in ps:
            L *= p
        s = math.ceil(math.log2(q))
        t = s - 1
        T = agg(ps, t, qm1)
        R_allowed = L - T
        lower, source = KNOWN_U_LOWER[L]
        decided = R_allowed < lower
        rows.append(dict(
            k=k, primes=ps, L=L, q=q, s=s, t=t,
            T=T, R_allowed=R_allowed,
            R_allowed_density=R_allowed / L,
            known_u_lower=lower, known_u_lower_source=source,
            decided_infeasible_by_counting=decided,
        ))
    return rows


def coprime_class_size(k):
    ps = PRIMES[:k]
    prod = 1
    for p in ps:
        prod *= (p - 1)
    return prod


def run():
    rows = table()
    coprime_sizes = {row["k"]: coprime_class_size(row["k"]) for row in rows}
    result = dict(
        scope=("Aligned construction test (G-residues zero) for the slice extension "
               "L -> Lq: exact counting only, no solver. Decides k <= 6 infeasible; "
               "k = 7, 8 remain ambiguous by counting alone."),
        rows=rows,
        coprime_class_size=coprime_sizes,
        note=("coprime_class_size[k] is |{x mod L : gcd(x,L)=1}|, the largest single "
              "pattern class; it is already far above R_allowed for k=8, and no "
              "residue-0 modulus can ever intersect it, showing the naive "
              "'residue 0 on small divisors' realization of this construction fails."),
    )
    for row in rows:
        print(json.dumps({k: v for k, v in row.items() if k != "primes"}))
    Path("data/aligned-construction-test.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    run()
