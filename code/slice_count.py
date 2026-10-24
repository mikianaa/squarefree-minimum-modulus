"""Trivial slice count for the extension L -> L q (q a new prime), least modulus 3.

In a cover of Z_{Lq} with distinct squarefree moduli >= 3, the q-free progressions cover the same set
in every slice x = c mod q; their uncovered set R has |R| >= u(L) := minimum uncovered count of a
family over L (all divisors >= 3 of L, one progression each).  The modulus q covers one slice; each of
the other q-1 slices must cover R with the shapes q d (d | L, d > 1) assigned to it, and a progression
mod d meets R in at most min(|R|, L/d) points.  Hence a cover at L q needs
      f(L, u) := Sum_{d | L, d > 1} min(u, L/d)  >=  (q - 1) u      for some feasible u >= u(L),
and f(L,u)/u is nonincreasing in u, so it suffices to check u = u(L).
Usage: python slice_count.py L q u_L
"""
import sys


def f(L, u): return sum(min(u, L // d) for d in range(2, L + 1) if L % d == 0)


if __name__ == "__main__":
    L, q, uL = map(int, sys.argv[1:4])
    val, need = f(L, uL), (q - 1) * uL
    bad = next((u for u in range(1, 10 * uL + 1) if f(L, u) < (q - 1) * u), None)
    print(f"L={L} q={q} u(L)={uL}: f={val} need={need} -> {'NO COVER at L*q' if val < need else 'not excluded'}; count fails for all u >= {bad}")
