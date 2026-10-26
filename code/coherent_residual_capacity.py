"""Closed-form cap(R)/|R| for a coherent (single full residue class) residual
R = {x in [0, L) : m | x}, for L the nine-prime primorial and various m | L.

See notes/coherence-capacity-tradeoff-2026-10-26.md.

For a divisor d > 1 of L: the points of R also lying in a fixed residue
class a mod d exist (by CRT) iff gcd(d, m) | a; whenever they exist there
are exactly L / lcm(d, m) of them, independent of which valid a (so a = 0,
always valid, already attains the maximum). Hence

    cap(R) = sum_{d | L, d > 1} L / lcm(d, m).

Validated against brute-force enumeration at L = 2310 in the note and in
test_coherent_residual_capacity.py.
"""
import itertools
import math

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]


def period(primes=PRIMES):
    L = 1
    for p in primes:
        L *= p
    return L


def divisors_gt1(primes=PRIMES):
    ds = []
    for size in range(1, len(primes) + 1):
        for combo in itertools.combinations(primes, size):
            d = 1
            for p in combo:
                d *= p
            if d > 1:
                ds.append(d)
    return ds


def cap_over_coherent(m, L=None, divs=None):
    L = L or period()
    divs = divs if divs is not None else divisors_gt1()
    cap = sum(L // (d * m // math.gcd(d, m)) for d in divs)
    Rsize = L // m
    return cap, Rsize


def table(primes=PRIMES):
    L = period(primes)
    divs = divisors_gt1(primes)
    rows = []
    for size in range(1, len(primes) + 1):
        for combo in itertools.combinations(primes, size):
            m = 1
            for p in combo:
                m *= p
            cap, Rsize = cap_over_coherent(m, L, divs)
            rows.append(dict(omega_m=size, m=m, primes_in_m=list(combo),
                              R_size=Rsize, cap=cap, ratio=cap / Rsize))
    return rows


if __name__ == "__main__":
    import json
    rows = table()
    rows.sort(key=lambda r: (r["omega_m"], r["m"]))
    print(json.dumps(rows[:5], indent=1))
    print("...")
    print(json.dumps(rows[-5:], indent=1))
