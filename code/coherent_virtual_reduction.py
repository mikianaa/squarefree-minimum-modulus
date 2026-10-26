"""Exact reduction of the static two-family model for a COHERENT target
R = {x in [0,L) : m | x} down to a tiny equivalent problem, and an exact
CP-SAT solve of that reduced problem.

Key fact (proved in the accompanying note): writing x = m*k for
k in [0, L/m), a real divisor d | L (d > 1) with g = gcd(d, m) can cover
R-points only via the compatible residue class a (a multiple of g), and
doing so restricts k to a SINGLE residue class modulo d' = d/g (since
m/g is invertible mod d'). So real divisor d acts exactly like a "virtual"
divisor d' of L/m on the k-domain. Divisors of m itself (d' = 1) cover
EVERY k unconditionally, for free: there are (number of divisors of m,
excluding 1 and excluding 2 if 2 | m) of them, each usable with residue 0.
Every other divisor of L maps to some d' > 1 dividing L/m, and distinct
real divisors of L mapping to the SAME d' are independent "copies" (since
they are different integers, each gets its own residue choice): the number
of real divisors of L mapping to a given d' > 1 is exactly 2^omega(m) (one
per subset of m's own prime factors to multiply in), independent of d'.

So the whole problem reduces to: on Z_{L/m} (period L/m), with
copies = 2^omega(m) independent "slots" per divisor v > 1 of L/m, does
some choice of residues give every k at least (q - 1 - free) additional
hits, where free = (number of divisors of m that are usable, i.e. >= 3)?

This is tiny compared to the original (L/m instead of L) whenever m is a
product of several of the largest primes, and is exactly solvable by
CP-SAT at that scale.
"""
import argparse
import itertools
import json
import time
from pathlib import Path

from ortools.sat.python import cp_model

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]


def divisors(n, primes=PRIMES):
    ds = {1}
    for p in primes:
        if n % p == 0:
            ds |= {d * p for d in ds if n % (d * p) == 0}
    return sorted(ds)


def solve(m, q, time_limit=120.0, workers=8, primes=PRIMES):
    L = 1
    for p in primes:
        L *= p
    assert L % m == 0
    Lm = L // m
    # G's divisors d | L need only d > 1 (not F's d >= 3: G's real moduli are the
    # shapes q*d, automatically >= 2q >= 3, so d=2 is a legitimate G choice too)
    m_divs = [d for d in divisors(m) if d > 1]
    free = len(m_divs)
    need_extra = (q - 1) - free
    virtual_divs = [d for d in divisors(Lm) if d > 1]
    omega_m = len([p for p in primes if m % p == 0])
    copies = 2 ** omega_m

    model = cp_model.CpModel()
    y = {}
    for v in virtual_divs:
        for c in range(copies):
            y[(v, c)] = [model.NewBoolVar(f"y_{v}_{c}_{a}") for a in range(v)]
            model.AddExactlyOne(y[(v, c)])
    for k in range(Lm):
        model.Add(sum(y[(v, c)][k % v] for v in virtual_divs for c in range(copies)) >= need_extra)

    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = workers
    t0 = time.monotonic()
    status = s.Solve(model)
    seconds = round(time.monotonic() - t0, 1)
    num_vars = sum(len(y[(v, c)]) for v in virtual_divs for c in range(copies))
    result = dict(
        L=L, m=m, q=q, Lm=Lm, free=free, need_extra=need_extra,
        omega_m=omega_m, copies=copies, num_virtual_divs=len(virtual_divs),
        num_vars=num_vars, status=s.StatusName(status), seconds=seconds,
    )
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        assignment = {(v, c): int(next(a for a in range(v) if s.Value(y[(v, c)][a])))
                      for v in virtual_divs for c in range(copies)}
        for k in range(Lm):
            mult = sum(1 for v in virtual_divs for c in range(copies) if k % v == assignment[(v, c)])
            assert mult >= need_extra, (k, mult, need_extra)
        result["status"] = "SAT"

        # Translate back to REAL (divisor, residue) pairs on Z_L and verify directly
        # on actual points x = m*k, not just virtually on k -- an independent check of
        # the reduction itself, not just of the reduced model's internal consistency.
        m_divs_all = divisors(m)  # includes 1; copies index over these, 0..copies-1
        assert len(m_divs_all) == copies
        real_system = []
        for v in virtual_divs:
            for c, s_ in enumerate(m_divs_all):
                b = assignment[(v, c)]
                d = v * s_
                m_prime = m // s_
                inv = pow(m_prime, -1, v) if v > 1 else 0
                a1 = (m_prime * b) % v if v > 1 else 0
                a_real = s_ * a1
                real_system.append((d, a_real))
        for d in m_divs:
            real_system.append((d, 0))  # the "free" baseline divisors of m
        assert len({d for d, _ in real_system}) == len(real_system), "duplicate real divisor!"
        # direct check on every actual multiple of m, using real arithmetic only
        for k in range(Lm):
            x = m * k
            mult = sum(1 for d, a in real_system if x % d == a)
            assert mult >= q - 1, (x, mult, q - 1)
        result["real_system_size"] = len(real_system)
        result["real_verification"] = "passed: every real multiple of m gets >= q-1 real divisors matching"
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--time", type=float, default=120.0)
    ap.add_argument("--primes", type=int, nargs="*", default=None,
                     help="override the 9-prime default, for small-scale validation")
    ap.add_argument("--output")
    args = ap.parse_args()
    primes = args.primes if args.primes else PRIMES
    out = solve(args.m, args.q, args.time, primes=primes)
    print(json.dumps(out, indent=2))
    if args.output:
        Path(args.output).write_text(json.dumps(out, indent=2) + "\n")
