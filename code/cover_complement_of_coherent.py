"""F-side probe for the coherent-residual route, take 2.

code/safe_residue_construction.py required every chosen residue to be
PROVABLY incapable of ever touching a multiple of m -- which excludes
every divisor coprime to m outright and throws away a lot of otherwise
useful covering power for no real benefit: it is perfectly fine (even
good) if the family also covers some multiples of m "by accident"; all
that is actually needed is that it covers every NON-multiple of m. This
script drops that restriction: target T = Z_L \\ {m | x}, and a plain
one-pass greedy (smallest divisors first, as in
code/greedy_cover_nine_primes.py) is scored only by how many T-points it
newly covers, with no penalty or restriction for also touching {m | x}.

Usage: python code/cover_complement_of_coherent.py --m 4199 --primes 2 3 5 7 11 13 17 19 --output data/cover-complement-4199.json
"""
import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]


def divisors_ge3(primes):
    L = 1
    for p in primes:
        L *= p
    ds = []
    for size in range(1, len(primes) + 1):
        for combo in itertools.combinations(primes, size):
            d = 1
            for p in combo:
                d *= p
            if d >= 3:
                ds.append(d)
    return L, sorted(ds)


def run(m, output=None, primes=PRIMES):
    L, seq = divisors_ge3(primes)
    assert L % m == 0
    need = np.ones(L, dtype=np.uint8)
    need.reshape(L // m, m)[:, 0] = 0  # multiples of m never need covering
    coherent_target = L // m
    t0 = time.monotonic()
    system = []
    for d in seq:
        view = need.reshape(L // d, d)
        counts = view.sum(axis=0)
        a = int(counts.argmax())
        gain = int(counts[a])
        if gain > 0:
            view[:, a] = 0
            system.append(dict(modulus=d, residue=a, gain=gain))
    seconds = round(time.monotonic() - t0, 1)
    uncovered_T = int(need.sum())  # non-multiples of m still uncovered: a genuine shortfall
    result = dict(
        L=L, m=m, moduli_available=len(seq), moduli_used=len(system),
        coherent_target_size=coherent_target, uncovered_T_points=uncovered_T,
        seconds=seconds,
        scope=("one-pass greedy targeting T = Z_L \\ {m|x} directly, with no restriction "
               "on also touching multiples of m. uncovered_T_points is the number of "
               "NON-multiples of m this family still fails to cover -- the true remaining "
               "gap to the coherent target, now measured correctly (the earlier "
               "safe_residue_construction.py's 'extra_beyond_coherent' conflated this with "
               "an unnecessary self-imposed restriction)."),
        system=system,
    )
    if output:
        Path(output).write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--primes", type=int, nargs="*", default=None)
    ap.add_argument("--output")
    args = ap.parse_args()
    primes = args.primes if args.primes else PRIMES
    out = run(args.m, args.output, primes=primes)
    print(json.dumps({k: v for k, v in out.items() if k != "system"}, indent=2))
