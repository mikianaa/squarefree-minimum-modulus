"""Follow-up to notes/coherence-capacity-tradeoff-2026-10-26.md: try to
actually realize a residual close to the coherent {m | x} case with a real
distinct-moduli family, instead of just asserting its theoretical capacity.

For a target divisor m, restrict every candidate modulus d to "safe"
residues a (those with a mod gcd(d, m) != 0): a progression at a safe
residue can never contain a multiple of m, by construction, so a family
built only from safe choices has its uncovered residual automatically
containing every multiple of m -- and the question becomes how much of the
REST (the non-multiples of m) it also fails to cover. If it covers all of
the rest, the realized residual is exactly {m | x} (the coherent case,
already known to need a duplicate to be exact -- but this safe-greedy
construction is NOT claiming to realize it exactly, just approximating it
with an actual family and measuring what's left).

Usage: python code/safe_residue_construction.py --m 390 --output data/safe-residue-390.json
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


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def run(m, output=None, primes=PRIMES):
    L, seq = divisors_ge3(primes)
    assert L % m == 0
    uncovered = np.ones(L, dtype=np.uint8)
    system = []
    t0 = time.monotonic()
    for d in seq:
        g = gcd(d, m)
        if g == 1:
            # x ≡ a (mod d) and x ≡ 0 (mod m) are simultaneously solvable whenever
            # a ≡ 0 (mod gcd(d, m)); with gcd=1 that's every a, so EVERY residue of a
            # divisor coprime to m risks containing a multiple of m. No safe choice
            # exists; skip this divisor entirely rather than silently use an unsafe one.
            continue
        view = uncovered.reshape(L // d, d)
        counts = view.sum(axis=0).astype(np.int64)
        # mask out unsafe residues (those ≡ 0 mod g): solvable jointly with x≡0 (mod m)
        counts[0::g] = -1
        a = int(counts.argmax())
        gain = int(counts[a])
        if gain > 0:
            view[:, a] = 0
            system.append(dict(modulus=d, residue=a, gain=gain))
    seconds = round(time.monotonic() - t0, 1)
    total_uncovered = int(uncovered.sum())
    coherent_target = L // m
    # invariant: every multiple of m must still be uncovered, by construction
    assert bool(uncovered.reshape(L // m, m)[:, 0].all())
    extra = total_uncovered - coherent_target
    assert extra >= 0, extra
    result = dict(
        L=L, m=m, moduli_available=len(seq), moduli_used=len(system),
        uncovered=total_uncovered, coherent_target_size=coherent_target,
        extra_beyond_coherent=extra, seconds=seconds,
        scope=("safe-residue-restricted greedy: an ACTUAL distinct family whose uncovered "
               "set automatically contains every multiple of m; 'extra' is how many "
               "additional (non-multiple-of-m) points it also failed to cover, i.e. how far "
               "short of the pure coherent target this achievable construction falls."),
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
