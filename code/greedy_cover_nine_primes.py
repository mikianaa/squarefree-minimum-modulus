"""Direction 4 of notes/direction-analysis-2026-10-25.md: take M_sf >= 3 seriously
and search constructively, with a heuristic rather than an exact solver, at the
first open primorial.

L = 223092870 = 2*3*5*7*11*13*17*19*23 (nine primes). The eight-prime
obstruction (notes/eight-prime-obstruction-2026-10-25.md) proves every
counterexample omitting modulus 2 needs >= 9 distinct prime factors; this L is
the smallest untested case. This script does NOT prove anything either way;
it is a one-pass greedy construction attempt, reported honestly as a
heuristic upper witness on u(L) (the minimum uncovered count over all
distinct families with moduli >= 3 dividing L), using every available
divisor >= 3 of L (510 of them) in one fixed order, each committed once to
its locally best residue against the state so far (no backtracking, no local
search afterwards -- this is cheaper than code/improve_cover_witness.py,
which needs a per-residue lookup table and is not memory-tractable at this L;
here each modulus's best residue is found in O(L) time via one reshape+sum,
independent of d, using numpy).

Usage: python code/greedy_cover_nine_primes.py [--order smalld|bigd|random] [--seed N]
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
                ds.append((len(combo), d))
    return L, ds


def order_divisors(ds, order, rng):
    if order == "smalld":
        return [d for _, d in sorted(ds, key=lambda t: (t[0], t[1]))]
    if order == "bigd":
        return [d for _, d in sorted(ds, key=lambda t: (-t[0], -t[1]))]
    if order == "random":
        ds2 = [d for _, d in ds]
        rng.shuffle(ds2)
        return ds2
    raise ValueError(order)


def run(order="smalld", seed=1, output=None):
    rng = __import__("random").Random(seed)
    L, ds = divisors_ge3(PRIMES)
    seq = order_divisors(ds, order, rng)
    assert len(seq) == 2 ** len(PRIMES) - 2
    uncovered = np.ones(L, dtype=np.uint8)
    total = L
    system = []
    t0 = time.monotonic()
    milestones = []
    for i, d in enumerate(seq):
        view = uncovered.reshape(L // d, d)
        counts = view.sum(axis=0)
        a = int(counts.argmax())
        gain = int(counts[a])
        if gain:
            view[:, a] = 0
            total -= gain
            system.append(dict(modulus=d, residue=a, gain=gain))
        if i % 50 == 0 or i == len(seq) - 1:
            milestones.append(dict(step=i, modulus=d, uncovered=total,
                                    seconds=round(time.monotonic() - t0, 1)))
            print(json.dumps(milestones[-1]), flush=True)
    result = dict(
        L=L, primes=PRIMES, order=order, seed=seed,
        moduli_available=len(seq), moduli_used=len(system),
        uncovered=total, uncovered_density=total / L,
        seconds=round(time.monotonic() - t0, 1),
        milestones=milestones,
        scope=("one-pass greedy heuristic upper witness on u(L); NOT a lower-bound "
               "proof, NOT a claim of optimality, NOT verified against a full-period "
               "pass yet (see the separate verification step)."),
        system=system,
    )
    if output:
        Path(output).write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", choices=["smalld", "bigd", "random"], default="smalld")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--output")
    args = ap.parse_args()
    out = run(args.order, args.seed, args.output)
    print(json.dumps({k: v for k, v in out.items() if k not in ("system", "milestones")}))
