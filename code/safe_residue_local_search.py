"""Coordinate-descent refinement of code/safe_residue_construction.py: does
local search get closer to the pure coherent target {m | x} than the
one-pass safe greedy did? Divisors coprime to m are still skipped entirely
(no safe residue exists for them, as established there); divisors sharing a
factor with m are revisited every sweep and allowed to change residue among
their safe options, Gauss-Seidel style (as in code/local_search_nine_primes.py).

Usage: python code/safe_residue_local_search.py --init data/safe-residue-390.json --sweeps 6
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def run(init_path, sweeps, output=None):
    init = json.loads(Path(init_path).read_text())
    L, m = init["L"], init["m"]
    residue = {e["modulus"]: e["residue"] for e in init["system"]}
    order = sorted(residue)
    count = np.zeros(L, dtype=np.int16)
    for d, a in residue.items():
        count.reshape(L // d, d)[:, a] += 1
    start_uncovered = int((count == 0).sum())
    print(json.dumps(dict(init_uncovered=start_uncovered)), flush=True)

    t0 = time.monotonic()
    uncovered_count = start_uncovered
    history = [dict(sweep=0, uncovered=uncovered_count, seconds=0.0)]
    for s in range(1, sweeps + 1):
        sweep_start = uncovered_count
        for d in order:
            g = gcd(d, m)
            a_old = residue[d]
            view = count.reshape(L // d, d)
            col = view[:, a_old]
            newly_uncovered = int((col == 1).sum())
            col -= 1
            uncovered_count += newly_uncovered
            zero_mask = (view == 0)
            gains = zero_mask.sum(axis=0).astype(np.int64)
            gains[0::g] = -1  # unsafe residues stay forbidden
            a_new = int(gains.argmax())
            gain = int(gains[a_new])
            view[:, a_new] += 1
            uncovered_count -= gain
            residue[d] = a_new
        history.append(dict(sweep=s, uncovered=uncovered_count,
                             seconds=round(time.monotonic() - t0, 1)))
        print(json.dumps(history[-1]), flush=True)
        if uncovered_count == sweep_start:
            print(json.dumps(dict(converged_at_sweep=s)), flush=True)
            break
    final_check = int((count == 0).sum())
    assert final_check == uncovered_count
    assert bool((count.reshape(L // m, m)[:, 0] == 0).all())  # multiples of m still untouched
    coherent_target = L // m
    result = dict(
        L=L, m=m, init_source=init_path, sweeps_requested=sweeps, sweeps_run=len(history) - 1,
        init_uncovered=start_uncovered, final_uncovered=uncovered_count,
        coherent_target_size=coherent_target, extra_beyond_coherent=uncovered_count - coherent_target,
        seconds=round(time.monotonic() - t0, 1), history=history,
        system=[dict(modulus=d, residue=a) for d, a in sorted(residue.items())],
    )
    if output:
        Path(output).write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--init", required=True)
    ap.add_argument("--sweeps", type=int, default=6)
    ap.add_argument("--output")
    args = ap.parse_args()
    out = run(args.init, args.sweeps, args.output)
    print(json.dumps({k: v for k, v in out.items() if k not in ("system", "history")}))
