"""Coordinate-descent local search refinement of the one-pass greedy in
code/greedy_cover_nine_primes.py, at the same scale (L = 223092870, nine
primes). Each sweep revisits every modulus already in use and lets it pick
a better residue given everyone else's current choice (Gauss-Seidel style):
temporarily remove its contribution from a coverage-count array, recompute
its best residue against what is uncovered without it, then re-add. This
can only keep or reduce the uncovered count on each sweep -- never increase
it -- unlike the one-pass greedy, which commits each modulus once and never
revisits it.

This is architecturally what code/improve_cover_witness.py does at smaller
L, but using an O(L) reshape+sum per modulus (independent of d) instead of a
per-residue dictionary (which needs about sum_{d|L} d ~ 8*10^8 entries here
and is not memory-tractable).

Usage: python code/local_search_nine_primes.py --init data/greedy-nine-prime-smalld.json --sweeps 8 --output data/local-search-nine-prime.json
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np


def run(init_path, sweeps, output=None):
    init = json.loads(Path(init_path).read_text())
    L = init["L"]
    order = [e["modulus"] for e in sorted(init["system"], key=lambda e: e["modulus"])]
    # use the same ascending-value order as the stored system for determinism
    residue = {e["modulus"]: e["residue"] for e in init["system"]}
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
            a_old = residue[d]
            view = count.reshape(L // d, d)
            # remove d's own contribution
            col = view[:, a_old]
            newly_uncovered = int((col == 1).sum())
            col -= 1
            uncovered_count += newly_uncovered
            # best residue given everyone else
            zero_mask = (view == 0)
            gains = zero_mask.sum(axis=0)
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
    assert final_check == uncovered_count, (final_check, uncovered_count)
    result = dict(
        L=L, init_source=init_path, sweeps_requested=sweeps, sweeps_run=len(history) - 1,
        init_uncovered=start_uncovered, final_uncovered=uncovered_count,
        uncovered_density=uncovered_count / L,
        seconds=round(time.monotonic() - t0, 1), history=history,
        scope=("local-search-refined heuristic upper witness on u(L); NOT a lower-bound "
               "proof, NOT a claim of optimality."),
        system=[dict(modulus=d, residue=a) for d, a in sorted(residue.items())],
    )
    if output:
        Path(output).write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--init", required=True)
    ap.add_argument("--sweeps", type=int, default=8)
    ap.add_argument("--output")
    args = ap.parse_args()
    out = run(args.init, args.sweeps, args.output)
    print(json.dumps({k: v for k, v in out.items() if k not in ("system", "history")}))
