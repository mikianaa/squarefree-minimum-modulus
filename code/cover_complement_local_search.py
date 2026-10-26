"""Coordinate-descent refinement of code/cover_complement_of_coherent.py:
same target T = Z_L \\ {m|x}, same "no restriction on touching m" freedom,
revisited over multiple sweeps.

Usage: python code/cover_complement_local_search.py --init data/cover-complement-4199.json --sweeps 8
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np


def run(init_path, sweeps, output=None):
    init = json.loads(Path(init_path).read_text())
    L, m = init["L"], init["m"]
    residue = {e["modulus"]: e["residue"] for e in init["system"]}
    order = sorted(residue)
    count = np.zeros(L, dtype=np.int16)
    for d, a in residue.items():
        count.reshape(L // d, d)[:, a] += 1
    need = np.ones(L, dtype=np.uint8)
    need.reshape(L // m, m)[:, 0] = 0
    start_uncovered_T = int(((count == 0) & (need == 1)).sum())
    print(json.dumps(dict(init_uncovered_T=start_uncovered_T)), flush=True)

    t0 = time.monotonic()
    uncovered_T = start_uncovered_T
    history = [dict(sweep=0, uncovered_T=uncovered_T, seconds=0.0)]
    for s in range(1, sweeps + 1):
        sweep_start = uncovered_T
        for d in order:
            a_old = residue[d]
            view = count.reshape(L // d, d)
            need_view = need.reshape(L // d, d)
            col = view[:, a_old]
            newly_uncovered_T = int(((col == 1) & (need_view[:, a_old] == 1)).sum())
            col -= 1
            uncovered_T += newly_uncovered_T
            # gain only counted on points that actually need covering (not multiples of m)
            zero_mask = (view == 0) & (need_view == 1)
            gains = zero_mask.sum(axis=0).astype(np.int64)
            a_new = int(gains.argmax())
            gain = int(gains[a_new])
            view[:, a_new] += 1
            uncovered_T -= gain
            residue[d] = a_new
        history.append(dict(sweep=s, uncovered_T=uncovered_T,
                             seconds=round(time.monotonic() - t0, 1)))
        print(json.dumps(history[-1]), flush=True)
        if uncovered_T == sweep_start:
            print(json.dumps(dict(converged_at_sweep=s)), flush=True)
            break
    final_check = int(((count == 0) & (need == 1)).sum())
    assert final_check == uncovered_T
    result = dict(
        L=L, m=m, init_source=init_path, sweeps_requested=sweeps, sweeps_run=len(history) - 1,
        init_uncovered_T=start_uncovered_T, final_uncovered_T=uncovered_T,
        coherent_target_size=L // m, seconds=round(time.monotonic() - t0, 1), history=history,
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
