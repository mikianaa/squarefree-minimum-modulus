"""Reproduce a solver-free obstruction to minimum 3 with at most six primes.

1. Enumerate axis lines to bound u(210), for each anchor parity pattern.
2. Propagate just these total bounds through 11, retaining the two parity counts
   during the calculation. Every possible allocation of even shapes is included.
3. Test every two-layer split of that lower bound against the capacity at 13.

See notes/six-prime-obstruction-2026-10-25.md for the proof of the reductions.
"""
import json
from pathlib import Path
from three_axis_bound import coupled_bound
from parity_capacity import capacities, extend


def final_deficits(pattern, lower):
    cap = capacities((3,5,7,11),pattern)
    return [min(12*lower - (a if b == 0 else lower-a)
                - sum(min(a,c[0])+min(lower-a,c[1])+max(min(a,c[0]),min(lower-a,c[1])) for c in cap)
                for a in range(lower+1)) for b in (0,1)]


def verify():
    out = []
    for row in coupled_bound()['cases']:
        u0 = row['uncovered_lower']
        for b11 in (0,1):
            pattern = list(row['pattern'])+[b11]
            lower_pairs = [(a,u0-a) for a in range(u0+1)]
            frontier = extend((3,5,7),row['pattern'],lower_pairs,11,b11)
            lower = min(map(sum,frontier))
            margins = final_deficits(pattern,lower)
            assert min(margins) > 0
            result = dict(pattern=pattern,uncovered_210_lower=u0,
                          uncovered_2310_lower=lower,extension_13_deficits=margins)
            out.append(result)
            print(json.dumps(result),flush=True)
    return out


if __name__ == '__main__':
    result = verify()
    path = Path(__file__).resolve().parents[1]/'data/six-prime-obstruction.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
