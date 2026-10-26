"""Weighted whole-block capacity test for a fixed optimistic core profile.

Retains the profile's per-divisor intersection capacities, but ignores positions
in the higher prime coordinates. Positive gap excludes the profile; a nonpositive
gap is only an inconclusive necessary condition, never a cover.
"""
import itertools
import argparse
import json
import math
from pathlib import Path

CORE_ODD = (1,3,5,7,15,21,35,105)


def weighted_obstruction(profile,high,high_pattern):
    residual = profile[:2]
    dims = [[p-1-(high_pattern[i]==e) for i,p in enumerate(high)] for e in (0,1)]
    demand = [residual[e]*math.prod(dims[e]) for e in (0,1)]
    odd_capacity = [0,0]
    even_capacities = []
    for k in range(1,len(high)+1):
        for subset in itertools.combinations(range(len(high)),k):
            volume = [math.prod(dims[e][i] for i in range(len(high)) if i not in subset) for e in (0,1)]
            for j,d in enumerate(CORE_ODD):
                if k==1 and d==1:
                    continue
                cap = residual if d==1 else profile[2+2*(j-1):4+2*(j-1)]
                pair = tuple(cap[e]*volume[e] for e in (0,1))
                for e in (0,1):
                    odd_capacity[e] += pair[e]
                even_capacities.append(pair)
    need = [demand[e]-odd_capacity[e] for e in (0,1)]
    weights = {(1,0),(0,1)}
    for c0,c1 in even_capacities:
        if c0 and c1:
            g = math.gcd(c0,c1)
            weights.add((c1//g,c0//g))
    best_gap = None
    best_weight = None
    for w0,w1 in weights:
        gap = w0*need[0]+w1*need[1]-sum(max(w0*c0,w1*c1) for c0,c1 in even_capacities)
        if best_gap is None or gap*sum(best_weight) > best_gap*(w0+w1):
            best_gap,best_weight = gap,(w0,w1)
    return dict(excluded=best_gap>0,gap=best_gap,weights=best_weight,
                demand=demand,odd_capacity=odd_capacity)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/core-geometry-profile-audit.json')
    ap.add_argument('--output',default='data/block-profile-capacity-510510.json')
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    results = []
    for row in data:
        for pattern in itertools.product((0,1),repeat=3):
            tests = [weighted_obstruction(p,(11,13,17),pattern) for p in row['easiest_profiles']]
            result = dict(core_pattern=row['pattern'],high_pattern=pattern,profiles=len(tests),
                          excluded=sum(t['excluded'] for t in tests),
                          surviving_indices=[i for i,t in enumerate(tests) if not t['excluded']],
                          tests=tests)
            results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('surviving_indices','tests')}),flush=True)
    Path(args.output).write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    main()
