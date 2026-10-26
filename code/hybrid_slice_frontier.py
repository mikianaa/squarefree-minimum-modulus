"""Scalable parity extension using compression and even-resource allocation.

Odd resources have total capacity at most the sharp F_e(r_e). The even
resources are assigned to one layer, with individual caps min(r_e,C_S(e))
and total cap F_e(r_e) in that layer. This is a necessary relaxation.
For any allocation, both cap sums and F_e are concave and vanish at zero;
their minimum is concave too. Clipped deficits are thus nondecreasing, so
componentwise lower residual pairs suffice for quantitative propagation.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from extremal_slice_capacity import maximum_capacity
from parity_capacity import pareto


def individual_caps(dimensions,counts):
    k = len(dimensions[0])
    return [tuple(min(counts[e],math.prod(dimensions[e][i] for i in range(k) if not s>>i&1))
                  for e in (0,1)) for s in range(1,1<<k)]


def pair_outputs(dimensions,counts,q,b):
    totals = [maximum_capacity(d,r) for d,r in zip(dimensions,counts)]
    caps = individual_caps(dimensions,counts)
    budget = sum(c[1] for c in caps)
    assert budget<10_000_000, 'explicit memory guard'
    impossible = -10**12
    dp = np.full(budget+1,impossible,dtype=np.int64)
    dp[0] = 0
    used = 0
    for gain,cost in caps:
        # Copy source before modifying a possibly overlapping destination.
        values = dp[:used+1].copy()+gain
        dp[cost:cost+used+1] = np.maximum(dp[cost:cost+used+1],values)
        used += cost
    costs = np.flatnonzero(dp>=0)
    even0 = np.minimum(totals[0],dp[costs])
    even1 = np.minimum(totals[1],budget-costs)
    left0 = np.maximum(0,(q-1-(b==0))*counts[0]-totals[0]-even0)
    left1 = np.maximum(0,(q-1-(b==1))*counts[1]-totals[1]-even1)
    return [list(p) for p in pareto(zip(left0.tolist(),left1.tolist()))]


def extend(dimensions,frontier,q,b):
    out = set()
    for pair in frontier:
        out.update(map(tuple,pair_outputs(dimensions,pair,q,b)))
    result = [list(p) for p in pareto(out)]
    return dict(frontier=result,uncovered_lower=min(map(sum,result)))


def run(source,output,q):
    start = time.monotonic()
    raw = json.loads(Path(source).read_text())
    rows = raw['at30030'] if 'at30030' in raw else raw['cases']
    primes = [3,5,7,11,13] if 'at30030' in raw else raw['odd_primes']
    assert q not in primes
    results = []
    for row in rows:
        if row.get('status')!='integer_only':
            continue
        pattern = row['pattern']
        dims = [[p-1-(c==e) for p,c in zip(primes,pattern)] for e in (0,1)]
        for b in (0,1):
            result = dict(pattern=pattern+[b],status='integer_only',**extend(dims,row['frontier'],q,b))
            results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k!='frontier'}),flush=True)
    result = dict(odd_primes=primes+[q],cases=results,seconds=round(time.monotonic()-start,3),
                  global_lower=min(r['uncovered_lower'] for r in results),
                  scope='compression plus individual even-cylinder capacity allocation; feasibility is not a cover')
    Path(output).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/core-capacity-frontiers.json')
    ap.add_argument('--prime',type=int,default=17)
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    run(a.source,a.output,a.prime)
