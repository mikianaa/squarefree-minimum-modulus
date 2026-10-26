"""Exact dual lower bounds, independent of the allocation knapsack producer.

For F_e bounding all cylinder capacities, A_e=D_e*r_e-2F_e and
B_e=F_e+min(A_e,0), every 0<=t_e<=1 gives the valid total deficit bound
sum max(A_e,0) + sum t_e*B_e - sum_S max(t_0*c_S0,t_1*c_S1).
Try integer-ratio multipliers; each reported fraction is a direct certificate.
"""
import itertools
import json
import math
from pathlib import Path
import time


def explicit_prefix_weights(dimensions):
    weights = sorted((2**sum(v==0 for v in point) for point in itertools.product(*(range(d) for d in dimensions))),reverse=True)
    prefix = [0]
    for w in weights:
        prefix.append(prefix[-1]+w)
    return prefix


def dual_bound(dimensions,counts,q,b,prefixes=None):
    if prefixes is None:
        prefixes = list(map(explicit_prefix_weights,dimensions))
    totals = [prefixes[e][counts[e]]-counts[e] for e in (0,1)]
    k = len(dimensions[0])
    caps = [tuple(min(counts[e],math.prod(dimensions[e][i] for i in range(k) if not mask>>i&1)) for e in (0,1))
            for mask in range(1,1<<k)]
    A = [(q-1-(b==e))*counts[e]-2*totals[e] for e in (0,1)]
    constant = sum(max(a,0) for a in A)
    B = [totals[e]+min(A[e],0) for e in (0,1)]
    multipliers = {(0,1),(1,0),(1,1)}|{(v,u) for u,v in caps if u or v}
    best = dict(numerator=constant,denominator=1,weights=[0,0])
    for a,c in multipliers:
        denominator = max(a,c)
        numerator = constant*denominator+B[0]*a+B[1]*c-sum(max(a*u,c*v) for u,v in caps)
        if numerator*best['denominator']>best['numerator']*denominator:
            best = dict(numerator=numerator,denominator=denominator,weights=[a,c])
    best['integer_lower'] = -(-best['numerator']//best['denominator'])
    return best


def main():
    start = time.monotonic()
    source = json.loads(Path('data/core-capacity-frontiers.json').read_text())['at30030']
    producer = json.loads(Path('data/hybrid-frontier-510510.json').read_text())
    results = []
    for row in source:
        if row['status']!='integer_only':
            continue
        p = row['pattern']
        dims = [[v-1-(a==e) for v,a in zip((3,5,7,11,13),p)] for e in (0,1)]
        prefixes = list(map(explicit_prefix_weights,dims))
        for b in (0,1):
            certificates = [dict(input_counts=r,**dual_bound(dims,r,17,b,prefixes)) for r in row['frontier']]
            lower = min(r['integer_lower'] for r in certificates)
            target = next(r for r in producer['cases'] if r['pattern']==p+[b])
            assert 0<=lower<=target['uncovered_lower']
            result = dict(pattern=p+[b],dual_lower=lower,knapsack_lower=target['uncovered_lower'],
                          certificates=certificates)
            results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k!='certificates'}),flush=True)
    result = dict(global_lower=min(r['dual_lower'] for r in results),verified=True,cases=results,
                  seconds=round(time.monotonic()-start,3),
                  scope='independent exact dual lower bounds; not verification of every knapsack output frontier')
    Path('data/audit-hybrid-slice-510510.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
