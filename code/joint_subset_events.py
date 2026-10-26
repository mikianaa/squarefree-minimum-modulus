"""Exact joint threshold-event separation on small core subsets.

A resource chooses one part of the residue partition of the selected points.
Boolean convolution enumerates attainable hit vectors; downward closure gives
all feasible lower thresholds. No solver is used for this finite check.
"""
import argparse
from collections import Counter
from functools import lru_cache
import itertools
import json
from pathlib import Path
import time
import numpy as np


def partitions(points,shapes):
    result = []
    for d in shapes:
        parts = tuple(sorted(sum(1<<i for i,x in enumerate(points) if x%d==a) for a in {x%d for x in points}))
        result.append(parts)
    return tuple(sorted(result))


@lru_cache(maxsize=256)
def attainable(partition_list,size):
    limit = len(partition_list)
    possible = np.zeros((limit+1,)*size,dtype=bool)
    possible[(0,)*size] = True
    for parts in partition_list:
        expanded = np.zeros_like(possible)
        for part in parts:
            before = tuple(slice(None,-1) if part>>i&1 else slice(None) for i in range(size))
            after = tuple(slice(1,None) if part>>i&1 else slice(None) for i in range(size))
            expanded[after] |= possible[before]
        possible = expanded
    for axis in range(size):
        possible = np.flip(np.maximum.accumulate(np.flip(possible,axis=axis),axis=axis),axis=axis)
    return possible


def feasible_integer(points,shapes,thresholds):
    """Independent set-of-tuples enumeration, truncated at required thresholds."""
    states = {(0,)*len(points)}
    target = tuple(thresholds)
    for d in shapes:
        parts = [tuple(int(x%d==a) for x in points) for a in {x%d for x in points}]
        states = {tuple(min(t,v+a) for t,v,a in zip(target,state,part)) for state in states for part in parts}
    return target in states


def separate(data,size=4):
    start = time.monotonic()
    source,witness = data['input'],data['result']['point_primal_witness']
    R,n,scale = source['core_residual'],len(source['high_primes']),witness['denominator']
    groups = (1<<n)-1
    hist = np.zeros((len(R),groups,17),dtype=np.int64)
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[row['point_index'],h,v] += row['weight']
    tails = np.cumsum(hist[:,:,::-1],axis=2)[:,:,::-1]
    # Only the largest threshold at each tail probability can be strongest.
    options = [[np.flatnonzero(hist[x,h,1:])+1 for h in range(groups)] for x in range(len(R))]
    core = [d for d in range(3,211) if 210%d==0]
    violations = []
    combos = candidates = profiles = 0
    for selected in itertools.combinations(range(len(R)),size):
        points = [R[x] for x in selected]
        for kind,shapes in ((0,core),(1,[1,2]+core)):
            possible = None
            for h in (range(n) if kind==0 else range(n,groups)):
                combos += 1
                choices = [options[x][h] for x in selected]
                if any(len(v)==0 for v in choices):
                    continue
                mesh = np.meshgrid(*choices,indexing='ij')
                score = sum(tails[x,h,ts] for x,ts in zip(selected,mesh))
                active = score>(size-1)*scale
                if not active.any():
                    continue
                indices = tuple(ts[active] for ts in mesh)
                candidates += len(indices[0])
                if possible is None:
                    possible = attainable(partitions(points,shapes),size)
                    profiles += 1
                bad = ~possible[indices]
                if bad.any():
                    gaps = score[active][bad]-(size-1)*scale
                    winner = int(gaps.argmax())
                    thresholds = [int(ts[bad][winner]) for ts in indices]
                    gap = int(gaps[winner])
                    assert not feasible_integer(points,shapes,thresholds)
                    cap = sum(max(Counter(x%d for x in points).values()) for d in shapes)
                    violations.append(dict(axis=h,points=list(selected),thresholds=thresholds,
                                           violation_numerator=gap,projection_sum_capacity=cap,
                                           detected_by_sum_capacity=sum(thresholds)>cap))
    violations.sort(key=lambda r:r['violation_numerator'],reverse=True)
    return dict(violations=violations,subset_size=size,subset_group_combinations_checked=combos,
                candidate_event_tuples_checked=candidates,projection_profiles_requested=profiles,
                denominator=scale,seconds=time.monotonic()-start,
                scope='exact joint event compatibility on fixed-size subsets; not a covering-system obstruction')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--size',type=int,default=4)
    args = ap.parse_args()
    result = separate(json.loads(Path(args.source).read_text()),args.size)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='violations'}|dict(violations=len(result['violations']),largest=result['violations'][:3])))
