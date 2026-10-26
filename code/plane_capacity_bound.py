"""Majorization bound retaining indivisible single-coordinate plane weights.

The seven odd core shapes and their assigned even doubles give loads on the
labels of each higher prime. The most balanced integer vector compatible with
the largest individual weights is majorized by every actual bin-load vector.
Convexity of (r-x-y-z)_+ then bounds the deficit after all plane classes.
Remaining line/point classes are granted their full total capacities.
All numerical arrays are exact int64; no floating solver result is used.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np

CORE_ODD = (3,5,7,15,21,35,105)
ALLOCATIONS = tuple(itertools.product((0,1),repeat=7))


def balanced_lower(weights,bins):
    """An integer load vector majorized by every packing of weights into bins."""
    assert bins>0 and all(isinstance(w,int) and w>=0 for w in weights)
    weights = sorted(weights,reverse=True)
    if len(weights)<=bins:
        return weights+[0]*(bins-len(weights))
    fixed = []
    total = sum(weights)
    for w in weights:
        left = bins-len(fixed)
        if w*left>total:
            fixed.append(w)
            total -= w
        else:
            q,t = divmod(total,left)
            return fixed+[q+1]*t+[q]*(left-t)
    raise AssertionError('unreachable')


def plane_deficits(profile,high,pattern,odd_prime3=None):
    assert len(high)==len(pattern)==3 and len(profile)==16
    assert all(isinstance(x,int) and 0<=x<=210 for x in profile)
    assert max(high)<=101, 'int64 safety restriction for this finite experiment'
    r = profile[:2]
    caps = [profile[2+2*j:4+2*j] for j in range(7)]
    dims = [[p-1-(pattern[i]==e) for i,p in enumerate(high)] for e in (0,1)]
    arrays = []
    for e in (0,1):
        if not r[e]:
            arrays.append(np.zeros((128,)*3,dtype=np.int64))
            continue
        histograms = []
        for axis,D in enumerate(dims[e]):
            rows = []
            for bits in ALLOCATIONS:
                odd_weights = [c[e] for c in caps]
                if odd_prime3 is not None:
                    assert 0<=odd_prime3[axis][e]<=caps[0][e]
                    odd_weights[0] = odd_prime3[axis][e]
                weights = odd_weights+[c[e] for bit,c in zip(bits,caps) if bit==e]
                loads = balanced_lower(weights,D)
                rows.append([loads.count(i) for i in range(r[e])])
            histograms.append(np.array(rows,dtype=np.int64))
        loss = np.fromfunction(lambda i,j,k:np.maximum(0,r[e]-i-j-k),(r[e],)*3,dtype=int).astype(np.int64)
        arrays.append(np.einsum('ai,bj,ck,ijk->abc',*histograms,loss,optimize=True))
    return arrays,dims,caps


def bound(profile,high,pattern,odd_prime3=None):
    start = time.monotonic()
    arrays,dims,caps = plane_deficits(profile,high,pattern,odd_prime3)
    totals = arrays[0]+arrays[1]
    index = np.unravel_index(int(totals.argmin()),totals.shape)
    minimum = int(totals[index])
    capacity = 0
    for k in (2,3):
        for subset in itertools.combinations(range(3),k):
            volume = [math.prod(dims[e][i] for i in range(3) if i not in subset) for e in (0,1)]
            for cap in [profile[:2]]+caps:
                pair = [cap[e]*volume[e] for e in (0,1)]
                capacity += sum(pair)+max(pair)
    gap = minimum-capacity
    return dict(plane_deficit_lower=minimum,line_point_capacity_upper=capacity,gap=gap,excluded=gap>0,
                minimizing_allocation_indices=[int(i) for i in index],
                minimizing_layer_deficits=[int(a[index]) for a in arrays],
                seconds=round(time.monotonic()-start,3),
                scope='necessary capacity bound; nonpositive gap does not give a covering')


def refined_bound(profile,high,pattern):
    """Keep the two plane deficits coupled; allocate later even classes exactly."""
    arrays,dims,caps = plane_deficits(profile,high,pattern)
    sentinel = 10**12
    least_second = np.full(int(arrays[0].max())+1,sentinel,dtype=np.int64)
    np.minimum.at(least_second,arrays[0].ravel(),arrays[1].ravel())
    frontier = []
    previous = sentinel
    for x,y in enumerate(least_second):
        if y<previous:
            frontier.append((x,int(y)))
            previous = y
    items = []
    for k in (2,3):
        for subset in itertools.combinations(range(3),k):
            volume = [math.prod(dims[e][i] for i in range(3) if i not in subset) for e in (0,1)]
            items.extend(tuple(c[e]*volume[e] for e in (0,1)) for c in [profile[:2]]+caps)
    odd = [sum(c[e] for c in items) for e in (0,1)]
    # Put every even resource in parity 1, then transfer a subset to parity 0.
    dp = np.full(odd[1]+1,-sentinel,dtype=np.int64)
    dp[0] = 0
    for gain,cost in items:
        if cost:
            dp[cost:] = np.maximum(dp[cost:],dp[:-cost]+gain)
        else:
            dp += gain
    dp = np.maximum.accumulate(dp)
    survivors = []
    margins = []
    for x,y in frontier:
        budget = odd[1]-max(0,y-odd[1])
        need = max(0,x-odd[0])
        if budget<0:
            margins.append(-budget)
        else:
            margin = need-int(dp[min(odd[1],budget)])
            margins.append(margin)
            if margin<=0:
                survivors.append((x,y))
    return dict(excluded=not survivors,frontier=frontier,surviving_lower_pairs=survivors,
                minimum_allocation_deficit=min(margins),odd_capacity=odd,even_capacity_items=items,
                scope='integer plane majorization plus exact two-parity line/point allocation')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/core-geometry-profile-audit.json')
    ap.add_argument('--core-pattern',help='optionally restrict to one of the four core patterns')
    ap.add_argument('--high-pattern',help='optionally restrict to one high-prime pattern')
    ap.add_argument('--output')
    ap.add_argument('--refine-survivors',action='store_true')
    ap.add_argument('--base',default='data/plane-capacity-510510.json')
    a = ap.parse_args()
    if a.output is None:
        a.output = 'data/refined-plane-capacity-510510.json' if a.refine_survivors else 'data/plane-capacity-510510.json'
    results = []
    sources = json.loads(Path(a.source).read_text())
    if a.refine_survivors:
        for row in json.loads(Path(a.base).read_text()):
            profiles = next(x['easiest_profiles'] for x in sources if x['pattern']==row['core_pattern'])
            for i in row['surviving_indices']:
                result = refined_bound(profiles[i],(11,13,17),row['high_pattern'])
                result.update(core_pattern=row['core_pattern'],high_pattern=row['high_pattern'],profile_index=i)
                results.append(result)
                print(json.dumps({k:v for k,v in result.items() if k not in ('frontier','surviving_lower_pairs','even_capacity_items')}),flush=True)
        Path(a.output).write_text(json.dumps(results,indent=2)+'\n')
        return
    for row in sources:
        if a.core_pattern and row['pattern']!=list(map(int,a.core_pattern)):
            continue
        for pattern in itertools.product((0,1),repeat=3):
            if a.high_pattern and tuple(map(int,a.high_pattern))!=pattern:
                continue
            start = time.monotonic()
            tests = [bound(p,(11,13,17),pattern) for p in row['easiest_profiles']]
            result = dict(core_pattern=row['pattern'],high_pattern=pattern,profiles=len(tests),
                          excluded=sum(t['excluded'] for t in tests),minimum_gap=min(t['gap'] for t in tests),
                          surviving_indices=[i for i,t in enumerate(tests) if not t['excluded']],
                          tests=tests,seconds=round(time.monotonic()-start,3))
            results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('tests','surviving_indices')}),flush=True)
            Path(a.output).write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    main()
