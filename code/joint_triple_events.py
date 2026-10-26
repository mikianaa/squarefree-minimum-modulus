"""Exact simultaneous threshold feasibility for three core points.

Every resource partitions three points into one, two, or three residue
classes. Enumerate how many two-class resources choose their singleton;
the three-class resources then fill the remaining point demands.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import time
import numpy as np


def partition_counts(points,shapes):
    common,separate,pairs = 0,0,[0,0,0]
    for d in shapes:
        labels = [x%d for x in points]
        if len(set(labels))==1:
            common += 1
        elif len(set(labels))==3:
            separate += 1
        else:
            i = next(i for i in range(3) if labels.count(labels[i])==1)
            pairs[i] += 1
    return common,separate,tuple(pairs)


def feasible(points,shapes,thresholds):
    common,separate,pairs = partition_counts(points,shapes)
    for singletons in itertools.product(*(range(v+1) for v in pairs)):
        hits = [common+singletons[j]+sum(pairs[i]-singletons[i] for i in range(3) if i!=j) for j in range(3)]
        if sum(max(0,t-v) for t,v in zip(thresholds,hits))<=separate:
            return True
    return False


def separate(data):
    start = time.monotonic()
    source,witness = data['input'],data['result']['point_primal_witness']
    R,n,scale = source['core_residual'],len(source['high_primes']),witness['denominator']
    groups = (1<<n)-1
    hist = np.zeros((len(R),groups,17),dtype=np.int64)
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[row['point_index'],h,v] += row['weight']
    tails = np.cumsum(hist[:,:,::-1],axis=2)[:,:,::-1]
    core = [d for d in range(3,211) if 210%d==0]
    cache,violations = {},[]
    checks = 0
    for kind,shapes in ((0,core),(1,[1,2]+core)):
        limit = len(shapes)
        thresholds = np.indices((limit,)*3)+1
        for selected in itertools.combinations(range(len(R)),3):
            points = [R[x] for x in selected]
            key = partition_counts(points,shapes)
            if key not in cache:
                common,singles,pairs = key
                possible = np.zeros((limit,)*3,dtype=bool)
                for choices in itertools.product(*(range(v+1) for v in pairs)):
                    hits = [common+choices[j]+sum(pairs[i]-choices[i] for i in range(3) if i!=j) for j in range(3)]
                    possible |= sum(np.maximum(0,thresholds[j]-hits[j]) for j in range(3))<=singles
                cache[key] = possible
            possible = cache[key]
            x,y,z = selected
            for h in (range(n) if kind==0 else range(n,groups)):
                scores = tails[x,h,1:limit+1,None,None]+tails[y,h,None,1:limit+1,None]+tails[z,h,None,None,1:limit+1]
                scores = np.where(possible,-1,scores)
                index = np.unravel_index(scores.argmax(),scores.shape)
                gap = int(scores[index])-2*scale
                checks += 1
                if gap>0:
                    ts = [int(i)+1 for i in index]
                    assert not feasible(points,shapes,ts)
                    used = sum(row['weight'] for row in witness['support'] for x,t in zip(selected,ts)
                               if row['point_index']==x and row['counts'][h]>=t)
                    assert used-2*scale==gap
                    cap = sum(max(Counter(x%d for x in points).values()) for d in shapes)
                    violations.append(dict(axis=h,points=list(selected),thresholds=ts,
                                           violation_numerator=gap,projection_sum_capacity=cap,
                                           detected_by_sum_capacity=sum(ts)>cap))
    violations.sort(key=lambda r:r['violation_numerator'],reverse=True)
    return dict(violations=violations,combinations_checked=checks,partition_profiles=len(cache),
                denominator=scale,seconds=time.monotonic()-start,
                scope='exact joint triple event cuts; not a covering-system obstruction')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    result = separate(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='violations'}|dict(violations=len(result['violations']),largest=result['violations'][:3])))
