"""Separate impossible simultaneous threshold events on three core points.

For a point set T, sum_{x in T} z(x)<=sum_d max_a |T intersect (a mod d)|.
If the proposed thresholds sum to more than this capacity, not all events
can hold, giving sum_{x in T} 1[z(x)>=t_x] <= |T|-1.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import time
import numpy as np


def separate(data,subset_size=3):
    start = time.monotonic()
    source = data['input']
    witness = data['result']['point_primal_witness']
    R = source['core_residual']
    n = len(source['high_primes'])
    groups = (1<<n)-1
    hist = np.zeros((len(R),groups,17),dtype=np.int64)
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[row['point_index'],h,v] += row['weight']
    tails = np.cumsum(hist[:,:,::-1],axis=2)[:,:,::-1]
    core = [d for d in range(3,211) if 210%d==0]
    violations = []
    checks = 0
    for kind,shapes in ((0,core),(1,[1,2]+core)):
        limit = len(shapes)
        t = np.arange(1,limit+1)[:,None]
        u = np.arange(1,limit+1)[None,:]
        for points in itertools.combinations(range(len(R)),subset_size):
            cap = sum(max(Counter(R[x]%d for x in points).values()) for d in shapes)
            if subset_size!=3:
                axes = list(range(n) if kind==0 else range(n,groups))
                width = subset_size*limit+1
                dp = np.full((len(axes),width),-10**18,dtype=np.int64)
                dp[:,0] = 0
                choices = []
                for x in points:
                    next_dp = np.full_like(dp,-10**18)
                    choice = np.zeros(dp.shape,dtype=np.int8)
                    for threshold in range(1,limit+1):
                        candidate = dp[:,:-threshold]+tails[x,axes,threshold][:,None]
                        improved = candidate>next_dp[:,threshold:]
                        np.copyto(next_dp[:,threshold:],candidate,where=improved)
                        np.copyto(choice[:,threshold:],threshold,where=improved)
                    dp = next_dp
                    choices.append(choice)
                for index,h in enumerate(axes):
                    checks += 1
                    total = cap+1+int(dp[index,cap+1:].argmax())
                    gap = int(dp[index,total])-(subset_size-1)*witness['denominator']
                    if gap<=0:
                        continue
                    thresholds = []
                    for choice in reversed(choices):
                        threshold = int(choice[index,total])
                        thresholds.append(threshold)
                        total -= threshold
                    thresholds.reverse()
                    assert total==0 and sum(thresholds)>cap
                    used = sum(row['weight'] for row in witness['support']
                               for point,threshold in zip(points,thresholds)
                               if row['point_index']==point and row['counts'][h]>=threshold)
                    assert used-(subset_size-1)*witness['denominator']==gap
                    violations.append(dict(axis=h,points=list(points),thresholds=thresholds,
                                           projection_sum_capacity=cap,violation_numerator=gap))
                continue
            v = np.maximum(1,cap+1-t-u)
            allowed = v<=limit
            clipped = np.minimum(v,limit)
            x,y,z = points
            for h in (range(n) if kind==0 else range(n,groups)):
                scores = tails[x,h,t]+tails[y,h,u]+tails[z,h,clipped]
                scores = np.where(allowed,scores,-1)
                i,j = np.unravel_index(scores.argmax(),scores.shape)
                gap = int(scores[i,j])-2*witness['denominator']
                checks += 1
                if gap>0:
                    thresholds = [int(t[i,0]),int(u[0,j]),int(v[i,j])]
                    # Independently repeat this one proposed violation with
                    # Python integers and direct support sums.
                    assert sum(thresholds)>cap
                    used = sum(row['weight'] for row in witness['support']
                               for point,threshold in zip(points,thresholds)
                               if row['point_index']==point and row['counts'][h]>=threshold)
                    assert used-2*witness['denominator']==gap
                    violations.append(dict(axis=h,points=list(points),thresholds=thresholds,
                                           projection_sum_capacity=cap,violation_numerator=gap))
    violations.sort(key=lambda row:row['violation_numerator'],reverse=True)
    return dict(subset_size=subset_size,subset_group_combinations_checked=checks,violations=violations,
                denominator=witness['denominator'],seconds=round(time.monotonic()-start,3),
                scope='valid cuts excluding the supplied fractional witness; not a covering obstruction')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--subset-size',type=int,default=3)
    args = ap.parse_args()
    result = separate(json.loads(Path(args.source).read_text()),args.subset_size)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='violations'}|dict(violations=len(result['violations']),largest=result['violations'][:3])))
