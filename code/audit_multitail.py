"""Independent Python-integer checks of multi-prime tail LP witnesses.

Capacity validity is a separate premise. A feasible negative-cost primal is
only a limitation of this relaxation, never a covering construction.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def local_minima(source,certificate,parity_only=None):
    high,pattern = source['high_primes'],source['high_pattern']
    n = len(high)
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    complements = [[i for i in range(n) if i not in S] for S in groups[n:]]
    scale = certificate['denominator']
    assert type(scale) is int and scale>0
    common = [[[0]*(15 if h<n else 17) for h in range(len(groups))] for e in (0,1)]
    for row in certificate['tail_multipliers']:
        h,y = row['axis'],row['multiplier']
        assert type(h) is int and 0<=h<len(groups) and type(y) is int and y<=0
        assert len(row['thresholds'])==2
        for e,t in enumerate(row['thresholds']):
            assert type(t) is int and 0<=t<len(common[e][h])
            for v in range(len(common[e][h])):
                common[e][h][v] += y*max(0,v-t)
    for row in certificate.get('rank_multipliers',[]):
        h,y = row['axis'],row['multiplier']
        assert type(h) is int and 0<=h<len(groups) and type(y) is int and y<=0
        assert len(row['thresholds'])==len(row['weights'])==2
        for e,(t,w) in enumerate(zip(row['thresholds'],row['weights'])):
            assert type(t) is int and 1<=t<len(common[e][h])
            assert type(w) is int and w>=0
            for v in range(t,len(common[e][h])):
                common[e][h][v] += y*w
    result,checked = [],0
    for e in ((0,1) if parity_only is None else (parity_only,)):
        D = [q-1-(b==e) for q,b in zip(high,pattern)]
        tables = [[min(-scale*v*z-common[e][h][z] for z in range(17))
                   for v in range(math.prod(D[i] for i in indices)+1)]
                  for h,indices in enumerate(complements,n)]
        minimum = None
        for counts in itertools.product(*(range(min(14,d)+1) for d in D)):
            f = [d-v for d,v in zip(D,counts)]
            value = scale*math.prod(f)-sum(common[e][i][v] for i,v in enumerate(counts))
            value += sum(table[math.prod(f[i] for i in indices)] for table,indices in zip(tables,complements))
            minimum = value if minimum is None else min(minimum,value)
            checked += 1
        result.append(minimum)
    return result,checked


def audit(data,primal=None):
    source,result = data['input'],data['result']
    high,pattern,counts,H = [source[k] for k in ('high_primes','high_pattern','counts','H_upper')]
    n = len(high)
    assert n>0 and len(pattern)==n and all(b in (0,1) for b in pattern)
    assert len(counts)==2 and all(type(v) is int and v>=0 for v in counts)
    assert len(H)==counts[0]+1 and all(len(row)==counts[1]+1 for row in H)
    assert all(type(v) is int and v>=0 for row in H for v in row)
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    limits = [14 if h<n else 16 for h in range(len(groups))]
    primal = bool(result.get('primal_relaxation_witness')) if primal is None else primal
    if primal:
        witness = result['primal_relaxation_witness']
        scale = witness['denominator']
        assert type(scale) is int and scale>0
        hist = [[[0]*(limit+1) for limit in limits] for e in (0,1)]
        masses,objective = [0,0],0
        for row in witness['support']:
            e,state,w = row['parity'],row['counts'],row['weight']
            assert e in (0,1) and type(w) is int and w>0 and len(state)==len(groups)
            D = [q-1-(b==e) for q,b in zip(high,pattern)]
            for h,v in enumerate(state):
                assert type(v) is int and 0<=v<=(min(14,D[h]) if h<n else 16)
                hist[e][h][v] += w
            masses[e] += w
            f = [d-v for d,v in zip(D,state[:n])]
            cost = math.prod(f)-sum(state[h]*math.prod(f[i] for i in range(n) if i not in S)
                                    for h,S in enumerate(groups[n:],n))
            objective += w*cost
        assert masses==[v*scale for v in counts]
        checked = 0
        for h,limit in enumerate(limits):
            for t0,t1 in itertools.product(range(limit+1),repeat=2):
                used = sum(w*max(0,v-t) for e,t in enumerate((t0,t1)) for v,w in enumerate(hist[e][h]))
                bound = max(v+(a+b+max(a,b) if h>=n else 0)-t0*a-t1*b for a,row in enumerate(H) for b,v in enumerate(row))
                assert used<=scale*bound,(h,t0,t1,used,bound)
                checked += 1
        rank_checked = 0
        if result.get('rank_cuts_enabled'):
            # Reconstruct all bounds from integer superlevel-set sizes,
            # independently of the generator's NumPy catalogue.
            for h,limit in enumerate(limits):
                for t0,t1 in itertools.product(range(1,limit+1),repeat=2):
                    feasible = [(a,b) for a,line in enumerate(H) for b,v in enumerate(line)
                                if t0*a+t1*b<=v+(a+b+max(a,b) if h>=n else 0)]
                    for w0,w1 in ((1,0),(0,1),(1,1),(2,1),(1,2)):
                        bound = max(w0*a+w1*b for a,b in feasible)
                        used = w0*sum(hist[0][h][t0:])+w1*sum(hist[1][h][t1:])
                        assert used<=scale*bound
                        rank_checked += 1
        assert objective==witness['objective_numerator'] and objective<0
        return dict(verified=True,high_primes=high,high_pattern=pattern,counts=counts,
                    objective_numerator=objective,denominator=scale,support_size=len(witness['support']),
                    concentration_inequalities_checked=checked,rank_inequalities_checked=rank_checked,
                    scope='feasible negative-cost aggregate-relaxation distribution; not an actual cover')
    certificate = result['certificate']
    points,checked = local_minima(source,certificate)
    assert points==certificate['point_multipliers']
    numerator = sum(v*p for v,p in zip(counts,points))
    for row in certificate['tail_multipliers']:
        t0,t1 = row['thresholds']
        bound = max(v+(a+b+max(a,b) if row['axis']>=n else 0)-t0*a-t1*b for a,line in enumerate(H) for b,v in enumerate(line))
        assert bound==row['bound']
        numerator += row['multiplier']*bound
    for row in certificate.get('rank_multipliers',[]):
        t0,t1 = row['thresholds']
        w0,w1 = row['weights']
        feasible = [(a,b) for a,line in enumerate(H) for b,v in enumerate(line)
                    if t0*a+t1*b<=v+(a+b+max(a,b) if row['axis']>=n else 0)]
        bound = max(w0*a+w1*b for a,b in feasible)
        assert bound==row['bound']
        numerator += row['multiplier']*bound
    assert numerator==certificate['numerator']
    lower = -(-numerator//certificate['denominator'])
    assert lower==certificate['integer_lower']
    return dict(verified=True,high_primes=high,high_pattern=pattern,counts=counts,
                integer_lower=lower,positive_obstruction=lower>0,plane_tuples_checked=checked,
                scope='dual bound conditional on separately verified input capacity; not a whole-branch result')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--dual',action='store_true')
    args = ap.parse_args()
    result = audit(json.loads(Path(args.source).read_text()),False if args.dual else None)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
