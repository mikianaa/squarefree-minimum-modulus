"""Exact size-only cylinder capacity in two anchor boxes (at most four axes).

Odd cylinders may hit both layers; their even doubles choose one. Coordinate
compression makes zero the maximizing label in both layers. For each assignment
of even shapes, selecting the largest point weights independently in each box
therefore gives the sharp capacity for arbitrary subsets with specified sizes.
This does not assert that such subsets are residuals of congruence families.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np


def weight_tables(dimensions):
    k = len(dimensions[0])
    assert len(dimensions)==2 and len(dimensions[1])==k and k<=4
    assert all(d>=1 for dims in dimensions for d in dims)
    masks = list(range(1,1<<k))
    assignments = np.arange(1<<len(masks),dtype=np.int64)
    first = np.zeros((len(assignments),1<<k),dtype=np.int64)
    base = np.array([(1<<z.bit_count())-1 for z in range(1<<k)],dtype=np.int64)
    for i,s in enumerate(masks):
        included = np.array([int(s&z==s) for z in range(1<<k)])
        first += ((assignments>>i)&1)[:,None]*included[None,:]
    tables = []
    for e in (0,1):
        weights = base[None,:]+(first if e==0 else base[None,:]-first)
        counts = np.array([math.prod(dimensions[e][i]-1 for i in range(k) if not z>>i&1)
                           for z in range(1<<k)],dtype=np.int64)
        order = np.argsort(-weights,axis=1,kind='stable')
        tables.append((np.take_along_axis(weights,order,axis=1),counts[order]))
    return tables


def top_sums(table,limit):
    weights,counts = table
    assert 0<=limit<=int(counts[0].sum())
    sizes = np.arange(limit+1,dtype=np.int64)[None,:]
    total = np.zeros((len(weights),limit+1),dtype=np.int64)
    used = np.zeros(len(weights),dtype=np.int64)
    for j in range(weights.shape[1]):
        take = np.minimum(np.maximum(sizes-used[:,None],0),counts[:,j,None])
        total += take*weights[:,j,None]
        used += counts[:,j]
    return total


def capacity(dimensions,counts):
    tables = weight_tables(dimensions)
    return int(np.max(sum(top_sums(t,r)[:,r] for t,r in zip(tables,counts))))


def extend_total(dimensions,lower,q,new_parity,tables=None):
    """Valid total lower bound from r0+r1>=lower in the indicated boxes.

For any fixed assignment and layer, f(r)=D*r-top_weight_sum(r) is convex
with f(0)=0. Thus max(0,f(r)) is nondecreasing. It suffices to minimize on
the boundary r0+r1=lower. Each layer's deficit is clipped separately.
"""
    if tables is None:
        tables = weight_tables(dimensions)
    limits = [min(lower,math.prod(d)) for d in dimensions]
    deficits = [np.maximum(0,(q-1-(new_parity==e))*np.arange(limits[e]+1)[None,:]
                           -top_sums(tables[e],limits[e])) for e in (0,1)]
    best = None
    minimizing = None
    for r0 in range(max(0,lower-limits[1]),limits[0]+1):
        values = deficits[0][:,r0]+deficits[1][:,lower-r0]
        choice = int(np.argmin(values))
        v = int(values[choice])
        if best is None or v<best:
            best,minimizing = v,[r0,lower-r0,choice]
    return dict(uncovered_lower=best,minimizing_sizes_and_assignment=minimizing)


def extend_frontier(dimensions,frontier,q,new_parity,tables=None):
    """Propagate componentwise lower pairs, retaining the output Pareto frontier."""
    if tables is None:
        tables = weight_tables(dimensions)
    limits = [max(r[e] for r in frontier) for e in (0,1)]
    deficits = [np.maximum(0,(q-1-(new_parity==e))*np.arange(limits[e]+1)[None,:]
                           -top_sums(tables[e],limits[e])) for e in (0,1)]
    sentinel = 10**15
    minima = np.full(int(deficits[0].max())+1,sentinel,dtype=np.int64)
    for a,b in frontier:
        np.minimum.at(minima,deficits[0][:,a],deficits[1][:,b])
    result = []
    smallest = sentinel
    for a,b in enumerate(minima):
        if b<smallest:
            result.append([a,int(b)])
            smallest = b
    return dict(uncovered_lower=min(map(sum,result)),frontier=result)


def require_total(frontier,lower):
    pairs = set()
    for a,b in frontier:
        for extra0 in range(max(0,lower-a-b)+1):
            pairs.add((a+extra0,max(b,lower-a-extra0)))
    best = math.inf
    result = []
    for a,b in sorted(pairs):
        if b<best:
            result.append((a,b))
            best = b
    return result


def run_frontiers():
    optima = {tuple(x['pattern']):x['uncovered'] for x in
              json.loads(Path('data/min-uncovered-2310-five-patterns-optimal.json').read_text())}
    base = json.loads(Path('data/three-axis-210-certificate.json').read_text())['cases']
    results = []
    first = []
    for row in base:
        p = row['pattern']
        dims = [[v-1-(b==e) for v,b in zip((3,5,7),p)] for e in (0,1)]
        tables = weight_tables(dims)
        for b in (0,1):
            out = dict(pattern=p+[b],**extend_frontier(dims,row['residual_pair_frontier'],11,b,tables))
            first.append(out)
            print(json.dumps({k:v for k,v in out.items() if k!='frontier'}),flush=True)
    for row in first:
        p = row['pattern']
        dims = [[v-1-(b==e) for v,b in zip((3,5,7,11),p)] for e in (0,1)]
        tables = weight_tables(dims)
        for status in ('integer_only','including_solver'):
            front = row['frontier']
            if status=='including_solver':
                front = require_total(front,optima.get(tuple(p),0))
            for b in (0,1):
                out = dict(pattern=p+[b],status=status,**extend_frontier(dims,front,13,b,tables))
                results.append(out)
                print(json.dumps({k:v for k,v in out.items() if k!='frontier'}),flush=True)
    Path('data/compressed-parity-frontiers.json').write_text(json.dumps(dict(at2310=first,at30030=results),indent=2)+'\n')


def run():
    rows = json.loads(Path('data/six-prime-obstruction.json').read_text())
    optima = {tuple(x['pattern']):x['uncovered'] for x in
              json.loads(Path('data/min-uncovered-2310-five-patterns-optimal.json').read_text())}
    results = []
    for row in rows:
        pattern = row['pattern']
        dimensions = [[p-1-(b==e) for p,b in zip((3,5,7,11),pattern)] for e in (0,1)]
        tables = weight_tables(dimensions)
        for status,lower in [('integer_only',row['uncovered_2310_lower']),
                             ('including_solver',max(row['uncovered_2310_lower'],optima.get(tuple(pattern),0)))]:
            for b in (0,1):
                result = dict(pattern=pattern+[b],input_lower=lower,status=status,
                              **extend_total(dimensions,lower,13,b,tables))
                results.append(result)
                print(json.dumps(result),flush=True)
    Path('data/compressed-parity-30030.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    import sys
    run_frontiers() if '--frontiers' in sys.argv else run()
