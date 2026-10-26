"""Quantitative parity propagation retaining audited 210-core intersections.

Uses ALL raw optimistic profiles, never normalized feasibility dominance.
For fixed even-resource allocation and layer, D*r-sum min(r,C_i) is convex
and zero at r=0; its positive part is nondecreasing. Thus lower core counts
and original pair-residual capacities give valid absolute deficit bounds.
"""
import itertools
import json
from pathlib import Path
import time
import numpy as np
from core_residual_geometry import pair_residuals,profile
from compressed_parity_capacity import weight_tables,extend_frontier,require_total
from parity_capacity import pareto


def intersect_frontiers(first,second):
    return [list(p) for p in pareto((max(a,c),max(b,d)) for a,b in first for c,d in second)]


def enhanced_profiles(points,masks):
    """Retain shared odd-3q choices as well as marginal class capacities."""
    assert all(a[0]==2 for e,a in points if e==0)
    types = [sum(1<<j for j,(e,a) in enumerate(points) if e==1 and a[0]==v) for v in (1,2)]
    rows = set()
    for mask in masks:
        p = profile(mask,points)
        counts3 = [(mask&t).bit_count() for t in types]
        for extra in (0,1):
            r = [max(0,p[e]-1-(e==extra)) for e in (0,1)]
            caps = [min(r[e],p[2+2*j+e]) for j in range(7) for e in (0,1)]
            rows.add(tuple(r+caps+[min(r[1],v) for v in counts3]))
    return sorted(rows)


def raw_profile_frontier(profiles,q,b):
    ar = np.array(profiles,dtype=np.int64)
    sizes = ar[:,:2]
    caps = ar[:,2:16].reshape(-1,7,2)
    pairs = set()
    for choice in (0,1):
        odd3 = np.stack((np.zeros(len(sizes),dtype=np.int64) if choice==0 else sizes[:,0],
                         ar[:,16+choice]),axis=1)
        odd = caps.sum(axis=1)-caps[:,0,:]+odd3
        initial = sizes*np.array([q-1-(b==e) for e in (0,1)])[None,:]-odd
        for bits in itertools.product((0,1),repeat=7):
            allocated = np.array(bits)[None,:,None]==np.arange(2)[None,None,:]
            deficit = np.maximum(0,initial-(caps*allocated).sum(axis=1))
            pairs.update(map(tuple,deficit.tolist()))
    return [list(p) for p in pareto(pairs)]


def run():
    start = time.monotonic()
    previous = json.loads(Path('data/compressed-parity-frontiers.json').read_text())
    optima = {tuple(r['pattern']):r['uncovered'] for r in
              json.loads(Path('data/min-uncovered-2310-five-patterns-optimal.json').read_text())}
    first = []
    second = []
    for tail in itertools.product((0,1),repeat=2):
        p = (0,)+tail
        points,masks,_ = pair_residuals(p)
        profiles = enhanced_profiles(points,masks)
        for b in (0,1):
            pattern = list(p)+[b]
            raw = raw_profile_frontier(profiles,11,b)
            old = next(r for r in previous['at2310'] if r['pattern']==pattern)['frontier']
            combined = intersect_frontiers(raw,old)
            row = dict(pattern=pattern,enhanced_profile_count=len(profiles),raw_frontier=raw,
                       frontier=combined,uncovered_lower=min(map(sum,combined)))
            first.append(row)
            print(json.dumps({k:v for k,v in row.items() if 'frontier' not in k}),flush=True)
    for row in first:
        p = row['pattern']
        dims = [[v-1-(b==e) for v,b in zip((3,5,7,11),p)] for e in (0,1)]
        tables = weight_tables(dims)
        for status in ('integer_only','including_solver'):
            front = row['frontier']
            if status=='including_solver':
                front = require_total(front,optima.get(tuple(p),0))
            for b in (0,1):
                result = dict(pattern=p+[b],status=status,**extend_frontier(dims,front,13,b,tables))
                second.append(result)
                print(json.dumps({k:v for k,v in result.items() if k!='frontier'}),flush=True)
    result = dict(at2310=first,at30030=second,seconds=round(time.monotonic()-start,3),
                  scope='quantitative core capacities with shared odd-3q labels and compressed parity propagation; solver-dependent rows labeled separately')
    Path('data/core-capacity-frontiers.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    run()
