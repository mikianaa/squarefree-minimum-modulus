"""Apply aggregate conditional-tail duals to arbitrary normalized core masks.

The dual point inequalities depend only on the four higher primes and their
anchor parities. Reuse them for other cores by recomputing valid concentration
right-hand sides. Here H is bounded by joint parity slice counts and by the
sharp compression bound for subsets of the anchor boxes. For pair residuals,
remove up to one point in each parity and a second in either parity to allow
for the final classes 105 and 210. Artificial extra removals are permitted:
any actual final residual contains a subset of one of the two target sizes.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from core_residual_geometry import points_for
from compressed_parity_capacity import weight_tables,top_sums


class CoreCapacityUpper:
    def __init__(self,pattern):
        self.points = points_for(pattern)
        self.even = sum(1<<j for j,(e,a) in enumerate(self.points) if e==0)
        crt = {(n%2,(n%3,n%5,n%7)):n for n in range(210)}
        self.numbers = [crt[e,tuple(a)] for e,a in self.points]
        self.slices = [[sum(1<<j for j,n in enumerate(self.numbers) if n%d==a) for a in range(d)]
                       for d in range(3,211) if 210%d==0]
        dimensions = [[p-1-(b==e) for p,b in zip((3,5,7),pattern)] for e in (0,1)]
        sizes = [sum(e==parity for e,a in self.points) for parity in (0,1)]
        a,b = [top_sums(t,n) for t,n in zip(weight_tables(dimensions),sizes)]
        self.compression = np.max(a[:,:,None]+b[:,None,:],axis=0)

    def upper(self,mask):
        r0 = (mask&self.even).bit_count()
        r1 = mask.bit_count()-r0
        a = np.arange(r0+1)[:,None]
        b = np.arange(r1+1)[None,:]
        value = np.zeros((r0+1,r1+1),dtype=np.int64)
        for slices in self.slices:
            capacities = set()
            for s in slices:
                c0 = (mask&s&self.even).bit_count()
                c1 = (mask&s).bit_count()-c0
                capacities.add((c0,c1))
            best = np.zeros_like(value)
            for c0,c1 in capacities:
                best = np.maximum(best,np.minimum(a,c0)+np.minimum(b,c1))
            value += best
        return np.minimum(value,self.compression[:r0+1,:r1+1])


def concentration_bounds(H,shapes=14):
    a = np.arange(H.shape[0])
    b = np.arange(H.shape[1])
    if shapes==16:
        H = H+a[:,None]+b[None,:]+np.maximum(a[:,None],b[None,:])
    assert shapes in (14,16)
    thresholds = np.arange(17)
    first = np.max(H[None,:,:]-thresholds[:,None,None]*a[None,:,None],axis=1)
    return np.max(first[:,None,:]-thresholds[None,:,None]*b[None,None,:],axis=2)


class PortableDuals:
    def __init__(self,data):
        self.patterns = [r['high_pattern'] for r in data['cases']]
        self.points = []
        self.weights = []
        self.scales = []
        for row in data['cases']:
            c = row['certificate']
            assert row['aggregate'] and not c['resource_multipliers']
            assert all(v==0 for line in c['link_multipliers'] for v in line)
            weights = np.zeros((2,17,17),dtype=np.int64)
            for item in c['tail_multipliers']:
                t0,t1 = item['thresholds']
                assert item['multiplier']<=0
                weights[int(item['axis']>=4),t0,t1] += item['multiplier']
            self.weights.append(weights)
            self.points.append(c['point_multipliers'])
            self.scales.append(c['denominator'])
        self.weights = np.array(self.weights)
        self.points = np.array(self.points,dtype=np.int64)
        self.scales = np.array(self.scales,dtype=np.int64)

    def bounds(self,H):
        counts = np.array([s-1 for s in H.shape],dtype=np.int64)
        K = np.stack([concentration_bounds(H,s) for s in (14,16)])
        numerator = self.points@counts+np.sum(self.weights*K[None,:,:,:],axis=(1,2,3))
        return -((-numerator)//self.scales)


def run(source):
    start = time.monotonic()
    duals = PortableDuals(json.loads(Path('data/conditional-tail-core-pattern-certificates.json').read_text()))
    catalogue = json.loads(Path(source).read_text())
    results = []
    for row in catalogue['cases']:
        model = CoreCapacityUpper(row['pattern'])
        masks = row['minimal_masks']
        values = []
        for i,mask in enumerate(masks):
            H = model.upper(mask)
            r0,r1 = [s-1 for s in H.shape]
            branches = [duals.bounds(H[:max(0,r0-1-(e==0))+1,:max(0,r1-1-(e==1))+1]) for e in (0,1)]
            values.append(np.minimum(*branches).tolist())
        counts = np.sum(np.array(values)>0,axis=0).tolist()
        out = dict(core_pattern=row['pattern'],pair_masks=masks,high_patterns=duals.patterns,
                   lower_bounds=values,excluded_pair_orbits_per_high_pattern=counts)
        results.append(out)
        print(json.dumps(dict(pattern=row['pattern'],pair_orbits=len(masks),excluded=counts,seconds=round(time.monotonic()-start,3))),flush=True)
    return dict(cases=results,seconds=round(time.monotonic()-start,3),
                scope='portable dual bounds with optimistic triple deletions; nonpositive values are inconclusive')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/core-orbit-reduction.json')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    Path(a.output).write_text(json.dumps(run(a.source),indent=2)+'\n')
