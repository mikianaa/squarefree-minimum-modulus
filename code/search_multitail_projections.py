"""Search integer core-projection choices for the conditional union bound.

For fixed plane choices, optimize every tail core projection independently.
A negative expression would limit this particular bound, not construct a cover.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np


def search(source,seconds=60,seed=1,identity_witness=None):
    start = time.monotonic()
    rng = np.random.default_rng(seed)
    R = source['core_residual']
    high,pattern = source['high_primes'],source['high_pattern']
    n = len(high)
    D = np.array([[q-1-(b==x%2) for x in R] for q,b in zip(high,pattern)],dtype=np.int64)
    core = [d for d in range(3,211) if 210%d==0]
    tails = [S for k in range(2,n+1) for S in itertools.combinations(range(n),k)]
    moduli = [1,2]+core
    residues = {d:sorted({x%d for x in R}) for d in moduli}
    matrices = {d:np.array([[int(x%d==a) for x in R] for a in residues[d]],dtype=np.int64) for d in moduli}
    offsets = np.cumsum([0]+[len(residues[d]) for d in moduli])
    matrix = np.concatenate([matrices[d] for d in moduli])
    complements = [[i for i in range(n) if i not in S] for S in tails]
    resources = [(i,d) for i in range(n) for d in core]
    def values(free):
        # free: candidates x axes x core points
        volume = np.stack([np.prod(free[:,indices,:],axis=1) for indices in complements],axis=2)
        scores = np.einsum('ar,crs->cas',matrix,volume,optimize=True)
        maximum = np.maximum.reduceat(scores,offsets[:-1],axis=1)
        return np.prod(free,axis=1).sum(axis=1)-maximum.sum(axis=(1,2)),scores
    choices = [int(rng.integers(len(residues[d]))) for i,d in resources]
    if identity_witness:
        by_key = {(r['axis'],r['modulus']):r for r in identity_witness['projections']}
        for j,(i,d) in enumerate(resources):
            a = max(by_key[i,d]['assignments'],key=lambda r:r['weight'])['residue']
            choices[j] = residues[d].index(a)
    counts = np.zeros_like(D)
    for (i,d),j in zip(resources,choices):
        counts[i] += matrices[d][j]
    best,best_choices = None,None
    passes = 0
    while passes==0 or time.monotonic()-start<seconds:
        improved = False
        for j in rng.permutation(len(resources)):
            i,d = resources[j]
            counts[i] -= matrices[d][choices[j]]
            free = np.broadcast_to(np.maximum(0,D-counts),(len(residues[d]),n,len(R))).copy()
            free[:,i,:] = np.maximum(0,D[i]-counts[i]-matrices[d])
            costs,_ = values(free)
            candidate = int(rng.choice(np.flatnonzero(costs==costs.min())))
            improved |= int(costs[candidate])<int(costs[choices[j]])
            choices[j] = candidate
            counts[i] += matrices[d][candidate]
        cost = int(values(np.maximum(0,D-counts)[None,:,:])[0][0])
        passes += 1
        if best is None or cost<best:
            best,best_choices = cost,list(choices)
            print(json.dumps(dict(passes=passes,integer_expression=best,seconds=round(time.monotonic()-start,3))),flush=True)
        if best<0:
            break
        if not improved:
            # Perturb several choices, then resume exact coordinate descent.
            for j in rng.choice(len(resources),size=min(10,len(resources)),replace=False):
                i,d = resources[j]
                counts[i] -= matrices[d][choices[j]]
                choices[j] = int(rng.integers(len(residues[d])))
                counts[i] += matrices[d][choices[j]]
    counts = np.zeros_like(D)
    choices_out = []
    for (i,d),j in zip(resources,best_choices):
        counts[i] += matrices[d][j]
        choices_out.append(dict(axis=i,modulus=d,residue=residues[d][j]))
    costs,scores = values(np.maximum(0,D-counts)[None,:,:])
    assert int(costs[0])==best
    for h in range(len(tails)):
        for k,d in enumerate(moduli):
            j = int(scores[0,offsets[k]:offsets[k+1],h].argmax())
            choices_out.append(dict(axis=n+h,modulus=d,residue=residues[d][j]))
    return dict(integer_expression=best,choices=choices_out,passes=passes,seconds=round(time.monotonic()-start,3),
                scope='integer core projections for a conditional union-bound expression; higher labels unspecified, not an actual cover')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=60)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    result = search(data['input'],args.time,identity_witness=data['result'].get('point_primal_witness'))
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='choices'}),flush=True)
