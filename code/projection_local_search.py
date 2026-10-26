"""Construct assignments in the high-coordinate capacity relaxation.

Even a zero-deficit assignment is NOT an actual covering system: separate
classes may each be credited with incompatible or overlapping core capacities.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np


def solve(high,pattern,residual,caps,seconds=30,seed=1):
    rng = np.random.default_rng(seed)
    alphabets = [[tuple(v for v in range(1,p) if v!=1 or pattern[i+3]!=e) for i,p in enumerate(high)] for e in (0,1)]
    cells = [(e,a) for e in (0,1) for a in itertools.product(*alphabets[e])]
    parity = np.array([e for e,a in cells])
    coordinates = np.array([a for e,a in cells])
    target = np.array([residual[e] for e,a in cells],dtype=np.int64)
    resources = []
    for k in range(1,len(high)+1):
        for subset in itertools.combinations(range(len(high)),k):
            size = math.prod(high[i]-1 for i in subset)
            code = np.zeros(len(cells),dtype=np.int64)
            for i in subset:
                code = code*(high[i]-1)+coordinates[:,i]-1
            for d in sorted(caps):
                if k==1 and d in (1,2):
                    continue
                cap = np.array([caps[d][e] for e,a in cells],dtype=np.int64)
                if not cap.any():
                    continue
                codes = code+parity*size if d%2==0 else code
                allowed = np.unique(codes[cap>0])
                resources.append(dict(modulus=d*math.prod(high[i] for i in subset),core_divisor=d,
                                      subset=subset,size=size,even=d%2==0,codes=codes,cap=cap,allowed=allowed))
    resources.sort(key=lambda r:r['modulus'])
    load = np.zeros(len(cells),dtype=np.int64)
    chosen = [-1]*len(resources)
    best = int(target.sum())
    best_choices = None
    deadline = time.monotonic()+seconds
    passes = 0
    while passes==0 or time.monotonic()<deadline:
        order = rng.permutation(len(resources)) if passes else np.arange(len(resources))
        for i in order:
            r = resources[i]
            if chosen[i]>=0:
                match = r['codes']==chosen[i]
                load[match] -= r['cap'][match]
            missing = np.maximum(0,target-load)
            gains = np.bincount(r['codes'],weights=np.minimum(missing,r['cap']),minlength=r['size']*(2 if r['even'] else 1))
            maximum = gains[r['allowed']].max()
            choice = int(rng.choice(r['allowed'][gains[r['allowed']]==maximum]))
            chosen[i] = choice
            match = r['codes']==choice
            load[match] += r['cap'][match]
        passes += 1
        deficit = int(np.maximum(0,target-load).sum())
        if deficit<best:
            best,best_choices = deficit,list(chosen)
            print(json.dumps(dict(pass_number=passes,relaxation_deficit=best)),flush=True)
        if best==0:
            break
        if passes%7==0:
            chosen = list(best_choices)
            load[:] = 0
            for r,choice in zip(resources,chosen):
                match = r['codes']==choice
                load[match] += r['cap'][match]
        pool = min(40,len(resources)) if passes%3==0 else len(resources)
        for i in rng.choice(pool,size=min(4,pool),replace=False):
            r = resources[i]
            match = r['codes']==chosen[i]
            load[match] -= r['cap'][match]
            chosen[i] = int(rng.choice(r['allowed']))
            match = r['codes']==chosen[i]
            load[match] += r['cap'][match]
    placements = []
    for r,choice in zip(resources,best_choices):
        e = choice//r['size'] if r['even'] else None
        code = choice%r['size']
        labels = []
        for i in reversed(r['subset']):
            labels.append(code%(high[i]-1)+1)
            code //= high[i]-1
        placements.append(dict(modulus=r['modulus'],core_divisor=r['core_divisor'],subset=r['subset'],
                               parity=e,labels=list(reversed(labels)),capacity=caps[r['core_divisor']]))
    # An independent pure-integer pass, without NumPy scores or stored loads.
    checked = sum(max(0,residual[e]-sum(v['capacity'][e] for v in placements
                    if (v['parity'] is None or v['parity']==e)
                    and all(a[i]==b for i,b in zip(v['subset'],v['labels'])))) for e,a in cells)
    assert checked==best
    return dict(high=high,pattern=pattern,residual=residual,relaxation_deficit=best,
                passes=passes,seed=seed,placements=placements,verified_cells=len(cells),
                scope='capacity relaxation only; zero deficit would not be an actual covering system')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/core-actual-profile-510510-unknown.json')
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=30)
    ap.add_argument('--seed',type=int,default=1)
    args = ap.parse_args()
    source = json.loads(Path(args.source).read_text())
    family = source['core_family']
    residuals = [[n for n in range(e,210,2) if not any(n%d==a for d,a in family)] for e in (0,1)]
    caps = {d:[max(sum(n%d==a for n in r) for a in range(d)) for r in residuals] for d in range(1,211) if 210%d==0}
    result = solve([11,13,17],source['pattern'],list(map(len,residuals)),caps,args.time,args.seed)
    result.update(core_family=family,core_capacities=caps)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('placements','core_family','core_capacities')}),flush=True)
