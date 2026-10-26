"""Construct partial covering families by local improvement of a fixed prefix.

This supplies upper witnesses for the number of uncovered residues, never a
nonexistence proof. Every reported family is checked on its complete period.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import random
import time


def extend(family, period, q, parity, seconds=10, seed=1):
    assert math.gcd(period,q)==1 and all(period%d==0 and d>=3 for d,a in family)
    rng = random.Random(seed)
    residual = [n for n in range(period) if not any(n%d==a for d,a in family)]
    points = [(n,b) for b in range(1,q) for n in residual if b!=1 or n%2!=parity]
    whole = (1 << len(points))-1
    resources = []
    for d in range(3,period+1):
        if period%d:
            continue
        lookup = {}
        for j,(n,b) in enumerate(points):
            key = (n%d,b)
            lookup[key] = lookup.get(key,0) | 1 << j
        # Equal masks are interchangeable here; retain one legal residue.
        masks = {}
        for (a,b),mask in lookup.items():
            residue = a+d*((b-a)*pow(d,-1,q)%q)
            masks[mask] = residue
        resources.append((q*d,list(masks.items())))
    best = len(points)+1
    best_choices = None
    deadline = time.monotonic()+seconds
    restarts = 0
    while restarts==0 or time.monotonic()<deadline:
        chosen = [(0,0)]*len(resources)
        order = list(range(len(resources)))
        # A large-first initial pass, followed by random coordinate passes.
        if restarts%2:
            rng.shuffle(order)
        changed = True
        passes = 0
        while changed and passes<12:
            changed = False
            for i in order:
                others = 0
                for j,(mask,a) in enumerate(chosen):
                    if i!=j:
                        others |= mask
                missing = whole ^ others
                candidates = resources[i][1]
                gain = max((mask&missing).bit_count() for mask,a in candidates)
                ties = [(mask,a) for mask,a in candidates if (mask&missing).bit_count()==gain]
                replacement = rng.choice(ties)
                if replacement!=chosen[i]:
                    changed = True
                chosen[i] = replacement
            passes += 1
            rng.shuffle(order)
            union = 0
            for mask,a in chosen:
                union |= mask
            uncovered = len(points)-union.bit_count()
            if uncovered<best:
                best,best_choices = uncovered,list(chosen)
        restarts += 1
    # q=0 and 2q has high-coordinate label 1 in the selected parity.
    anchor = 1 if parity else q+1
    system = family+[(q,0),(2*q,anchor)]+[(m,a) for (m,_),(mask,a) in zip(resources,best_choices)]
    assert len({d for d,a in system})==len(system)
    covered = bytearray(period*q)
    for d,a in system:
        assert period*q%d==0 and d>=3 and 0<=a<d
        covered[a::d] = b'\1'*len(covered[a::d])
    checked = covered.count(0)
    assert checked==best
    return dict(N=period*q,q=q,new_parity=parity,uncovered=best,restarts=restarts,
                verified_full_period=True,seed=seed,
                scope='upper witness for uncovered count; no optimality claim',
                system=[dict(modulus=d,residue=a) for d,a in sorted(system)])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/min-uncovered-2310-five-patterns-optimal.json')
    ap.add_argument('--pattern',default='0111')
    ap.add_argument('--q',type=int,default=13)
    ap.add_argument('--parity',type=int,choices=(0,1),default=1)
    ap.add_argument('--time',type=float,default=10)
    ap.add_argument('--seed',type=int,default=1)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    source = next(x for x in data if x['pattern']==list(map(int,args.pattern))) if isinstance(data,list) else data
    family = [(x['modulus'],x['residue']) for x in source['system']]
    result = extend(family,source['N'],args.q,args.parity,args.time,args.seed)
    result['pattern'] = source['pattern']+[args.parity]
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='system'}),flush=True)
