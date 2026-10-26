"""Local coverage search with prime and twice-prime classes held fixed.

Outputs are full-period-verified partial families, not lower-bound proofs.
"""
import argparse
import json
from pathlib import Path
import random
import time


def improve(source,seconds=30,seed=1):
    period = source['N']
    family = {v['modulus']:v['residue'] for v in source['system']}
    primes = [p for p in range(3,period+1) if period%p==0 and all(p%d for d in range(2,int(p**0.5)+1))]
    anchors = {d:a for d,a in family.items() if d in primes or d in [2*p for p in primes]}
    points = [n for n in range(period) if not any(n%d==a for d,a in anchors.items())]
    whole = (1<<len(points))-1
    resources = []
    initial = []
    for d in sorted(family):
        if d in anchors:
            continue
        lookup = {}
        for j,n in enumerate(points):
            a = n%d
            lookup[a] = lookup.get(a,0) | (1<<j)
        choices = [(mask,a) for a,mask in lookup.items()]
        resources.append((d,choices))
        initial.append((lookup.get(family[d],0),family[d]))
    def missed(chosen):
        union = 0
        for mask,a in chosen:
            union |= mask
        return len(points)-union.bit_count()
    rng = random.Random(seed)
    chosen = list(initial)
    best = missed(chosen)
    best_choices = list(chosen)
    deadline = time.monotonic()+seconds
    passes = 0
    while passes==0 or time.monotonic()<deadline:
        order = list(range(len(resources)))
        rng.shuffle(order)
        strict = False
        for i in order:
            others = 0
            for j,(mask,a) in enumerate(chosen):
                if i!=j:
                    others |= mask
            missing = whole ^ others
            previous = (chosen[i][0]&missing).bit_count()
            maximum = previous
            ties = [chosen[i]]
            for mask,a in resources[i][1]:
                gain = (mask&missing).bit_count()
                if gain>maximum:
                    maximum,ties = gain,[(mask,a)]
                elif gain==maximum:
                    ties.append((mask,a))
            chosen[i] = rng.choice(ties)
            strict |= maximum>previous
        passes += 1
        uncovered = missed(chosen)
        if uncovered<best:
            best,best_choices = uncovered,list(chosen)
            print(json.dumps(dict(pass_number=passes,uncovered=best)),flush=True)
        if not strict:
            # Escape flat local optima, periodically returning to the best family.
            if passes%7==0:
                chosen = list(best_choices)
            # Include smaller moduli in some kicks: uniformly selected classes
            # are mostly tiny cylinders and leave the large-scale layout fixed.
            pool = list(range(min(16,len(resources)))) if passes%3==0 else list(range(len(resources)))
            for i in rng.sample(pool,min(4,len(pool))):
                chosen[i] = rng.choice(resources[i][1])
    system = sorted(list(anchors.items())+[(d,a) for (d,_),(mask,a) in zip(resources,best_choices)])
    coverage = bytearray(period)
    for d,a in system:
        assert period%d==0 and d>=3 and 0<=a<d
        coverage[a::d] = b'\1'*len(coverage[a::d])
    assert coverage.count(0)==best and len({d for d,a in system})==len(system)
    return dict(N=period,pattern=source['pattern'],uncovered=best,initial_uncovered=missed(initial),
                passes=passes,verified_full_period=True,seed=seed,
                scope='upper witness for uncovered count; no optimality claim',
                system=[dict(modulus=d,residue=a) for d,a in system])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--time',type=float,default=30)
    ap.add_argument('--seed',type=int,default=1)
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    out = improve(json.loads(Path(a.source).read_text()),a.time,a.seed)
    Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='system'}),flush=True)
