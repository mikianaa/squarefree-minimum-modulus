"""Actual squarefree partial-cover search beyond a fixed 210 core.

Scores ordinary congruence intersections, not per-class capacity surrogates.
Every emitted family is checked on the full integer period. A positive residual
count is only an upper witness, never a nonexistence or optimality proof.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np


def crt_extend(a,modulus,labels,primes):
    for b,p in zip(labels,primes):
        a += modulus*((b-a)*pow(modulus,-1,p)%p)
        modulus *= p
    return a,modulus


def search(source,high,high_pattern,seconds=60,seed=1):
    start = time.monotonic()
    assert len(high)==len(high_pattern) and all(b in (0,1) for b in high_pattern)
    assert list(high)==sorted(set(high)) and all(p>7 and all(p%d for d in range(2,math.isqrt(p)+1)) for p in high)
    family = {v['modulus']:v['residue'] for v in source['system']}
    core = {d:a for d,a in family.items() if 210%d==0}
    assert len(core)==14 and all(d>=3 for d in core)
    period = 210*math.prod(high)
    assert period<=10_000_000, 'explicit memory/time guard for this experiment'
    residual = [[n for n in range(e,210,2) if not any(n%d==a for d,a in core.items())] for e in (0,1)]
    arrays = []
    core_arrays = []
    for e in (0,1):
        alphabets = [tuple(v for v in range(1,p) if v!=1 or b!=e) for p,b in zip(high,high_pattern)]
        grid = np.array(list(itertools.product(*alphabets)),dtype=np.int16)
        arrays.append(np.tile(grid,(len(residual[e]),1)))
        core_arrays.append(np.repeat(np.array(residual[e],dtype=np.int16),len(grid)))
    coordinates = np.concatenate(arrays)
    core_values = np.concatenate(core_arrays)
    del arrays,core_arrays
    high_codes = {}
    for k in range(1,len(high)+1):
        for subset in itertools.combinations(range(len(high)),k):
            code = np.zeros(len(core_values),dtype=np.int32)
            for i in subset:
                code = code*(high[i]-1)+coordinates[:,i]-1
            high_codes[subset] = (code,math.prod(high[i]-1 for i in subset))
    core_codes = {}
    for d in range(1,211):
        if 210%d:
            continue
        residues = sorted({n%d for r in residual for n in r})
        code = np.searchsorted(residues,core_values%d).astype(np.int16)
        core_codes[d] = (code,residues)
    resources = []
    for subset,(code,size) in high_codes.items():
        for d,(cc,residues) in core_codes.items():
            if len(subset)==1 and d in (1,2):
                continue
            modulus = d*math.prod(high[i] for i in subset)
            resources.append(dict(modulus=modulus,d=d,subset=subset,size=size,options=size*len(residues)))
    resources.sort(key=lambda r:r['modulus'])
    load = np.zeros(len(core_values),dtype=np.int16)
    choices = [-1]*len(resources)
    def codes(r):
        return high_codes[r['subset']][0]+core_codes[r['d']][0].astype(np.int32)*r['size']
    for i,r in enumerate(resources):
        if r['modulus'] not in family:
            continue
        a = family[r['modulus']]
        residues = core_codes[r['d']][1]
        labels = [a%high[j] for j in r['subset']]
        if a%r['d'] not in residues or any(b==0 for b in labels):
            continue
        code = 0
        for b,j in zip(labels,r['subset']):
            code = code*(high[j]-1)+b-1
        choices[i] = residues.index(a%r['d'])*r['size']+code
        load[codes(r)==choices[i]] += 1
    print(json.dumps(dict(stage='built',period=period,core_residual=list(map(len,residual)),
                          searched_points=len(load),resources=len(resources),seconds=round(time.monotonic()-start,3))),flush=True)
    rng = np.random.default_rng(seed)
    best = len(load)+1
    best_choices = None
    passes = 0
    deadline = time.monotonic()+seconds
    while passes==0 or time.monotonic()<deadline:
        order = np.arange(len(resources)) if passes==0 else rng.permutation(len(resources))
        for i in order:
            r = resources[i]
            code = codes(r)
            if choices[i]>=0:
                load[code==choices[i]] -= 1
            gain = np.bincount(code,weights=(load==0),minlength=r['options'])
            candidates = np.flatnonzero(gain==gain.max())
            choices[i] = int(rng.choice(candidates))
            load[code==choices[i]] += 1
        passes += 1
        uncovered = int(np.count_nonzero(load==0))
        if uncovered<best:
            best,best_choices = uncovered,list(choices)
            print(json.dumps(dict(pass_number=passes,candidate_uncovered=best,elapsed=round(time.monotonic()-start,3))),flush=True)
        if best==0:
            break
        # Occasional perturbations of substantial classes, with periodic resets
        # to the best actual family. No lower bound is inferred from stagnation.
        if passes%7==0:
            choices = list(best_choices)
            load[:] = 0
            for r,choice in zip(resources,choices):
                load[codes(r)==choice] += 1
        pool = min(50,len(resources)) if passes%3==0 else len(resources)
        for i in rng.choice(pool,size=min(4,pool),replace=False):
            r = resources[i]
            code = codes(r)
            load[code==choices[i]] -= 1
            choices[i] = int(rng.integers(r['options']))
            load[code==choices[i]] += 1
    system = list(core.items())
    system += [(p,0) for p in high]
    system += [(2*p,1 if b else p+1) for p,b in zip(high,high_pattern)]
    for r,choice in zip(resources,best_choices):
        a = core_codes[r['d']][1][choice//r['size']]
        code = choice%r['size']
        labels = []
        for i in reversed(r['subset']):
            labels.append(code%(high[i]-1)+1)
            code //= high[i]-1
        labels.reverse()
        a,d = crt_extend(a,r['d'],labels,[high[i] for i in r['subset']])
        assert d==r['modulus']
        system.append((d,a))
    assert len({d for d,a in system})==len(system)==2**(4+len(high))-2
    coverage = bytearray(period)
    for d,a in system:
        assert d>=3 and period%d==0 and 0<=a<d
        coverage[a::d] = b'\1'*len(coverage[a::d])
    verified = coverage.count(0)
    assert verified==best
    return dict(N=period,pattern=[core[2*p]%2 for p in (3,5,7)]+list(high_pattern),core_residual=list(map(len,residual)),
                uncovered=verified,verified_full_period=True,passes=passes,seed=seed,
                seconds=round(time.monotonic()-start,3),
                scope='actual partial family; upper witness only, no optimality claim',
                system=[dict(modulus=d,residue=a) for d,a in sorted(system)])


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/min-uncovered-30030-witness.json')
    ap.add_argument('--high',type=int,nargs='+',default=[11,13,17,19])
    ap.add_argument('--high-pattern',default='1111')
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--seed',type=int,default=1)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    out = search(json.loads(Path(args.source).read_text()),args.high,list(map(int,args.high_pattern)),args.time,args.seed)
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='system'}),flush=True)
