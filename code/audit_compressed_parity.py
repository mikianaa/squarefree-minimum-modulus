"""Standard-library audit of compressed parity certificates.

Builds point-type multiplicities by ordinary residue enumeration. Each assignment
is evaluated using explicit integer weight histograms and running partial sums;
does not import the NumPy producer or its weight/propagation routines.
"""
from collections import Counter
import itertools
import json
import math
from pathlib import Path
import time


def boundary(pairs):
    best = math.inf
    out = []
    for a,b in sorted(set(map(tuple,pairs))):
        if b<best:
            out.append([a,b])
            best = b
    return out


def audit_extension(primes,pattern,frontier,q,b):
    period = 2*math.prod(primes)
    types = [Counter(),Counter()]
    for n in range(period):
        e = n%2
        if any(n%p==0 or (e==parity and n%p==1) for p,parity in zip(primes,pattern)):
            continue
        mask = sum(1<<i for i,p in enumerate(primes) if n%p==2)
        types[e][mask] += 1
    subset_bits = [sum(1<<(s-1) for s in range(1,1<<len(primes)) if s&z==s)
                   for z in range(1<<len(primes))]
    limits = [max(r[e] for r in frontier) for e in (0,1)]
    caches = [{},{}]
    minimal = {}
    for assignment in range(1<<((1<<len(primes))-1)):
        tables = []
        for e in (0,1):
            histogram = Counter()
            for z,count in types[e].items():
                all_odd = subset_bits[z].bit_count()
                assigned_first = (assignment&subset_bits[z]).bit_count()
                weight = all_odd+(assigned_first if e==0 else all_odd-assigned_first)
                histogram[weight] += count
            key = tuple(sorted(histogram.items(),reverse=True))
            if key not in caches[e]:
                values = [0]
                used = mass = 0
                for weight,count in key:
                    for _ in range(min(count,limits[e]-used)):
                        used += 1
                        mass += weight
                        values.append(max(0,(q-1-(b==e))*used-mass))
                    if used==limits[e]:
                        break
                assert len(values)==limits[e]+1
                caches[e][key] = values
            tables.append(caches[e][key])
        for r0,r1 in frontier:
            a,c = tables[0][r0],tables[1][r1]
            minimal[a] = min(minimal.get(a,math.inf),c)
    return boundary(minimal.items())


def main():
    start = time.monotonic()
    source = json.loads(Path('data/three-axis-210-certificate.json').read_text())['cases']
    certificate = json.loads(Path('data/compressed-parity-frontiers.json').read_text())
    optima = {tuple(r['pattern']):r['uncovered'] for r in
              json.loads(Path('data/min-uncovered-2310-five-patterns-optimal.json').read_text())}
    results = []
    for row in certificate['at2310']:
        p = row['pattern']
        base = next(r for r in source if r['pattern']==p[:3])
        computed = audit_extension((3,5,7),p[:3],base['residual_pair_frontier'],11,p[3])
        assert computed==row['frontier'] and min(map(sum,computed))==row['uncovered_lower']
        result = dict(period=2310,pattern=p,verified=True,lower=row['uncovered_lower'])
        results.append(result)
        print(json.dumps(result),flush=True)
    for row in certificate['at30030']:
        p = row['pattern']
        base = next(r for r in certificate['at2310'] if r['pattern']==p[:4])['frontier']
        if row['status']=='including_solver':
            lower = optima.get(tuple(p[:4]),0)
            base = boundary((a+t,max(b,lower-a-t)) for a,b in base for t in range(max(0,lower-a-b)+1))
        computed = audit_extension((3,5,7,11),p[:4],base,13,p[4])
        assert computed==row['frontier'] and min(map(sum,computed))==row['uncovered_lower']
        result = dict(period=30030,pattern=p,status=row['status'],verified=True,lower=row['uncovered_lower'])
        results.append(result)
        print(json.dumps(result),flush=True)
    out = dict(verified=True,seconds=round(time.monotonic()-start,3),cases=results,
               scope='arithmetic audit; solver-dependent inputs remain solver-dependent')
    Path('data/audit-compressed-parity.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
