"""Independent integer check of an aggregate-relaxation primal witness."""
import argparse
import itertools
import json
import math
from pathlib import Path


def audit(data):
    source = data['input']
    witness = data['result']['primal_relaxation_witness']
    assert witness is not None
    scale = witness['denominator']
    assert type(scale) is int and scale>0
    H = source['H_upper']
    counts = source['counts']
    assert len(H)==counts[0]+1 and all(len(row)==counts[1]+1 for row in H)
    pattern = source['high_pattern']
    assert len(pattern)==4 and all(b in (0,1) for b in pattern)
    hist = [[[0]*(15 if h<4 else 17) for h in range(15)] for e in (0,1)]
    mass = [0,0]
    objective = 0
    groups = [S for k in range(2,5) for S in itertools.combinations(range(4),k)]
    for row in witness['support']:
        e,ns,w = row['parity'],row['counts'],row['weight']
        assert e in (0,1) and type(w) is int and w>0 and len(ns)==15
        D = [q-1-(b==e) for q,b in zip((11,13,17,19),pattern)]
        for h,n in enumerate(ns):
            assert type(n) is int and 0<=n<=(min(14,D[h]) if h<4 else 16)
            hist[e][h][n] += w
        mass[e] += w
        free = [d-n for d,n in zip(D,ns[:4])]
        objective += w*(math.prod(free)-sum(ns[h]*math.prod(free[i] for i in range(4) if i not in S) for h,S in enumerate(groups,4)))
    assert mass==[n*scale for n in counts]
    inequalities = 0
    for h in range(15):
        maximum = 14 if h<4 else 16
        excess = [[sum(w*max(0,n-t) for n,w in enumerate(hist[e][h])) for t in range(maximum+1)] for e in (0,1)]
        for t0,t1 in itertools.product(range(maximum+1),repeat=2):
            bound = max(v+(a+b+max(a,b) if h>=4 else 0)-t0*a-t1*b for a,row in enumerate(H) for b,v in enumerate(row))
            assert excess[0][t0]+excess[1][t1]<=scale*bound
            inequalities += 1
    assert objective==witness['objective_numerator'] and objective<0
    return dict(verified=True,counts=counts,high_pattern=pattern,denominator=scale,
                objective_numerator=objective,concentration_inequalities_checked=inequalities,
                support_size=len(witness['support']),
                scope='feasible negative-cost aggregate-relaxation distribution; not an actual covering system')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/joint-tail-010-primal-probe.json')
    ap.add_argument('--output',default='data/audit-joint-tail-010-primal.json')
    a = ap.parse_args()
    result = audit(json.loads(Path(a.source).read_text()))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
