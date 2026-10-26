"""Independent integer audit of full per-group configuration mixtures."""
import argparse
import itertools
import json
import math
from pathlib import Path


def audit(data):
    source,witness = data['input'],data['result']['coherent_primal_witness']
    R,high,pattern = [source[k] for k in ('core_residual','high_primes','high_pattern')]
    n = len(high)
    assert len(pattern)==n and all(b in (0,1) for b in pattern)
    assert len(R)==len(set(R)) and all(type(v) is int and 0<=v<210 for v in R)
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    shapes = [[d for d in range(1 if h>=n else 3,211) if 210%d==0] for h in range(len(groups))]
    scale = witness['denominator']
    assert type(scale) is int and scale>0
    supply = [[[0]*(len(shape)+1) for x in R] for shape in shapes]
    seen,assignments = set(),0
    for row in witness['groups']:
        h = row['axis']
        assert type(h) is int and 0<=h<len(groups) and h not in seen
        seen.add(h)
        assert sum(a['weight'] for a in row['assignments'])==scale
        for assignment in row['assignments']:
            residues,w = assignment['residues'],assignment['weight']
            assert type(w) is int and w>0 and len(residues)==len(shapes[h])
            assert all(type(a) is int and 0<=a<d for a,d in zip(residues,shapes[h]))
            for x,v in enumerate(R):
                hits = sum(v%d==a for d,a in zip(shapes[h],residues))
                for t in range(1,hits+1):
                    supply[h][x][t] += w
            assignments += 1
    assert seen==set(range(len(groups)))
    demand = [[[0]*(len(shape)+1) for x in R] for shape in shapes]
    mass,objective = [0]*len(R),0
    point_costs = [0]*len(R)
    for row in witness['support']:
        x,state,w = [row[k] for k in ('point_index','counts','weight')]
        assert type(x) is int and 0<=x<len(R) and type(w) is int and w>0
        assert len(state)==len(groups)
        D = [q-1-(b==R[x]%2) for q,b in zip(high,pattern)]
        mass[x] += w
        for h,v in enumerate(state):
            assert type(v) is int and 0<=v<=(min(14,D[h]) if h<n else 16)
            for t in range(1,v+1):
                demand[h][x][t] += w
        free = [d-v for d,v in zip(D,state[:n])]
        cost = math.prod(free)-sum(state[h]*math.prod(free[i] for i in range(n) if i not in S)
                                   for h,S in enumerate(groups[n:],n))
        objective += w*cost
        point_costs[x] += w*cost
    assert mass==[scale]*len(R)
    checks = 0
    for h,shape in enumerate(shapes):
        for x in range(len(R)):
            for t in range(1,len(shape)+1):
                assert demand[h][x][t]<=supply[h][x][t]
                checks += 1
    assert objective==witness['objective_numerator']
    positive_means = sum(max(0,v) for v in point_costs)
    if 'positive_point_means_numerator' in witness:
        assert positive_means==witness['positive_point_means_numerator']
    return dict(verified=True,negative_cost=objective<0,objective_numerator=objective,denominator=scale,
                positive_point_means_numerator=positive_means,point_cost_numerators=point_costs,
                all_point_means_nonpositive=positive_means==0,
                groups_checked=len(groups),resource_assignments_checked=assignments,threshold_links_checked=checks,
                support_size=len(witness['support']),
                scope='consistent assignment mixture for each group separately; not a joint covering system')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    result = audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
