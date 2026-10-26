"""Independent standard-library verification of a compressed-tail branch.

Enumerates ordinary core classes, explicit point-weight compression tables,
and every local plane-count tuple. Tail counts are minimized independently
over all 17 possibilities. No orbit reduction, optimization solver, or
floating-point arithmetic is used.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
from audit_eight_prime_anchor_branch import ordinary_frontier,explicit_capacity


def local_minima(certificate,pattern):
    scale = certificate['denominator']
    assert type(scale) is int and scale>0
    assert not certificate['resource_multipliers']
    assert all(v==0 for row in certificate['link_multipliers'] for v in row)
    common = [[[0]*(15 if h<4 else 17) for h in range(15)] for e in (0,1)]
    for row in certificate['tail_multipliers']:
        h = row['axis']
        y = row['multiplier']
        assert type(h) is int and 0<=h<15 and type(y) is int and y<=0
        assert len(row['thresholds'])==2
        for e,t in enumerate(row['thresholds']):
            assert type(t) is int and 0<=t<len(common[e][h])
            for n in range(len(common[e][h])):
                common[e][h][n] += y*max(0,n-t)
    subsets = [S for k in range(2,5) for S in itertools.combinations(range(4),k)]
    free_indices = [[i for i in range(4) if i not in S] for S in subsets]
    result = []
    checked = 0
    for e in (0,1):
        D = [q-1-(b==e) for q,b in zip((11,13,17,19),pattern)]
        tables = [[min(-scale*v*n-common[e][h][n] for n in range(17))
                   for v in range(math.prod(D[i] for i in indices)+1)] for h,indices in enumerate(free_indices,4)]
        minimum = None
        for ns in itertools.product(*(range(min(14,d)+1) for d in D)):
            free = [d-n for d,n in zip(D,ns)]
            value = scale*math.prod(free)-sum(common[e][i][n] for i,n in enumerate(ns))
            value += sum(table[math.prod(free[i] for i in indices)] for table,indices in zip(tables,free_indices))
            minimum = value if minimum is None else min(minimum,value)
            checked += 1
        result.append(minimum)
    return result,checked


def audit(data):
    start = time.monotonic()
    pattern = data['core_pattern']
    assert pattern in ([0,0,0],[0,0,1],[0,1,0],[0,1,1]) and data['period']==9699690
    assert data['complete_positive_bounds']
    expected_patterns = [list(p) for p in itertools.product((0,1),repeat=4)]
    assert data['high_patterns']==expected_patterns
    frontier,combinations = ordinary_frontier(pattern)
    assert frontier==data['frontier']
    assert len(data['selection'])==len(data['bounds'])==len(frontier)
    assert all(len(row)==16 for row in data['selection']+data['bounds'])
    print(json.dumps(dict(stage='ordinary_frontier',pattern=pattern,pairs=frontier,combinations=combinations,seconds=round(time.monotonic()-start,3))),flush=True)
    H = explicit_capacity(pattern)
    points = {}
    K_cache = {}
    checked = 0
    result_bounds = []
    for i,(a,b) in enumerate(frontier):
        line = []
        for j,high in enumerate(expected_patterns):
            index = data['selection'][i][j]
            assert type(index) is int and 0<=index<len(data['pool'])
            item = data['pool'][index]
            c = item['certificate']
            key = index,j
            if key not in points:
                values,count = local_minima(c,high)
                assert values==item['point_multipliers_by_pattern'][j]
                points[key] = values
                checked += count
            numerator = a*points[key][0]+b*points[key][1]
            for row in c['tail_multipliers']:
                t0,t1 = row['thresholds']
                tail = row['axis']>=4
                k = a,b,tail,t0,t1
                if k not in K_cache:
                    K_cache[k] = max(H[m0][m1]+(m0+m1+max(m0,m1) if tail else 0)-t0*m0-t1*m1
                                     for m0 in range(a+1) for m1 in range(b+1))
                numerator += row['multiplier']*K_cache[k]
            lower = -(-numerator//c['denominator'])
            assert lower==data['bounds'][i][j] and lower>0
            line.append(lower)
        result_bounds.append(line)
        print(json.dumps(dict(stage='frontier_pair_verified',counts=[a,b],minimum=min(line),seconds=round(time.monotonic()-start,3))),flush=True)
    minimum = min(map(min,result_bounds))
    assert minimum==data['uniform_lower']
    return dict(verified=True,period=9699690,core_pattern=pattern,frontier=frontier,high_patterns=expected_patterns,
                bounds=result_bounds,uniform_uncovered_lower=minimum,ordinary_pair_checks=combinations,
                local_inequalities_checked=len(points),plane_tuples_checked=checked,seconds=round(time.monotonic()-start,3),
                scope='complete normalized low-anchor branch, all cores and all higher anchor patterns; not a global conjecture proof')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='001')
    a = ap.parse_args()
    data = json.loads(Path(f'data/compressed-tail-branch-{a.pattern}.json').read_text())
    result = audit(data)
    Path(f'data/audit-compressed-tail-branch-{a.pattern}.json').write_text(json.dumps(result,indent=2)+'\n')
    data['independent_audit'] = True
    data['audit_report'] = f'data/audit-compressed-tail-branch-{a.pattern}.json'
    Path(f'data/compressed-tail-branch-{a.pattern}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('frontier','high_patterns','bounds')}))
