"""Independent ordinary-residue / Python-integer conditional-tail dual audit."""
import itertools
import json
import math
from pathlib import Path
import time
from audit_fibre_concentration import gray_capacities


def audit(data,concentration,incumbent,check_subsets=True):
    start = time.monotonic()
    core = {row['modulus']:row['residue'] for row in incumbent['core_system']}
    assert sorted(core)==[d for d in range(3,211) if 210%d==0]
    R = [n for n in range(210) if not any(n%d==a for d,a in core.items())]
    assert R==data['core_residual']==concentration['core_residual']
    assert concentration['moduli']==sorted(core)
    H = concentration['two_layer_capacities']
    if check_subsets:
        assert gray_capacities(R,sorted(core))==H
        print(json.dumps(dict(stage='subset_capacity',verified=True,count=1<<len(R))),flush=True)
    high,pattern = data['high_primes'],data['high_pattern']
    assert high==[11,13,17,19] and len(pattern)==4 and all(b in (0,1) for b in pattern)
    groups = [[i] for i in range(4)]+[list(S) for k in range(2,5) for S in itertools.combinations(range(4),k)]
    assert groups==data['groups']
    shapes = [sorted(core)]*4+[[1,2]+sorted(core)]*(len(groups)-4)
    c = data['certificate']
    scale = c['denominator']
    assert type(scale) is int and scale>0
    common = [[[0]*(len(s)+1) for s in shapes] for e in (0,1)]
    tail_term = 0
    for row in c['tail_multipliers']:
        h = row['axis']
        t0,t1 = row['thresholds']
        assert 0<=h<len(groups) and 0<=t0<=len(shapes[h]) and 0<=t1<=len(shapes[h])
        coefficient = row['multiplier']
        assert type(coefficient) is int and coefficient<=0
        bound = max(v+(a+b+max(a,b) if h>=4 else 0)-t0*a-t1*b for a,line in enumerate(H) for b,v in enumerate(line))
        assert bound==row['bound']
        tail_term += coefficient*bound
        for e,t in enumerate((t0,t1)):
            for n in range(len(shapes[h])+1):
                common[e][h][n] += coefficient*max(0,n-t)
    aggregate = data.get('aggregate',False)
    point_residues = [0,1] if aggregate else R
    masses = [sum(n%2==e for n in R) for e in (0,1)] if aggregate else [1]*len(R)
    if 'point_masses' in data:
        assert masses==data['point_masses']
    links = c['link_multipliers']
    assert len(links)==len(point_residues) and all(len(row)==len(groups) and all(type(v) is int and v<=0 for v in row) for row in links)
    if aggregate:
        assert all(v==0 for row in links for v in row)
    checked = 0
    point_minima = []
    for x,value in enumerate(point_residues):
        e = value%2
        D = [q-1-(b==e) for q,b in zip(high,pattern)]
        free_indices = [[i for i in range(4) if i not in S] for S in groups[4:]]
        tables = []
        for h,indices in enumerate(free_indices,4):
            maximum_volume = math.prod(D[i] for i in indices)
            table = [min(-scale*v*n-links[x][h]*n-common[e][h][n] for n in range(17)) for v in range(maximum_volume+1)]
            tables.append(table)
        minimum = None
        for ns in itertools.product(*(range(min(14,d)+1) for d in D)):
            free = [d-n for d,n in zip(D,ns)]
            actual = scale*math.prod(free)-sum(links[x][i]*n+common[e][i][n] for i,n in enumerate(ns))
            actual += sum(table[math.prod(free[i] for i in indices)] for table,indices in zip(tables,free_indices))
            if minimum is None or actual<minimum:
                minimum = actual
            checked += 1
        assert minimum==c['point_multipliers'][x],(x,minimum,c['point_multipliers'][x])
        point_minima.append(minimum)
    resources = []
    if not aggregate:
        for h,moduli in enumerate(shapes):
            for d in moduli:
                resources.append(min(sum(links[x][h] for x,n in enumerate(R) if n%d==a) for a in {n%d for n in R}))
    assert resources==c['resource_multipliers']
    numerator = sum(m*v for m,v in zip(masses,point_minima))+sum(resources)+tail_term
    assert numerator==c['numerator']
    lower = -(-numerator//scale)
    assert lower==c['integer_lower']
    return dict(verified=True,exact_numerator=numerator,exact_denominator=scale,
                conditional_residual_lower=lower,excludes_this_core_and_pattern=lower>0,
                core_subsets_checked=(1<<len(R)) if check_subsets else None,
                plane_count_tuples_checked=checked,tail_counts_per_group_checked=17,
                core_system=incumbent['core_system'],high_primes=high,high_pattern=pattern,
                seconds=round(time.monotonic()-start,3),
                scope='one specified core and high-prime anchor pattern only; not all eight-prime covers')


if __name__=='__main__':
    read = lambda name: json.loads(Path('data',name).read_text())
    result = audit(read('conditional-tail-eight-lp.json'),read('core-hit-concentration-7-15.json'),read('core-fibre-eight-plane.json'))
    Path('data/audit-conditional-tail-eight.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='core_system'}))
