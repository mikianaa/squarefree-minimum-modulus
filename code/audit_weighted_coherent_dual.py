"""Solver-free integer audit of weighted point duals and group capacities.

Can certify a bound over a stated finite configuration pool and a possibly
weaker unrestricted bound using elementary capacity majorants. CP-SAT bounds
are deliberately not used as proof inputs.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import numpy as np


def minimum(D,coefficient,potentials,groups):
    # Enumerate free-coordinate values, independently of the generator's
    # count-grid construction. Tail minimization uses 17 explicit choices.
    ranges = [np.arange(max(0,d-14),d+1,dtype=np.int64) for d in D]
    free = [v.reshape(-1) for v in np.meshgrid(*ranges,indexing='ij')]
    n = len(D)
    risk = coefficient*math.prod(D)+sum(max(map(abs,row)) for row in potentials)
    risk += sum(16*coefficient*math.prod(D[i] for i in range(n) if i not in S) for S in groups[n:])
    assert risk<2**63
    values = coefficient*np.prod(np.stack(free),axis=0)
    for h in range(n):values -= np.array(potentials[h],dtype=np.int64)[D[h]-free[h]]
    for h,S in enumerate(groups[n:],n):
        indices = [i for i in range(n) if i not in S]
        max_volume = math.prod(D[i] for i in indices)
        volumes = np.arange(max_volume+1,dtype=np.int64)
        table = np.full(max_volume+1,-potentials[h][0],dtype=np.int64)
        for count in range(1,17):
            table = np.minimum(table,-coefficient*count*volumes-potentials[h][count])
        actual = np.ones(len(values),dtype=np.int64)
        for i in indices:actual *= free[i]
        values += table[actual]
    return int(values.min()),len(values)


def capacity(points,shapes,events,certificate):
    if certificate.get('method')=='clique_affine':
        from audit_clique_event_majorant import audit_bound
        return audit_bound(points,shapes,events,certificate)
    if certificate.get('method')=='pair_affine':
        from audit_pair_event_majorant import audit_bound
        return audit_bound(points,shapes,events,certificate)
    if certificate.get('method')=='branch_affine':
        from audit_branched_event_bound import audit_bound
        return audit_bound(points,shapes,events,certificate)
    base = int(1 in shapes);limit=len(shapes)
    assert certificate['base']==base
    slopes = [Fraction(*v) for v in certificate['slopes']]
    assert len(slopes)==len(points) and all(v>=0 for v in slopes)
    rewards = [[sum(e['weight'] for e in events if e['point']==x and e['threshold']<=k)
                for k in range(limit+1)] for x in range(len(points))]
    offsets=[Fraction(*v) for v in certificate['point_intercepts']] if 'point_intercepts' in certificate else [Fraction(row[base]) for row in rewards]
    assert len(offsets)==len(points)
    for row,slope,offset in zip(rewards,slopes,offsets):
        assert all(row[k]<=offset+(k-base)*slope for k in range(base,limit+1))
    intercept = sum(offsets)
    if 'intercept_fraction' in certificate:assert intercept==Fraction(*certificate['intercept_fraction'])
    else:assert intercept==certificate['intercept']
    capacities = []
    for d in shapes:
        if d==1:continue
        capacities.append(max(sum(w for x,w in zip(points,slopes) if x%d==a) for a in range(d)))
    assert capacities==[Fraction(*v) for v in certificate['resource_bounds']]
    total = intercept+sum(capacities)
    bound = min(sum(row[-1] for row in rewards),total.numerator//total.denominator)
    assert bound==certificate['bound']
    return bound


def audit(data,index=-1):
    source,result = data['input'],data['result']
    c = result['pricing_history'][index]
    R,high,pattern = [source[k] for k in ('core_residual','high_primes','high_pattern')]
    n = len(high)
    assert n==len(pattern) and all(b in (0,1) for b in pattern)
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    shapes = [[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(len(groups))]
    weight_den,event_scale = c['point_weight_denominator'],c['event_scale']
    assert type(weight_den) is int and weight_den>0 and type(event_scale) is int and event_scale>0
    assert c['denominator']==weight_den*event_scale
    theta,multipliers = c['point_weights'],c['event_multipliers']
    assert len(theta)==len(R) and all(type(v) is int and 0<=v<=weight_den for v in theta)
    assert len(multipliers)==len(groups)
    for h,group in enumerate(multipliers):
        assert len(group)==len(R)
        assert all(len(row)==len(shapes[h])+1 and row[0]==0 and all(type(v) is int and v<=0 for v in row) for row in group)
    minima,checked = [],0
    for x,point in enumerate(R):
        D = [q-1-(b==point%2) for q,b in zip(high,pattern)]
        potentials = [[weight_den*sum(group[x][:t+1]) for t in range(len(group[x]))] for group in multipliers]
        value,count = minimum(D,event_scale*theta[x],potentials,groups)
        minima.append(value);checked += count
    assert minima==c['point_multipliers']
    capacity_bounds,pool_bounds = [],[]
    configurations_checked = 0
    for h,(shape,row) in enumerate(zip(shapes,c['group_pricing'])):
        expected = [dict(point=x,threshold=t,weight=-multipliers[h][x][t]) for x in range(len(R))
                    for t in range(1,len(shape)+1) if multipliers[h][x][t]]
        assert row['events']==expected
        capacity_bounds.append(capacity(R,shape,expected,row['capacity_certificate']))
        if c.get('pool_sizes') is not None:
            configs = result['group_columns'][h][:c['pool_sizes'][h]]
            assert len(configs)==c['pool_sizes'][h] and configs
            maximum = 0
            for residues in configs:
                assert len(residues)==len(shape) and all(type(a) is int and 0<=a<d for a,d in zip(residues,shape))
                counts = [sum(x%d==a for a,d in zip(residues,shape)) for x in R]
                maximum = max(maximum,sum(e['weight'] for e in expected if counts[e['point']]>=e['threshold']))
                configurations_checked += 1
            pool_bounds.append(maximum)
    numerator = sum(minima)-weight_den*sum(capacity_bounds)
    assert numerator==c['capacity_numerator']
    pool_numerator = None
    if pool_bounds:
        assert pool_bounds==c['pool_upper_bounds']
        pool_numerator = sum(minima)-weight_den*sum(pool_bounds)
        assert pool_numerator==c['pool_numerator']
    return dict(verified=True,denominator=c['denominator'],unrestricted_capacity_numerator=numerator,
                pool_numerator=pool_numerator,plane_states_checked=checked,groups_checked=len(groups),
                configurations_checked=configurations_checked,pricing_index=index,
                excludes_unrestricted_core_pattern=numerator>0,
                excludes_stated_pool=pool_numerator is not None and pool_numerator>0,
                scope='weighted positive-point-expectation bound for this core/pattern; CP-SAT bounds not used')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--index',type=int,default=-1)
    args = ap.parse_args()
    report = audit(json.loads(Path(args.source).read_text()),args.index)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
