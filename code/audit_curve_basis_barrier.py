"""Independent integer/rational audit of a specified curve-cone barrier.

Only the standard library is used. Numerical objectives and generator repair
claims are ignored. Fractional group states never constitute a covering.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path


def positive_int(v):return type(v) is int and v>0


def rational(row):
    assert len(row)==2 and type(row[0]) is int and positive_int(row[1])
    return Fraction(*row)


def audit(data,comparison=None):
    source=data['input'];points=source['core_residual'];high=source['high_primes'];pattern=source['high_pattern']
    n=len(high);assert 1<=n<=5 and len(pattern)==n and all(b in (0,1) for b in pattern)
    assert high==sorted(set(high)) and all(type(q) is int and q>2 and all(q%d for d in range(2,math.isqrt(q)+1)) for q in high)
    assert points and len(set(points))==len(points) and all(type(p) is int and 0<=p<210 for p in points)
    G=2**n-1;P=len(points);subsets=[s for k in range(1,n+1) for s in itertools.combinations(range(n),k)]
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    events=data['basis_event_multipliers'];assert events
    def check_events(axes):
        assert len(axes)==G
        for h,axis in enumerate(axes):
            assert len(axis)==P
            for row in axis:assert len(row)==len(shapes[h])+1 and row[0]==0 and all(type(v) is int and v<=0 for v in row)
    for axes in events:check_events(axes)
    D=data['probability_denominator'];assert positive_int(D)
    point_mix=data['point_mixtures'];assert sorted(row['point_index'] for row in point_mix)==list(range(P))
    demand=[[[0]*(len(shape)+1) for p in points] for shape in shapes];point_den={};local_count=0
    for row in point_mix:
        x=row['point_index'];den=row['denominator'];assert positive_int(den)
        assert row['states'] and all(positive_int(s['weight']) for s in row['states'])
        assert sum(s['weight'] for s in row['states'])==den;point_den[x]=den
        available=[q-1-int(b==points[x]%2) for q,b in zip(high,pattern)]
        for s in row['states']:
            k=s['counts'];assert len(k)==G and all(type(v) is int and 0<=v<=(min(14,available[h]) if h<n else 16) for h,v in enumerate(k))
            free=[a-v for a,v in zip(available,k)]
            assert math.prod(free)<=sum(k[h]*math.prod(free[i] for i in range(n) if i not in S) for h,S in enumerate(subsets[n:],n))
            for h,v in enumerate(k):
                for t in range(v+1):demand[h][x][t]+=s['weight']
            local_count+=1
    group_mix=data['group_mixtures'];assert sorted(row['axis'] for row in group_mix)==list(range(G))
    supply=[[[0]*(len(shape)+1) for p in points] for shape in shapes];group_den={}
    actual_count=fractional_count=pair_checks=resource_checks=0
    for group in group_mix:
        h=group['axis'];shape=shapes[h];N=len(shape);den=group['denominator'];assert positive_int(den)
        assert group['states'] and all(positive_int(s['weight']) for s in group['states'])
        assert sum(s['weight'] for s in group['states'])==den;group_den[h]=den*D
        for s in group['states']:
            hist=[[0]*(N+1) for p in points]
            if s['kind']=='actual':
                assignments=s['assignment'];assert len(assignments)==N
                assert all(type(a) is int and 0<=a<d for a,d in zip(assignments,shape))
                for x,p in enumerate(points):hist[x][sum(p%d==a for a,d in zip(assignments,shape))]=D
                actual_count+=1
            else:
                assert s['kind']=='fractional'
                fixed=s['fixed'];assert all(len(r)==2 and all(type(v) is int for v in r) for r in fixed)
                assert len(dict(fixed))==len(fixed) and all(d in shape and 0<=a<d for d,a in fixed)
                remaining=[d for d in shape if d not in dict(fixed)]
                base=int(1 in remaining);hits=[sum(p%d==a for d,a in fixed) for p in points]
                resources=s['resource_distributions'];assert sorted(r['resource'] for r in resources)==[d for d in remaining if d!=1]
                supply_mean=[0]*P
                for row in resources:
                    d=row['resource'];residues=row['residues'];units=row['probability_numerators']
                    assert residues and len(set(residues))==len(residues) and all(type(a) is int and 0<=a<d for a in residues)
                    assert len(units)==len(residues) and all(type(v) is int and v>=0 for v in units) and sum(units)==D
                    for x,p in enumerate(points):supply_mean[x]+=sum(v for a,v in zip(residues,units) if p%d==a)
                    resource_checks+=1
                counts=s['count_numerators'];assert len(counts)==P
                for x,row in enumerate(counts):
                    assert len(row)==len(remaining)-base+1 and all(type(v) is int and v>=0 for v in row) and sum(row)==D
                    assert sum(j*v for j,v in enumerate(row))<=supply_mean[x]
                    for j,v in enumerate(row):hist[x][hits[x]+base+j]=v
                fractional_count+=1
            tail=[]
            for row in hist:
                running=0;line=[]
                for v in reversed(row):running+=v;line.append(running)
                tail.append(list(reversed(line)))
            if s['kind']=='fractional':
                for x,y in itertools.combinations(range(P),2):
                    capacity=hits[x]+hits[y]+sum(max(int(points[x]%d==a)+int(points[y]%d==a)
                        for a in {points[x]%d,points[y]%d}) for d in remaining)
                    for t in range(1,N+1):
                        u=capacity+1-t
                        if 1<=u<=N:
                            assert tail[x][t]+tail[y][u]<=D
                            pair_checks+=1
            for x,row in enumerate(tail):
                for t,v in enumerate(row):supply[h][x][t]+=s['weight']*v
    claimed={(r['basis'],r['axis'],r['point_index']):rational(r['surplus']) for r in data['coordinate_surpluses']}
    assert len(claimed)==len(data['coordinate_surpluses']);expected=set();normalized=[]
    for b,axes in enumerate(events):
        for h,axis in enumerate(axes):
            for x,row in enumerate(axis):
                mass=-sum(row)
                if not mass or not -sum(events[0][h][x]):continue
                value=Fraction(sum(-v*supply[h][x][t] for t,v in enumerate(row)),group_den[h])-Fraction(
                    sum(-v*demand[h][x][t] for t,v in enumerate(row)),point_den[x])
                key=(b,h,x);expected.add(key);assert value>=0 and claimed[key]==value
                normalized.append(value/mass)
    assert expected==set(claimed) and normalized
    margin=min(normalized);assert margin==rational(data['minimum_normalized_surplus'])
    deficits=[]
    for h,shape in enumerate(shapes):
        for x in range(P):
            if not -sum(events[0][h][x]):continue
            for t in range(1,len(shape)+1):
                gap=Fraction(supply[h][x][t],group_den[h])-Fraction(demand[h][x][t],point_den[x])
                if gap<0:deficits.append(dict(axis=h,point_index=x,threshold=t,supply_minus_demand=[gap.numerator,gap.denominator]))
    report=dict(verified=True,curve_cone_barrier_verified=True,basis_count=len(events),coordinate_inequalities=len(expected),
        active_local_states=local_count,actual_group_states=actual_count,fractional_group_states=fractional_count,
        conditional_pair_inequalities_checked=pair_checks,resource_distributions_checked=resource_checks,
        minimum_normalized_surplus=[margin.numerator,margin.denominator],
        all_allowed_nonnegative_threshold_weights_barrier_verified=not deficits,
        deficient_thresholds=len(deficits),threshold_deficits=deficits,
        actual_assignment_only_barrier=fractional_count==0,covering_constructed=False,global_problem_settled=False,
        scope='Only the specified curve cone with positive baseline caps and group bounds admitting the displayed actual and fractional case states. A deeper case refinement, stronger compatibility, or additional curves can escape this barrier.')
    if comparison is not None:
        check_events(comparison);value=Fraction(0)
        for h,axis in enumerate(comparison):
            for x,row in enumerate(axis):
                value+=Fraction(sum(-v*demand[h][x][t] for t,v in enumerate(row)),point_den[x])-Fraction(
                    sum(-v*supply[h][x][t] for t,v in enumerate(row)),group_den[h])
        report.update(comparison_local_minus_group=[value.numerator,value.denominator],comparison_defeats_this_witness=value>0)
    return report


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True);ap.add_argument('--output',required=True);ap.add_argument('--compare')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());comparison=None
    if args.compare:
        other=json.loads(Path(args.compare).read_text());assert other['input']==data['input'];comparison=other['result']['pricing_history'][-1]['event_multipliers']
    report=audit(data,comparison);Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='threshold_deficits'}))
