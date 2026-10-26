"""Independent rational audit of a fixed-case, fixed-event-curve barrier.

Uses only the standard library and reconstructs every reward and local cost.
No solver objective, numerical scenario payoff, or generator is trusted.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path


def rational(value):
    assert len(value)==2 and all(type(v) is int for v in value) and value[1]>0
    return Fraction(*value)


def audit_tree(points,shapes,nodes):
    assert nodes and nodes[0]['fixed']==[]
    visited=set();leaves=set()
    def visit(j,expected):
        assert type(j) is int and 0<=j<len(nodes) and j not in visited
        visited.add(j);node=nodes[j];assert node['fixed']==expected
        if 'children' not in node:
            assert 'resource' not in node
            leaves.add(tuple(map(tuple,expected)));return
        d=node['resource'];assert type(d) is int and d in shapes and d!=1 and d not in dict(expected)
        residues=sorted({p%d for p in points})
        assert len(node['children'])==len(residues)
        for child,a in zip(node['children'],residues):visit(child,expected+[[d,a]])
    visit(0,[]);assert visited==set(range(len(nodes)))
    return leaves


def audit(data,comparison_multipliers=None):
    source=data['input'];R=source['core_residual'];high=source['high_primes'];pattern=source['high_pattern']
    n=len(high);G=2**n-1
    assert n and len(pattern)==n and all(b in (0,1) for b in pattern)
    assert R and len(set(R))==len(R) and all(type(p) is int and 0<=p<210 for p in R)
    assert all(type(q) is int and q>2 and all(q%d for d in range(2,math.isqrt(q)+1)) for q in high)
    assert high==sorted(set(high))
    subsets=[s for k in range(1,n+1) for s in itertools.combinations(range(n),k)]
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    events=data['baseline_event_multipliers'];assert len(events)==G
    rewards=[]
    for h,axis in enumerate(events):
        assert len(axis)==len(R);curves=[]
        for row in axis:
            assert len(row)==len(shapes[h])+1 and row[0]==0 and all(type(v) is int and v<=0 for v in row)
            curves.append([-sum(row[:k+1]) for k in range(len(row))])
        rewards.append(curves)
    trees=data['case_trees'];assert sorted(row['axis'] for row in trees)==list(range(G))
    leaves={row['axis']:audit_tree(R,shapes[row['axis']],row['nodes']) for row in trees}
    demand=[[Fraction(0) for p in R] for h in range(G)];supply=[[Fraction(0) for p in R] for h in range(G)]
    demand_tail=[[[Fraction(0) for k in curve] for curve in axis] for axis in rewards]
    supply_tail=[[[Fraction(0) for k in curve] for curve in axis] for axis in rewards]
    point_mix=data['point_mixtures'];assert sorted(row['point_index'] for row in point_mix)==list(range(len(R)))
    local_count=0
    for row in point_mix:
        x=row['point_index'];den=row['denominator'];assert type(den) is int and den>0
        assert sum(s['weight'] for s in row['states'])==den
        D=[q-1-int(b==R[x]%2) for q,b in zip(high,pattern)]
        for s in row['states']:
            weight=s['weight'];counts=s['counts'];assert type(weight) is int and weight>0 and len(counts)==G
            assert all(type(k) is int and 0<=k<=(min(14,D[h]) if h<n else 16) for h,k in enumerate(counts))
            free=[d-k for d,k in zip(D,counts)]
            cost=math.prod(free)-sum(counts[h]*math.prod(free[i] for i in range(n) if i not in S) for h,S in enumerate(subsets[n:],n))
            assert cost<=0
            for h in range(G):
                demand[h][x]+=Fraction(weight,den)*rewards[h][x][counts[h]]
                for t in range(1,counts[h]+1):demand_tail[h][x][t]+=Fraction(weight,den)
            local_count+=1
    group_mix=data['group_mixtures'];assert sorted(row['axis'] for row in group_mix)==list(range(G))
    scenario_count=0
    for group in group_mix:
        h=group['axis'];den=group['denominator'];assert type(den) is int and den>0
        assert sum(s['weight'] for s in group['scenarios'])==den
        for s in group['scenarios']:
            w=s['weight'];assert type(w) is int and w>0
            fixed=s['fixed'];assert tuple(map(tuple,fixed)) in leaves[h]
            used=dict(fixed);remaining=[d for d in shapes[h] if d not in used];base=int(1 in remaining)
            hits=[sum(p%d==a for d,a in fixed) for p in R]
            resource_rows=s['resource_distributions']
            assert sorted(row['resource'] for row in resource_rows)==[d for d in remaining if d!=1]
            capacities=[Fraction(0) for p in R]
            for row in resource_rows:
                d=row['resource'];residues=row['residues'];ps=list(map(rational,row['probabilities']))
                assert residues==sorted({p%d for p in R}) and len(ps)==len(residues)
                assert min(ps)>=0 and sum(ps)==1
                for x,p in enumerate(R):capacities[x]+=ps[residues.index(p%d)]
            assert len(s['count_distributions'])==len(R)
            for x,row in enumerate(s['count_distributions']):
                ps=list(map(rational,row));assert len(ps)==len(remaining)-base+1
                assert min(ps)>=0 and sum(ps)==1
                assert sum(j*p for j,p in enumerate(ps))<=capacities[x]
                payoff=sum(p*rewards[h][x][hits[x]+base+j] for j,p in enumerate(ps))
                supply[h][x]+=Fraction(w,den)*payoff
                for t in range(1,len(rewards[h][x])):
                    supply_tail[h][x][t]+=Fraction(w,den)*sum(p for j,p in enumerate(ps) if hits[x]+base+j>=t)
            scenario_count+=1
    excess=[[s-d for s,d in zip(sr,dr)] for sr,dr in zip(supply,demand)]
    assert len(data['coordinate_surpluses'])==G
    for h,row in enumerate(excess):
        assert len(data['coordinate_surpluses'][h])==len(R)
        assert row==list(map(rational,data['coordinate_surpluses'][h])) and min(row)>=0
    nonzero=[excess[h][x]/rewards[h][x][-1] for h in range(G) for x in range(len(R)) if rewards[h][x][-1]]
    margin=min(nonzero,default=Fraction(0));assert margin==rational(data['minimum_normalized_surplus'])
    absolute=min((excess[h][x] for h in range(G) for x in range(len(R)) if rewards[h][x][-1]),default=Fraction(0))
    deficits=[]
    for h in range(G):
        for x in range(len(R)):
            for t in range(1,len(rewards[h][x])):
                delta=supply_tail[h][x][t]-demand_tail[h][x][t]
                if delta<0: deficits.append((delta,h,x,t))
    deficits.sort()
    report=dict(verified=True,case_relaxation_barrier_verified=True,all_nonnegative_group_point_coefficients=True,
                coordinate_inequalities=G*len(R),active_local_states=local_count,active_fractional_scenarios=scenario_count,
                case_tree_nodes=sum(len(row['nodes']) for row in trees),case_tree_leaves=sum(map(len,leaves.values())),
                minimum_normalized_surplus=[margin.numerator,margin.denominator],
                minimum_active_coordinate_surplus=[absolute.numerator,absolute.denominator],
                integer_rounding_cannot_restore_positive_margin=absolute>=1,
                all_nonnegative_event_weights_barrier_verified=not deficits,
                deficient_thresholds=len(deficits),
                threshold_deficits=[dict(axis=h,point_index=x,threshold=t,
                    supply_minus_demand=[delta.numerator,delta.denominator],baseline_event_weight=-events[h][x][t])
                    for delta,h,x,t in deficits],
                actual_assignment_barrier_verified=False,covering_constructed=False,global_problem_settled=False,
                scope='Fixed event reward curves and the supplied residue-case count relaxations only; stronger compatibility or different reward curves are not excluded')
    if comparison_multipliers is not None:
        assert len(comparison_multipliers)==G
        comparison=Fraction(0)
        for h,axis in enumerate(comparison_multipliers):
            assert len(axis)==len(R)
            for x,row in enumerate(axis):
                assert len(row)==len(shapes[h])+1 and row[0]==0 and all(type(v) is int and v<=0 for v in row)
                comparison+=sum((-v)*(demand_tail[h][x][t]-supply_tail[h][x][t]) for t,v in enumerate(row) if t)
        report['comparison_local_minus_group']=[comparison.numerator,comparison.denominator]
        report['comparison_defeats_this_witness']=comparison>0
        report['comparison_scope']='Comparison of the same two mixtures under new curves; not an obstruction or a general event-weight barrier'
    return report


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True);ap.add_argument('--output',required=True);ap.add_argument('--compare')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());comparison=None
    if args.compare:
        other=json.loads(Path(args.compare).read_text());assert other['input']==data['input']
        comparison=other['result']['pricing_history'][-1]['event_multipliers']
    report=audit(data,comparison)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
