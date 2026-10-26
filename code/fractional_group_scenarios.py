"""Exploratory resource/count LP scenarios for fixed residue cases.

Scenarios are fractional relaxation states, never actual residue assignments.
They constrain coefficient search; final certificates must use rational bounds.
"""
import argparse
import itertools
import json
from pathlib import Path
from ortools.linear_solver import pywraplp


def maximize(points, shapes, rewards, coefficients=None, fixed=(),pair_compatibility=False,return_dual=False,clique_compatibility=False):
    assert points and len(set(shapes)) == len(shapes)
    assert len(rewards) == len(points)
    assert all(len(row) == len(shapes)+1 and all(a <= b for a,b in zip(row,row[1:]))
               and min(row) >= 0 for row in rewards)
    assert len(dict(fixed)) == len(fixed)
    assert all(d in shapes and 0 <= a < d for d,a in fixed)
    coefficients = [1.0]*len(points) if coefficients is None else coefficients
    assert len(coefficients) == len(points) and min(coefficients) >= 0
    remaining = [d for d in shapes if d not in dict(fixed)]
    base = int(1 in remaining)
    hits = [sum(p%d == a for d,a in fixed) for p in points]
    solver = pywraplp.Solver.CreateSolver('GLOP')
    z = {}
    for d in remaining:
        if d == 1: continue
        z[d] = {a: solver.NumVar(0,solver.infinity(),f'z{d}_{a}') for a in sorted({p%d for p in points})}
        solver.Add(sum(z[d].values()) == 1)
    u = [[solver.NumVar(0,solver.infinity(),f'u{x}_{k}') for k in range(base,len(remaining)+1)]
         for x in range(len(points))]
    norm = max(1,max(c*row[-1] for c,row in zip(coefficients,rewards)))
    objective = solver.Objective(); objective.SetMaximization()
    mean_rows=[];pair_rows={}
    for x,p in enumerate(points):
        solver.Add(sum(u[x]) == 1)
        mean_rows.append(solver.Add(sum(j*v for j,v in enumerate(u[x])) <= sum(row[p%d] for d,row in z.items())))
        for j,v in enumerate(u[x]):
            objective.SetCoefficient(v,coefficients[x]*rewards[x][hits[x]+base+j]/norm)
    pairs=[];added=set();rounds=0;clique_rows={}
    if pair_compatibility:
        for x,y in itertools.combinations(range(len(points)),2):
            capacity=hits[x]+hits[y]+len(remaining)+sum(points[x]%d==points[y]%d for d in remaining)
            for t in range(hits[x]+base+1,hits[x]+len(remaining)+1):
                other=capacity+1-t
                if hits[y]+base<other<=hits[y]+len(remaining):pairs.append((x,y,t,other))
    while True:
        if solver.Solve() != pywraplp.Solver.OPTIMAL: raise RuntimeError('Fractional scenario LP failed')
        rounds+=1
        if not pairs and not clique_compatibility:break
        ps=[[v.solution_value() for v in row] for row in u]
        tails=[[sum(p for j,p in enumerate(row) if hits[x]+base+j>=t) for t in range(len(shapes)+1)] for x,row in enumerate(ps)]
        violations=[]
        for i,(x,y,t,other) in enumerate(pairs):
            gap=tails[x][t]+tails[y][other]-1
            if gap>1e-8:
                if i in added:raise RuntimeError('Numerical pair constraint violation')
                violations.append((gap,i))
        clique=None
        if clique_compatibility:
            from clique_case_constraints import strongest
            cleaned=[[v if v>1e-10 else 0. for v in row] for row in ps]
            score,selected=strongest(points,remaining,cleaned)
            if score>1+1e-8:
                clique=tuple(sorted((e['point'],e['threshold']) for e in selected))
                if clique in clique_rows:raise RuntimeError('Numerical clique constraint violation')
        if not violations and clique is None:break
        for gap,i in sorted(violations,reverse=True)[:128]:
            x,y,t,other=pairs[i];added.add(i)
            pair_rows[i]=solver.Add(sum(v for j,v in enumerate(u[x]) if hits[x]+base+j>=t)
                                   +sum(v for j,v in enumerate(u[y]) if hits[y]+base+j>=other)<=1)
        if clique is not None:
            clique_rows[clique]=solver.Add(sum(v for x,t in clique for j,v in enumerate(u[x]) if base+j>=t)<=1)
    distributions = [[v.solution_value() for v in row] for row in u]
    payoffs = [sum(prob*rewards[x][hits[x]+base+j] for j,prob in enumerate(row))
               for x,row in enumerate(distributions)]
    result=dict(fixed=list(fixed),payoffs=payoffs,
                pair_compatibility=pair_compatibility,pair_constraints_added=len(added),pair_constraints_checked=len(pairs),lp_solves=rounds,
                clique_compatibility=clique_compatibility,clique_constraints_added=len(clique_rows),
                objective=sum(c*s for c,s in zip(coefficients,payoffs)),
                count_base=base,fixed_hits=hits,count_distributions=distributions,
                resource_distributions=[dict(resource=d,residues=list(row),
                    probabilities=[v.solution_value() for v in row.values()]) for d,row in z.items()],
                scope='Numerical fractional relaxation state; not an actual assignment or proof')
    if return_dual:
        result['numerical_dual']=dict(slopes=[c.dual_value()*norm for c in mean_rows],
            pair_cuts=[dict(points=list(pairs[i][:2]),thresholds=list(pairs[i][2:]),multiplier=c.dual_value()*norm)
                       for i,c in sorted(pair_rows.items()) if c.dual_value()>0],
            clique_cuts=[dict(events=[dict(point=x,threshold=t+hits[x]) for x,t in clique],multiplier=c.dual_value()*norm)
                         for clique,c in sorted(clique_rows.items()) if c.dual_value()>0])
    return result


def generate(data, weighted, trees, per_group=8,pair_compatibility=False):
    assert data['input'] == weighted['input'] == trees['input']
    source = data['input']; base = data['result']['pricing_history'][-1]
    current = weighted['result']['pricing_history'][-1]
    n = len(source['high_primes']); result=[]
    for entry in trees['group_bounds']:
        h=entry['axis']; shapes=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        rewards=[[-sum(row[:k+1]) for k in range(len(shapes)+1)] for row in base['event_multipliers'][h]]
        coefficients=[]
        for old,new in zip(base['event_multipliers'][h],current['event_multipliers'][h]):
            pairs=[(a,b) for a,b in zip(old,new) if a]
            ratio=pairs[0][1]/pairs[0][0] if pairs else 0
            assert all(b == a*ratio for a,b in zip(old,new))
            coefficients.append(ratio)
        leaves=[(j,row) for j,row in enumerate(entry['certificate']['nodes']) if 'children' not in row]
        leaves.sort(key=lambda pair:pair[1]['effective_bound'],reverse=True)
        selected=leaves if per_group==0 else leaves[:per_group]
        rows=[]
        for j,node in selected:
            scenario=maximize(source['core_residual'],shapes,rewards,coefficients,node['fixed'],pair_compatibility)
            scenario['node_index']=j;rows.append(scenario)
        result.append(dict(axis=h,scenarios=rows))
        print(json.dumps(dict(axis=h,scenarios=len(rows),maximum=max(row['objective'] for row in rows))),flush=True)
    return dict(input=source,baseline_event_multipliers=base['event_multipliers'],group_scenarios=result,
                scope='Optimization cuts from selected case relaxations; not unrestricted score certificates')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--weights',required=True)
    ap.add_argument('--trees',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--per-group',type=int,default=8)
    ap.add_argument('--pairs',action='store_true')
    args=ap.parse_args()
    result=generate(*[json.loads(Path(p).read_text()) for p in (args.source,args.weights,args.trees)],args.per_group,args.pairs)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
