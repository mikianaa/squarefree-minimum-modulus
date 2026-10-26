"""Find an exact mixture barrier for fixed event curves and case relaxations.

This cannot prove that the covering conjecture or the full method fails.
The group mixtures are fractional case states, not ordinary assignments.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from ortools.linear_solver import pywraplp


def probabilities(values):
    values=[Fraction(max(0,float(v))).limit_denominator(10**6) for v in values]
    total=sum(values);assert total>0
    return [v/total for v in values]


def pair(v):return [v.numerator,v.denominator]


def rational_scenario(points,shapes,scenario,rewards):
    fixed=scenario['fixed'];used=dict(fixed)
    remaining=[d for d in shapes if d not in used];base=int(1 in remaining)
    hits=[sum(p%d==a for d,a in fixed) for p in points]
    resources=[];supply=[Fraction(0)]*len(points)
    for row in scenario['resource_distributions']:
        d=row['resource'];ps=probabilities(row['probabilities'])
        table=dict(zip(row['residues'],ps))
        for x,p in enumerate(points):supply[x]+=table.get(p%d,0)
        resources.append(dict(resource=d,residues=row['residues'],probabilities=list(map(pair,ps))))
    counts=[];payoffs=[]
    for x,row in enumerate(scenario['count_distributions']):
        ps=probabilities(row);demand=sum(k*v for k,v in enumerate(ps))
        # Rounding can introduce tiny demand excess. Transfer non-base mass
        # to the base count, which restores feasibility and only lowers reward.
        if demand>supply[x]:
            factor=supply[x]/demand
            ps=[1-sum(v*factor for v in ps[1:])]+[v*factor for v in ps[1:]]
        counts.append(list(map(pair,ps)))
        payoffs.append(sum(v*rewards[x][hits[x]+base+k] for k,v in enumerate(ps)))
    return dict(fixed=fixed,resource_distributions=resources,count_distributions=counts),payoffs


def rounded_mixture(variables):
    rows=[(j,round(max(0,var.solution_value())*10**10)) for j,var in enumerate(variables)]
    rows=[(j,w) for j,w in rows if w]
    return rows,sum(w for j,w in rows)


def generate(data,seed,scenarios,trees,seconds=120):
    source=data['input'];assert source==seed['input']==scenarios['input']==trees['input']
    events=data['result']['pricing_history'][-1]['event_multipliers']
    assert events==scenarios['baseline_event_multipliers']
    R=source['core_residual'];n=len(source['high_primes']);G=2**n-1
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    rewards=[[[-sum(row[:k+1]) for k in range(len(row))] for row in axis] for axis in events]
    states=[[] for p in R]
    for row in seed['result']['support_seed']:states[row['point_index']].append(row['counts'])
    groups={row['axis']:row['scenarios'] for row in scenarios['group_scenarios']}
    solver=pywraplp.Solver.CreateSolver('GLOP');solver.SetTimeLimit(int(seconds*1000))
    margin=solver.NumVar(-solver.infinity(),1,'margin');solver.Maximize(margin)
    local=[[solver.NumVar(0,1,f'l{x}_{j}') for j in range(len(rows))] for x,rows in enumerate(states)]
    global_rows=[[solver.NumVar(0,1,f'g{h}_{j}') for j in range(len(groups[h]))] for h in range(G)]
    for row in local+global_rows:solver.Add(sum(row)==1)
    constraints={}
    for h in range(G):
        for x in range(len(R)):
            scale=rewards[h][x][-1]
            if not scale:continue
            c=solver.Constraint(0,solver.infinity());c.SetCoefficient(margin,-1);constraints[h,x]=c
            for state,var in zip(states[x],local[x]):c.SetCoefficient(var,-rewards[h][x][state[h]]/scale)
    for h in range(G):
        for row,var in zip(groups[h],global_rows[h]):
            for x,score in enumerate(row['payoffs']):
                if (h,x) in constraints:constraints[h,x].SetCoefficient(var,score/rewards[h][x][-1])
    status=solver.Solve();report=dict(status={solver.OPTIMAL:'OPTIMAL',solver.FEASIBLE:'FEASIBLE',solver.INFEASIBLE:'INFEASIBLE'}.get(status,str(status)))
    if status not in (pywraplp.Solver.OPTIMAL,pywraplp.Solver.FEASIBLE):return report
    report['numerical_margin']=margin.solution_value();print(json.dumps(report),flush=True)
    if margin.solution_value()<=0:return report
    demands=[[Fraction(0) for p in R] for h in range(G)];supplies=[[Fraction(0) for p in R] for h in range(G)]
    point_mix=[];group_mix=[]
    for x,variables in enumerate(local):
        rows,den=rounded_mixture(variables);assert den>0
        point_mix.append(dict(point_index=x,denominator=den,states=[dict(weight=w,counts=states[x][j]) for j,w in rows]))
        for j,w in rows:
            for h in range(G):demands[h][x]+=Fraction(w,den)*rewards[h][x][states[x][j][h]]
    for h,variables in enumerate(global_rows):
        rows,den=rounded_mixture(variables);assert den>0;entries=[]
        for j,w in rows:
            row,payoffs=rational_scenario(R,shapes[h],groups[h][j],rewards[h]);row['weight']=w;entries.append(row)
            for x,v in enumerate(payoffs):supplies[h][x]+=Fraction(w,den)*v
        group_mix.append(dict(axis=h,denominator=den,scenarios=entries))
    surpluses=[[s-d for s,d in zip(supply,demand)] for supply,demand in zip(supplies,demands)]
    assert all(v>=0 for row in surpluses for v in row),'Rational reconstruction did not retain dominance'
    normalized=min(surpluses[h][x]/rewards[h][x][-1] for h,x in constraints)
    topology=[dict(axis=row['axis'],nodes=[{k:v for k,v in node.items() if k in ('fixed','children','resource')}
                                          for node in row['certificate']['nodes']]) for row in trees['group_bounds']]
    return dict(input=source,baseline_event_multipliers=events,point_mixtures=point_mix,group_mixtures=group_mix,
                case_trees=topology,coordinate_surpluses=[[pair(v) for v in row] for row in surpluses],
                minimum_normalized_surplus=pair(normalized),numerical_search=report,
                scope='All nonnegative group/point coefficients of these fixed event curves with these case relaxations; not actual group assignments or a covering')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','seed','scenarios','trees','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--time',type=float,default=120)
    args=ap.parse_args();result=generate(*[json.loads(Path(getattr(args,k)).read_text()) for k in ('source','seed','scenarios','trees')],args.time)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
