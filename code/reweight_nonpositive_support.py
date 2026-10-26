"""Optimize nonnegative combinations of audited point support inequalities.

The numerical LP concerns the explicitly supplied finite configuration pool.
Its rounded output is evaluated exactly; unrestricted bounds use separately
checkable elementary event majorants, never finite-pool maxima.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from event_score_majorant import upper_bound


def point_demands(result):
    refined=result.get('integer_support_refinement') or result.get('refined_support_slopes')
    if refined:
        assert refined['complete']
        assert [r['point_index'] for r in refined['points']]==list(range(len(result['positive_state_pricing']['point_multipliers'])))
        return [r['lower_numerator'] for r in refined['points']]
    return result['positive_state_pricing']['point_multipliers']


def barrier_witness(data,matrices,indices,normalizer):
    """Find supplied event scores dominating every fixed point demand.

    This certifies a limitation of nonnegative point reweighting only; it
    does not satisfy all count-marginal constraints of the covering model.
    """
    result=data['result'];base=result['pricing_history'][-1]
    beta=point_demands(result)
    lp=pywraplp.Solver.CreateSolver('GLOP')
    margin=lp.NumVar(-lp.infinity(),lp.infinity(),'margin')
    constraints=[lp.Constraint(b/base['point_weight_denominator']/normalizer,lp.infinity()) for b in beta]
    for c in constraints:c.SetCoefficient(margin,-1)
    variables=[]
    for h,selected in enumerate(indices):
        selected=sorted(selected or {0})
        row=[lp.NumVar(0,lp.infinity(),f'm{h}_{j}') for j in selected]
        lp.Add(sum(row)==1)
        for j,var in zip(selected,row):
            for x,c in enumerate(constraints):c.SetCoefficient(var,float(matrices[h][j,x])/normalizer)
        variables.append(list(zip(selected,row)))
    lp.Maximize(margin);lp.SetTimeLimit(10000)
    if lp.Solve()!=pywraplp.Solver.OPTIMAL:return None
    scale=10**9;groups=[];scores=[0]*len(beta)
    for h,row in enumerate(variables):
        values=[max(Fraction(0),Fraction(str(v.solution_value()))) for j,v in row];total=sum(values)
        weights=[int(v*scale//total) for v in values];weights[0]+=scale-sum(weights)
        assignments=[]
        for (j,var),w in zip(row,weights):
            if not w:continue
            assignments.append(dict(residues=result['group_columns'][h][j],weight=w))
            for x in range(len(beta)):scores[x]+=w*int(matrices[h][j,x])
        groups.append(dict(axis=h,assignments=assignments))
    excess=[v*base['point_weight_denominator']-b*scale for v,b in zip(scores,beta)]
    if min(excess)<0:return None
    return dict(denominator=scale,groups=groups,point_score_excess_numerators=excess,
                score_denominator=scale*base['denominator'],
                scope='Barrier only to nonnegative combinations of these fixed point event inequalities; not a covering or a full coherent witness')


def optimize(data,seconds=60):
    source,result = data['input'],data['result']
    base = result['pricing_history'][-1];support = result['positive_state_pricing']
    assert support['mode']=='nonpositive_support'
    beta=point_demands(result)
    points,n = source['core_residual'],len(source['high_primes'])
    event_scale,weight_den = base['event_scale'],base['point_weight_denominator']
    shapes = [[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(2**n-1)]
    matrices = []
    for h,(shape,configs) in enumerate(zip(shapes,result['group_columns'])):
        rewards = [[sum(e['weight'] for e in base['group_pricing'][h]['events']
                        if e['point']==x and e['threshold']<=k) for k in range(len(shape)+1)] for x in range(len(points))]
        values=[[rewards[x][sum(v%d==a for a,d in zip(residues,shape))]
                 for x,v in enumerate(points)] for residues in configs]
        assert all(0<=v<2**53 for row in values for v in row)
        matrices.append(np.asarray(values,dtype=np.float64))
    normalizer = max(1,max(beta)/weight_den)
    solver = pywraplp.Solver.CreateSolver('GLOP')
    variables = [solver.NumVar(0,1,f'p{x}') for x in points]
    bounds = [solver.NumVar(0,solver.infinity(),f'b{h}') for h in shapes]
    objective = solver.Objective();objective.SetMaximization()
    for var,b in zip(variables,beta):objective.SetCoefficient(var,b/weight_den/normalizer)
    for var in bounds:objective.SetCoefficient(var,-1)
    seen = [set() for h in shapes];history=[];start=time.monotonic();best=None
    while time.monotonic()-start<seconds:
        solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
        if solver.Solve()!=pywraplp.Solver.OPTIMAL:break
        weights = np.array([max(0,v.solution_value()) for v in variables])
        scores = [matrix@weights/normalizer for matrix in matrices]
        achieved = sum(b*w/weight_den/normalizer for b,w in zip(beta,weights))-sum(float(s.max()) for s in scores)
        if best is None or achieved>best[0]:best=(achieved,weights.copy())
        upper = objective.Value();bound_values=[v.solution_value() for v in bounds]
        added=0
        for h,values in enumerate(scores):
            for j in np.argsort(values)[-4:]:
                j=int(j)
                if j not in seen[h] and values[j]>bound_values[h]+1e-9:
                    seen[h].add(j);added+=1
                    c=solver.Constraint(-solver.infinity(),0)
                    c.SetCoefficient(bounds[h],-1)
                    for var,v in zip(variables,matrices[h][j]):
                        if v:c.SetCoefficient(var,float(v)/normalizer)
        history.append(dict(iteration=len(history)+1,upper=upper*normalizer/event_scale,
                            evaluated_value=achieved*normalizer/event_scale,added=added))
        if not added:break
    assert best is not None
    scale=10**6
    weights=[max(0,round(float(v)*scale)) for v in best[1]]
    groups=[]
    for h,(shape,configs) in enumerate(zip(shapes,result['group_columns'])):
        events=[dict(point=e['point'],threshold=e['threshold'],weight=e['weight']*weights[e['point']])
                for e in base['group_pricing'][h]['events'] if weights[e['point']]*e['weight']]
        maximum=0
        for residues in configs:
            counts=[sum(v%d==a for a,d in zip(residues,shape)) for v in points]
            maximum=max(maximum,sum(e['weight'] for e in events if counts[e['point']]>=e['threshold']))
        cap=upper_bound(points,shape,[(e['point'],e['threshold']) for e in events],[e['weight'] for e in events])
        groups.append(dict(events=events,pool_bound=maximum,capacity_certificate=cap))
    demand=sum(w*b for w,b in zip(weights,beta))
    barrier=barrier_witness(data,matrices,seen,normalizer) if best[0]<=1e-9 else None
    return dict(point_weights=weights,point_weight_denominator=scale,denominator=base['denominator']*scale,
                group_pricing=groups,pool_numerator=demand-weight_den*sum(g['pool_bound'] for g in groups),
                capacity_numerator=demand-weight_den*sum(g['capacity_certificate']['bound'] for g in groups),
                history=history,seconds=time.monotonic()-start,barrier_witness=barrier,
                scope='Nonnegative combination of point support inequalities; pool and unrestricted capacity bounds are separate')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=60)
    args=ap.parse_args();data=json.loads(Path(args.source).read_text())
    result=optimize(data,args.time);data['result']['support_reweighting']=result
    data['support_reweighting_source']=args.source
    Path(args.output).write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('group_pricing','history','barrier_witness')}))
    if result['barrier_witness']:
        b=result['barrier_witness']
        print(json.dumps(dict(barrier=True,minimum_score_excess_numerator=min(b['point_score_excess_numerators']),
                              denominator=b['score_denominator'])))
