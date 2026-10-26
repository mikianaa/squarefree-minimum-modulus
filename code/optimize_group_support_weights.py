"""Optimize separate group penalties at each point, with integer state pricing.

Numerical master bounds are exploratory. Each saved candidate includes a
rounded integer certificate, full-domain point pricing, and separate finite
pool and elementary unrestricted group bounds for independent auditing.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from multitail_lp import LocalStates
from refine_integer_support import refine_point
from event_score_majorant import upper_bound


def optimize(data,seconds=180,min_weight=0,seed=None,checkpoint=None,scenarios=None):
    assert 0<=min_weight<=1
    source=data['input'];previous=data['result'];base=previous['pricing_history'][-1]
    R,n=source['core_residual'],len(source['high_primes']);G=2**n-1
    model=LocalStates(source['high_primes'],source['high_pattern'])
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    rewards=[[np.asarray([-sum(base['event_multipliers'][h][x][:k+1]) for k in range(len(shapes[h])+1)],dtype=np.int64)
              for x in range(len(R))] for h in range(G)]
    matrices=[]
    for h,(shape,configs) in enumerate(zip(shapes,previous['group_columns'])):
        matrix=np.asarray([[rewards[h][x][sum(v%d==a for a,d in zip(residues,shape))]
                            for x,v in enumerate(R)] for residues in configs],dtype=np.int64)
        assert int(matrix.max())<2**53;matrices.append(matrix)
    separation_matrices=list(matrices)
    if scenarios:
        assert scenarios['input']==source and scenarios['baseline_event_multipliers']==base['event_multipliers']
        assert sorted(row['axis'] for row in scenarios['group_scenarios'])==list(range(G))
        for row in scenarios['group_scenarios']:
            h=row['axis'];extra=np.asarray([s['payoffs'] for s in row['scenarios']],dtype=float)
            assert extra.ndim==2 and extra.shape[1]==len(R) and np.all(np.isfinite(extra)) and np.all(extra>=0)
            separation_matrices[h]=np.vstack((matrices[h],extra))
    normalizer=max(1,max(sum(int(rewards[h][x][-1]) for h in range(G)) for x in range(len(R))))
    solver=pywraplp.Solver.CreateSolver('GLOP')
    alpha=[[solver.NumVar(min_weight,1,f'a{h}_{x}') for x in range(len(R))] for h in range(G)]
    beta=[solver.NumVar(0,solver.infinity(),f'p{x}') for x in range(len(R))]
    bound=[solver.NumVar(0,solver.infinity(),f'b{h}') for h in range(G)]
    objective=solver.Objective();objective.SetMaximization()
    for var in beta:objective.SetCoefficient(var,1)
    for var in bound:objective.SetCoefficient(var,-1)
    point_seen=[set() for x in R];group_seen=[set() for h in range(G)]
    def add_state(x,state):
        state=tuple(state)
        limits=[min(14,q-1-(b==R[x]%2)) for q,b in zip(model.high,model.pattern)]+[16]*(G-n)
        assert len(state)==len(limits) and all(type(k) is int and 0<=k<=limit for k,limit in zip(state,limits))
        assert len(state)==G and model.cost(R[x]%2,state)<=0
        if state in point_seen[x]:return 0
        point_seen[x].add(state)
        c=solver.Constraint(-solver.infinity(),0);c.SetCoefficient(beta[x],1)
        for h,k in enumerate(state):c.SetCoefficient(alpha[h][x],-int(rewards[h][x][k])/normalizer)
        return 1
    for row in previous['integer_support_refinement']['points']:add_state(row['point_index'],row['attaining_state'])
    for x,point in enumerate(R):
        add_state(x,[min(14,q-1-(b==point%2)) for q,b in zip(model.high,model.pattern)]+[16]*(G-n))
    if seed:
        for row in seed['support']:
            x=row['point_index']
            if model.cost(R[x]%2,row['counts'])<=0:add_state(x,row['counts'])
    start=time.monotonic();history=[];best=None;best_value=None;scale=1000
    event_scale=base['event_scale']*scale
    coefficients=sorted({k*event_scale//8 for k in range(9)}|{k*event_scale for k in (2,4,8,16)})
    while time.monotonic()-start<seconds:
        # Close all catalogue constraints before spending time on point pricing.
        while time.monotonic()-start<seconds:
            solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
            if solver.Solve()!=pywraplp.Solver.OPTIMAL:return best
            values=np.asarray([[v.solution_value() for v in row] for row in alpha])
            upper=objective.Value()*normalizer/base['event_scale']
            group_bounds=[v.solution_value() for v in bound];added=0
            for h,matrix in enumerate(separation_matrices):
                scores=matrix@values[h]/normalizer
                for j in np.argsort(scores)[-4:]:
                    j=int(j)
                    if j not in group_seen[h] and scores[j]>group_bounds[h]+1e-9:
                        group_seen[h].add(j);added+=1
                        c=solver.Constraint(-solver.infinity(),0);c.SetCoefficient(bound[h],-1)
                        for x,v in enumerate(matrix[j]):c.SetCoefficient(alpha[h][x],float(v)/normalizer)
            if not added:break
        else:break
        integers=[[max(0,round(float(v)*scale)) for v in row] for row in values]
        multipliers=[[[v*integers[h][x] for v in base['event_multipliers'][h][x]] for x in range(len(R))] for h in range(G)]
        raw_minima=[];coarse_minima=[];points=[];new_states=0
        for x,point in enumerate(R):
            potentials=[np.cumsum(np.asarray(g[x],dtype=np.int64)) for g in multipliers]
            priced=refine_point(model,point%2,[list(map(int,row)) for row in potentials],coefficients,event_scale)
            priced['point_index']=x;points.append(priced)
            raw_minima.append(model.price(point%2,potentials,event_scale)[0][0])
            coarse_minima.append(priced['initial_envelope_minimum'])
            new_states+=add_state(x,priced['attaining_state'])
            for special in priced['refined_planes']:
                if 'tail_counts' in special:new_states+=add_state(x,special['plane']+special['tail_counts'])
        groups=[];pool=[];capacities=[]
        for h,shape in enumerate(shapes):
            events=[dict(point=x,threshold=t,weight=-v) for x,row in enumerate(multipliers[h]) for t,v in enumerate(row) if t and v]
            cap=upper_bound(R,shape,[(e['point'],e['threshold']) for e in events],[e['weight'] for e in events])
            # The integer dot product evaluates every saved assignment exactly.
            assert int(matrices[h].max())*sum(integers[h])<2**63
            maximum=int((matrices[h]@np.asarray(integers[h],dtype=np.int64)).max())
            pool.append(maximum);capacities.append(cap['bound'])
            groups.append(dict(events=events,capacity_certificate=cap))
        c=dict(point_weight_denominator=1,event_scale=event_scale,denominator=event_scale,
               point_weights=[1]*len(R),event_multipliers=multipliers,point_multipliers=raw_minima,
               group_pricing=groups,pool_sizes=list(map(len,previous['group_columns'])),pool_upper_bounds=pool,
               pool_numerator=sum(raw_minima)-sum(pool),capacity_numerator=sum(raw_minima)-sum(capacities))
        coarse=dict(mode='nonpositive_support',denominator=event_scale,point_multipliers=coarse_minima,
                    point_slope_numerators=[coefficients for x in R],pool_numerator=sum(coarse_minima)-sum(pool),
                    unrestricted_capacity_numerator=sum(coarse_minima)-sum(capacities))
        exact=sum(row['lower_numerator'] for row in points)-sum(pool)
        search_margin=(sum(row['lower_numerator'] for row in points)-sum(
            float((matrix@np.asarray(integers[h])).max()) for h,matrix in enumerate(separation_matrices)))/event_scale
        history.append(dict(iteration=len(history)+1,numerical_master_upper=upper,
                            exact_pool_lower=exact/event_scale,new_states=new_states,
                            numerical_scenario_margin=search_margin,
                            point_states=sum(map(len,point_seen)),group_constraints=sum(map(len,group_seen)),seconds=time.monotonic()-start))
        candidate=dict(input=source,result=dict(group_columns=previous['group_columns'],pricing_history=[c],
                 positive_state_pricing=coarse,integer_support_refinement=dict(points=points,complete=True),
                 group_support_optimization=dict(weights=integers,weight_denominator=scale,minimum_weight=min_weight,
                    exact_pool_numerator=exact,denominator=event_scale,history=list(history),
                    exact_capacity_numerator=sum(row['lower_numerator'] for row in points)-sum(capacities),
                    scenario_count=sum(len(row['scenarios']) for row in scenarios['group_scenarios']) if scenarios else 0,
                    numerical_scenario_margin=search_margin,
                    scope='Master may include numerical fractional scenarios; proof pool contains only actual assignments; unrestricted bounds are separate'),
                 support_seed=[dict(point_index=x,counts=list(state)) for x,states in enumerate(point_seen) for state in sorted(states)]))
        selection_value=search_margin if scenarios else exact
        if best_value is None or selection_value>best_value:best,best_value=candidate,selection_value
        best['result']['group_support_optimization']['all_iterations']=list(history)
        best['result']['support_seed']=candidate['result']['support_seed']
        if checkpoint:Path(checkpoint).write_text(json.dumps(best,separators=(',',':'))+'\n')
        print(json.dumps(history[-1]),flush=True)
        if (exact>0 and not scenarios) or not new_states:break
    return best


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=180);ap.add_argument('--minimum-weight',type=float,default=0)
    ap.add_argument('--seed')
    ap.add_argument('--scenarios')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text())
    seed=None
    if args.seed:
        saved=json.loads(Path(args.seed).read_text());assert saved['input']==data['input']
        previous=saved['result']
        seed=dict(support=previous['support_seed']) if 'support_seed' in previous else previous['coherent_primal_witness']
    scenarios=json.loads(Path(args.scenarios).read_text()) if args.scenarios else None
    result=optimize(data,args.time,args.minimum_weight,seed,args.output,scenarios)
    if result is None:raise RuntimeError('No complete pricing iteration; no certificate saved')
