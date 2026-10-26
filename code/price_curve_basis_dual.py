"""Round a curve-cone LP dual direction and price all local integer states.

LP dual weights are search suggestions. Only the resulting integer certificate
and its independent audit can establish an obstruction for the stated scope.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from optimize_threshold_support import curve_directions,scenario_tails
from cached_support_slopes import coefficients,safe_maximum
from refine_integer_support import refine_point
from multitail_lp import LocalStates
from event_score_majorant import optimized_upper_bound


def price(data,seed,scenarios,basis,probe,precision=10**12,use_cliques=False):
    start=time.monotonic();source=data['input'];assert source==seed['input']==scenarios['input']==probe['input']
    assert all(entry['input']==source for entry in basis)
    assert positive_precision(precision) and probe['status'] in ('OPTIMAL','FEASIBLE')
    assert probe.get('coordinate_dual_weights'),'Numerical coefficient weights are required; mixture status alone is insufficient'
    R=source['core_residual'];n=len(source['high_primes']);G=2**n-1;scale=10**6
    model=LocalStates(source['high_primes'],source['high_pattern'])
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    directions=curve_directions(data['result']['pricing_history'][-1],basis,shapes,len(R))
    weights=[np.zeros(axis.shape[:2]) for axis in directions]
    for row in probe['coordinate_dual_weights']:
        h,x,b=row['axis'],row['point_index'],row['basis'];w=row['weight']
        assert w>=0 and np.isfinite(w);weights[h][x,b]+=w
    multipliers=[]
    for axis,w in zip(directions,weights):
        events=np.einsum('pb,pbt->pt',w,axis)
        multipliers.append([[0]+[-max(0,round(float(v)*precision)) for v in row] for row in events])
    old_result=seed['result'];history=old_result.get('pricing_history',[])
    old_points={r['point_index']:r for r in old_result.get('integer_support_refinement',{}).get('points',[])}
    old_pricing=old_result.get('positive_state_pricing',{})
    defaults=sorted({k*scale//8 for k in range(9)}|{k*scale for k in (2,4,8,16)})
    retained=[set() for p in R]
    for row in old_result['support_seed']:retained[row['point_index']].add(tuple(row['counts']))
    points=[];raw=[];coarse=[];point_slopes=[]
    for x,p in enumerate(R):
        potentials=[list(map(int,np.cumsum(axis[x],dtype=np.int64))) for axis in multipliers]
        old_mass=(history[-1].get('point_weight_denominator',1)*sum(-sum(axis[x]) for axis in history[-1]['event_multipliers'])) if history else 0
        ratio=(sum(-row[-1] for row in potentials),old_mass) if old_mass else None
        chosen=coefficients(defaults,old_points.get(x,{}).get('refined_planes',[]),old_pricing.get('denominator',scale),scale,
            12,safe_maximum(model,p%2,potentials),old_pricing.get('point_slope_numerators',[[] for p in R])[x],ratio)
        row=refine_point(model,p%2,potentials,chosen,scale,1);row['point_index']=x;points.append(row);point_slopes.append(chosen)
        coarse.append(row['initial_envelope_minimum'])
        raw.append(model.price(p%2,[np.asarray(a,dtype=np.int64) for a in potentials],scale)[0][0])
        retained[x].add(tuple(row['attaining_state']))
        for special in row['refined_planes']:
            if 'tail_counts' in special:retained[x].add(tuple(special['plane']+special['tail_counts']))
        for extra in row.get('diverse_support',[]):retained[x].add(tuple(extra['counts']))
        print(json.dumps(dict(stage='exact_local_pricing',point=x,minimum=row['lower_numerator']/scale,
            exceptional_planes=len(row['refined_planes']),integer_dp_calls=row['integer_dp_calls'])),flush=True)
    groups={row['axis']:row['scenarios'] for row in scenarios['group_scenarios']}
    prior_rows=old_result.get('threshold_support_optimization',{}).get('group_constraint_rows',[[] for h in range(G)])
    group_rows=[set(rows) for rows in prior_rows];pool=[];capacity=[];prices=[];scenario_bounds=[]
    for h,shape in enumerate(shapes):
        columns=data['result']['group_columns'][h]
        counts=np.asarray([[sum(p%d==a for a,d in zip(column,shape)) for p in R] for column in columns],dtype=np.int64)
        rewards=-np.cumsum(np.asarray(multipliers[h],dtype=np.int64),axis=1)
        scores=sum(rewards[x,counts[:,x]] for x in range(len(R)));pool.append(int(scores.max()))
        group_rows[h].update(map(int,np.argsort(scores)[-4:]))
        fractional=np.asarray([scenario_tails(R,shape,s) for s in groups[h]]).reshape(-1,len(R),len(shape))
        sampled=float(pool[-1])
        if len(fractional):
            scores=np.einsum('npt,pt->n',fractional,-np.asarray(multipliers[h])[:,1:],optimize=True)
            sampled=max(sampled,float(scores.max()));group_rows[h].update(len(columns)+int(j) for j in np.argsort(scores)[-4:])
        scenario_bounds.append(sampled)
        events=[dict(point=x,threshold=t,weight=-v) for x,line in enumerate(multipliers[h]) for t,v in enumerate(line) if t and v]
        majorant=optimized_upper_bound
        if use_cliques:
            from clique_event_majorant import upper_bound as majorant
        bound=majorant(R,shape,[(e['point'],e['threshold']) for e in events],[e['weight'] for e in events])
        capacity.append(bound['bound']);prices.append(dict(events=events,capacity_certificate=bound))
    local=sum(row['lower_numerator'] for row in points);mass=sum(-sum(row) for axis in multipliers for row in axis)
    certificate=dict(point_weight_denominator=1,event_scale=scale,denominator=scale,point_weights=[1]*len(R),
        event_multipliers=multipliers,point_multipliers=raw,group_pricing=prices,pool_sizes=list(map(len,data['result']['group_columns'])),
        pool_upper_bounds=pool,pool_numerator=sum(raw)-sum(pool),capacity_numerator=sum(raw)-sum(capacity))
    summary=dict(exact_pool_numerator=local-sum(pool),exact_capacity_numerator=local-sum(capacity),denominator=scale,
        normalized_pool_gap=(local-sum(pool))/mass if mass else None,numerical_scenario_margin=(local-sum(scenario_bounds))/scale,
        total_event_weight=mass,point_states=sum(map(len,retained)),group_constraints=sum(map(len,group_rows)),seconds=time.monotonic()-start)
    summary['group_bound_method']='clique_affine_or_fallback' if use_cliques else 'ordinary_affine'
    print(json.dumps(summary),flush=True)
    return dict(input=source,result=dict(group_columns=data['result']['group_columns'],pricing_history=[certificate],
        positive_state_pricing=dict(mode='nonpositive_support',denominator=scale,point_multipliers=coarse,point_slope_numerators=point_slopes,
            pool_numerator=sum(coarse)-sum(pool),unrestricted_capacity_numerator=sum(coarse)-sum(capacity)),
        integer_support_refinement=dict(points=points,complete=True),curve_basis_dual_pricing=summary,
        threshold_support_optimization=dict(group_constraint_rows=[sorted(rows) for rows in group_rows]),
        support_seed=[dict(point_index=x,counts=list(state)) for x,states in enumerate(retained) for state in sorted(states)]))


def positive_precision(v):return type(v) is int and 0<v<10**14


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','seed','scenarios','probe','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--basis',action='append',default=[]);ap.add_argument('--precision',type=int,default=10**12)
    ap.add_argument('--cliques',action='store_true')
    args=ap.parse_args();numerical=json.loads(Path(args.probe).read_text())
    result=price(*[json.loads(Path(getattr(args,k)).read_text()) for k in ('source','seed','scenarios')],
        [json.loads(Path(p).read_text()) for p in args.basis],numerical,args.precision,args.cliques)
    result['curve_dual_provenance']=dict(source=args.source,seed=args.seed,scenarios=args.scenarios,probe=args.probe,basis=args.basis,
        cliques=args.cliques,probe_status=numerical['status'],
        scope='Rounded numerical coefficient suggestion (optimality is not assumed), followed by complete exact local pricing and exact actual-pool maxima. Fractional scenario scores are search diagnostics only.')
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
