"""Search arbitrary nonnegative threshold curves with exact local support pricing.

Each curve's total mass is capped by its original total. Fractional scenarios
are search constraints only; exact certificates keep actual assignments apart.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from multitail_lp import LocalStates
from refine_integer_support import refine_point,diverse_support_states
from event_score_majorant import optimized_upper_bound
from cached_support_slopes import coefficients as cached_coefficients,safe_maximum


def scenario_tails(points,shape,row):
    fixed=row['fixed'];remaining=len(shape)-len(fixed)
    base=int(1 in shape and 1 not in dict(fixed))
    hits=[sum(p%d==a for d,a in fixed) for p in points]
    assert len(row['count_distributions'])==len(points)
    tails=[]
    for x,probs in enumerate(row['count_distributions']):
        assert len(probs)==remaining-base+1 and min(probs)>=-1e-8 and abs(sum(probs)-1)<1e-8
        tails.append([sum(p for j,p in enumerate(probs) if hits[x]+base+j>=t) for t in range(1,len(shape)+1)])
    return np.asarray(tails,dtype=float)


def valid_state_pool(source,rows,share_parity=False):
    n=len(source['high_primes']);G=2**n-1;R=source['core_residual']
    import itertools
    groups=[S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    separated=[[] for p in R]
    for row in rows:
        x=row['point_index'];state=row['counts']
        assert type(x) is int and 0<=x<len(R) and len(state)==G and all(type(k) is int for k in state)
        separated[x].append(state)
    result=[]
    for x,states in enumerate(separated):
        array=np.asarray(states,dtype=np.int64).reshape(-1,G)
        D=np.asarray([q-1-int(b==R[x]%2) for q,b in zip(source['high_primes'],source['high_pattern'])])
        limits=np.asarray([min(14,int(d)) for d in D]+[16]*(G-n))
        valid=np.all((array>=0)&(array<=limits),axis=1);array=array[valid]
        free=D-array[:,:n];cost=np.prod(free,axis=1)
        for h,S in enumerate(groups[n:],n):cost-=array[:,h]*np.prod(free[:,[i for i in range(n) if i not in S]],axis=1)
        result.append(np.unique(array[cost<=0],axis=0))
    if share_parity:
        # Local domain and F depend on the core point only through parity.
        # Never share across parities: both plane caps and free counts differ.
        shared=[]
        for parity in (0,1):
            arrays=[array for x,array in enumerate(result) if R[x]%2==parity]
            shared.append(np.unique(np.concatenate(arrays),axis=0) if arrays else np.empty((0,G),dtype=np.int64))
        return [shared[p%2] for p in R]
    return result


def curve_directions(base,basis,shapes,points):
    if basis is None:return [np.broadcast_to(np.eye(len(shape)),(points,len(shape),len(shape))).copy() for shape in shapes]
    sources=[base]+[entry['result']['pricing_history'][-1] for entry in basis]
    result=[]
    for h,shape in enumerate(shapes):
        axis=[]
        for x in range(points):
            rows=[]
            for source in sources:
                row=source['event_multipliers'][h][x]
                assert len(row)==len(shape)+1 and row[0]==0 and all(type(v) is int and v<=0 for v in row)
                total=-sum(row)
                rows.append([-v/total if total else 0 for v in row[1:]])
            axis.append(rows)
        result.append(np.asarray(axis,dtype=float))
    return result


def optimize(data,seed,scenarios=None,seconds=180,checkpoint=None,iterations=0,point_pool=None,basis=None,diverse=0,share_parity=False,cached_slopes=0,dual_simplex=False,compact_master=False):
    source=data['input'];assert seed['input']==source
    previous=data['result'];base=previous['pricing_history'][-1]
    R=source['core_residual'];n=len(source['high_primes']);G=2**n-1;model=LocalStates(source['high_primes'],source['high_pattern'])
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    caps=np.asarray([[-sum(row) for row in axis] for axis in base['event_multipliers']],dtype=np.int64)
    assert caps.shape==(G,len(R)) and np.all(caps>=0)
    normalizer=max(1,int(caps.sum(axis=0).max()))
    if basis is not None:assert all(entry['input']==source for entry in basis)
    directions=curve_directions(base,basis,shapes,len(R))
    cumulative=[np.concatenate((np.zeros((len(R),axis.shape[1],1)),np.cumsum(axis,axis=2)),axis=2) for axis in directions]
    counts=[];matrices=[];fractional_tails=[]
    extra={} if scenarios is None else {row['axis']:row['scenarios'] for row in scenarios['group_scenarios']}
    if scenarios is not None:assert scenarios['input']==source and set(extra)==set(range(G))
    for h,shape in enumerate(shapes):
        actual=np.asarray([[sum(p%d==a for a,d in zip(row,shape)) for p in R] for row in previous['group_columns'][h]],dtype=np.int64)
        counts.append(actual)
        matrix=(actual[:,:,None]>=np.arange(1,len(shape)+1)[None,None,:]).astype(float)
        additional=np.asarray([scenario_tails(R,shape,s) for s in extra.get(h,[])]).reshape(-1,len(R),len(shape))
        fractional_tails.append(additional)
        if len(additional):matrix=np.concatenate((matrix,additional))
        matrix*=caps[h][None,:,None]/normalizer
        projected=np.einsum('npt,pjt->npj',matrix,directions[h],optimize=True)
        matrices.append(projected.reshape(len(projected),-1))
    matrix_hashes=[hashlib.sha256(str(matrix.shape).encode()+matrix.tobytes()).hexdigest() for matrix in matrices]
    solver=pywraplp.Solver.CreateSolver('GLOP')
    if dual_simplex:assert solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    variables=[[[solver.NumVar(0,1,f'w{h}_{x}_{j}') for j in range(axis.shape[1])] for x in range(len(R))] for h,axis in enumerate(directions)]
    for axis in variables:
        for row in axis:solver.Add(sum(row)<=1)
    flat=[[v for row in axis for v in row] for axis in variables]
    beta=[solver.NumVar(0,solver.infinity(),f'p{x}') for x in range(len(R))]
    bounds=[solver.NumVar(0,solver.infinity(),f'b{h}') for h in range(G)]
    solver.Maximize(sum(beta)-sum(bounds))
    seen=[set() for p in R];retained=[set() for p in R];group_seen=[set() for h in range(G)]
    def add_group(h,j):
        assert type(j) is int and 0<=j<len(matrices[h])
        if j in group_seen[h]:return 0
        group_seen[h].add(j);c=solver.Constraint(-solver.infinity(),0);c.SetCoefficient(bounds[h],-1)
        for k in np.flatnonzero(matrices[h][j]):c.SetCoefficient(flat[h][int(k)],float(matrices[h][j,k]))
        return 1
    def remember_state(x,state):
        state=tuple(state)
        limits=[min(14,q-1-int(b==R[x]%2)) for q,b in zip(model.high,model.pattern)]+[16]*(G-n)
        assert len(state)==G and all(type(k) is int and 0<=k<=limit for k,limit in zip(state,limits)) and model.cost(R[x]%2,state)<=0
        retained[x].add(state)
        return state
    def add_state(x,state):
        state=remember_state(x,state)
        if state in seen[x]:return 0
        seen[x].add(state);c=solver.Constraint(-solver.infinity(),0);c.SetCoefficient(beta[x],1)
        for h,k in enumerate(state):
            for j,var in enumerate(variables[h][x]):
                value=cumulative[h][x,j,k]
                if value:c.SetCoefficient(var,-float(value)*int(caps[h,x])/normalizer)
        return 1
    for row in seed['result']['support_seed']:
        (remember_state if compact_master else add_state)(row['point_index'],row['counts'])
    if compact_master:
        for row in seed['result'].get('integer_support_refinement',{}).get('points',[]):
            add_state(row['point_index'],row['attaining_state'])
    for x,p in enumerate(R):add_state(x,[min(14,q-1-int(b==p%2)) for q,b in zip(model.high,model.pattern)]+[16]*(G-n))
    known=valid_state_pool(source,point_pool) if point_pool is not None and not share_parity else []
    warm=seed['result'].get('threshold_support_optimization',{});reused=0
    if warm.get('group_matrix_hashes')==matrix_hashes:
        assert len(warm['group_constraint_rows'])==G
        for h,rows in enumerate(warm['group_constraint_rows']):reused+=sum(add_group(h,j) for j in rows)
    enrichment=0
    if diverse:
        seed_base=seed['result']['pricing_history'][-1];seed_slopes=seed['result']['positive_state_pricing']
        for x,p in enumerate(R):
            potentials=[[seed_base['point_weight_denominator']*int(v) for v in np.cumsum(axis[x],dtype=np.int64)] for axis in seed_base['event_multipliers']]
            for row in diverse_support_states(model,p%2,potentials,seed_slopes['point_slope_numerators'][x],seed_base['denominator'],diverse):
                enrichment+=add_state(x,row['counts'])
        print(json.dumps(dict(stage='seed_enrichment',new_active_constraints=enrichment,reused_group_constraints=reused,
            point_states=sum(map(len,retained)),active_point_constraints=sum(map(len,seen)))),flush=True)
    def refresh_shared_states():
        return valid_state_pool(source,(point_pool or [])+[dict(point_index=x,counts=list(state))
            for x,states in enumerate(retained) for state in states],share_parity)
    if share_parity or compact_master:known=refresh_shared_states()
    start=time.monotonic();history=[];best=None;best_value=None;scale=1000;event_scale=base['event_scale']*scale
    slopes=sorted({k*event_scale//8 for k in range(9)}|{k*event_scale for k in (2,4,8,16)})
    slope_seed=seed
    while time.monotonic()-start<seconds and (not iterations or len(history)<iterations):
        values=None;closed=False;separation_round=0
        while time.monotonic()-start<seconds:
            solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
            status=solver.Solve()
            if status!=solver.OPTIMAL:break
            values=[np.asarray([v.solution_value() for v in axis]) for axis in flat]
            numerical=solver.Objective().Value()*normalizer/base['event_scale']
            added=0;bound_values=[v.solution_value() for v in bounds];beta_values=[v.solution_value() for v in beta]
            for h,matrix in enumerate(matrices):
                scores=matrix@values[h]
                for j in np.argsort(scores)[-4:]:
                    j=int(j)
                    if j not in group_seen[h] and scores[j]>bound_values[h]+1e-9:
                        added+=add_group(h,j)
            if known:
                curves=[np.einsum('pj,pjk->pk',v.reshape(len(R),directions[h].shape[1]),cumulative[h])
                        *caps[h,:,None]/normalizer for h,v in enumerate(values)]
                for x,states in enumerate(known):
                    costs=sum(curves[h][x,states[:,h]] for h in range(G))
                    for j in np.argsort(costs)[:4]:
                        if costs[j]<beta_values[x]-1e-9:added+=add_state(x,list(map(int,states[j])))
            separation_round+=1
            if separation_round==1 or separation_round%10==0 or not added:
                print(json.dumps(dict(stage='group_separation',round=separation_round,numerical_upper=numerical,
                    added=added,constraints=sum(map(len,group_seen)),seconds=time.monotonic()-start)),flush=True)
            if not added:closed=True;break
        if values is None:return best
        multipliers=[]
        for h,shape in enumerate(shapes):
            rows=np.einsum('pj,pjt->pt',values[h].reshape(len(R),directions[h].shape[1]),directions[h])
            multipliers.append([[0]+[-max(0,round(float(v)*int(caps[h,x])*scale)) for v in row] for x,row in enumerate(rows)])
        points=[];raw=[];coarse=[];point_slopes=[];added_states=0
        old_points={row['point_index']:row for row in slope_seed['result'].get('integer_support_refinement',{}).get('points',[])}
        old_pricing=slope_seed['result'].get('positive_state_pricing',{})
        old_history=slope_seed['result'].get('pricing_history',[])
        for x,p in enumerate(R):
            potentials=[list(map(int,np.cumsum(axis[x],dtype=np.int64))) for axis in multipliers]
            old_mass=(old_history[-1].get('point_weight_denominator',1)*sum(-sum(axis[x])
                for axis in old_history[-1]['event_multipliers'])) if old_history else 0
            ratio=(sum(-row[-1] for row in potentials),old_mass) if old_mass else None
            selected=cached_coefficients(slopes,old_points.get(x,{}).get('refined_planes',[]),
                old_pricing.get('denominator',event_scale),event_scale,cached_slopes,
                safe_maximum(model,p%2,potentials),old_pricing.get('point_slope_numerators',[[] for p in R])[x],ratio)
            point_slopes.append(selected)
            priced=refine_point(model,p%2,potentials,selected,event_scale,diverse);priced['point_index']=x;points.append(priced)
            raw.append(model.price(p%2,[np.asarray(row,dtype=np.int64) for row in potentials],event_scale)[0][0]);coarse.append(priced['initial_envelope_minimum'])
            added_states+=add_state(x,priced['attaining_state'])
            for row in priced['refined_planes']:
                if 'tail_counts' in row:added_states+=add_state(x,row['plane']+row['tail_counts'])
            for row in priced.get('diverse_support',[]):added_states+=add_state(x,row['counts'])
            print(json.dumps(dict(iteration=len(history)+1,point=x,minimum=priced['lower_numerator']/event_scale,
                                 exceptions=len(priced['refined_planes']),integer_dp_calls=priced['integer_dp_calls'])),flush=True)
        group_prices=[];pool=[];capacity=[];scenario_bounds=[]
        for h,shape in enumerate(shapes):
            events=[dict(point=x,threshold=t,weight=-v) for x,row in enumerate(multipliers[h]) for t,v in enumerate(row) if t and v]
            cap=optimized_upper_bound(R,shape,[(e['point'],e['threshold']) for e in events],[e['weight'] for e in events])
            rewards=-np.cumsum(np.asarray(multipliers[h],dtype=np.int64),axis=1)
            actual_scores=sum(rewards[x,counts[h][:,x]] for x in range(len(R)))
            maximum=int(actual_scores.max());pool.append(maximum);capacity.append(cap['bound'])
            group_prices.append(dict(events=events,capacity_certificate=cap))
            relaxed=maximum
            if len(fractional_tails[h]):
                relaxed=max(relaxed,float(np.einsum('npt,pt->n',fractional_tails[h],-np.asarray(multipliers[h])[:,1:],optimize=True).max()))
            scenario_bounds.append(relaxed)
        local=sum(row['lower_numerator'] for row in points)
        score=(local-sum(scenario_bounds))/event_scale
        history.append(dict(iteration=len(history)+1,numerical_master_upper=numerical,
                            sampled_group_constraints_closed=closed,
                            exact_pool_lower=(local-sum(pool))/event_scale,numerical_scenario_margin=score,
                            point_states=sum(map(len,retained)),active_point_constraints=sum(map(len,seen)),new_states=added_states,
                            group_constraints=sum(map(len,group_seen)),seconds=time.monotonic()-start))
        certificate=dict(point_weight_denominator=1,event_scale=event_scale,denominator=event_scale,
                         point_weights=[1]*len(R),event_multipliers=multipliers,point_multipliers=raw,group_pricing=group_prices,
                         pool_sizes=list(map(len,previous['group_columns'])),pool_upper_bounds=pool,
                         pool_numerator=sum(raw)-sum(pool),capacity_numerator=sum(raw)-sum(capacity))
        result=dict(input=source,result=dict(group_columns=previous['group_columns'],pricing_history=[certificate],
            positive_state_pricing=dict(mode='nonpositive_support',denominator=event_scale,point_multipliers=coarse,
                point_slope_numerators=point_slopes,pool_numerator=sum(coarse)-sum(pool),unrestricted_capacity_numerator=sum(coarse)-sum(capacity)),
            integer_support_refinement=dict(points=points,complete=True),
            threshold_support_optimization=dict(history=list(history),curve_mass_caps=caps.tolist(),
                curve_basis_size=None if basis is None else 1+len(basis),
                diverse_per_slice=diverse,initial_enrichment_states=enrichment,reused_group_constraints=reused,
                shared_parity_states=share_parity,
                cached_slope_limit=cached_slopes,
                cached_slope_scaling='point_reward_mass',
                dual_simplex=dual_simplex,
                compact_master=compact_master,
                known_nonpositive_states=sum(map(len,known)),
                exact_pool_numerator=local-sum(pool),exact_capacity_numerator=local-sum(capacity),denominator=event_scale,
                numerical_scenario_margin=score,scope='Arbitrary threshold curves; fractional scenarios guide search but are excluded from the actual-assignment proof pool'),
            support_seed=[dict(point_index=x,counts=list(state)) for x,states in enumerate(retained) for state in sorted(states)]))
        if best_value is None or score>best_value:best,best_value=result,score
        best['result']['threshold_support_optimization']['all_iterations']=list(history)
        best['result']['threshold_support_optimization']['group_matrix_hashes']=matrix_hashes
        best['result']['threshold_support_optimization']['group_constraint_rows']=[sorted(rows) for rows in group_seen]
        best['result']['support_seed']=result['result']['support_seed']
        if checkpoint:Path(checkpoint).write_text(json.dumps(best,separators=(',',':'))+'\n')
        print(json.dumps(history[-1]),flush=True)
        slope_seed=result
        if not added_states:break
        if share_parity or compact_master:known=refresh_shared_states()
    return best


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','seed','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--scenarios');ap.add_argument('--time',type=float,default=180);ap.add_argument('--iterations',type=int,default=0);ap.add_argument('--point-pool');ap.add_argument('--basis',action='append');ap.add_argument('--diverse',type=int,default=0)
    ap.add_argument('--share-parity',action='store_true')
    ap.add_argument('--cached-slopes',type=int,default=0)
    ap.add_argument('--dual-simplex',action='store_true')
    ap.add_argument('--compact-master',action='store_true')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());seed=json.loads(Path(args.seed).read_text())
    scenarios=json.loads(Path(args.scenarios).read_text()) if args.scenarios else None
    point_pool=None
    if args.point_pool:
        saved=json.loads(Path(args.point_pool).read_text());assert saved['input']==data['input']
        point_pool=saved['result']['point_columns']
    basis=[json.loads(Path(p).read_text()) for p in args.basis] if args.basis else None
    answer=optimize(data,seed,scenarios,args.time,args.output,args.iterations,point_pool,basis,args.diverse,args.share_parity,args.cached_slopes,args.dual_simplex,args.compact_master)
    if answer is None:raise RuntimeError('No complete pricing iteration; no certificate saved')
