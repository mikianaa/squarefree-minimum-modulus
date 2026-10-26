"""Numerical mixture test for a specified cone of threshold reward curves.

A positive numerical margin is only a candidate for rational reconstruction.
Fractional group scenarios are not residue assignments or covering systems.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from optimize_threshold_support import curve_directions,scenario_tails,valid_state_pool


def projection_representatives(matrix):
    """Keep one original row for each identical floating LP column, without rounding."""
    array=np.asarray(matrix,dtype=float)
    assert array.ndim==2 and np.all(np.isfinite(array))
    _,indices=np.unique(array,axis=0,return_index=True)
    return np.sort(indices)


def prepare(data,seed,scenarios,basis,group_pool='active',extra_points=None,deduplicate=False):
    source=data['input'];assert source==seed['input']==scenarios['input']
    assert all(entry['input']==source for entry in basis)
    R=source['core_residual'];n=len(source['high_primes']);G=2**n-1
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    base=data['result']['pricing_history'][-1]
    directions=curve_directions(base,basis,shapes,len(R))
    curves=[np.concatenate((np.zeros((len(R),a.shape[1],1)),np.cumsum(a,axis=2)),axis=2) for a in directions]
    caps=[[-sum(row) for row in axis] for axis in base['event_multipliers']]
    states=valid_state_pool(source,seed['result']['support_seed']+(extra_points or []))
    input_local_count=sum(map(len,states))
    if deduplicate:
        for x,array in enumerate(states):
            projection=np.concatenate([curves[h][x][:,array[:,h]].T for h in range(G)],axis=1)
            states[x]=array[projection_representatives(projection)]
    groups={row['axis']:row['scenarios'] for row in scenarios['group_scenarios']}
    selected=[];scores=[];input_group_count=0
    cached=seed['result'].get('threshold_support_optimization',{}).get('group_constraint_rows')
    for h,shape in enumerate(shapes):
        actual=data['result']['group_columns'][h];M=len(actual);total=M+len(groups[h])
        indices=range(total) if group_pool=='all' else (range(M,total) if group_pool=='scenarios' else cached[h])
        rows=[];values=[]
        for j in indices:
            assert type(j) is int and 0<=j<total
            if j<M:
                counts=[sum(p%d==a for a,d in zip(actual[j],shape)) for p in R]
                value=np.asarray([curves[h][x,:,k] for x,k in enumerate(counts)])
                row=dict(kind='actual',index=j)
            else:
                value=np.einsum('pt,pbt->pb',scenario_tails(R,shape,groups[h][j-M]),directions[h])
                row=dict(kind='scenario',index=j-M)
            rows.append(row);values.append(value)
        assert rows
        input_group_count+=len(rows);values=np.asarray(values)
        if deduplicate:
            indices=projection_representatives(values.reshape(len(values),-1))
            rows=[rows[j] for j in indices];values=values[indices]
        selected.append(rows);scores.append(values)
    return dict(source=source,points=R,groups=G,curves=curves,caps=caps,states=states,
                selected=selected,scores=scores,input_local_count=input_local_count,input_group_count=input_group_count)


def probe(data,seed,scenarios,basis,seconds=120,group_pool='active',extra_points=None,dual_simplex=True,
          solver_name='GLOP',deduplicate=False):
    start=time.monotonic()
    p=prepare(data,seed,scenarios,basis,group_pool,extra_points,deduplicate)
    source=p['source'];R=p['points'];G=p['groups'];curves=p['curves'];caps=p['caps'];states=p['states']
    selected=p['selected'];scores=p['scores'];input_local_count=p['input_local_count'];input_group_count=p['input_group_count']
    assert solver_name in ('GLOP','CLP')
    solver=pywraplp.Solver.CreateSolver(solver_name)
    assert solver is not None,'Requested LP backend is unavailable'
    solver.SetTimeLimit(int(seconds*1000))
    if solver_name=='GLOP' and dual_simplex:solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    margin=solver.NumVar(-solver.infinity(),1,'margin');solver.Maximize(margin)
    local=[[solver.NumVar(0,solver.infinity(),f'l{x}_{j}') for j in range(len(rows))] for x,rows in enumerate(states)]
    global_rows=[[solver.NumVar(0,solver.infinity(),f'g{h}_{j}') for j in range(len(rows))] for h,rows in enumerate(selected)]
    for variables in local+global_rows:solver.Add(sum(variables)==1)
    coordinates=0;coordinate_rows=[]
    for h in range(G):
        for x in range(len(R)):
            for b in range(curves[h].shape[1]):
                if not caps[h][x] or not curves[h][x,b,-1]:continue
                c=solver.Constraint(0,solver.infinity());c.SetCoefficient(margin,-1);coordinates+=1
                coordinate_rows.append((h,x,b,c))
                for j in np.flatnonzero(curves[h][x,b,states[x][:,h]]):
                    c.SetCoefficient(local[x][int(j)],-float(curves[h][x,b,states[x][j,h]]))
                for j in np.flatnonzero(scores[h][:,x,b]):c.SetCoefficient(global_rows[h][int(j)],float(scores[h][j,x,b]))
    print(json.dumps(dict(stage='built',coordinates=coordinates,local_states=sum(map(len,states)),
        group_states=sum(map(len,selected)),input_local_states=input_local_count,input_group_states=input_group_count,
        solver=solver_name,deduplicated=deduplicate,seconds=time.monotonic()-start)),flush=True)
    parameters=pywraplp.MPSolverParameters()
    parameters.SetIntegerParam(parameters.LP_ALGORITHM,parameters.DUAL if dual_simplex else parameters.PRIMAL)
    status=solver.Solve(parameters)
    statuses={solver.OPTIMAL:'OPTIMAL',solver.FEASIBLE:'FEASIBLE',solver.INFEASIBLE:'INFEASIBLE',
              solver.UNBOUNDED:'UNBOUNDED',solver.ABNORMAL:'ABNORMAL',solver.MODEL_INVALID:'MODEL_INVALID',solver.NOT_SOLVED:'NOT_SOLVED'}
    result=dict(input=source,status=statuses.get(status,str(status)),dual_simplex=dual_simplex,solver=solver_name,
        coordinates=coordinates,local_state_count=sum(map(len,states)),group_state_count=sum(map(len,selected)),
        input_local_state_count=input_local_count,input_group_state_count=input_group_count,deduplicated=deduplicate,
        seconds=time.monotonic()-start,verified=False,
        scope='Numerical mixture test for the supplied curve cone and group-state sample; no covering or global method limitation.')
    if status in (solver.OPTIMAL,solver.FEASIBLE):
        result['numerical_margin']=margin.solution_value()
        result['point_mixtures']=[dict(point_index=x,states=[dict(counts=list(map(int,state)),weight=v.solution_value())
            for state,v in zip(states[x],variables) if v.solution_value()>1e-12]) for x,variables in enumerate(local)]
        result['group_mixtures']=[dict(axis=h,states=[dict(**row,weight=v.solution_value())
            for row,v in zip(selected[h],variables) if v.solution_value()>1e-12]) for h,variables in enumerate(global_rows)]
        if status==solver.OPTIMAL:
            result['coordinate_dual_weights']=[dict(axis=h,point_index=x,basis=b,weight=-c.dual_value())
                for h,x,b,c in coordinate_rows if -c.dual_value()>1e-12]
            result['coordinate_dual_weight_sum']=sum(row['weight'] for row in result['coordinate_dual_weights'])
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','seed','scenarios','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--basis',action='append',default=[]);ap.add_argument('--time',type=float,default=120)
    ap.add_argument('--group-pool',choices=['active','scenarios','all'],default='active');ap.add_argument('--point-pool')
    ap.add_argument('--primal-simplex',action='store_true')
    ap.add_argument('--solver',choices=['GLOP','CLP'],default='GLOP')
    ap.add_argument('--deduplicate',action='store_true')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());seed=json.loads(Path(args.seed).read_text())
    scenarios=json.loads(Path(args.scenarios).read_text());basis=[json.loads(Path(p).read_text()) for p in args.basis]
    extra=None
    if args.point_pool:
        points=json.loads(Path(args.point_pool).read_text());assert points['input']==data['input'];extra=points['result']['point_columns']
    result=probe(data,seed,scenarios,basis,args.time,args.group_pool,extra,not args.primal_simplex,args.solver,args.deduplicate)
    result['sources']=dict(source=args.source,seed=args.seed,scenarios=args.scenarios,basis=args.basis,
                           group_pool=args.group_pool,point_pool=args.point_pool)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('point_mixtures','group_mixtures','coordinate_dual_weights','input','sources')}),flush=True)
