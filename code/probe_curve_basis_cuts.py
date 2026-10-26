"""Close a supplied finite curve-cone LP by pricing all known rows each round.

This is the coefficient dual of the mixture probe. The complete finite pools
remain in memory, while only violated inequalities enter the solver. An
OPTIMAL status requires all supplied rows to pass the final numerical check;
neither that status nor a numerical margin is an arithmetic proof.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from probe_curve_basis_barrier import prepare


def probe(data,seed,scenarios,basis,seconds=180,group_pool='all',extra_points=None,
          solver_name='GLOP',deduplicate=True,batch=4,tolerance=1e-8,warm_start=None,checkpoint=None):
    start=time.monotonic()
    p=prepare(data,seed,scenarios,basis,group_pool,extra_points,deduplicate)
    R,G,curves,caps,states=p['points'],p['groups'],p['curves'],p['caps'],p['states']
    selected,scores=p['selected'],p['scores']
    coordinates=[(h,x,b) for h in range(G) for x in range(len(R)) for b in range(curves[h].shape[1])
                 if caps[h][x] and curves[h][x,b,-1]]
    point_columns=[[i for i,(_,y,_) in enumerate(coordinates) if y==x] for x in range(len(R))]
    group_columns=[[i for i,(j,_,_) in enumerate(coordinates) if j==h] for h in range(G)]
    point_matrices=[np.column_stack([curves[h][x,b,rows[:,h]] for h,_,b in [coordinates[i] for i in point_columns[x]]])
                    if point_columns[x] else np.zeros((len(rows),0)) for x,rows in enumerate(states)]
    group_matrices=[np.column_stack([scores[h][:,x,b] for _,x,b in [coordinates[i] for i in group_columns[h]]])
                    if group_columns[h] else np.zeros((len(rows),0)) for h,rows in enumerate(selected)]
    assert solver_name in ('GLOP','CLP') and batch>0 and tolerance>0
    solver=pywraplp.Solver.CreateSolver(solver_name);assert solver is not None
    if solver_name=='GLOP':solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    w=[solver.NumVar(0,solver.infinity(),f'w{i}') for i in range(len(coordinates))]
    solver.Add(sum(w)==1)
    alpha=[solver.NumVar(-solver.infinity(),solver.infinity(),f'a{x}') for x in R]
    beta=[solver.NumVar(-solver.infinity(),solver.infinity(),f'b{h}') for h in range(G)]
    objective=solver.Objective();objective.SetMinimization()
    for v in alpha:objective.SetCoefficient(v,-1)
    for v in beta:objective.SetCoefficient(v,1)
    point_rows=[{} for x in R];group_rows=[{} for h in range(G)]

    def add(kind,index,j):
        local=kind=='point'
        rows=point_rows if local else group_rows
        if j in rows[index]:return
        matrix=point_matrices[index] if local else group_matrices[index]
        columns=point_columns[index] if local else group_columns[index]
        row=solver.Constraint(-solver.infinity(),0) if local else solver.Constraint(0,solver.infinity())
        row.SetCoefficient(alpha[index] if local else beta[index],1)
        for k in np.flatnonzero(matrix[j]):row.SetCoefficient(w[columns[int(k)]],-float(matrix[j,k]))
        rows[index][j]=row

    for x in range(len(R)):add('point',x,0)
    for h in range(G):add('group',h,0)
    if warm_start:
        assert warm_start['input']==p['source']
        warm_points=warm_start.get('active_point_counts')
        if warm_points is None:
            by_point={r['point_index']:[s['counts'] for s in r['states']] for r in warm_start.get('point_mixtures',[])}
            warm_points=[by_point.get(x,[]) for x in range(len(R))]
        warm_groups=warm_start.get('active_group_states')
        if warm_groups is None:
            by_group={r['axis']:r['states'] for r in warm_start.get('group_mixtures',[])}
            warm_groups=[by_group.get(h,[]) for h in range(G)]
        for x,rows in enumerate(warm_points):
            indices={tuple(map(int,state)):j for j,state in enumerate(states[x])}
            for state in rows:
                if tuple(state) in indices:add('point',x,indices[tuple(state)])
        for h,rows in enumerate(warm_groups):
            indices={(s['kind'],s['index']):j for j,s in enumerate(selected[h])}
            for state in rows:
                key=(state['kind'],state['index'])
                if key in indices:add('group',h,indices[key])
    parameters=pywraplp.MPSolverParameters();parameters.SetIntegerParam(parameters.LP_ALGORITHM,parameters.DUAL)
    parameters.SetDoubleParam(parameters.PRIMAL_TOLERANCE,min(1e-10,tolerance/10))
    parameters.SetDoubleParam(parameters.DUAL_TOLERANCE,min(1e-10,tolerance/10))
    deadline=time.monotonic()+seconds
    history=[];snapshot=None;closed=False;terminal='TIME_LIMIT'
    print(json.dumps(dict(stage='prepared',coordinates=len(coordinates),local_states=sum(map(len,states)),
                          group_states=sum(map(len,selected)),solver=solver_name,seconds=time.monotonic()-start)),flush=True)
    while time.monotonic()<deadline:
        solver.SetTimeLimit(max(1,int((deadline-time.monotonic())*1000)))
        status=solver.Solve(parameters)
        if status!=solver.OPTIMAL:
            terminal={solver.FEASIBLE:'FEASIBLE',solver.INFEASIBLE:'INFEASIBLE',solver.ABNORMAL:'ABNORMAL',
                      solver.NOT_SOLVED:'NOT_SOLVED',solver.UNBOUNDED:'UNBOUNDED'}.get(status,str(status))
            break
        values=np.asarray([v.solution_value() for v in w]);a=[v.solution_value() for v in alpha];b=[v.solution_value() for v in beta]
        point_scores=[matrix@values[cols] for matrix,cols in zip(point_matrices,point_columns)]
        group_scores=[matrix@values[cols] for matrix,cols in zip(group_matrices,group_columns)]
        point_violation=max([0.]+[a[x]-float(cost.min()) for x,cost in enumerate(point_scores)])
        group_violation=max([0.]+[float(cost.max())-b[h] for h,cost in enumerate(group_scores)])
        priced_margin=sum(float(cost.max()) for cost in group_scores)-sum(float(cost.min()) for cost in point_scores)
        snapshot=dict(input=p['source'],numerical_margin=objective.Value(),
            finite_pool_priced_margin=priced_margin,finite_pool_numerical_duality_gap=priced_margin-objective.Value(),
            point_mixtures=[dict(point_index=x,states=[dict(counts=list(map(int,states[x][j])),weight=-row.dual_value())
                            for j,row in rows.items() if -row.dual_value()>1e-12]) for x,rows in enumerate(point_rows)],
            group_mixtures=[dict(axis=h,states=[dict(**selected[h][j],weight=row.dual_value())
                            for j,row in rows.items() if row.dual_value()>1e-12]) for h,rows in enumerate(group_rows)],
            coordinate_dual_weights=[dict(axis=h,point_index=x,basis=c,weight=float(v))
                                     for (h,x,c),v in zip(coordinates,values) if v>1e-12],
            coordinate_dual_weight_sum=float(values.sum()))
        closed=max(point_violation,group_violation)<=tolerance
        row=dict(iteration=len(history)+1,margin=objective.Value(),point_violation=point_violation,
                 group_violation=group_violation,priced_margin=priced_margin,point_constraints=sum(map(len,point_rows)),
                 group_constraints=sum(map(len,group_rows)),seconds=time.monotonic()-start)
        history.append(row)
        retained_violation=max([0.]+[a[x]-float(point_scores[x][j]) for x,rows in enumerate(point_rows) for j in rows]
                               +[float(group_scores[h][j])-b[h] for h,rows in enumerate(group_rows) for j in rows])
        if retained_violation>tolerance:
            terminal='NUMERICAL_CONSTRAINT_VIOLATION';closed=False
            snapshot['retained_constraint_violation']=retained_violation
            print(json.dumps(row|dict(retained_constraint_violation=retained_violation)),flush=True)
            break
        if closed:
            terminal='OPTIMAL';print(json.dumps(row|dict(closed=True)),flush=True);break
        added=0
        for x,cost in enumerate(point_scores):
            indices=np.flatnonzero(cost<a[x]-tolerance)
            for j in indices[np.argsort(cost[indices])[:batch]]:
                assert int(j) not in point_rows[x], 'Retained point constraint violated'
                add('point',x,int(j));added+=1
        for h,cost in enumerate(group_scores):
            indices=np.flatnonzero(cost>b[h]+tolerance)
            for j in indices[np.argsort(-cost[indices])[:batch]]:
                assert int(j) not in group_rows[h], 'Retained group constraint violated'
                add('group',h,int(j));added+=1
        assert added
        if len(history)%10==0 or len(history)==1:
            print(json.dumps(row|dict(added=added)),flush=True)
            if checkpoint:
                saved=dict(**snapshot,status='FEASIBLE',finite_pool_closed=False,terminal_solver_status='RUNNING',
                           history=history,active_point_counts=[[list(map(int,states[x][j])) for j in rows] for x,rows in enumerate(point_rows)],
                           active_group_states=[[selected[h][j] for j in rows] for h,rows in enumerate(group_rows)],verified=False)
                temporary=Path(checkpoint).with_suffix('.tmp');temporary.write_text(json.dumps(saved,separators=(',',':'))+'\n')
                temporary.replace(checkpoint)
    result=snapshot or dict(input=p['source'])
    result.update(status='OPTIMAL' if closed else ('FEASIBLE' if snapshot else terminal),terminal_solver_status=terminal,
                  finite_pool_closed=closed,history=history,solver=solver_name,method='coefficient_constraint_generation',
                  coordinates=len(coordinates),local_state_count=sum(map(len,states)),group_state_count=sum(map(len,selected)),
                  input_local_state_count=p['input_local_count'],input_group_state_count=p['input_group_count'],
                  deduplicated=deduplicate,maximum_point_violation=history[-1]['point_violation'] if history else None,
                  maximum_group_violation=history[-1]['group_violation'] if history else None,
                  seconds=time.monotonic()-start,verified=False,
                  scope='Numerical closure of the supplied finite state pools only. Exact local pricing, actual-group bounds, and independent arithmetic audits remain necessary.')
    result['active_point_counts']=[[list(map(int,states[x][j])) for j in rows] for x,rows in enumerate(point_rows)]
    result['active_group_states']=[[selected[h][j] for j in rows] for h,rows in enumerate(group_rows)]
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','seed','scenarios','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--basis',action='append',default=[]);ap.add_argument('--time',type=float,default=180)
    ap.add_argument('--group-pool',choices=['active','scenarios','all'],default='all');ap.add_argument('--point-pool')
    ap.add_argument('--solver',choices=['GLOP','CLP'],default='GLOP');ap.add_argument('--batch',type=int,default=4)
    ap.add_argument('--warm-start')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());seed=json.loads(Path(args.seed).read_text())
    scenarios=json.loads(Path(args.scenarios).read_text());basis=[json.loads(Path(p).read_text()) for p in args.basis]
    extra=None
    if args.point_pool:
        points=json.loads(Path(args.point_pool).read_text());assert points['input']==data['input'];extra=points['result']['point_columns']
    warm=json.loads(Path(args.warm_start).read_text()) if args.warm_start else None
    result=probe(data,seed,scenarios,basis,args.time,args.group_pool,extra,args.solver,True,args.batch,warm_start=warm,checkpoint=args.output)
    result['sources']=dict(source=args.source,seed=args.seed,scenarios=args.scenarios,basis=args.basis,
                           group_pool=args.group_pool,point_pool=args.point_pool)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('point_mixtures','group_mixtures','coordinate_dual_weights','input','sources','history','active_point_counts','active_group_states')}),flush=True)
