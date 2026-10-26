"""Conditional-tail LP retaining actual core projection identities.

Each core modulus in each higher-prime group chooses a probability
distribution over its nonempty core slices. Local hit-count distributions
must be dominated pointwise by those choices. This is still a relaxation.
"""
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from multitail_lp import LocalStates


def rational_primal(source,model,columns,resources,rows,rank_rows=(),pair_rows=(),subset_rows=(),group_rows=()):
    R = source['core_residual']
    scale = 10**9
    projections = []
    rhs = [[0]*len(model.groups) for x in R]
    for resource in resources:
        h,d,options = resource
        values = [max(Fraction(0),Fraction(str(v.solution_value()))) for a,indices,v in options]
        total = sum(values)
        assert total>0
        weights = [int(value*scale//total) for value in values]
        weights[0] += scale-sum(weights)
        assignments = []
        for (a,indices,v),w in zip(options,weights):
            if w:
                assignments.append(dict(residue=a,weight=w))
                for x in indices:
                    rhs[x][h] += w
        projections.append(dict(axis=h,modulus=d,assignments=assignments))
    support = [(x,state,max(0,math.floor(v.solution_value()*scale))) for x,state,v in columns if v.solution_value()>0]
    support = [(x,state,w) for x,state,w in support if w]
    ratio = Fraction(1)
    for x in range(len(R)):
        used = sum(w for y,state,w in support if y==x)
        if used:
            ratio = min(ratio,Fraction(scale,used))
        for h in range(len(model.groups)):
            used = sum(w*state[h] for y,state,w in support if y==x)
            if used:
                ratio = min(ratio,Fraction(rhs[x][h],used))
    for h,t0,t1,bound,c in rows:
        used = sum(w*max(0,state[h]-(t0 if R[x]%2==0 else t1)) for x,state,w in support)
        if used:
            ratio = min(ratio,Fraction(scale*bound,used))
    for h,t0,t1,w0,w1,bound,c in rank_rows:
        used = sum(w*(w0 if R[x]%2==0 else w1)*int(state[h]>=(t0 if R[x]%2==0 else t1)) for x,state,w in support)
        if used:
            ratio = min(ratio,Fraction(scale*bound,used))
    for h,x,y,t,u,c in pair_rows:
        used = sum(w for z,state,w in support if (z==x and state[h]>=t) or (z==y and state[h]>=u))
        if used:
            ratio = min(ratio,Fraction(scale,used))
    for h,points,thresholds,c in subset_rows:
        selected = dict(zip(points,thresholds))
        used = sum(w for x,state,w in support if x in selected and state[h]>=selected[x])
        if used:
            ratio = min(ratio,Fraction(scale*(len(points)-1),used))
    for cut,c in group_rows:
        used = sum(w*event['weight'] for x,state,w in support for event in cut['events']
                   if x==event['point'] and state[cut['axis']]>=event['threshold'])
        if used:
            ratio = min(ratio,Fraction(scale*cut['bound'],used))
    support = [(x,state,w*ratio.numerator//ratio.denominator) for x,state,w in support]
    support = [(x,state,w) for x,state,w in support if w]
    for x in range(len(R)):
        missing = scale-sum(w for y,state,w in support if y==x)
        assert missing>=0
        if missing:
            support.append((x,(0,)*len(model.groups),missing))
    objective = sum(w*model.cost(R[x]%2,state) for x,state,w in support)
    return dict(denominator=scale,objective_numerator=objective,projections=projections,
                support=[dict(point_index=x,counts=list(state),weight=w) for x,state,w in support],
                scope='rational point-identity relaxation witness, not an actual covering family')


def solve(source,seconds=180,batch=2,target_lower=1,seed=None,rank_cuts=False,pair_cuts=False,warm_start=None,triple_cuts=False,subset_cuts=False,group_cuts=()):
    start = time.monotonic()
    R = source['core_residual']
    assert len(R)==len(set(R)) and all(0<=v<210 for v in R)
    assert source['counts']==[sum(v%2==e for v in R) for e in (0,1)]
    model = LocalStates(source['high_primes'],source['high_pattern'])
    n,groups = model.n,model.groups
    core_moduli = [d for d in range(3,211) if 210%d==0]
    shapes = [core_moduli if h<n else [1,2]+core_moduli for h in range(len(groups))]
    H = np.array(source['H_upper'],dtype=np.int64)
    a,b = np.indices(H.shape)
    Hs = [H,H+a+b+np.maximum(a,b)]
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    objective = solver.Objective()
    objective.SetMinimization()
    mass = [solver.Constraint(1,1) for x in R]
    links = [[solver.Constraint(-solver.infinity(),0) for h in groups] for x in R]
    resources = []
    for h,moduli in enumerate(shapes):
        for d in moduli:
            constraint = solver.Constraint(1,1)
            options = []
            for residue in sorted({v%d for v in R}):
                indices = [i for i,v in enumerate(R) if v%d==residue]
                variable = solver.NumVar(0,solver.infinity(),f'p{h}_{d}_{residue}')
                constraint.SetCoefficient(variable,1)
                for x in indices:
                    links[x][h].SetCoefficient(variable,-1)
                options.append((residue,indices,variable))
            resources.append((h,d,options))
    rows = []
    lookup = [[[[] for v in range(len(shape)+1)] for shape in shapes] for e in (0,1)]
    for h,shape in enumerate(shapes):
        limit = len(shape)
        for t0,t1 in itertools.product(range(limit+1),repeat=2):
            if t0==t1==limit:
                continue
            bound = int(np.max(Hs[int(h>=n)]-t0*a-t1*b))
            constraint = solver.Constraint(-solver.infinity(),bound)
            rows.append((h,t0,t1,bound,constraint))
            for e,t in enumerate((t0,t1)):
                for v in range(t+1,limit+1):
                    lookup[e][h][v].append((constraint,v-t))
    seen = [set() for x in R]
    columns = []
    point_lookup = [[[[] for v in range(len(shape)+1)] for shape in shapes] for x in R]
    def add(x,state):
        if state in seen[x]:
            return 0
        seen[x].add(state)
        variable = solver.NumVar(0,solver.infinity(),f'z{x}_{len(seen[x])}')
        columns.append((x,state,variable))
        mass[x].SetCoefficient(variable,1)
        for h,v in enumerate(state):
            if v:
                links[x][h].SetCoefficient(variable,v)
                for c,coefficient in lookup[R[x]%2][h][v]:
                    c.SetCoefficient(variable,coefficient)
                for c,coefficient in point_lookup[x][h][v]:
                    c.SetCoefficient(variable,coefficient)
        objective.SetCoefficient(variable,model.cost(R[x]%2,state))
        return 1
    for x in range(len(R)):
        add(x,(0,)*len(groups))
    if seed:
        for row in seed['support']:
            if 'point_index' in row:
                add(row['point_index'],tuple(row['counts']))
            else:
                for x,v in enumerate(R):
                    if v%2==row['parity']:
                        add(x,tuple(row['counts']))
    rank_rows,rank_seen = [],set()
    pair_rows,pair_seen = [],set()
    subset_rows,subset_seen = [],set()
    group_rows,group_seen = [],set()
    rank_catalogue = []
    if rank_cuts:
        from rank_capacity import catalogue,violations
        rank_catalogue = catalogue(H,n)
    if pair_cuts:
        from point_pair_capacity import catalogue as pair_catalogue_builder,violations as pair_violations
        pair_catalogue = pair_catalogue_builder(R,n)
    if triple_cuts:
        from core_subset_thresholds import separate as separate_triples
    if subset_cuts:
        from separate_subset_events import separate as separate_subsets
    def add_rank(row):
        if row in rank_seen:
            return
        rank_seen.add(row)
        h,t0,t1,w0,w1,bound = row
        c = solver.Constraint(-solver.infinity(),bound)
        rank_rows.append((*row,c))
        for e,t,w in ((0,t0,w0),(1,t1,w1)):
            if w:
                for v in range(t,len(shapes[h])+1):
                    lookup[e][h][v].append((c,w))
        for x,state,variable in columns:
            coefficient = (w0 if R[x]%2==0 else w1)*int(state[h]>=(t0 if R[x]%2==0 else t1))
            if coefficient:
                c.SetCoefficient(variable,coefficient)
    def add_pair(row):
        if row in pair_seen:
            return
        h,x,y,t,u = row
        assert type(h) is int and 0<=h<len(groups)
        assert 0<=x<y<len(R) and 1<=t<=len(shapes[h]) and 1<=u<=len(shapes[h])
        assert t+u>len(shapes[h])+sum(R[x]%d==R[y]%d for d in shapes[h])
        pair_seen.add(row)
        c = solver.Constraint(-solver.infinity(),1)
        pair_rows.append((*row,c))
        for z,threshold in ((x,t),(y,u)):
            for v in range(threshold,len(shapes[h])+1):
                point_lookup[z][h][v].append((c,1))
        for z,state,variable in columns:
            if (z==x and state[h]>=t) or (z==y and state[h]>=u):
                c.SetCoefficient(variable,1)
    def add_subset(row):
        h,points,thresholds = row
        points,thresholds = tuple(points),tuple(thresholds)
        key = h,points,thresholds
        if key in subset_seen:
            return
        assert len(points)==len(set(points))==len(thresholds) and len(points)>=3
        assert all(0<=x<len(R) for x in points) and all(1<=t<=len(shapes[h]) for t in thresholds)
        cap = sum(max(Counter(R[x]%d for x in points).values()) for d in shapes[h])
        assert sum(thresholds)>cap
        subset_seen.add(key)
        c = solver.Constraint(-solver.infinity(),len(points)-1)
        subset_rows.append((h,points,thresholds,c))
        selected = dict(zip(points,thresholds))
        for x,t in selected.items():
            for v in range(t,len(shapes[h])+1):
                point_lookup[x][h][v].append((c,1))
        for x,state,variable in columns:
            if x in selected and state[h]>=selected[x]:
                c.SetCoefficient(variable,1)
    def add_group(cut):
        h = cut['axis']
        key = (h,cut['bound'],tuple((r['point'],r['threshold'],r['weight']) for r in cut['events']))
        if key in group_seen:
            return
        group_seen.add(key)
        c = solver.Constraint(-solver.infinity(),cut['bound'])
        group_rows.append((cut,c))
        for x in range(len(R)):
            for v in range(len(shapes[h])+1):
                score = sum(r['weight'] for r in cut['events'] if r['point']==x and v>=r['threshold'])
                if score:
                    point_lookup[x][h][v].append((c,score))
        for x,state,variable in columns:
            score = sum(r['weight'] for r in cut['events'] if r['point']==x and state[h]>=r['threshold'])
            if score:
                c.SetCoefficient(variable,score)
    for cut in group_cuts:
        add_group(cut)
    if warm_start:
        for cut in warm_start.get('group_event_constraints',[]):
            add_group(cut)
        if rank_cuts:
            allowed = set(rank_catalogue)
            for row in warm_start.get('rank_constraints',[]):
                assert tuple(row) in allowed
                add_rank(tuple(row))
        if pair_cuts:
            for row in warm_start.get('pair_constraints',[]):
                add_pair(tuple(row))
        if triple_cuts or subset_cuts:
            for row in warm_start.get('subset_constraints',[]):
                add_subset(row)
    scale = 10**6
    iterations,status,best,witness = 0,None,(warm_start or {}).get('certificate'),None
    while time.monotonic()-start<seconds:
        solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
        status = solver.Solve()
        if status!=pywraplp.Solver.OPTIMAL:
            break
        iterations += 1
        master = objective.Value()
        if master<0:
            candidate = rational_primal(source,model,columns,resources,rows,rank_rows,pair_rows,subset_rows,group_rows)
            if candidate['objective_numerator']<0:
                aggregated = dict(denominator=candidate['denominator'],support=[dict(parity=R[row['point_index']]%2,counts=row['counts'],weight=row['weight']) for row in candidate['support']])
                violated = violations(aggregated,n,rank_catalogue) if rank_cuts else []
                bad_pairs = pair_violations(candidate,R,n,pair_catalogue) if pair_cuts else []
                bad_triples = separate_triples(dict(input=source,result=dict(point_primal_witness=candidate)))['violations'] if triple_cuts else []
                if subset_cuts:
                    subset_report = separate_subsets(dict(input=source,result=dict(point_primal_witness=candidate)))
                    bad_triples += subset_report['violations']
                if violated or bad_pairs or bad_triples:
                    for gap,row in violated:
                        add_rank(row)
                    for gap,row in bad_pairs:
                        add_pair(row)
                    for row in bad_triples:
                        if len(row['points'])==2:
                            add_pair((row['axis'],*row['points'],*row['thresholds']))
                        else:
                            add_subset((row['axis'],row['points'],row['thresholds']))
                    print(json.dumps(dict(stage='cuts',rank_added=len(violated),rank_total=len(rank_rows),pair_added=len(bad_pairs),pair_total=len(pair_rows),triple_added=len(bad_triples),subset_total=len(subset_rows),seconds=round(time.monotonic()-start,3))),flush=True)
                    continue
                witness = candidate
                break
        ys = np.minimum(0,np.rint(np.array([r[4].dual_value() for r in rows])*scale).astype(np.int64))
        lambdas = np.minimum(0,np.rint(np.array([[c.dual_value() for c in line] for line in links])*scale).astype(np.int64))
        common = [[np.zeros(len(shape)+1,dtype=np.int64) for shape in shapes] for e in (0,1)]
        for y,(h,t0,t1,bound,c) in zip(ys,rows):
            if y:
                for e,t in enumerate((t0,t1)):
                    common[e][h] += y*np.maximum(0,np.arange(len(shapes[h])+1)-t)
        rank_multipliers = [min(0,round(row[-1].dual_value()*scale)) for row in rank_rows]
        for y,(h,t0,t1,w0,w1,bound,c) in zip(rank_multipliers,rank_rows):
            for e,t,w in ((0,t0,w0),(1,t1,w1)):
                common[e][h] += y*w*(np.arange(len(shapes[h])+1)>=t)
        pair_multipliers = [min(0,round(row[-1].dual_value()*scale)) for row in pair_rows]
        point_penalties = [[np.zeros(len(shape)+1,dtype=np.int64) for shape in shapes] for x in R]
        for coefficient,(h,x,y,t,u,c) in zip(pair_multipliers,pair_rows):
            point_penalties[x][h] += coefficient*(np.arange(len(shapes[h])+1)>=t)
            point_penalties[y][h] += coefficient*(np.arange(len(shapes[h])+1)>=u)
        subset_multipliers = [min(0,round(row[-1].dual_value()*scale)) for row in subset_rows]
        for coefficient,(h,selected,thresholds,c) in zip(subset_multipliers,subset_rows):
            for x,t in zip(selected,thresholds):
                point_penalties[x][h] += coefficient*(np.arange(len(shapes[h])+1)>=t)
        group_multipliers = [min(0,round(c.dual_value()*scale)) for cut,c in group_rows]
        for coefficient,(cut,c) in zip(group_multipliers,group_rows):
            h = cut['axis']
            for event in cut['events']:
                point_penalties[event['point']][h] += coefficient*event['weight']*(np.arange(len(shapes[h])+1)>=event['threshold'])
        prices = []
        for x,v in enumerate(R):
            penalties = [common[v%2][h]+lambdas[x,h]*np.arange(len(shape)+1)+point_penalties[x][h] for h,shape in enumerate(shapes)]
            prices.append(model.price(v%2,penalties,scale,batch))
        points = [row[0][0] for row in prices]
        resource_terms = [min(sum(int(lambdas[x,h]) for x in indices) for residue,indices,var in options) for h,d,options in resources]
        numerator = sum(points)+sum(resource_terms)+sum(int(y)*r[3] for y,r in zip(ys,rows))
        numerator += sum(y*r[5] for y,r in zip(rank_multipliers,rank_rows))
        numerator += sum(pair_multipliers)
        numerator += sum(y*(len(r[1])-1) for y,r in zip(subset_multipliers,subset_rows))
        numerator += sum(y*cut['bound'] for y,(cut,c) in zip(group_multipliers,group_rows))
        certificate = dict(numerator=numerator,denominator=scale,integer_lower=-(-numerator//scale),point_multipliers=points,
                           link_multipliers=lambdas.tolist(),resource_multipliers=resource_terms,
                           tail_multipliers=[dict(axis=r[0],thresholds=list(r[1:3]),bound=r[3],multiplier=int(y)) for y,r in zip(ys,rows) if y],
                           rank_multipliers=[dict(axis=r[0],thresholds=list(r[1:3]),weights=list(r[3:5]),bound=r[5],multiplier=y) for y,r in zip(rank_multipliers,rank_rows) if y],
                           pair_multipliers=[dict(axis=r[0],points=list(r[1:3]),thresholds=list(r[3:5]),multiplier=y) for y,r in zip(pair_multipliers,pair_rows) if y],
                           subset_multipliers=[dict(axis=r[0],points=list(r[1]),thresholds=list(r[2]),multiplier=y) for y,r in zip(subset_multipliers,subset_rows) if y],
                           group_event_multipliers=[dict(cut,multiplier=y) for y,(cut,c) in zip(group_multipliers,group_rows) if y])
        if best is None or numerator>best['numerator']:
            best = certificate
        print(json.dumps(dict(iteration=iterations,columns=len(columns),master=master,lower=best['integer_lower'],seconds=round(time.monotonic()-start,3))),flush=True)
        if best['integer_lower']>=target_lower:
            break
        mass_duals = [c.dual_value() for c in mass]
        additions = sum(add(x,state) for x,row in enumerate(prices) for value,state in row if mass_duals[x]-value/scale>0.001)
        if not additions:
            break
    return dict(certificate=best,point_primal_witness=witness,iterations=iterations,master_status=status,
                columns=len(columns),seconds=round(time.monotonic()-start,3),groups=groups,
                rank_cuts_enabled=rank_cuts,rank_constraints=[list(row[:-1]) for row in rank_rows],
                pair_cuts_enabled=pair_cuts,pair_constraints=[list(row[:-1]) for row in pair_rows],
                 triple_cuts_enabled=triple_cuts,subset_cuts_enabled=subset_cuts,subset_constraints=[list(row[:-1]) for row in subset_rows],
                column_states=[dict(point_index=x,counts=list(state)) for x,state,var in columns],
                group_event_constraints=[cut for cut,c in group_rows],
                scope='point-identity conditional-tail relaxation for one specified core and pattern only')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=180)
    ap.add_argument('--batch',type=int,default=2)
    ap.add_argument('--rank-cuts',action='store_true')
    ap.add_argument('--pair-cuts',action='store_true')
    ap.add_argument('--triple-cuts',action='store_true')
    ap.add_argument('--subset-cuts',action='store_true')
    ap.add_argument('--group-cuts',help='JSON file containing core_residual and cuts, each with independently audited proof metadata')
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    seed = data['result'].get('point_primal_witness') or data['result'].get('primal_relaxation_witness')
    if not seed and 'column_states' in data['result']:
        seed = dict(support=data['result']['column_states'])
    cuts = []
    if args.group_cuts:
        supplied = json.loads(Path(args.group_cuts).read_text())
        assert supplied['core_residual']==data['input']['core_residual']
        assert all(c.get('proof_audit',{}).get('verified') for c in supplied['cuts'])
        cuts = supplied['cuts']
    result = solve(data['input'],args.time,args.batch,seed=seed,rank_cuts=args.rank_cuts,pair_cuts=args.pair_cuts,warm_start=data['result'],triple_cuts=args.triple_cuts,subset_cuts=args.subset_cuts,group_cuts=cuts)
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps(dict(iterations=result['iterations'],seconds=result['seconds'],
                          lower=result['certificate']['integer_lower'] if result['certificate'] else None,
                          negative_primal=bool(result['point_primal_witness']))),flush=True)
