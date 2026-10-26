"""Conditional-tail relaxation with full per-group configuration mixtures.

Local point count distributions are stochastically dominated by explicit
mixtures of consistent resource assignments, separately for every higher
group. Groups need not have a common joint distribution. A negative rational
witness is checkable directly; positive bounds still require pricing proofs.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import random
import time
import numpy as np
from ortools.linear_solver import pywraplp
from group_event_coherence import maximize
from multitail_lp import LocalStates


def rational_witness(source,model,columns,configurations,shapes):
    R,scale = source['core_residual'],10**9
    supply = [[[0]*(len(shape)+1) for x in R] for shape in shapes]
    groups = []
    for h,configs in enumerate(configurations):
        values = [max(Fraction(0),Fraction(str(var.solution_value()))) for residues,counts,var in configs]
        total = sum(values)
        assert total>0
        weights = [int(value*scale//total) for value in values]
        weights[0] += scale-sum(weights)
        saved = []
        for (residues,counts,var),weight in zip(configs,weights):
            if not weight:
                continue
            saved.append(dict(residues=list(residues),weight=weight))
            for x,v in enumerate(counts):
                for t in range(1,v+1):
                    supply[h][x][t] += weight
        groups.append(dict(axis=h,assignments=saved))
    support = [(x,state,max(0,math.floor(var.solution_value()*scale))) for x,state,var in columns]
    support = [(x,state,w) for x,state,w in support if w]
    ratio = Fraction(1)
    for x in range(len(R)):
        mass = sum(w for y,state,w in support if y==x)
        if mass:
            ratio = min(ratio,Fraction(scale,mass))
        for h,shape in enumerate(shapes):
            for t in range(1,len(shape)+1):
                used = sum(w for y,state,w in support if y==x and state[h]>=t)
                if used:
                    ratio = min(ratio,Fraction(supply[h][x][t],used))
    support = [(x,state,w*ratio.numerator//ratio.denominator) for x,state,w in support]
    support = [(x,state,w) for x,state,w in support if w]
    for x in range(len(R)):
        missing = scale-sum(w for y,state,w in support if y==x)
        if missing:
            support.append((x,(0,)*len(shapes),missing))
    objective = sum(w*model.cost(R[x]%2,state) for x,state,w in support)
    return dict(denominator=scale,objective_numerator=objective,groups=groups,
                support=[dict(point_index=x,counts=list(state),weight=w) for x,state,w in support])


def solve(source,seconds=180,seed=None,pricing_seconds=.3,warm_start=None,integer_seed=None,target_lower=1,stabilization=0,dual_seed=None):
    start = time.monotonic()
    R = source['core_residual']
    model = LocalStates(source['high_primes'],source['high_pattern'])
    n,G = model.n,len(model.groups)
    shapes = [[d for d in range(1 if h>=n else 3,211) if 210%d==0] for h in range(G)]
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    objective = solver.Objective()
    objective.SetMinimization()
    mass = [solver.Constraint(1,1) for x in R]
    group_mass = [solver.Constraint(1,1) for h in range(G)]
    links = [[[None]+[solver.Constraint(-solver.infinity(),0) for t in shape] for x in R] for shape in shapes]
    columns,seen = [],[set() for x in R]
    configurations = [[] for h in range(G)]
    config_seen = [set() for h in range(G)]
    def add_point(x,state):
        state = tuple(state)
        if state in seen[x]:
            return 0
        seen[x].add(state)
        var = solver.NumVar(0,solver.infinity(),f'z{x}_{len(seen[x])}')
        mass[x].SetCoefficient(var,1)
        objective.SetCoefficient(var,model.cost(R[x]%2,state))
        for h,v in enumerate(state):
            for t in range(1,v+1):
                links[h][x][t].SetCoefficient(var,1)
        columns.append((x,state,var))
        return 1
    def add_config(h,residues):
        assert len(residues)==len(shapes[h]) and all(0<=a<d for a,d in zip(residues,shapes[h]))
        counts = tuple(sum(v%d==a for d,a in zip(shapes[h],residues)) for v in R)
        if counts in config_seen[h]:
            return 0
        config_seen[h].add(counts)
        var = solver.NumVar(0,solver.infinity(),f'g{h}_{len(config_seen[h])}')
        group_mass[h].SetCoefficient(var,1)
        for x,v in enumerate(counts):
            for t in range(1,v+1):
                links[h][x][t].SetCoefficient(var,-1)
        configurations[h].append((tuple(residues),counts,var))
        return 1
    for x in range(len(R)):
        add_point(x,(0,)*G)
    if seed:
        for row in seed['support']:
            add_point(row['point_index'],row['counts'])
    if warm_start:
        for h,configs in enumerate(warm_start.get('group_columns',[])):
            for residues in configs:
                add_config(h,residues)
    if integer_seed:
        choices = {(r['axis'],r['modulus']):r['residue'] for r in integer_seed['choices']}
        state = [[0]*G for x in R]
        for h,shape in enumerate(shapes):
            residues = [choices[h,d] for d in shape]
            add_config(h,residues)
            for x,v in enumerate(R):
                hits = sum(v%d==a for d,a in zip(shape,residues))
                state[x][h] = min(hits,source['high_primes'][h]-1-(source['high_pattern'][h]==v%2)) if h<n else hits
        assert sum(model.cost(v%2,state[x]) for x,v in enumerate(R))==integer_seed['integer_expression']
        for x,row in enumerate(state):
            add_point(x,row)
    rng = random.Random(0)
    for h,shape in enumerate(shapes):
        for x in R:
            add_config(h,[x%d for d in shape])
        for repeat in range(12):
            add_config(h,[rng.choice(sorted({x%d for x in R})) for d in shape])
    scale,best,witness = 10**6,(warm_start or {}).get('certificate'),None
    assert 0<=stabilization<1
    center = dual_seed or ((warm_start or {}).get('stabilization_center'))
    if center is None and best:
        center = best['event_multipliers']
    center_value = None
    mixing = stabilization
    iterations,history = 0,[]
    while time.monotonic()-start<seconds:
        solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
        if solver.Solve()!=pywraplp.Solver.OPTIMAL:
            break
        iterations += 1
        master = objective.Value()
        if master<0:
            candidate = rational_witness(source,model,columns,configurations,shapes)
            if candidate['objective_numerator']<0:
                witness = candidate
                break
        # Snapshot every dual before adding any column.
        point_duals = [c.dual_value() for c in mass]
        group_duals = [c.dual_value() for c in group_mass]
        raw = [[[0]+[min(0,round(c.dual_value()*scale)) for c in line[1:]] for line in group] for group in links]
        if center is not None and mixing:
            fraction = 1 if dual_seed is not None and iterations==1 else mixing
            lambdas = [[[round(fraction*a+(1-fraction)*b) for a,b in zip(old,new)] for old,new in zip(oldgroup,newgroup)] for oldgroup,newgroup in zip(center,raw)]
        else:
            lambdas = raw
        priced_points = []
        for x,v in enumerate(R):
            penalties = [np.cumsum(np.array(lambdas[h][x],dtype=np.int64)) for h in range(G)]
            priced_points.append(model.price(v%2,penalties,scale,batch=16))
        priced_groups = []
        for h,shape in enumerate(shapes):
            events = [(x,t) for x in range(len(R)) for t in range(1,len(shape)+1) if lambdas[h][x][t]]
            weights = [-lambdas[h][x][t] for x,t in events]
            priced_groups.append(maximize(R,shape,events,weights,pricing_seconds,collect=True))
        if all('upper_bound' in row for row in priced_groups):
            numerator = sum(row[0][0] for row in priced_points)-sum(row['upper_bound'] for row in priced_groups)
            if center_value is None or numerator>center_value:
                center,center_value = lambdas,numerator
            if best is None or numerator>best['numerator']:
                best = dict(numerator=numerator,denominator=scale,integer_lower=-(-numerator//scale),
                            point_multipliers=[row[0][0] for row in priced_points],event_multipliers=lambdas,
                            group_pricing=priced_groups,scope='CP-SAT upper bounds not independently proof-checked')
        additions = sum(add_point(x,state) for x,row in enumerate(priced_points) for value,state in row if mixing or value/scale<point_duals[x]-.001)
        for h,row in enumerate(priced_groups):
            candidates = row.get('candidate_assignments',[])+([row] if 'objective' in row else [])
            for candidate in candidates:
                if mixing or -candidate['objective']/scale<group_duals[h]-.001:
                    additions += add_config(h,candidate['residues'])
        row = dict(iteration=iterations,master=master,point_columns=len(columns),group_columns=sum(map(len,configurations)),
                   additions=additions,lower=best['integer_lower'] if best else None,stabilization=mixing,seconds=time.monotonic()-start)
        history.append(row)
        print(json.dumps(row),flush=True)
        if not additions and mixing:
            mixing = mixing/2 if mixing>.01 else 0
        elif not additions:
            break
        if best and best['integer_lower']>=target_lower:
            break
    return dict(coherent_primal_witness=witness,certificate=best,iterations=iterations,history=history,
                point_columns=[dict(point_index=x,counts=list(state)) for x,state,var in columns],
                group_columns=[[list(residues) for residues,counts,var in configs] for configs in configurations],
                stabilization_center=center,
                seconds=time.monotonic()-start,scope='full per-group consistency; cross-group coupling remains relaxed')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=180)
    ap.add_argument('--integer-source')
    ap.add_argument('--stabilization',type=float,default=0)
    ap.add_argument('--dual-seed')
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    seed = data['result'].get('point_primal_witness') or data['result'].get('coherent_primal_witness')
    if not seed:
        seed = dict(support=data['result'].get('point_columns',data['result'].get('column_states',[])))
    integer_seed = None
    if args.integer_source:
        integer = json.loads(Path(args.integer_source).read_text())
        assert integer['input']==data['input']
        integer_seed = integer['result']
    dual_seed = None
    if args.dual_seed:
        from event_dual_seed import convert
        old = json.loads(Path(args.dual_seed).read_text())
        assert old['input']==data['input']
        dual_seed = convert(data['input'],old['result']['certificate'])
    result = solve(data['input'],args.time,seed,warm_start=data['result'],integer_seed=integer_seed,stabilization=args.stabilization,dual_seed=dual_seed)
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps(dict(iterations=result['iterations'],negative_primal=bool(result['coherent_primal_witness']),
                         lower=result['certificate']['integer_lower'] if result['certificate'] else None)),flush=True)
