"""Search coherent group mixtures via small, independent point couplings.

Each higher group mixes two actual residue assignments. Plane counts can be
decreased; their three relevant values are 0 and the two available hit counts.
There are at most 3**5 plane states per point. Tail count marginals are coupled
to those plane states by a transport LP, then realized by nested intervals.
Every saved rational witness is checked by audit_coherent_groups.
"""
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import random
import time
from ortools.linear_solver import pywraplp


def structure(source):
    n = len(source['high_primes'])
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    shapes = [[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(len(groups))]
    return groups,shapes


def counts_for(source,assignments):
    groups,shapes = structure(source)
    assert len(assignments)==len(groups)
    for row,shape in zip(assignments,shapes):
        assert len(row)==len(shape) and all(type(a) is int and 0<=a<d for a,d in zip(row,shape))
    return [[sum(x%d==a for a,d in zip(row,shape)) for row,shape in zip(assignments,shapes)]
            for x in source['core_residual']]


def cost(source,x,state):
    groups,_ = structure(source)
    n = len(source['high_primes'])
    free = [q-1-(b==x%2)-state[i] for i,(q,b) in enumerate(zip(source['high_primes'],source['high_pattern']))]
    return math.prod(free)-sum(state[h]*math.prod(free[i] for i in range(n) if i not in S)
                               for h,S in enumerate(groups[n:],n))


def probabilities(value,group_count):
    row = list(value) if isinstance(value,(list,tuple)) else [value]*group_count
    assert len(row)==group_count and all(0<=p<=1 for p in row)
    return row


def point_coupling(source,point,a,b,probability):
    groups,_ = structure(source)
    n = len(source['high_primes'])
    D = [q-1-(v==point%2) for q,v in zip(source['high_primes'],source['high_pattern'])]
    a = [min(v,D[h]) if h<n else v for h,v in enumerate(a)]
    b = [min(v,D[h]) if h<n else v for h,v in enumerate(b)]
    low = [min(v,w) for v,w in zip(a,b)]
    high = [max(v,w) for v,w in zip(a,b)]
    theta = probabilities(probability,len(groups))
    high_prob = [float(p if w>v else 1-p) for v,w,p in zip(a,b,theta)]
    states = list(itertools.product(*(sorted({0,a[h],b[h]}) for h in range(n))))
    solver = pywraplp.Solver.CreateSolver('GLOP')
    objective = solver.Objective()
    objective.SetMinimization()
    mass = [solver.NumVar(0,1,f'p{j}') for j in range(len(states))]
    solver.Add(sum(mass)==1)
    marginal_constraints = {}
    for h in range(n):
        if low[h]<high[h]:
            marginal_constraints[h] = solver.Add(sum(v for state,v in zip(states,mass) if state[h]==high[h])<=high_prob[h])
    transport = [[] for state in states]
    tail_mass = {h:[] for h in range(n,len(groups)) if low[h]<high[h]}
    for j,(state,var) in enumerate(zip(states,mass)):
        free = [d-v for d,v in zip(D,state)]
        volumes = [math.prod(free[i] for i in range(n) if i not in S) for S in groups[n:]]
        objective.SetCoefficient(var,math.prod(free)-sum(low[h]*vol for h,vol in enumerate(volumes,n)))
        for h,vol in enumerate(volumes,n):
            if low[h]==high[h]:
                continue
            y = solver.NumVar(0,1,f'y{j}_{h}')
            solver.Add(y<=var)
            objective.SetCoefficient(y,-(high[h]-low[h])*vol)
            tail_mass[h].append(y)
            transport[j].append((h,y))
    for h,line in tail_mass.items():
        marginal_constraints[h] = solver.Add(sum(line)==high_prob[h])
    assert solver.Solve()==pywraplp.Solver.OPTIMAL
    support = []
    for state,var,line in zip(states,mass,transport):
        p = max(0,var.solution_value())
        if p<1e-12:
            continue
        heights = [(h,min(p,max(0,y.solution_value()))) for h,y in line]
        endpoints = sorted({0,p,*[v for h,v in heights]})
        for left,right in zip(endpoints,endpoints[1:]):
            if right-left<1e-12:
                continue
            row = list(state)+low[n:]
            for h,height in heights:
                if height>left:
                    row[h] = high[h]
            support.append((row,right-left))
    direct = sum(weight*cost(source,point,row) for row,weight in support)
    assert abs(direct-objective.Value())<1e-5*max(1,abs(direct))
    gradient = [((1 if b[h]>a[h] else -1)*marginal_constraints[h].dual_value())
                if h in marginal_constraints else 0 for h in range(len(groups))]
    return dict(objective=direct,support=support,plane_states=len(states),gradient=gradient,
                supporting_intercept=objective.Value()-sum(g*float(p) for g,p in zip(gradient,theta)))


def rational_witness(source,A,B,probability,solutions,scale=10**9):
    groups,shapes = structure(source)
    theta = probabilities(probability,len(groups))
    assert all(isinstance(p,Fraction) and scale%p.denominator==0 for p in theta)
    wb = [int(scale*p) for p in theta]
    wa = [scale-w for w in wb]
    saved_groups = [dict(axis=h,assignments=[dict(residues=list(row),weight=w)
                    for row,w in ((A[h],wa[h]),(B[h],wb[h])) if w]) for h in range(len(groups))]
    ca,cb = counts_for(source,A),counts_for(source,B)
    support = []
    for x,solution in enumerate(solutions):
        merged = Counter()
        for state,p in solution['support']:
            merged[tuple(state)] += max(0,math.floor(p*scale))
        ratio = min(Fraction(1),Fraction(scale,sum(merged.values())))
        for h,shape in enumerate(shapes):
            for t in range(1,len(shape)+1):
                used = sum(w for state,w in merged.items() if state[h]>=t)
                available = wa[h]*(ca[x][h]>=t)+wb[h]*(cb[x][h]>=t)
                if used:
                    ratio = min(ratio,Fraction(available,used))
        repaired = {state:w*ratio.numerator//ratio.denominator for state,w in merged.items()}
        repaired[(0,)*len(groups)] = repaired.get((0,)*len(groups),0)+scale-sum(repaired.values())
        support.extend(dict(point_index=x,counts=list(state),weight=w) for state,w in repaired.items() if w)
    objective = sum(r['weight']*cost(source,source['core_residual'][r['point_index']],r['counts']) for r in support)
    return dict(denominator=scale,objective_numerator=objective,groups=saved_groups,support=support)


def evaluate(source,A,B,probability=Fraction(1,2),save_witness=False):
    ca,cb = counts_for(source,A),counts_for(source,B)
    solutions = [point_coupling(source,x,a,b,probability) for x,a,b in zip(source['core_residual'],ca,cb)]
    result = dict(objective=sum(r['objective'] for r in solutions),plane_states=sum(r['plane_states'] for r in solutions))
    if save_witness:
        from audit_coherent_groups import audit
        result['coherent_primal_witness'] = rational_witness(source,A,B,probability,solutions)
        result['audit'] = audit(dict(input=source,result=result))
    return result


def assignments_from_integer(source,result):
    groups,shapes = structure(source)
    choices = {(r['axis'],r['modulus']):r['residue'] for r in result['choices']}
    return [[choices[h,d] for d in shape] for h,shape in enumerate(shapes)]


def search(source,integer,pool,seconds=120,seed=0):
    start = time.monotonic()
    rng = random.Random(seed)
    A = assignments_from_integer(source,integer)
    n = len(source['high_primes'])
    configurations = pool['group_columns']
    prices = pool['certificate']['event_multipliers']
    _,shapes = structure(source)
    ranked = []
    for h,rows in enumerate(configurations):
        def score(residues):
            hits = [sum(x%d==a for d,a in zip(shapes[h],residues)) for x in source['core_residual']]
            return -sum(sum(prices[h][x][1:v+1]) for x,v in enumerate(hits))
        ranked.append(sorted(rows,key=score,reverse=True))
    baseline = evaluate(source,A,A,save_witness=True)
    best = dict(baseline,assignments_A=A,assignments_B=A,probability=[1,2])
    history = [dict(kind='baseline',objective=baseline['objective'])]
    iteration = 0
    while time.monotonic()-start<seconds:
        if iteration==0:
            B = [rows[0] for rows in ranked]
        elif iteration%3==1:
            B = [rng.choice(rows[:min(8,len(rows))]) for rows in ranked]
        else:
            B = [list(row) for row in A]
            # Small actual residue changes probe a neighborhood of the
            # strongest known integer assignment, with coherent choices.
            for h in rng.sample(range(len(B)),k=min(len(B),5+iteration%8)):
                for j in rng.sample(range(len(B[h])),k=min(3,len(B[h]))):
                    B[h][j] = rng.choice([x%shapes[h][j] for x in source['core_residual']])
        probability = (Fraction(1,2),Fraction(1,10),Fraction(9,10))[iteration%3]
        result = evaluate(source,A,B,probability)
        iteration += 1
        row = dict(iteration=iteration,objective=result['objective'],seconds=time.monotonic()-start)
        history.append(row)
        if result['objective']<best['objective']-1e-7:
            exact = evaluate(source,A,B,probability,save_witness=True)
            best = dict(exact,assignments_A=A,assignments_B=B,probability=[probability.numerator,probability.denominator])
            print(json.dumps(dict(row,best_rational=best['coherent_primal_witness']['objective_numerator'])),flush=True)
        if best['coherent_primal_witness']['objective_numerator']<0:
            break
    best.update(history=history,iterations=iteration,seconds=time.monotonic()-start,
                scope='two actual assignments per group with arbitrary local coupling; not a covering system')
    return best


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--integer',required=True)
    ap.add_argument('--pool',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=120)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    integer = json.loads(Path(args.integer).read_text())
    pool = json.loads(Path(args.pool).read_text())
    assert data['input']['core_residual']==integer['input']['core_residual']==pool['input']['core_residual']
    result = search(data['input'],integer['result'],pool['result'],args.time)
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps(dict(iterations=result['iterations'],seconds=result['seconds'],objective=result['objective'],audit=result['audit'])))
