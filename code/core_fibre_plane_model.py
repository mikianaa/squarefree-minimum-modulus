"""Keep actual core projections when bounding higher-coordinate planes.

For each core point x and higher prime q, let n_q(x) be the number of plane
resources whose core projections contain x. Distinct excluded q labels can
number at most n_q(x). Hence the residual after planes is at least
sum_x product_q max(0,D_q(x)-n_q(x)). Higher-codimension resources are then
granted their maximum capacity on the original core. A model solution is
only a relaxation witness, not an actual covering system.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import random
import time


def prepare(source,high,pattern):
    family = {r['modulus']:r['residue'] for r in source['system']}
    core = {d:a for d,a in family.items() if 210%d==0}
    assert set(core)=={d for d in range(3,211) if 210%d==0}
    residual = [x for x in range(210) if not any(x%d==a for d,a in core.items())]
    dimensions = [[q-1-(b==x%2) for x in residual] for q,b in zip(high,pattern)]
    resources = []
    for i,q in enumerate(high):
        for d in sorted(core):
            options = [(a,tuple(j for j,x in enumerate(residual) if x%d==a))
                       for a in sorted({x%d for x in residual})]
            resources.append(dict(axis=i,core_modulus=d,options=options))
    tail = 0
    for k in range(2,len(high)+1):
        for fixed in itertools.combinations(range(len(high)),k):
            volumes = [math.prod(q-1-(pattern[i]==e) for i,q in enumerate(high) if i not in fixed) for e in (0,1)]
            for d in range(1,106):
                if 105%d:
                    continue
                caps = [[sum(x%2==e and x%d==a for x in residual) for e in (0,1)] for a in range(d)]
                odd = max(sum(c[e]*volumes[e] for e in (0,1)) for c in caps)
                even = max(c[e]*volumes[e] for c in caps for e in (0,1))
                tail += odd+even
    return core,residual,dimensions,resources,tail


def search(source,high,pattern,seconds=3,seed=1):
    core,residual,dimensions,resources,tail = prepare(source,high,pattern)
    family = {r['modulus']:r['residue'] for r in source['system']}
    rng = random.Random(seed)
    choices = []
    counts = [[0]*len(residual) for _ in high]
    for resource in resources:
        i,d = resource['axis'],resource['core_modulus']
        a = family.get(high[i]*d,0)%d
        choice = next((j for j,(b,_) in enumerate(resource['options']) if a==b),rng.randrange(len(resource['options'])))
        choices.append(choice)
        for x in resource['options'][choice][1]:
            counts[i][x] += 1
    best = None
    best_choices = None
    start = time.monotonic()
    passes = 0
    while passes==0 or time.monotonic()-start<seconds:
        order = list(range(len(resources)))
        rng.shuffle(order)
        for j in order:
            resource = resources[j]
            i = resource['axis']
            for x in resource['options'][choices[j]][1]:
                counts[i][x] -= 1
            gains = [math.prod(max(0,dimensions[h][x]-counts[h][x]) for h in range(len(high)) if h!=i)
                     if counts[i][x]<dimensions[i][x] else 0 for x in range(len(residual))]
            scores = [sum(gains[x] for x in points) for _,points in resource['options']]
            maximum = max(scores)
            choices[j] = rng.choice([i for i,s in enumerate(scores) if s==maximum])
            for x in resource['options'][choices[j]][1]:
                counts[i][x] += 1
        value = sum(math.prod(max(0,dimensions[i][x]-counts[i][x]) for i in range(len(high))) for x in range(len(residual)))
        passes += 1
        if best is None or value<best:
            best,best_choices = value,list(choices)
        for j in rng.sample(range(len(resources)),4):
            resource = resources[j]
            i = resource['axis']
            for x in resource['options'][choices[j]][1]:
                counts[i][x] -= 1
            choices[j] = rng.randrange(len(resource['options']))
            for x in resource['options'][choices[j]][1]:
                counts[i][x] += 1
    chosen = [dict(high_prime=high[r['axis']],core_modulus=r['core_modulus'],core_residue=r['options'][c][0])
              for r,c in zip(resources,best_choices)]
    # Independently reconstruct projection counts from ordinary residues.
    actual = sum(math.prod(max(0,q-1-(b==x%2)-sum(x%r['core_modulus']==r['core_residue'] for r in chosen if r['high_prime']==q))
                           for q,b in zip(high,pattern)) for x in residual)
    assert actual==best
    return dict(core_system=[dict(modulus=d,residue=a) for d,a in sorted(core.items())],
                core_residual=residual,high_primes=list(high),high_pattern=list(pattern),
                minimum_found_plane_expression=best,tail_capacity=tail,passes=passes,
                choices=chosen,seconds=round(time.monotonic()-start,3),
                scope='relaxed projection assignment; no optimality or actual covering claim')


def solve_relaxation(incumbent,seconds=60,workers=4):
    from ortools.sat.python import cp_model
    high,pattern = incumbent['high_primes'],incumbent['high_pattern']
    _,residual,dimensions,resources,tail = prepare(dict(system=incumbent['core_system']),high,pattern)
    model = cp_model.CpModel()
    choices = []
    by_point = [[[] for _ in residual] for _ in high]
    hints = {(r['high_prime'],r['core_modulus']):r['core_residue'] for r in incumbent['choices']}
    for j,r in enumerate(resources):
        variables = [model.NewBoolVar(f'choice{j}_{a}') for a,points in r['options']]
        model.AddExactlyOne(variables)
        for (a,points),v in zip(r['options'],variables):
            model.AddHint(v,int(a==hints[high[r['axis']],r['core_modulus']]))
            for x in points:
                by_point[r['axis']][x].append(v)
        choices.append(variables)
    products = []
    for x in range(len(residual)):
        factors = []
        for i in range(len(high)):
            D = dimensions[i][x]
            z = model.NewIntVar(max(0,D-14),D,f'left_{i}_{x}')
            model.AddMaxEquality(z,[0,D-sum(by_point[i][x])])
            factors.append(z)
        value = model.NewIntVar(0,math.prod(dimensions[i][x] for i in range(len(high))),f'product{x}')
        model.AddMultiplicationEquality(value,factors)
        products.append(value)
    model.Minimize(sum(products))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.linearization_level = 2
    start = time.monotonic()
    status = solver.Solve(model)
    result = dict(status=solver.StatusName(status),seconds=round(time.monotonic()-start,3),
                  tail_capacity=tail,scope='solver-bound relaxation for this fixed actual core, not a global proof')
    if status in (cp_model.UNKNOWN,cp_model.FEASIBLE,cp_model.OPTIMAL):
        result['plane_expression_lower'] = max(0,math.ceil(solver.BestObjectiveBound()-1e-7))
    if status in (cp_model.FEASIBLE,cp_model.OPTIMAL):
        result['choices'] = [dict(high_prime=high[r['axis']],core_modulus=r['core_modulus'],
                                  core_residue=next(a for (a,pts),v in zip(r['options'],variables) if solver.Value(v)))
                             for r,variables in zip(resources,choices)]
        value = sum(math.prod(max(0,q-1-(b==x%2)-sum(x%r['core_modulus']==r['core_residue'] for r in result['choices'] if r['high_prime']==q))
                               for q,b in zip(high,pattern)) for x in residual)
        assert value==round(solver.ObjectiveValue()) and result['plane_expression_lower']<=value
        result['plane_expression_upper'] = value
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/min-uncovered-9699690-witness.json')
    ap.add_argument('--time',type=float,default=3)
    ap.add_argument('--seed',type=int,default=1)
    ap.add_argument('--output',required=True)
    ap.add_argument('--solver-time',type=float,default=0)
    a = ap.parse_args()
    result = search(json.loads(Path(a.source).read_text()),(11,13,17,19),(1,1,1,1),a.time,a.seed)
    if a.solver_time:
        result['solver'] = solve_relaxation(result,a.solver_time)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('core_system','core_residual','choices')}))
