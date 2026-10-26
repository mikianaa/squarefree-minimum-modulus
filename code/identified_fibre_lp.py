"""Column generation retaining individual core-point projection identities.

The master links each point's distribution of higher-axis hit counts to
fractional choices of actual core congruences. Exact rational dual checks price
every hit-count tuple and every core projection, so solver convergence is not
needed to validate a reported lower bound. This remains a relaxation for one
fixed actual core, not a proof about all covering systems.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from core_fibre_plane_model import prepare


def solve(incumbent,concentration,seconds=60,unit_mode='none'):
    high,pattern = incumbent['high_primes'],incumbent['high_pattern']
    _,residual,dimensions,resources,tail = prepare(dict(system=incumbent['core_system']),high,pattern)
    # After independently augmenting each point's plane deletions to the
    # maximum hit count, tail classes still have to cover the remaining set.
    # For core modulus 1 their capacity is exactly the sum of free-coordinate
    # volumes. For modulus 2 it is at most that same sum across both parities.
    assert unit_mode in ('none','odd','both')
    unit_coefficient = {'none':0,'odd':1,'both':2}[unit_mode]
    high_subsets = [S for k in range(2,len(high)+1) for S in itertools.combinations(range(len(high)),k)]
    if unit_coefficient:
        for S in high_subsets:
            volumes = [math.prod(dimensions[i][x] for i in range(len(high)) if i not in S) for x in range(len(residual))]
            tail -= sum(volumes)
            if unit_mode=='both':
                tail -= max(sum(v for v,n in zip(volumes,residual) if n%2==e) for e in (0,1))
    assert residual==concentration['core_residual']
    H = np.array(concentration['two_layer_capacities'],dtype=np.int64)
    m0,m1 = np.indices(H.shape)
    maxhits = len(concentration['moduli'])
    start = time.monotonic()
    solver = pywraplp.Solver.CreateSolver('GLOP')
    objective = solver.Objective()
    objective.SetMinimization()
    mass = [solver.Constraint(1,1) for _ in residual]
    links = [[solver.Constraint(-solver.infinity(),0) for _ in high] for _ in residual]
    resource_constraints = []
    for j,r in enumerate(resources):
        constraint = solver.Constraint(1,1)
        resource_constraints.append(constraint)
        for a,points in r['options']:
            v = solver.NumVar(0,solver.infinity(),f'r{j}_{a}')
            constraint.SetCoefficient(v,1)
            for x in points:
                links[x][r['axis']].SetCoefficient(v,-1)
    # All mixed thresholds are inexpensive in the small generated master.
    tail_rows = []
    for axis in range(len(high)):
        for t0,t1 in itertools.product(range(maxhits+1),repeat=2):
            if t0==maxhits and t1==maxhits:
                continue
            bound = int(np.max(H-t0*m0-t1*m1))
            tail_rows.append((axis,t0,t1,bound,solver.Constraint(-solver.infinity(),bound)))
    tuples,costs = [],[]
    for e in (0,1):
        D = np.array([q-1-(b==e) for q,b in zip(high,pattern)],dtype=np.int64)
        grid = np.array(list(itertools.product(*(range(min(maxhits,int(d))+1) for d in D))),dtype=np.int16)
        tuples.append(grid)
        free = D[None,:]-grid
        cost = np.prod(free,axis=1)
        if unit_coefficient:
            for S in high_subsets:
                cost -= unit_coefficient*np.prod(free[:,[i for i in range(len(high)) if i not in S]],axis=1)
        costs.append(cost)
    seen = [set() for _ in residual]
    columns = []
    def add_column(x,values):
        values = tuple(map(int,values))
        if values in seen[x]:
            return False
        seen[x].add(values)
        v = solver.NumVar(0,solver.infinity(),f'z{x}_{len(seen[x])}')
        columns.append((x,values,v))
        mass[x].SetCoefficient(v,1)
        for i,n in enumerate(values):
            if n:
                links[x][i].SetCoefficient(v,n)
        for axis,t0,t1,bound,constraint in tail_rows:
            excess = max(0,values[axis]-(t0 if residual[x]%2==0 else t1))
            if excess:
                constraint.SetCoefficient(v,excess)
        cost = math.prod(dimensions[i][x]-n for i,n in enumerate(values))
        if unit_coefficient:
            cost -= unit_coefficient*sum(math.prod(dimensions[i][x]-values[i] for i in range(len(high)) if i not in S) for S in high_subsets)
        objective.SetCoefficient(v,cost)
        return True
    # Zero-hit columns guarantee initial feasibility without a heuristic input.
    for x in range(len(residual)):
        add_column(x,[0]*len(high))
        values = [min(dimensions[i][x],sum(residual[x]%r['core_modulus']==r['core_residue']
                  for r in incumbent['choices'] if r['high_prime']==q)) for i,q in enumerate(high)]
        add_column(x,values)
    iterations = 0
    best_certificate = None
    complete = False
    while time.monotonic()-start<seconds:
        solver.SetTimeLimit(max(1,int(1000*(seconds-(time.monotonic()-start)))))
        status = solver.Solve()
        if status!=pywraplp.Solver.OPTIMAL:
            break
        iterations += 1
        link_duals = np.array([[c.dual_value() for c in row] for row in links])
        tail_duals = np.array([row[4].dual_value() for row in tail_rows])
        common = [np.zeros(len(g),dtype=float) for g in tuples]
        scale = 10**6
        link_integer = np.minimum(0,np.rint(link_duals*scale).astype(np.int64))
        tail_integer = np.minimum(0,np.rint(tail_duals*scale).astype(np.int64))
        exact_common = [np.zeros(len(g),dtype=np.int64) for g in tuples]
        for e in (0,1):
            for k,(axis,t0,t1,bound,_) in enumerate(tail_rows):
                if tail_duals[k] or tail_integer[k]:
                    excess = np.maximum(0,tuples[e][:,axis]-(t0 if e==0 else t1))
                    common[e] += tail_duals[k]*excess
                    exact_common[e] += tail_integer[k]*excess
        prices = []
        point_integer = []
        max_violation = 0
        for x,n in enumerate(residual):
            e = n%2
            reduced = costs[e]-common[e]-tuples[e]@link_duals[x]-mass[x].dual_value()
            j = int(np.argmin(reduced))
            if reduced[j]<-1e-6:
                prices.append((x,tuples[e][j]))
                max_violation = max(max_violation,float(-reduced[j]))
            slack = costs[e]*scale-exact_common[e]-tuples[e].astype(np.int64)@link_integer[x]
            point_integer.append(int(slack.min()))
        resource_integer = [min(sum(int(link_integer[x,r['axis']]) for x in points) for a,points in r['options']) for r in resources]
        numerator = sum(point_integer)+sum(resource_integer)+sum(int(v)*row[3] for v,row in zip(tail_integer,tail_rows))
        certificate = dict(numerator=numerator,denominator=scale,integer_lower=-(-numerator//scale),
                           point_multipliers=point_integer,resource_multipliers=resource_integer,
                           link_multipliers=link_integer.tolist(),
                           tail_multipliers=[dict(axis=row[0],thresholds=list(row[1:3]),bound=row[3],multiplier=int(v))
                                             for v,row in zip(tail_integer,tail_rows) if v])
        if best_certificate is None or numerator>best_certificate['numerator']:
            best_certificate = certificate
        print(json.dumps(dict(iteration=iterations,columns=len(columns),master_value=objective.Value(),priced_lower=certificate['integer_lower'],best_lower=best_certificate['integer_lower'],maximum_reduced_cost_violation=max_violation,seconds=round(time.monotonic()-start,3))),flush=True)
        if not prices:
            complete = True
            break
        added = sum(add_column(x,v) for x,v in prices)
        if not added:
            break
    return dict(core_residual=residual,high_primes=high,high_pattern=pattern,tail_capacity=tail,
                complete_pricing=complete,iterations=iterations,seconds=round(time.monotonic()-start,3),
                certificate=best_certificate,unit_mode=unit_mode,
                scope='exact dual lower for the projection/hit relaxation of one fixed core; not a global covering proof')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--incumbent',required=True)
    ap.add_argument('--concentration',required=True)
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--unit-mode',choices=('none','odd','both'),default='none')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    result = solve(json.loads(Path(a.incumbent).read_text()),json.loads(Path(a.concentration).read_text()),a.time,a.unit_mode)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
