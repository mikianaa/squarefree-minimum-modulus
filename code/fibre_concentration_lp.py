"""Transport relaxation for core-fibre plane residuals with hit-tail bounds.

Each column is a possible four-axis hit-count tuple for one core point.
Subset-capacity bounds constrain excess hits jointly across the parity layers.
Fractional core points are allowed, so a feasible solution is not a covering.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp


def solve(data,all_thresholds=False):
    start = time.monotonic()
    high = (11,13,17,19)
    pattern = (1,1,1,1)
    sizes = data['layer_sizes']
    H = np.array(data['two_layer_capacities'],dtype=np.int64)
    maximum = len(data['moduli'])
    tuples = []
    costs = []
    parity = []
    for e in (0,1):
        D = np.array([q-1-(b==e) for q,b in zip(high,pattern)],dtype=np.int64)
        rows = np.array(list(itertools.product(*(range(min(maximum,int(d))+1) for d in D))),dtype=np.int16)
        tuples.append(rows)
        costs.extend(np.prod(D[None,:]-rows,axis=1).tolist())
        parity.extend([e]*len(rows))
    counts = np.concatenate(tuples)
    parity = np.array(parity,dtype=np.int8)
    costs = np.array(costs,dtype=np.int64)
    pairs = sorted({(t,t) for t in range(maximum)}|{(t,maximum) for t in range(maximum)}|{(maximum,t) for t in range(maximum)})
    coefficients = []
    rhs,descriptions = [],[]
    m0,m1 = np.indices(H.shape)
    for axis in range(len(high)):
        for t0,t1 in pairs:
            excess = np.maximum(0,counts[:,axis]-np.where(parity==0,t0,t1))
            indices = np.flatnonzero(excess)
            if not len(indices):
                continue
            coefficients.append(excess.astype(np.int16))
            rhs.append(int(np.max(H-t0*m0-t1*m1)))
            descriptions.append([axis,t0,t1])
    matrix = np.stack(coefficients)
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetTimeLimit(60000)
    variables = [solver.NumVar(0,solver.infinity(),f'z{j}') for j in range(len(costs))]
    constraints = []
    for row,b in zip(matrix,rhs):
        constraint = solver.Constraint(-solver.infinity(),b)
        for j in np.flatnonzero(row):
            constraint.SetCoefficient(variables[j],int(row[j]))
        constraints.append(constraint)
    for e in (0,1):
        constraint = solver.Constraint(sizes[e],sizes[e])
        for j in np.flatnonzero(parity==e):
            constraint.SetCoefficient(variables[j],1)
    objective = solver.Objective()
    for variable,c in zip(variables,costs):
        objective.SetCoefficient(variable,int(c))
    objective.SetMinimization()
    print(json.dumps(dict(stage='built',variables=len(costs),constraints=len(rhs),nonzeros=int(np.count_nonzero(matrix)),seconds=round(time.monotonic()-start,3))),flush=True)
    status = solver.Solve()
    cut_iterations = 0
    complete_threshold_check = not all_thresholds
    existing = set(map(tuple,descriptions))
    while all_thresholds and status==pywraplp.Solver.OPTIMAL and cut_iterations<20 and time.monotonic()-start<60:
        support = [(j,v.solution_value()) for j,v in enumerate(variables) if v.solution_value()>1e-9]
        violating = []
        for axis in range(len(high)):
            for t0,t1 in itertools.product(range(maximum+1),repeat=2):
                if (axis,t0,t1) in existing:
                    continue
                b = int(np.max(H-t0*m0-t1*m1))
                value = sum(mass*max(0,int(counts[j,axis])-(t0 if parity[j]==0 else t1)) for j,mass in support)
                if value>b+1e-6:
                    violating.append((value-b,axis,t0,t1,b))
        if not violating:
            complete_threshold_check = True
            break
        for _,axis,t0,t1,b in sorted(violating,reverse=True)[:16]:
            row = np.maximum(0,counts[:,axis]-np.where(parity==0,t0,t1)).astype(np.int16)
            constraint = solver.Constraint(-solver.infinity(),b)
            for j in np.flatnonzero(row):
                constraint.SetCoefficient(variables[j],int(row[j]))
            constraints.append(constraint)
            coefficients.append(row)
            rhs.append(b)
            descriptions.append([axis,t0,t1])
            existing.add((axis,t0,t1))
        cut_iterations += 1
        solver.SetTimeLimit(max(1,int(1000*(60-(time.monotonic()-start)))))
        status = solver.Solve()
        print(json.dumps(dict(stage='mixed_threshold_cuts',iteration=cut_iterations,constraints=len(rhs),status=int(status),value=objective.Value() if status==pywraplp.Solver.OPTIMAL else None)),flush=True)
    matrix = np.stack(coefficients)
    success = status==pywraplp.Solver.OPTIMAL
    out = dict(success=success,status=int(status),
               variables=len(costs),constraints=len(rhs),seconds=round(time.monotonic()-start,3),
               mixed_thresholds=all_thresholds,complete_threshold_check=complete_threshold_check,
               scope='fractional hit-distribution relaxation; does not represent an actual covering')
    if success:
        out['float_lower'] = float(objective.Value())
        # Round to integer multipliers, then recompute the largest valid equality
        # multiplier in each layer over ALL columns. This enforces exact dual
        # feasibility without relying on the floating solver's tolerances.
        scale = 10**6
        multipliers = np.minimum(0,np.rint(np.array([c.dual_value() for c in constraints])*scale).astype(np.int64))
        column_terms = np.asarray(matrix.astype(np.int64).T@multipliers).ravel()
        slacks = costs*scale-column_terms
        equal_multipliers = [int(slacks[parity==e].min()) for e in (0,1)]
        numerator = sum(int(v)*int(b) for v,b in zip(multipliers,rhs))+sum(a*b for a,b in zip(equal_multipliers,sizes))
        out.update(exact_dual_numerator=numerator,exact_dual_denominator=scale,
                   integer_lower=-(-numerator//scale),equality_multipliers=equal_multipliers,
                   active_multipliers=[dict(axis=d[0],thresholds=d[1:],bound=b,multiplier=int(v))
                                       for d,b,v in zip(descriptions,rhs,multipliers) if v],
                   primal_support=[dict(parity=int(parity[j]),counts=counts[j].tolist(),mass=variable.solution_value())
                                   for j,variable in enumerate(variables) if variable.solution_value()>1e-8])
    return out


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--all-thresholds',action='store_true')
    a = ap.parse_args()
    result = solve(json.loads(Path(a.source).read_text()),a.all_thresholds)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('active_multipliers','primal_support')}))
