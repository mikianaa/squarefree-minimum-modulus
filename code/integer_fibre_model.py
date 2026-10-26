"""Integer projection model with audited concentration and pointwise dual cuts.

All choices are core projections only. A feasible solution is not a covering
system. Even an infeasibility result applies only to the specified core.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
from ortools.sat.python import cp_model
from core_fibre_plane_model import prepare


def solve(incumbent, concentration, dual, seconds=60, workers=4):
    high, pattern = incumbent['high_primes'], incumbent['high_pattern']
    _, R, dims, resources, tail = prepare(dict(system=incumbent['core_system']), high, pattern)
    assert R == concentration['core_residual'] == dual['core_residual']
    certificate = dual['certificate']
    model = cp_model.CpModel()
    by_point = [[[] for _ in R] for _ in high]
    choices = []
    hints = {(r['high_prime'], r['core_modulus']): r['core_residue'] for r in incumbent['choices']}
    for j, r in enumerate(resources):
        vs = [model.NewBoolVar(f'w{j}_{a}') for a, pts in r['options']]
        model.AddExactlyOne(vs)
        for v, (a, pts) in zip(vs, r['options']):
            model.AddHint(v, int(a == hints[high[r['axis']], r['core_modulus']]))
            for x in pts:
                by_point[r['axis']][x].append(v)
        choices.append(vs)
    counts, factors, excess = [], [], []
    for x, value in enumerate(R):
        nc, fc, ec = [], [], []
        for i in range(len(high)):
            D = dims[i][x]
            maximum = min(14, D)
            n = model.NewIntVar(0, maximum, f'n{x}_{i}')
            model.AddMinEquality(n, [D, sum(by_point[i][x])])
            f = model.NewIntVar(D-maximum, D, f'f{x}_{i}')
            model.Add(f == D-n)
            tails = []
            for t in range(15):
                if t == 0:
                    z = n
                elif t >= maximum:
                    z = 0
                else:
                    z = model.NewIntVar(0, maximum-t, f'e{x}_{i}_{t}')
                    model.AddMaxEquality(z, [0, n-t])
                tails.append(z)
            nc.append(n)
            fc.append(f)
            ec.append(tails)
        counts.append(nc)
        factors.append(fc)
        excess.append(ec)
    H = concentration['two_layer_capacities']
    for i in range(len(high)):
        for t0, t1 in itertools.product(range(15), repeat=2):
            bound = max(v-t0*a-t1*b for a, row in enumerate(H) for b, v in enumerate(row))
            model.Add(sum(excess[x][i][t0 if value % 2 == 0 else t1] for x, value in enumerate(R)) <= bound)
    products = []
    for x, value in enumerate(R):
        pair_products = []
        for i, j in ((0, 1), (2, 3)):
            p = model.NewIntVar(0, dims[i][x]*dims[j][x], f'pair{x}_{i}')
            rows = [(a, b, a*b) for a in range(max(0, dims[i][x]-14), dims[i][x]+1)
                    for b in range(max(0, dims[j][x]-14), dims[j][x]+1)]
            model.AddAllowedAssignments([factors[x][i], factors[x][j], p], rows)
            pair_products.append(p)
        p = model.NewIntVar(0, math.prod(d[x] for d in dims), f'product{x}')
        model.AddMultiplicationEquality(p, pair_products)
        rhs = certificate['point_multipliers'][x] + sum(c*n for c,n in zip(certificate['link_multipliers'][x],counts[x]))
        for row in certificate['tail_multipliers']:
            rhs += row['multiplier'] * excess[x][row['axis']][row['thresholds'][value % 2]]
        model.Add(certificate['denominator']*p >= rhs)
        products.append(p)
    model.Add(sum(products) >= certificate['integer_lower'])
    model.Minimize(sum(products))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.linearization_level = 2
    start = time.monotonic()
    status = solver.Solve(model)
    out = dict(status=solver.StatusName(status), seconds=round(time.monotonic()-start,3),
               tail_capacity=tail, scope='integer projections for one specified core; no covering claim')
    if status in (cp_model.UNKNOWN, cp_model.FEASIBLE, cp_model.OPTIMAL):
        out['solver_lower'] = math.ceil(solver.BestObjectiveBound()-1e-7)
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        out['choices'] = [dict(high_prime=high[r['axis']], core_modulus=r['core_modulus'],
                              core_residue=next(a for (a,pts),v in zip(r['options'],vs) if solver.Value(v)))
                          for r,vs in zip(resources,choices)]
        direct = sum(math.prod(max(0,q-1-(b==x%2)-sum(x%r['core_modulus']==r['core_residue']
                        for r in out['choices'] if r['high_prime']==q)) for q,b in zip(high,pattern)) for x in R)
        assert direct == round(solver.ObjectiveValue())
        out['plane_expression_upper'] = direct
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--time', type=float, default=60)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    read = lambda name: json.loads(Path('data',name).read_text())
    result = solve(read('core-fibre-eight-plane.json'), read('core-hit-concentration-7-15.json'),
                   read('identified-fibre-eight-lp.json'), args.time)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='choices'}))
