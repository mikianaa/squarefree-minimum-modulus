"""Coupled integer optimization of actual core projections.

The model includes all local conditional-tail expressions, with shared
residue choices across points. Solver bounds are diagnostics, never proof
certificates. Returned assignments can be checked without the solver.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
from ortools.sat.python import cp_model
from search_joint_core_projections import JointProjections


def build(initial, axes=None, tables=True):
    source = initial['input']
    problem = JointProjections(source)
    assignments = initial['group_assignments']
    counts0 = problem.counts(assignments)
    initial_value = sum(max(int(v), 0) for v in problem.evaluate(counts0))
    assert initial_value == initial['positive_volume']
    active = set(range(len(problem.groups))) if axes is None else set(axes)
    assert active and active <= set(range(len(problem.groups)))
    model = cp_model.CpModel()
    choices = {}
    counts = []
    for h, shape in enumerate(problem.shapes):
        if h not in active:
            counts.append([(int(k), int(k), int(k), int(k)) for k in counts0[h]])
            continue
        rows = []
        for j, d in enumerate(shape):
            # The incumbent may use a residue missing every core point.
            # Such a choice is dominated for positive F, but is retained in
            # the model so its complete hint remains feasible and exact.
            old = assignments[h][j]
            residues = sorted({int(p % d) for p in problem.points} | {old})
            row = [(a, model.NewBoolVar(f'c{h}_{d}_{a}')) for a in residues]
            model.AddExactlyOne(v for a, v in row)
            for a, v in row:
                model.AddHint(v, int(a == old))
            rows.append(row)
        choices[h] = rows
        line = []
        for x, p in enumerate(problem.points):
            k = model.NewIntVar(0, len(shape), f'k{h}_{x}')
            model.Add(k == sum(v for d, row in zip(shape, rows) for a, v in row if p % d == a))
            model.AddHint(k, int(counts0[h, x]))
            line.append((k, 0, len(shape), int(counts0[h, x])))
        counts.append(line)

    def multiply(left, right, name):
        a, alo, ahi, ahint = left
        b, blo, bhi, bhint = right
        if ahi == 0 or bhi == 0:
            return (0, 0, 0, 0)
        if alo == ahi == 1:
            return right
        if blo == bhi == 1:
            return left
        lo, hi, hint = alo*blo, ahi*bhi, ahint*bhint
        if lo == hi:
            return (lo, lo, lo, lo)
        v = model.NewIntVar(lo, hi, name)
        model.AddMultiplicationEquality(v, [a, b])
        model.AddHint(v, hint)
        return (v, lo, hi, hint)

    objectives = []
    use_tables = tables and len(active) <= 2
    if use_tables:
        ordered = sorted(active)
        states = []
        for ks in itertools.product(*(range(len(problem.shapes[h])+1) for h in ordered)):
            state = counts0.copy()
            for h, k in zip(ordered, ks):
                state[h] = k
            states.append((ks, [max(int(v), 0) for v in problem.evaluate(state)]))
        for x in range(len(problem.points)):
            upper = math.prod(map(int, problem.D[:, x]))
            z = model.NewIntVar(0, upper, f'z{x}')
            model.AddAllowedAssignments([counts[h][x][0] for h in ordered]+[z],
                                        [list(ks)+[values[x]] for ks, values in states])
            model.AddHint(z, max(int(problem.evaluate(counts0)[x]), 0))
            objectives.append(z)
    else:
        for x in range(len(problem.points)):
            free = []
            for h in range(problem.n):
                D = int(problem.D[h, x])
                k, lo, hi, hint = counts[h][x]
                f_hint = max(0, D-hint)
                if lo == hi:
                    free.append((f_hint, f_hint, f_hint, f_hint))
                else:
                    f = model.NewIntVar(max(0, D-hi), D, f'f{h}_{x}')
                    model.AddMaxEquality(f, [0, D-k])
                    model.AddHint(f, f_hint)
                    free.append((f, max(0, D-hi), D, f_hint))
            products = [(1, 1, 1, 1)]
            for mask in range(1, 1 << problem.n):
                bit = mask & -mask
                products.append(multiply(products[mask ^ bit], free[bit.bit_length()-1], f'p{x}_{mask}'))
            terms = [multiply(counts[h][x], products[problem.full ^ mask], f't{x}_{h}')
                     for h, mask in enumerate(problem.masks[problem.n:], problem.n)]
            expression = products[-1][0]-sum(t[0] for t in terms)
            z = model.NewIntVar(0, products[-1][2], f'z{x}')
            model.AddMaxEquality(z, [0, expression])
            model.AddHint(z, max(products[-1][3]-sum(t[3] for t in terms), 0))
            objectives.append(z)
    model.Add(sum(objectives) <= initial_value)
    model.Minimize(sum(objectives))
    assert not model.Validate(), model.Validate()
    return problem, model, choices, objectives, sorted(active), use_tables


def solve(initial, seconds=120, axes=None, workers=1, tables=True):
    start = time.monotonic()
    problem, model, choices, objectives, active, used_tables = build(initial, axes, tables)
    print(json.dumps(dict(stage='prepared', active_groups=active, tables=used_tables,
                          variables=len(model.Proto().variables), constraints=len(model.Proto().constraints),
                          seconds=time.monotonic()-start)), flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 1
    class Progress(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.history = []
        def on_solution_callback(self):
            value = sum(self.Value(z) for z in objectives)
            self.history.append(dict(positive_volume=value, seconds=time.monotonic()-start))
            print(json.dumps(self.history[-1]), flush=True)
    progress = Progress()
    status = solver.Solve(model, progress)
    assignments = [row[:] for row in initial['group_assignments']]
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        for h, rows in choices.items():
            assignments[h] = [next(a for a, v in row if solver.Value(v)) for row in rows]
    counts = problem.counts(assignments)
    values = list(map(int, problem.evaluate(counts)))
    value = sum(max(v, 0) for v in values)
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        assert value == round(solver.ObjectiveValue())
    assert value <= initial['positive_volume']
    return dict(input=initial['input'], group_assignments=assignments, local_F=values, positive_volume=value,
                verified=False, optimization=dict(status=solver.StatusName(status), active_groups=active,
                tables=used_tables, initial_positive_volume=initial['positive_volume'],
                numerical_lower_bound=solver.BestObjectiveBound(), seconds=time.monotonic()-start,
                workers=workers, history=progress.history, solver_bounds_used_as_proofs=False),
                scope='Coupled actual core projections for the stated groups; high-prime labels omitted. An incumbent is not a covering; solver bounds are not independently verified impossibility proofs.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--time', type=float, default=120)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--axes', default='all', help='all, tail, or comma-separated zero-based groups')
    ap.add_argument('--no-tables', action='store_true')
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    n = len(data['input']['high_primes'])
    axes = None if args.axes == 'all' else (list(range(n, 2**n-1)) if args.axes == 'tail' else list(map(int, args.axes.split(','))))
    result = solve(data, args.time, axes, args.workers, not args.no_tables)
    result['provenance'] = dict(source=args.source)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('input', 'group_assignments')}))
