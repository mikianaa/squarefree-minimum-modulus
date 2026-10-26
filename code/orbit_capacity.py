"""Exact dual obstructions from a symmetrized fractional capacity model.

Higher-coordinate labels 2,...,p-1 are interchangeable; label 1 is distinguished
by the twice-prime anchor. Average the fractional placement model over these
permutations. Its cell orbits have at most nine types for three high primes.
LP output suggests dual weights only. Every exclusion is recomputed with exact
Fractions and integer orbit sizes, without trusting the solver's optimality.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import time
from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model

CORE_ODD = (1,3,5,7,15,21,35,105)


def model_data(profile,high,pattern):
    assert len(profile)==16 and len(high)==len(pattern)
    assert all(p>2 for p in high) and all(b in (0,1) for b in pattern)
    # Type 0 means label 1; type 1 means any of the p-2 other surviving labels.
    cells = [(e,t) for e in (0,1) if profile[e]
             for t in itertools.product((0,1),repeat=len(high))
             if all(t[i] or pattern[i]!=e for i in range(len(high)))]
    demand = [profile[e]*math.prod(high[i]-2 for i,v in enumerate(t) if v) for e,t in cells]
    resources = []
    for k in range(1,len(high)+1):
        for subset in itertools.combinations(range(len(high)),k):
            for j,d in enumerate(CORE_ODD):
                if k==1 and d==1:
                    continue
                cap = profile[:2] if d==1 else profile[2+2*(j-1):4+2*(j-1)]
                for even in (False,True):
                    choices = []
                    for parity in ((0,1) if even else (None,)):
                        for labels in itertools.product((0,1),repeat=k):
                            column = []
                            for e,t in cells:
                                allowed = (parity is None or parity==e) and all(t[i]==v for i,v in zip(subset,labels))
                                column.append(cap[e]*math.prod(high[i]-2 for i,v in enumerate(t) if v and i not in subset) if allowed else 0)
                            if any(column):
                                choices.append(dict(parity=parity,types=labels,column=tuple(column)))
                    if choices:
                        resources.append(dict(modulus=d*(2 if even else 1)*math.prod(high[i] for i in subset),
                                              core_odd=d,subset=subset,capacity=cap,even=even,choices=choices))
    return cells,demand,resources


def exact_gap(demand,resources,weights):
    """Weak duality certificate: weighted demand minus maximum resource capacity."""
    assert len(weights)==len(demand) and all(w>=0 for w in weights)
    rhs = sum(max(sum(w*c for w,c in zip(weights,v['column'])) for v in r['choices']) for r in resources)
    return sum(w*d for w,d in zip(weights,demand))-rhs


def solve(profile,high,pattern,keep_primal=False):
    cells,demand,resources = model_data(profile,high,pattern)
    solver = pywraplp.Solver.CreateSolver('GLOP')
    rows = [solver.Constraint(d,solver.infinity()) for d in demand]
    deficit = []
    for i,row in enumerate(rows):
        u = solver.NumVar(0,solver.infinity(),f'u{i}')
        row.SetCoefficient(u,1)
        deficit.append(u)
    all_variables = []
    for j,r in enumerate(resources):
        budget = solver.Constraint(0,1)
        variables = []
        for k,v in enumerate(r['choices']):
            x = solver.NumVar(0,solver.infinity(),f'x{j}_{k}')
            budget.SetCoefficient(x,1)
            for row,c in zip(rows,v['column']):
                if c:
                    row.SetCoefficient(x,c)
            variables.append(x)
        all_variables.append(variables)
    objective = solver.Objective()
    for u in deficit:
        objective.SetCoefficient(u,1)
    objective.SetMinimization()
    status = solver.Solve()
    weights = [Fraction(max(0,row.dual_value())).limit_denominator(10**7) for row in rows]
    gap = exact_gap(demand,resources,weights)
    out = dict(status=int(status),float_objective=objective.Value(),excluded=gap>0,
               exact_gap=[gap.numerator,gap.denominator],
               weights=[[w.numerator,w.denominator] for w in weights],
               cells=cells,resources=len(resources),variables=sum(map(len,all_variables)),
               scope='fractional optimistic-capacity model; only positive exact dual gaps exclude')
    if keep_primal and status==solver.OPTIMAL:
        primal = [[Fraction(x.solution_value()).limit_denominator(10**9) for x in vs] for vs in all_variables]
        capacity = [sum(sum(x*v['column'][i] for x,v in zip(xs,r['choices'])) for xs,r in zip(primal,resources))
                    for i in range(len(demand))]
        verified = all(all(x>=0 for x in xs) and sum(xs)<=1 for xs in primal) and all(c>=d for c,d in zip(capacity,demand))
        out['exact_fractional_feasibility_verified'] = verified
        if verified:
            out['primal'] = [[[x.numerator,x.denominator] for x in xs] for xs in primal]
    return out


def solve_integer(profile,high,pattern,seconds=1):
    """Necessary integral orbit totals, still ignoring positions within orbits."""
    cells,demand,resources = model_data(profile,high,pattern)
    mdl = cp_model.CpModel()
    terms = [[] for _ in cells]
    variables = []
    for j,r in enumerate(resources):
        xs = [mdl.NewBoolVar(f'x{j}_{k}') for k in range(len(r['choices']))]
        mdl.AddExactlyOne(xs)
        for x,v in zip(xs,r['choices']):
            for i,c in enumerate(v['column']):
                if c:
                    terms[i].append(c*x)
        variables.append(xs)
    for row,d in zip(terms,demand):
        mdl.Add(sum(row)>=d)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = 1
    status = solver.Solve(mdl)
    out = dict(status=solver.StatusName(status),seconds=solver.WallTime(),
               scope='integer orbit totals; feasibility is not a covering system')
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        selected = [next(k for k,x in enumerate(xs) if solver.Value(x)) for xs in variables]
        coverage = [sum(r['choices'][k]['column'][i] for r,k in zip(resources,selected)) for i in range(len(cells))]
        assert all(c>=d for c,d in zip(coverage,demand))
        out.update(selected_choices=selected,exact_aggregate_feasibility_verified=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/core-geometry-profile-audit.json')
    ap.add_argument('--output',default='data/orbit-capacity-510510.json')
    ap.add_argument('--integer',action='store_true')
    ap.add_argument('--fractional-source',default='data/orbit-capacity-510510.json')
    ap.add_argument('--time',type=float,default=1)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    fractional = json.loads(Path(args.fractional_source).read_text()) if args.integer else None
    results = []
    for row in data:
        for pattern in itertools.product((0,1),repeat=3):
            start = time.monotonic()
            if args.integer:
                prior = next(x for x in fractional if x['core_pattern']==row['pattern'] and tuple(x['high_pattern'])==pattern)
                tests = [dict(profile_index=i,**solve_integer(row['easiest_profiles'][i],(11,13,17),pattern,args.time))
                         for i in prior['surviving_indices']]
                result = dict(core_pattern=row['pattern'],high_pattern=pattern,profiles=len(tests),
                              statuses={s:sum(t['status']==s for t in tests) for s in ('OPTIMAL','FEASIBLE','INFEASIBLE','UNKNOWN')},
                              seconds=round(time.monotonic()-start,3),tests=tests)
                results.append(result)
                print(json.dumps({k:v for k,v in result.items() if k!='tests'}),flush=True)
                Path(args.output).write_text(json.dumps(results,indent=2)+'\n')
                continue
            tests = [solve(p,(11,13,17),pattern) for p in row['easiest_profiles']]
            result = dict(core_pattern=row['pattern'],high_pattern=pattern,profiles=len(tests),
                          excluded=sum(t['excluded'] for t in tests),
                          surviving_indices=[i for i,t in enumerate(tests) if not t['excluded']],
                          seconds=round(time.monotonic()-start,3),tests=tests)
            results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('surviving_indices','tests')}),flush=True)
            Path(args.output).write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    main()
