"""Separate point event marginals from mixtures of complete core assignments.

For nonnegative event weights w, every assignment obeys sum w*event <= B,
where B is the maximum over shared resource residue choices. Column generation
finds w; CP-SAT prices assignments. Reported upper bounds require independent
verification before use in a rigorous obstruction certificate.
"""
import argparse
import json
from pathlib import Path
import time
from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model


def maximize(R,shapes,events,weights,seconds=2,collect=False,hint=None):
    assert len(events)==len(weights) and all(type(w) is int and w>=0 for w in weights)
    model = cp_model.CpModel()
    resources = []
    for d in shapes:
        row = [(a,model.NewBoolVar(f'r{d}_{a}')) for a in sorted({v%d for v in R})]
        model.AddExactlyOne(v for a,v in row)
        resources.append((d,row))
    hint_row = None
    if hint is not None:
        assert len(hint)==len(shapes) and all(type(a) is int and 0<=a<d for a,d in zip(hint,shapes))
        # A residue missing every point can be replaced by an occupied one:
        # all event weights are nonnegative, so this cannot lose score.
        normalized = [a if any(a==b for b,var in row) else row[0][0]
                      for a,(d,row) in zip(hint,resources)]
        counts = [sum(v%d==a for d,a in zip(shapes,normalized)) for v in R]
        hint_row = dict(residues=normalized,counts=counts,
                        objective=sum(w for (x,t),w in zip(events,weights) if counts[x]>=t))
        for a,(d,row) in zip(normalized,resources):
            for b,var in row:model.AddHint(var,int(a==b))
    hits = []
    for x,v in enumerate(R):
        count = model.NewIntVar(0,len(shapes),f'h{x}')
        model.Add(count==sum(var for d,row in resources for a,var in row if v%d==a))
        hits.append(count)
    terms = []
    for j,((x,t),w) in enumerate(zip(events,weights)):
        if w:
            event = model.NewBoolVar(f'e{j}')
            model.Add(hits[x]>=t).OnlyEnforceIf(event)
            model.Add(hits[x]<t).OnlyEnforceIf(event.Not())
            terms.append(w*event)
    model.Maximize(sum(terms))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.max_time_in_seconds = seconds
    class Collector(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.rows = []
            self.seen = set()
        def on_solution_callback(self):
            residues = tuple(next(a for a,var in row if self.Value(var)) for d,row in resources)
            counts = tuple(sum(v%d==a for d,a in zip(shapes,residues)) for v in R)
            if counts not in self.seen:
                self.seen.add(counts)
                self.rows.append(dict(residues=list(residues),counts=list(counts),
                                      objective=sum(w for (x,t),w in zip(events,weights) if counts[x]>=t)))
    collector = Collector() if collect else None
    status = solver.Solve(model,collector) if collector else solver.Solve(model)
    if status not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        return dict(status=solver.StatusName(status),candidate_assignments=[hint_row] if hint_row else [])
    residues = [next(a for a,var in row if solver.Value(var)) for d,row in resources]
    counts = [sum(v%d==a for d,a in zip(shapes,residues)) for v in R]
    value = sum(w for (x,t),w in zip(events,weights) if counts[x]>=t)
    assert value==round(solver.ObjectiveValue())
    result = dict(status=solver.StatusName(status),counts=counts,residues=residues,
                  objective=value,upper_bound=int(solver.BestObjectiveBound()))
    if collector:
        result['candidate_assignments'] = collector.rows
    if hint_row:
        assert hint_row['objective']<=result['upper_bound']
        result.setdefault('candidate_assignments',[]).append(hint_row)
        if hint_row['objective']>value:result.update(hint_row)
    return result


def separate_axis(data,h,seconds=15):
    start = time.monotonic()
    source,witness = data['input'],data['result']['point_primal_witness']
    R,n,den = source['core_residual'],len(source['high_primes']),witness['denominator']
    shapes = [d for d in range(1 if h>=n else 3,211) if 210%d==0]
    hist = [[0]*(len(shapes)+1) for x in R]
    for row in witness['support']:
        hist[row['point_index']][row['counts'][h]] += row['weight']
    events = [(x,t) for x,line in enumerate(hist) for t,w in enumerate(line) if t>0 and w>0]
    targets = [sum(hist[x][t:]) for x,t in events]
    lp = pywraplp.Solver.CreateSolver('GLOP')
    weights = [lp.NumVar(0,1,f'w{j}') for j in range(len(events))]
    beta = lp.NumVar(0,1,'beta')
    lp.Add(sum(weights)<=1)
    objective = lp.Objective()
    objective.SetCoefficient(beta,-1)
    for variable,target in zip(weights,targets):
        objective.SetCoefficient(variable,target/den)
    objective.SetMaximization()
    profiles = []
    best = None
    iterations = 0
    scale = 10**6
    while time.monotonic()-start<seconds:
        if lp.Solve()!=pywraplp.Solver.OPTIMAL:
            break
        iterations += 1
        if objective.Value()<1e-9:
            break
        integers = [max(0,round(w.solution_value()*scale)) for w in weights]
        priced = maximize(R,shapes,events,integers,min(2,max(.1,seconds-time.monotonic()+start)))
        if 'counts' not in priced:
            break
        gap = sum(w*t for w,t in zip(integers,targets))-priced['upper_bound']*den
        if gap>0:
            best = dict(axis=h,events=[dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,integers) if w],
                        bound=priced['upper_bound'],violation_numerator=gap,
                        witness_denominator=den,weight_scale=scale,pricing=priced)
            break
        counts = priced['counts']
        if counts in profiles:
            break
        profiles.append(counts)
        lp.Add(sum(w for w,(x,t) in zip(weights,events) if counts[x]>=t)<=beta)
    return dict(axis=h,iterations=iterations,seconds=time.monotonic()-start,cut=best,
                final_restricted_gap=objective.Value(),profiles=profiles,
                scope='CP-SAT upper bound requires independent audit; no covering obstruction established')


def simplify_candidate(data,cut,scales=(16,24,32,48,64,96,128,192,256,384,512),seconds=3):
    source,witness = data['input'],data['result']['point_primal_witness']
    h,R = cut['axis'],source['core_residual']
    shapes = [d for d in range(1 if h>=len(source['high_primes']) else 3,211) if 210%d==0]
    events = [(r['point'],r['threshold']) for r in cut['events']]
    targets = [sum(r['weight'] for r in witness['support'] if r['point_index']==x and r['counts'][h]>=t) for x,t in events]
    for scale in scales:
        weights = [round(r['weight']*scale/cut['weight_scale']) for r in cut['events']]
        if not any(weights):
            continue
        priced = maximize(R,shapes,events,weights,seconds)
        if 'upper_bound' not in priced:
            continue
        gap = sum(w*t for w,t in zip(weights,targets))-witness['denominator']*priced['upper_bound']
        if gap>0:
            return dict(axis=h,events=[dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,weights) if w],
                        bound=priced['upper_bound'],violation_numerator=gap,denominator=witness['denominator'],
                        simplifying_scale=scale,pricing=priced)
    return None


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--axes',default='8')
    ap.add_argument('--time',type=float,default=15)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    results = []
    for h in map(int,args.axes.split(',')):
        row = separate_axis(data,h,args.time)
        results.append(row)
        Path(args.output).write_text(json.dumps(dict(source=args.source,results=results),indent=2)+'\n')
        print(json.dumps({k:v for k,v in row.items() if k!='profiles'}),flush=True)
