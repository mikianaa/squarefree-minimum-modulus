"""Linked capacity relaxation with residual variables at every prime prefix.

Every cell variable represents the number of actual uncovered core points.
It lies below its parent count and above parent minus the capacities of newly
available classes. Known prefix lower bounds constrain sums of these counts.
The constraints are necessary, not sufficient, for a covering system.
Every supplied prefix floor must have a separate proof; the model does not
establish its input assumptions.
"""
import argparse
import itertools
import json
import math
import time
from pathlib import Path
from ortools.sat.python import cp_model


def solve(high,pattern,frontier,floors,seconds=60,workers=8,minimize=False):
    assert high and high == sorted(set(high))
    assert all(p>7 and all(p%d for d in range(2,math.isqrt(p)+1)) for p in high)
    assert len(pattern)==3+len(high) and pattern[0]==0 and all(b in (0,1) for b in pattern)
    assert frontier and len(floors)<=len(high) and all(v>=0 for v in floors)
    start = time.monotonic()
    mdl = cp_model.CpModel()
    dims = [[p-1-(pattern[i] == e) for i,p in enumerate((3,5,7))] for e in (0,1)]
    sizes = [math.prod(row) for row in dims]
    root = [mdl.NewIntVar(0,sizes[e],f'r{e}') for e in (0,1)]
    choose = [mdl.NewBoolVar(f'frontier{i}') for i in range(len(frontier))]
    mdl.AddExactlyOne(choose)
    for e in (0,1):
        mdl.Add(root[e] >= sum(row[e]*v for row,v in zip(frontier,choose)))
    alphabet = [[tuple(a for a in range(1,p) if a != 1 or pattern[i+3] != e)
                 for i,p in enumerate(high)] for e in (0,1)]
    resources = []
    for k in range(1,len(high)+1):
        for subset in itertools.combinations(range(len(high)),k):
            for d in range(1,211):
                if 210 % d or (k == 1 and d in (1,2)):
                    continue
                resources.append((d*math.prod(high[i] for i in subset),d,subset))
    resources.sort()
    choices = []
    label_vars = [[] for _ in high]
    for modulus,d,subset in resources:
        caps = [math.prod(dims[e][i] for i,p in enumerate((3,5,7)) if d % p) for e in (0,1)]
        lookup = {}
        entries = []
        for parity in ((0,1) if d % 2 == 0 else (None,)):
            labels_space = [alphabet[parity][i] if parity is not None else range(1,high[i]) for i in subset]
            for labels in itertools.product(*labels_space):
                if parity is None and not any(all(a in alphabet[e][i] for i,a in zip(subset,labels)) for e in (0,1)):
                    continue
                x = mdl.NewBoolVar(f'x{modulus}_{parity}_{labels}')
                lookup[(parity,labels)] = x
                entries.append((parity,labels,x))
        mdl.AddExactlyOne(x for _,_,x in entries)
        for j,i in enumerate(subset):
            label = mdl.NewIntVar(1,high[i]-1,f'label{modulus}_{i}')
            mdl.Add(label == sum(labels[j]*x for _,labels,x in entries))
            label_vars[i].append(label)
        choices.append((modulus,d,subset,caps,lookup,entries))
    for i,variables in enumerate(label_vars):
        previous = 1
        for j,label in enumerate(variables):
            mdl.Add(label <= previous+1)
            new = mdl.NewIntVar(1,high[i]-1,f'max{i}_{j}')
            mdl.AddMaxEquality(new,[previous,label])
            previous = new
    previous = {(e,()):root[e] for e in (0,1)}
    stages = []
    for j,q in enumerate(high):
        stage = {}
        active = [entry for entry in choices if max(entry[2]) == j]
        for e in (0,1):
            for point in itertools.product(*alphabet[e][:j+1]):
                u = mdl.NewIntVar(0,sizes[e],f'u{j}_{e}_{point}')
                parent = previous[(e,point[:-1])]
                terms = []
                for modulus,d,subset,caps,lookup,entries in active:
                    key = (e if d % 2 == 0 else None,tuple(point[i] for i in subset))
                    terms.append(caps[e]*lookup[key])
                mdl.Add(u <= parent)
                mdl.Add(u + sum(terms) >= parent)
                stage[(e,point)] = u
        if j < len(floors):
            mdl.Add(sum(stage.values()) >= floors[j])
        stages.append(stage)
        previous = stage
    if minimize:
        mdl.Minimize(sum(previous.values()))
    else:
        mdl.Add(sum(previous.values()) == 0)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.linearization_level = 2
    built = time.monotonic()
    status = solver.Solve(mdl)
    result = dict(N=210*math.prod(high),high_primes=high,pattern=pattern,prefix_floors=floors,
                  status=solver.StatusName(status),minimize=minimize,
                  build_seconds=round(built-start,3),seconds=round(time.monotonic()-start,3),
                  stage_cell_counts=[len(stage) for stage in stages],
                  scope='necessary capacity constraints; a feasible result is not a covering')
    if minimize:
        result['uncovered_lower_bound'] = math.ceil(solver.BestObjectiveBound()-1e-7)
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        selected_root = [solver.Value(x) for x in root]
        assert any(all(a<=b for a,b in zip(row,selected_root)) for row in frontier)
        chosen = []
        for modulus,d,subset,caps,lookup,entries in choices:
            parity,labels,_ = next(v for v in entries if solver.Value(v[2]))
            chosen.append(dict(modulus=modulus,core_divisor=d,high_subset=list(subset),
                               parity=parity,labels=list(labels),capacity=caps))
        counts = [{key:solver.Value(v) for key,v in stage.items()} for stage in stages]
        parent_counts = {(e,()):selected_root[e] for e in (0,1)}
        for j,stage in enumerate(counts):
            active = [v for v in chosen if max(v['high_subset']) == j]
            for (e,point),value in stage.items():
                parent = parent_counts[(e,point[:-1])]
                capacity = sum(v['capacity'][e] for v in active
                               if (v['parity'] is None or v['parity']==e)
                               and all(point[i]==a for i,a in zip(v['high_subset'],v['labels'])))
                assert 0 <= value <= parent and value+capacity >= parent
            if j < len(floors):
                assert sum(stage.values()) >= floors[j]
            parent_counts = stage
        result['selected_core_residual'] = selected_root
        result['stage_residual_totals'] = [sum(stage.values()) for stage in counts]
        result['placements'] = chosen
        result['stage_residuals'] = [[dict(parity=e,point=list(point),count=value) for (e,point),value in stage.items()] for stage in counts]
        if not minimize:
            assert sum(counts[-1].values()) == 0
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--high',type=int,nargs='+',required=True)
    ap.add_argument('--pattern',required=True)
    ap.add_argument('--floors',type=int,nargs='*',default=[])
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--minimize',action='store_true')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    pattern = list(map(int,a.pattern))
    cert = json.loads((Path(__file__).resolve().parents[1]/'data/three-axis-210-certificate.json').read_text())
    frontier = next(row['residual_pair_frontier'] for row in cert['cases'] if row['pattern']==pattern[:3])
    result = solve(a.high,pattern,frontier,a.floors,a.time,a.workers,a.minimize)
    print(json.dumps({k:v for k,v in result.items() if k not in ('placements','stage_residuals')}),flush=True)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
