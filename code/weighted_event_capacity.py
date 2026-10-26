"""Find threshold incompatibilities with point-weighted residue capacities.

Every assignment of one class per modulus obeys
sum_x w_x z(x) <= sum_d max_a sum_{x == a (mod d)} w_x.
Thus thresholds exceeding this capacity cannot all occur. Optimization only
discovers weights; an integer-only auditor verifies each resulting cut.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import time


def capacities(points, weights, moduli):
    result = []
    for d in moduli:
        row = Counter()
        for x, w in zip(points, weights):
            row[x % d] += w
        result.append(max(row.values(), default=0))
    return result


def separate(data, seconds_per_axis=2, weight_limit=3, axes=None):
    from ortools.sat.python import cp_model
    from audit_multitail_identity import audit_group_cut
    start = time.monotonic()
    source, witness = data['input'], data['result']['point_primal_witness']
    R, n, scale = source['core_residual'], len(source['high_primes']), witness['denominator']
    group_count = (1 << n)-1
    hist = [[[0]*17 for h in range(group_count)] for x in R]
    for row in witness['support']:
        for h, v in enumerate(row['counts']):
            hist[row['point_index']][h][v] += row['weight']
    tails = [[[sum(line[t:]) for t in range(17)] for line in point] for point in hist]
    cuts, reports = [], []
    for h in (range(group_count) if axes is None else axes):
        moduli = [d for d in range(3 if h<n else 1, 211) if 210 % d == 0]
        model = cp_model.CpModel()
        event_rows, point_weights, deficits, selected = [], [], [], []
        for x in range(len(R)):
            # Only distribution breakpoints need be considered: raising a
            # threshold within a zero-mass interval strengthens the capacity
            # contradiction without changing the event probability.
            row = [(t, w, model.NewBoolVar(f'e{x}_{t}_{w}'))
                   for t in range(1, len(moduli)+1) if hist[x][h][t]>0
                   for w in range(1, weight_limit+1)]
            model.Add(sum(v for t,w,v in row)<=1)
            event_rows.append(row)
            point_weights.append(sum(w*v for t,w,v in row))
            selected.extend(v for t,w,v in row)
            deficits.extend((scale-tails[x][h][t])*v for t,w,v in row)
        resource_bounds = []
        for d in moduli:
            bound = model.NewIntVar(0, weight_limit*len(R), f'c{d}')
            for a in {v % d for v in R}:
                model.Add(bound>=sum(point_weights[x] for x,v in enumerate(R) if v%d==a))
            resource_bounds.append(bound)
        model.Add(sum(selected)>=2)
        model.Add(sum(t*w*v for row in event_rows for t,w,v in row)>sum(resource_bounds))
        model.Add(sum(deficits)<scale)
        model.Minimize(sum(deficits))
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 1
        solver.parameters.max_time_in_seconds = seconds_per_axis
        status = solver.Solve(model)
        reports.append(dict(axis=h,status=solver.StatusName(status),seconds=solver.WallTime()))
        if status not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            continue
        events, weights = [], []
        for x,row in enumerate(event_rows):
            for t,w,v in row:
                if solver.Value(v):
                    events.append(dict(point=x,threshold=t,weight=1))
                    weights.append(w)
        bounds = capacities([R[e['point']] for e in events],weights,moduli)
        cut = dict(axis=h,events=events,bound=len(events)-1,
                   capacity_certificate=dict(point_weights=weights,resource_capacities=bounds),
                   proof_audit=dict(kind='weighted_capacity',verified=True,axis=h,bound=len(events)-1))
        audit_group_cut(source,cut)
        cut['denominator'] = scale
        cut['violation_numerator'] = sum(tails[e['point']][h][e['threshold']] for e in events)-scale*cut['bound']
        assert cut['violation_numerator']>0
        cuts.append(cut)
    return dict(core_residual=R,cuts=cuts,axis_reports=reports,weight_limit=weight_limit,
                seconds=time.monotonic()-start,
                scope='arithmetic-certified cuts; solver absence of further cuts is exploratory')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=2)
    ap.add_argument('--weight-limit',type=int,default=3)
    ap.add_argument('--axes')
    args = ap.parse_args()
    result = separate(json.loads(Path(args.source).read_text()),args.time,args.weight_limit,
                      list(map(int,args.axes.split(','))) if args.axes else None)
    result['source'] = args.source
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(cuts=len(result['cuts']),seconds=result['seconds'],
                         statuses=dict(Counter(r['status'] for r in result['axis_reports'])),
                         gaps=[c['violation_numerator'] for c in result['cuts']])),flush=True)
