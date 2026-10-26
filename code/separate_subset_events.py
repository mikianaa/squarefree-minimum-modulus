"""Separate threshold-event cuts on arbitrary core subsets, not covers."""
import argparse
from collections import Counter
import json
from pathlib import Path
import time
from ortools.sat.python import cp_model


def separate(data,seconds_per_axis=2):
    start = time.monotonic()
    source,witness = data['input'],data['result']['point_primal_witness']
    R,n,scale = source['core_residual'],len(source['high_primes']),witness['denominator']
    groups = (1<<n)-1
    hist = [[[0]*17 for h in range(groups)] for x in R]
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[row['point_index']][h][v] += row['weight']
    tails = [[[sum(line[t:]) for t in range(17)] for line in point] for point in hist]
    core = [d for d in range(3,211) if 210%d==0]
    reports,violations = [],[]
    for h in range(groups):
        shapes = core if h<n else [1,2]+core
        limit = len(shapes)
        model = cp_model.CpModel()
        selected,events,deficits = [],[],[]
        for x in range(len(R)):
            choose = model.NewBoolVar(f's{x}')
            row = [(t,model.NewBoolVar(f'e{x}_{t}')) for t in range(1,limit+1)
                   if tails[x][h][t]>0]
            model.Add(sum(v for t,v in row)==choose)
            selected.append(choose)
            events.append(row)
            deficits.extend((scale-tails[x][h][t])*v for t,v in row)
        capacities = []
        for d in shapes:
            cap = model.NewIntVar(0,len(R),f'c{d}')
            for a in {x%d for x in R}:
                model.Add(cap>=sum(selected[x] for x,v in enumerate(R) if v%d==a))
            capacities.append(cap)
        model.Add(sum(selected)>=2)
        model.Add(sum(t*v for row in events for t,v in row)>=sum(capacities)+1)
        model.Add(sum(deficits)<=scale-1)
        model.Minimize(sum(deficits))
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = seconds_per_axis
        solver.parameters.num_search_workers = 1
        status = solver.Solve(model)
        reports.append(dict(axis=h,status=solver.StatusName(status),seconds=solver.WallTime()))
        if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            points,thresholds = [],[]
            for x,row in enumerate(events):
                for t,v in row:
                    if solver.Value(v):
                        points.append(x)
                        thresholds.append(t)
            cap = sum(max(Counter(R[x]%d for x in points).values()) for d in shapes)
            used = sum(tails[x][h][t] for x,t in zip(points,thresholds))
            gap = used-(len(points)-1)*scale
            assert len(points)>=2 and sum(thresholds)>cap and gap>0
            violations.append(dict(axis=h,points=points,thresholds=thresholds,
                                   projection_sum_capacity=cap,violation_numerator=gap))
    return dict(violations=violations,axis_reports=reports,
                all_axes_proved_no_violation=all(r['status']=='INFEASIBLE' for r in reports),
                denominator=scale,seconds=time.monotonic()-start,
                scope='threshold-event separation only; solver infeasibility is not an independently audited proof')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=2)
    args = ap.parse_args()
    result = separate(json.loads(Path(args.source).read_text()),args.time)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(violations=len(result['violations']),seconds=result['seconds'],
                         statuses=dict(Counter(r['status'] for r in result['axis_reports'])),
                         largest=sorted(result['violations'],key=lambda r:r['violation_numerator'],reverse=True)[:3])))
