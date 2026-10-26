"""Search for whole-assignment event cuts against displayed fractional cases.

CP-SAT bounds are search suggestions; independently certify any positive cut.
All thresholds here count only resources not fixed in the case.
"""
import argparse
import json
from pathlib import Path
import time
from ortools.linear_solver import pywraplp
from group_event_coherence import maximize


def case_events(points, shapes, state, denominator):
    fixed = dict(state['fixed'])
    assert len(fixed) == len(state['fixed'])
    assert all(d in shapes and 0 <= a < d for d, a in fixed.items())
    remaining = [d for d in shapes if d not in fixed]
    base = int(1 in remaining)
    events, targets = [], []
    for x, row in enumerate(state['count_numerators']):
        assert len(row) == len(remaining) - base + 1 and sum(row) == denominator
        for j, v in enumerate(row):
            if j and v:
                events.append((x, base + j))
                targets.append(sum(row[j:]))
    return remaining, events, targets


def separate(points, shapes, state, denominator, seconds=10):
    start = time.monotonic()
    remaining, events, targets = case_events(points, shapes, state, denominator)
    if not events:
        return dict(status='NO_NONCONSTANT_EVENTS', cut=None, assignments=[])
    lp = pywraplp.Solver.CreateSolver('GLOP')
    lp.SetTimeLimit(max(1, int(seconds * 1000)))
    weights = [lp.NumVar(0, 1, f'w{j}') for j in range(len(events))]
    beta = lp.NumVar(0, 1, 'beta')
    lp.Add(sum(weights) <= 1)
    objective = lp.Objective()
    objective.SetCoefficient(beta, -1)
    for v, t in zip(weights, targets):
        objective.SetCoefficient(v, t / denominator)
    objective.SetMaximization()
    profiles = set()
    assignments = []
    iterations = 0
    cut = None
    last_gap = None
    status = 'TIME_LIMIT'
    while time.monotonic() - start < seconds:
        if lp.Solve() != pywraplp.Solver.OPTIMAL:
            status = 'MASTER_UNRESOLVED'
            break
        iterations += 1
        last_gap = objective.Value()
        if last_gap < 1e-9:
            status = 'NO_NUMERICAL_SEPARATION'
            break
        # Small integer weights first: these can be independently certified
        # without a very large cardinality encoding.
        suggestions = [max(0, w.solution_value()) for w in weights]
        for scale in (16, 64, 256, 1000000):
            if time.monotonic() - start >= seconds:
                break
            integers = [max(0, round(w * scale)) for w in suggestions]
            if not any(integers):
                continue
            priced = maximize(points, remaining, events, integers,
                              min(1, max(.05, seconds - time.monotonic() + start)))
            if 'upper_bound' not in priced:
                continue
            gap = sum(w * t for w, t in zip(integers, targets)) - denominator * priced['upper_bound']
            if gap > 0:
                cut = dict(events=[dict(point=x, threshold=t, weight=w)
                                   for (x, t), w in zip(events, integers) if w],
                           bound=priced['upper_bound'], violation_numerator=gap,
                           denominator=denominator, weight_scale=scale, pricing=priced)
                status = 'CANDIDATE_CUT'
                break
            profile = tuple(int(priced['counts'][x] >= t) for x, t in events)
            if profile not in profiles:
                profiles.add(profile)
                assignments.append(priced['residues'])
                lp.Add(sum(w for w, hit in zip(weights, profile) if hit) <= beta)
        if cut:
            break
    return dict(status=status, cut=cut, remaining=remaining, events=events,
                targets=targets, iterations=iterations, numerical_master_gap=last_gap,
                assignments=assignments, seconds=time.monotonic()-start, verified=False)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--count', type=int, default=12)
    ap.add_argument('--time', type=float, default=10)
    ap.add_argument('--offset', type=int, default=0)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    points = data['input']['core_residual']
    n = len(data['input']['high_primes'])
    states = sorted([(s['weight']/g['denominator'], g['axis'], j, s)
                     for g in data['group_mixtures'] for j, s in enumerate(g['states'])
                     if s['kind'] == 'fractional'], key=lambda row: (-row[0], row[1], row[2]))
    report = dict(source=args.source, denominator=data['probability_denominator'],
                  results=[], verified=False, covering_obstruction=False,
                  scope='Whole-assignment cut search against selected fractional case states; every positive cut still needs independent proof.')
    for weight, h, j, s in states[args.offset:args.offset + args.count]:
        shape = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        row = separate(points, shape, s, data['probability_denominator'], args.time)
        row.update(axis=h, state_index=j, fixed=s['fixed'], mixture_weight=weight)
        report['results'].append(row)
        Path(args.output).write_text(json.dumps(report, separators=(',', ':')) + '\n')
        print(json.dumps({k: v for k, v in row.items() if k not in ('assignments', 'events', 'targets', 'cut')}
                         | dict(cut=row['cut'])), flush=True)
