"""Independent ordinary-residue check of conditional triple capacity cuts."""
import argparse
import json
from pathlib import Path


def audit(barrier, report):
    points = barrier['input']['core_residual']
    n = len(barrier['input']['high_primes'])
    D = barrier['probability_denominator']
    assert D == report['denominator']
    groups = {g['axis']: g for g in barrier['group_mixtures']}
    seen = set()
    largest = 0
    for cut in report['cuts']:
        h, j = cut['axis'], cut['state_index']
        assert (h, j) not in seen
        seen.add((h, j))
        state = groups[h]['states'][j]
        assert state['kind'] == 'fractional'
        shape = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        fixed = dict(state['fixed'])
        assert len(fixed) == len(state['fixed']) and all(d in shape and 0 <= a < d for d, a in fixed.items())
        indices, ts = cut['point_indices'], cut['thresholds']
        assert len(indices) == len(ts) == 3 and len(set(indices)) == 3
        assert all(type(x) is int and 0 <= x < len(points) for x in indices)
        assert all(type(t) is int and 1 <= t <= len(shape) for t in ts)
        chosen = [points[x] for x in indices]
        capacity = sum(sum(p % d == fixed[d] for p in chosen) if d in fixed else
                       max(sum(p % d == a for p in chosen) for a in range(d)) for d in shape)
        assert capacity == cut['capacity']
        # Independent capped-count dynamic program over ordinary residues.
        # Saturating at the thresholds preserves reachability of all events.
        reachable = {(0, 0, 0)}
        for d in shape:
            choices = [fixed[d]] if d in fixed else range(d)
            increments = {tuple(int(p % d == a) for p in chosen) for a in choices}
            reachable = {tuple(min(t, v + inc) for t, v, inc in zip(ts, old, delta))
                         for old in reachable for delta in increments}
        assert tuple(ts) not in reachable
        score = 0
        for x, t in zip(indices, ts):
            row = state['count_numerators'][x]
            assert all(type(v) is int and v >= 0 for v in row) and sum(row) == D
            offset = sum(points[x] % d == a for d, a in fixed.items()) + int(1 in shape and 1 not in fixed)
            score += sum(v for k, v in enumerate(row) if offset + k >= t)
        gap = score - 2 * D
        assert gap == cut['violation_numerator'] and gap > 0
        largest = max(largest, gap)
    return dict(verified=True, conditional_cuts_checked=len(seen),
                largest_violation_numerator=largest, denominator=D,
                displayed_barrier_decomposition_excluded=bool(seen),
                covering_obstruction=False, global_problem_settled=False,
                scope='Every checked cut holds for every actual assignment in its fixed-residue case and excludes the cited fractional state. Other barrier decompositions remain possible.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ('source', 'cuts', 'output'):
        ap.add_argument('--' + key, required=True)
    args = ap.parse_args()
    result = audit(json.loads(Path(args.source).read_text()), json.loads(Path(args.cuts).read_text()))
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
