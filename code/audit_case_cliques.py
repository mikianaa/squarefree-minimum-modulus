"""Independent ordinary-residue and integer audit of conditional clique cuts."""
import argparse
import itertools
import json
from pathlib import Path


def audit(barrier, report):
    source = barrier['input']
    points, n = source['core_residual'], len(source['high_primes'])
    D = barrier['probability_denominator']
    assert type(D) is int and D > 0 and report['denominator'] == D
    groups = {g['axis']: g for g in barrier['group_mixtures']}
    seen = set()
    pairs = 0
    worst = 0
    for cut in report['cuts']:
        h, j = cut['axis'], cut['state_index']
        assert (h, j) not in seen
        seen.add((h, j))
        state = groups[h]['states'][j]
        assert state['kind'] == 'fractional'
        shape = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        fixed = dict(state['fixed'])
        assert len(fixed) == len(state['fixed']) and all(d in shape and 0 <= a < d for d, a in fixed.items())
        remaining = [d for d in shape if d not in fixed]
        base = int(1 in remaining)
        events = [(e['point'], e['threshold']) for e in cut['events']]
        assert events and len(set(events)) == len(events)
        for x, t in events:
            assert type(x) is int and 0 <= x < len(points)
            assert type(t) is int and base < t <= len(remaining)
        for (x, t), (y, u) in itertools.combinations(events, 2):
            cap = sum(max(int(points[x] % d == a) + int(points[y] % d == a)
                          for a in range(d)) for d in remaining)
            assert t + u > cap
            pairs += 1
        score = 0
        for x, t in events:
            row = state['count_numerators'][x]
            assert len(row) == len(remaining) - base + 1
            assert all(type(v) is int and v >= 0 for v in row) and sum(row) == D
            score += sum(v for k, v in enumerate(row) if base + k >= t)
        assert score == cut['probability_sum_numerator'] and score-D == cut['violation_numerator'] > 0
        worst = max(worst, score-D)
    assert len(seen) == report['violated_states']
    return dict(verified=True, conditional_clique_cuts_checked=len(seen),
                ordinary_residue_pair_checks=pairs, largest_violation_numerator=worst,
                denominator=D, displayed_barrier_decomposition_excluded=bool(seen),
                covering_obstruction=False, global_problem_settled=False,
                scope='Within each stated residue case the events are pairwise disjoint, so their probabilities sum to at most one. This excludes the displayed fractional states, not every possible barrier mixture or a covering system.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ('source', 'cuts', 'output'):
        ap.add_argument('--' + key, required=True)
    args = ap.parse_args()
    report = audit(json.loads(Path(args.source).read_text()), json.loads(Path(args.cuts).read_text()))
    Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
