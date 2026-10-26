"""Mutually exclusive threshold-event families in fixed residue cases.

Pair incompatibility gives edges, not merely separate probability inequalities:
every clique of incompatible events has total probability at most one.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path


@lru_cache(maxsize=512)
def pair_capacities(points, remaining):
    return tuple(tuple(len(remaining) + sum(x % d == y % d for d in remaining)
                       for y in points) for x in points)


def maximum_clique(events, weights, capacities):
    assert len(events) == len(weights) and all(w >= 0 for w in weights)
    keep = sorted((j for j, w in enumerate(weights) if w > 0), key=lambda j: (-weights[j], j))
    vertices = [events[j] for j in keep]
    values = [weights[j] for j in keep]
    adj = [0] * len(keep)
    for i, (x, t) in enumerate(vertices):
        for j in range(i):
            y, u = vertices[j]
            if x != y and t + u > capacities[x][y]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    best, chosen = 0, []

    def visit(mask, score, clique):
        nonlocal best, chosen
        if score > best:
            best, chosen = score, clique
        total = sum(values[j] for j in range(len(values)) if mask >> j & 1)
        while mask and score + total > best:
            bit = mask & -mask
            j = bit.bit_length() - 1
            mask ^= bit
            total -= values[j]
            visit(mask & adj[j], score + values[j], clique + [j])

    visit((1 << len(values)) - 1, 0, [])
    return best, [keep[j] for j in chosen]


def strongest(points, remaining, distributions):
    base = int(1 in remaining)
    events, weights = [], []
    for x, row in enumerate(distributions):
        assert len(row) == len(remaining) - base + 1 and min(row) >= 0
        for j, value in enumerate(row):
            if j and value:
                events.append((x, base + j))
                weights.append(sum(row[j:]))
    total, indices = maximum_clique(events, weights, pair_capacities(tuple(points), tuple(remaining)))
    return total, [dict(point=events[j][0], threshold=events[j][1]) for j in indices]


def screen(data):
    points = data['input']['core_residual']
    n = len(data['input']['high_primes'])
    D = data['probability_denominator']
    cuts = []
    checked = 0
    for group in data['group_mixtures']:
        h = group['axis']
        shape = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        for j, state in enumerate(group['states']):
            if state['kind'] != 'fractional':
                continue
            remaining = [d for d in shape if d not in dict(state['fixed'])]
            assert all(sum(row) == D for row in state['count_numerators'])
            score, events = strongest(points, remaining, state['count_numerators'])
            checked += 1
            if score > D:
                cuts.append(dict(axis=h, state_index=j, events=events,
                                 probability_sum_numerator=score, violation_numerator=score-D))
    cuts.sort(key=lambda c: (-c['violation_numerator'], c['axis'], c['state_index']))
    return dict(fractional_states_checked=checked, violated_states=len(cuts),
                denominator=D, cuts=cuts, verified=False, covering_obstruction=False)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    report = screen(json.loads(Path(args.source).read_text()))
    report['source'] = args.source
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'cuts'} | dict(largest=report['cuts'][:1])))
