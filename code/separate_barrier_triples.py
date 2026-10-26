"""Find conditional three-point capacity cuts violated by exact barrier states.

For thresholds whose sum exceeds the conditional maximum count sum, at
most two events can hold. Rejecting a displayed state is not an obstruction.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import numpy as np


@lru_cache(maxsize=None)
def triple_frontier(N, common, separate, pairs):
    """Largest attainable third count for each pair of required counts."""
    ts = np.arange(N + 1)[:, None]
    us = np.arange(N + 1)[None, :]
    frontier = np.full((N + 1, N + 1), -1, dtype=np.int64)
    for choices in itertools.product(*(range(v + 1) for v in pairs)):
        hits = [common + choices[j] + sum(pairs[i] - choices[i] for i in range(3) if i != j) for j in range(3)]
        debt = np.maximum(0, ts - hits[0]) + np.maximum(0, us - hits[1])
        frontier = np.maximum(frontier, np.where(debt <= separate, hits[2] + separate - debt, -1))
    return frontier


def strongest(points, shapes, fixed, counts, denominator, exact=False):
    remaining = [d for d in shapes if d not in dict(fixed)]
    base = int(1 in remaining)
    hits = np.array([sum(p % d == a for d, a in fixed) for p in points])
    N = len(shapes)
    hist = np.zeros((len(points), N + 1), dtype=np.int64)
    for x, row in enumerate(counts):
        assert sum(row) == denominator and len(row) == len(remaining) - base + 1
        hist[x, hits[x] + base:hits[x] + len(remaining) + 1] = row
    tails = np.cumsum(hist[:, ::-1], axis=1)[:, ::-1]
    triples = np.array(list(itertools.combinations(range(len(points)), 3)))
    if not len(triples):
        return None
    capacity = hits[triples].sum(axis=1)
    common = np.zeros(len(triples), dtype=np.int64)
    separate = np.zeros(len(triples), dtype=np.int64)
    pairs = np.zeros((len(triples), 3), dtype=np.int64)
    for d in remaining:
        residues = np.array(points)[triples] % d
        a, b, c = residues.T
        capacity += 1 + ((a == b) | (a == c) | (b == c)) + ((a == b) & (b == c))
        if exact:
            common += (a == b) & (b == c)
            separate += (a != b) & (a != c) & (b != c)
            pairs[:, 0] += (b == c) & (a != b)
            pairs[:, 1] += (a == c) & (b != a)
            pairs[:, 2] += (a == b) & (c != a)
    ts = np.arange(1, N + 1)[None, :, None]
    us = np.arange(1, N + 1)[None, None, :]
    vs = np.maximum(1, capacity[:, None, None] + 1 - ts - us)
    if exact:
        for j, (x, y, z) in enumerate(triples):
            frontier = triple_frontier(N, int(common[j]), int(separate[j]), tuple(map(int, pairs[j])))
            a = np.maximum(0, ts[0] - hits[x])
            b = np.maximum(0, us[0] - hits[y])
            vs[j] = np.maximum(1, frontier[a, b] + hits[z] + 1)
    valid = (vs >= 1) & (vs <= N)
    x, y, z = triples.T
    scores = (tails[x[:, None, None], ts] + tails[y[:, None, None], us]
              + tails[z[:, None, None], np.clip(vs, 0, N)])
    scores = np.where(valid, scores, -1)
    j, t, u = np.unravel_index(scores.argmax(), scores.shape)
    excess = int(scores[j, t, u]) - 2 * denominator
    if excess <= 0:
        return None
    return dict(point_indices=list(map(int, triples[j])),
                thresholds=[int(t + 1), int(u + 1), int(vs[j, t, u])],
                capacity=int(capacity[j]), violation_numerator=excess,
                sum_capacity_detects=int(t + u + 2 + vs[j, t, u]) > int(capacity[j]))


def separate(data, exact=False):
    points = data['input']['core_residual']
    n = len(data['input']['high_primes'])
    D = data['probability_denominator']
    cuts = []
    checked = 0
    for group in data['group_mixtures']:
        h = group['axis']
        shapes = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        for j, state in enumerate(group['states']):
            if state['kind'] != 'fractional':
                continue
            cut = strongest(points, shapes, state['fixed'], state['count_numerators'], D, exact)
            checked += 1
            if cut:
                cuts.append(dict(axis=h, state_index=j, **cut))
    cuts.sort(key=lambda c: c['violation_numerator'], reverse=True)
    return dict(denominator=D, fractional_states_checked=checked, cuts=cuts,
                verified=False, covering_obstruction=False, exact_triples=exact,
                scope='At most one violated conditional triple event cut per displayed fractional state. No claim of compatibility of four or more points.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--exact', action='store_true')
    args = ap.parse_args()
    result = separate(json.loads(Path(args.source).read_text()), args.exact)
    result['source'] = args.source
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'cuts'} | dict(violated_states=len(result['cuts']), largest=result['cuts'][:1])))
