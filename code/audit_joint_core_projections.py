"""Independently check simultaneous projection counts using integer residues."""
import argparse
import itertools
import json
import math
from pathlib import Path


def audit(data):
    source = data['input']
    points, primes, pattern = [source[k] for k in ('core_residual', 'high_primes', 'high_pattern')]
    n = len(primes)
    assert n == len(pattern) and all(type(b) is int and b in (0, 1) for b in pattern)
    assert len(set(points)) == len(points) and all(type(x) is int and 0 <= x < 210 for x in points)
    assert primes == sorted(set(primes)) and all(type(q) is int and q > 7 and
        all(q % d for d in range(2, math.isqrt(q)+1)) for q in primes)
    groups = [s for k in range(1, n+1) for s in itertools.combinations(range(n), k)]
    assignments = data['group_assignments']
    assert len(assignments) == len(groups)
    counts = [[] for x in points]
    for h, residues in enumerate(assignments):
        shapes = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        assert len(residues) == len(shapes) and all(type(a) is int and 0 <= a < d for a, d in zip(residues, shapes))
        for x, row in zip(points, counts):
            row.append(sum(x % d == a for a, d in zip(residues, shapes)))
    values = []
    for x, row in zip(points, counts):
        free = [max(0, q-1-int(b == x % 2)-row[i]) for i, (q, b) in enumerate(zip(primes, pattern))]
        values.append(math.prod(free)-sum(row[h]*math.prod(free[i] for i in range(n) if i not in group)
                                         for h, group in enumerate(groups[n:], n)))
    assert values == data['local_F']
    positive = sum(max(v, 0) for v in values)
    assert positive == data['positive_volume']
    feasible = positive == 0
    return dict(verified=True, points_checked=len(points), groups_checked=len(groups), local_F=values,
                positive_volume=positive, violated_points=sum(v > 0 for v in values),
                simultaneous_projection_relaxation_feasible=feasible,
                covering_constructed=False, global_problem_settled=False,
                scope='Explicit group core residues and exact local F only. A positive objective is not a lower bound on other assignments.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    result = audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
