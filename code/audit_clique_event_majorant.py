"""Independent rational affine-majorant audit with ordinary-residue clique checks."""
from fractions import Fraction
import itertools


def number(v):
    assert len(v) == 2 and all(type(k) is int for k in v) and v[1] > 0
    return Fraction(*v)


def audit_bound(points, shapes, events, certificate):
    assert points and len(set(points)) == len(points)
    assert len(set(shapes)) == len(shapes) and all(type(d) is int and d > 0 for d in shapes)
    assert certificate['method'] == 'clique_affine'
    N = len(shapes)
    base = int(1 in shapes)
    assert certificate['base'] == base
    rewards = [[0]*(N+1) for p in points]
    for e in events:
        x, t, w = e['point'], e['threshold'], e['weight']
        assert type(x) is int and 0 <= x < len(points) and type(t) is int and 1 <= t <= N
        assert type(w) is int and w >= 0
        for k in range(t, N+1):
            rewards[x][k] += w
    slopes = list(map(number, certificate['slopes']))
    offsets = list(map(number, certificate['point_intercepts']))
    assert len(slopes) == len(offsets) == len(points) and min(slopes) >= 0
    corrections = [[Fraction(0) for k in range(N+1)] for p in points]
    cut_total = Fraction(0)
    for cut in certificate['clique_cuts']:
        clique = [(e['point'], e['threshold']) for e in cut['events']]
        assert clique and len(set(clique)) == len(clique)
        for x, t in clique:
            assert type(x) is int and 0 <= x < len(points) and type(t) is int and 1 <= t <= N
        for (x, t), (y, u) in itertools.combinations(clique, 2):
            cap = sum(max(int(points[x] % d == a)+int(points[y] % d == a) for a in range(d)) for d in shapes)
            assert t+u > cap
        mu = number(cut['multiplier'])
        assert mu >= 0
        cut_total += mu
        for x, t in clique:
            for k in range(t, N+1):
                corrections[x][k] += mu
    for x in range(len(points)):
        for k in range(base, N+1):
            assert rewards[x][k] <= offsets[x]+slopes[x]*(k-base)+corrections[x][k]
    assert sum(offsets) == number(certificate['intercept_fraction'])
    capacities = [max(sum(slopes[x] for x, p in enumerate(points) if p % d == a)
                      for a in range(d)) for d in shapes if d != 1]
    assert capacities == list(map(number, certificate['resource_bounds']))
    total = sum(offsets)+sum(capacities)+cut_total
    bound = min(sum(row[-1] for row in rewards), total.numerator//total.denominator)
    assert type(certificate['bound']) is int and certificate['bound'] == bound
    return bound
