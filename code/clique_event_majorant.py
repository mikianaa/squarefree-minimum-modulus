"""Rational event upper bounds using mutually exclusive event cliques."""
from fractions import Fraction
from fractional_group_scenarios import maximize
from pair_event_majorant import upper_bound as pair_bound


def pair(v):
    return [v.numerator, v.denominator]


def upper_bound(points, shapes, events, weights):
    assert len(events) == len(weights) and all(type(w) is int and w >= 0 for w in weights)
    N = len(shapes)
    base = int(1 in shapes)
    assert all(type(x) is int and 0 <= x < len(points) and type(t) is int and 1 <= t <= N for x, t in events)
    rewards = [[sum(w for (x, t), w in zip(events, weights) if x == j and t <= k)
                for k in range(N+1)] for j in range(len(points))]
    primal = maximize(points, shapes, rewards, pair_compatibility=True,
                      return_dual=True, clique_compatibility=True)
    dual = primal['numerical_dual']
    quantum = 10**6
    slopes = [Fraction(max(0, round(v*quantum)), quantum) for v in dual['slopes']]
    cuts = []
    adjustments = [[Fraction(0) for k in range(N+1)] for p in points]
    source_cuts = dual['clique_cuts'] + [dict(events=[dict(point=x, threshold=t)
                      for x, t in zip(r['points'], r['thresholds'])], multiplier=r['multiplier']) for r in dual['pair_cuts']]
    for row in source_cuts:
        mu = Fraction(max(0, round(row['multiplier']*quantum)), quantum)
        if not mu:
            continue
        cuts.append(dict(events=row['events'], multiplier=pair(mu)))
        for e in row['events']:
            for k in range(e['threshold'], N+1):
                adjustments[e['point']][k] += mu
    offsets = [max(Fraction(rewards[x][k])-(k-base)*slopes[x]-adjustments[x][k]
                   for k in range(base, N+1)) for x in range(len(points))]
    capacities = [max(sum(a for p, a in zip(points, slopes) if p % d == r)
                      for r in {p % d for p in points}) for d in shapes if d != 1]
    total = sum(offsets)+sum(capacities)+sum(Fraction(*r['multiplier']) for r in cuts)
    bound = min(sum(row[-1] for row in rewards), total.numerator//total.denominator)
    old = pair_bound(points, shapes, events, weights)
    if bound >= old['bound']:
        return old
    return dict(method='clique_affine', base=base, bound=bound, slopes=list(map(pair, slopes)),
                point_intercepts=list(map(pair, offsets)), intercept_fraction=pair(sum(offsets)),
                resource_bounds=list(map(pair, capacities)), clique_cuts=cuts)
