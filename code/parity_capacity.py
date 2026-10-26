"""Integer-only propagation of residual-count bounds in the two parity layers.

Each frontier pair is a componentwise lower bound: every realizable residual
dominates some listed pair. Bounds permit overlaps to be absent and are only
necessary conditions. In particular, passing a capacity test is not a cover.
"""
import itertools
import json
import math
from functools import lru_cache
from three_axis_bound import coupled_bound


def pareto(pairs):
    out = []
    smallest = math.inf
    for a,b in sorted(set(pairs)):
        if b < smallest:
            out.append((a,b))
            smallest = b
    return out


def capacities(primes, pattern):
    dims = [[p-1-(pattern[i] == e) for i,p in enumerate(primes)] for e in (0,1)]
    return [tuple(math.prod(dims[e][i] for i in range(len(primes)) if not mask >> i & 1)
                  for e in (0,1)) for mask in range(1,1 << len(primes))]


@lru_cache(maxsize=250000)
def minimum_deficit(lower, multiplier, caps):
    """min_{integer r >= lower} multiplier*r - sum min(r,c).

    The function is convex piecewise linear; its first nonnegative slope is
    reached at the (len(caps)-multiplier)-th smallest capacity, if positive.
    """
    excess = len(caps)-multiplier
    r = max(lower, caps[excess-1] if excess > 0 else lower)
    return max(0, multiplier*r-sum(min(r,c) for c in caps))


def extend(primes, pattern, frontier, q, new_parity):
    cap = capacities(primes,pattern)
    odd = [[c[e] for c in cap] for e in (0,1)]
    results = set()
    profiles = set()
    for bits in itertools.product((0,1),repeat=len(cap)):
        both = tuple(tuple(sorted(odd[e] + [c[e] for b,c in zip(bits,cap) if b == e])) for e in (0,1))
        profiles.add(both)
    for both in profiles:
        for lower in frontier:
            results.add(tuple(minimum_deficit(lower[e], q-1-(new_parity == e), both[e]) for e in (0,1)))
    return pareto(results)


def can_cover_counts(primes, pattern, lower, q, new_parity):
    """Can the capacity relaxation cover both layers from this residual pair?

    For a larger residual, each class covers a smaller fraction of that layer.
    Thus failure at every minimal frontier point excludes every larger pair.
    A knapsack DP allocates even shapes; odd shapes are granted capacity in both.
    """
    cap = [tuple(min(lower[e],c[e]) for e in (0,1)) for c in capacities(primes,pattern)]
    need = [max(0,(q-1-(new_parity == e))*lower[e]-sum(c[e] for c in cap)) for e in (0,1)]
    # Starting with all even shapes in layer 1, transfer a subset to layer 0.
    budget = sum(c[1] for c in cap)-need[1]
    if budget < 0:
        return False
    dp = {0:0}
    for gain,cost in cap:
        updated = dict(dp)
        for used,value in dp.items():
            if used+cost <= budget:
                updated[used+cost] = max(updated.get(used+cost,-1),value+gain)
        dp = updated
    return max(dp.values()) >= need[0]


def run():
    out = []
    for row in coupled_bound()['cases']:
        for b11 in (0,1):
            pattern = tuple(row['pattern'])+(b11,)
            frontier = extend((3,5,7), row['pattern'], row['residual_pair_frontier'], 11,b11)
            outcomes = []
            for b13 in (0,1):
                survivors = [r for r in frontier if can_cover_counts((3,5,7,11),pattern,r,13,b13)]
                outcomes.append(dict(b13=b13,excluded=not survivors,surviving_lower_pairs=survivors))
            result = dict(pattern=pattern, residual_lower=min(sum(r) for r in frontier),
                          frontier=frontier, extension_13=outcomes)
            out.append(result)
            print(json.dumps(result),flush=True)
    return out


if __name__ == '__main__':
    run()
