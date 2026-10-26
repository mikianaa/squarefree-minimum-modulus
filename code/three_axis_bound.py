"""Solver-free bound for squarefree minimum 3 on N = 2*p*q*r.

After prime and twice-prime anchors, each parity layer is a 3D box.
The odd pair moduli give at most one axis line in each direction in each
layer. Each even pair modulus adds one line in one of the two layers.
The odd triple modulus and its double cover at most three further points.
Ignore consistency of the odd lines between layers: this is a relaxation.
Enumerate the maximum union of axis lines exactly, using coordinate
permutations to fix the first direction and disjointness in the last.
"""
import itertools
import json
import math
from functools import lru_cache


@lru_cache(None)
def max_axis_union(dims, counts):
    x, y, z = dims
    a, b, c = counts
    # Counts are 1 or 2, capped when there are fewer distinct lines.
    a, b, c = min(a, y*z), min(b, x*z), min(c, x*y)
    first = [((0, 0),)] if a == 1 else []
    if a == 2:
        if z > 1:
            first.append(((0, 0), (0, 1)))
        if y > 1:
            first.append(((0, 0), (1, 0)))
        if y > 1 and z > 1:
            first.append(((0, 0), (1, 1)))
    def mask(points):
        return sum(1 << ((i*y+j)*z+k) for i,j,k in points)
    y_lines = [mask((i,j,k) for j in range(y)) for i in range(x) for k in range(z)]
    z_lines = [mask((i,j,k) for k in range(z)) for i in range(x) for j in range(y)]
    best = 0
    for yz in first:
        base = mask((i,j,k) for j,k in yz for i in range(x))
        for chosen in itertools.combinations(y_lines, b):
            covered = base
            for line in chosen:
                covered |= line
            # Parallel z-lines are disjoint, so their marginal gains add.
            gains = sorted(((line & ~covered).bit_count() for line in z_lines), reverse=True)
            best = max(best, covered.bit_count() + sum(gains[:c]))
    return best


def bound(primes=(3,5,7)):
    rows = []
    for tail in itertools.product((0,1), repeat=2):
        pattern = (0,) + tail
        dims = [tuple(p-1-(pattern[i] == e) for i,p in enumerate(primes)) for e in (0,1)]
        best = 0
        allocation = None
        for parity in itertools.product((0,1), repeat=3):
            total = sum(max_axis_union(dims[e], tuple(1+(b == e) for b in parity)) for e in (0,1))
            if total > best:
                best, allocation = total, parity
        size = sum(math.prod(d) for d in dims)
        rows.append(dict(pattern=pattern, boxes=dims, anchor_residual=size,
                         pair_union_upper=best, triple_capacity=3,
                         uncovered_lower=max(0,size-best-3), maximizing_allocation=allocation))
    return dict(primes=primes, N=2*math.prod(primes), uncovered_lower=min(r['uncovered_lower'] for r in rows), cases=rows)


def coupled_bound(primes=(3,5,7), canonical=True):
    """Keep the three odd pair lines identical across the two parity layers.

    Each coordinate appears in exactly two odd pair lines. Up to permutations
    fixing label 1 their label pair is (1,1), (1,2), (2,1), (2,2), or (2,3).
    The last option requires p >= 4. After fixing those lines, enumerate all
    choices of the three even lines. Triple classes add at most three points.
    """
    rows = []
    label_pairs = [([(1,1),(1,2),(2,1),(2,2)] + ([(2,3)] if p > 3 else [])) for p in primes]
    if not canonical:
        label_pairs = [list(itertools.product(range(1,p), repeat=2)) for p in primes]
    for tail in itertools.product((0,1), repeat=2):
        pattern = (0,)+tail
        points = [(e,a) for e in (0,1) for a in itertools.product(*(range(1,p) for p in primes))
                  if all(a[i] != 1 or pattern[i] != e for i in range(3))]
        masks = []
        even = []
        for direction in range(3):
            fixed = [i for i in range(3) if i != direction]
            lookup = {}
            for labels in itertools.product(*(range(1,primes[i]) for i in fixed)):
                pair = [sum(1 << j for j,(e,a) in enumerate(points)
                            if e == parity and all(a[i] == v for i,v in zip(fixed,labels)))
                        for parity in (0,1)]
                lookup[labels] = pair[0] | pair[1]
            masks.append(lookup)
            # Even lines may choose either parity and any labels.
            even.append(sorted(set(sum(1 << j for j,(e,a) in enumerate(points)
                                       if e == parity and all(a[i] == v for i,v in zip(fixed,labels)))
                                   for parity in (0,1)
                                   for labels in itertools.product(*(range(1,primes[i]) for i in fixed)))))
        best = 0
        residual_pairs = set()
        size0 = sum(e == 0 for e,a in points)
        mask0 = (1 << size0) - 1
        for pairs in itertools.product(*label_pairs):
            base = 0
            for direction in range(3):
                fixed = [i for i in range(3) if i != direction]
                labels = tuple(pairs[i][[j for j in range(3) if j != i].index(direction)] for i in fixed)
                base |= masks[direction][labels]
            candidates = [list(set(line & ~base for line in lines)) for lines in even]
            for first in candidates[0]:
                first_union = base | first
                for second in candidates[1]:
                    second_union = first_union | second
                    for third in candidates[2]:
                        union = second_union | third
                        covered = union.bit_count()
                        best = max(best, covered)
                        left0 = size0 - (union & mask0).bit_count()
                        left1 = len(points) - size0 - (covered - (union & mask0).bit_count())
                        # Optimistically let the odd triple remove one point
                        # in each layer and the even triple remove one more.
                        residual_pairs.add((max(0,left0-2),max(0,left1-1)))
                        residual_pairs.add((max(0,left0-1),max(0,left1-2)))
        frontier = []
        best_second = math.inf
        for pair in sorted(residual_pairs):
            if pair[1] < best_second:
                frontier.append(pair)
                best_second = pair[1]
        rows.append(dict(pattern=pattern, anchor_residual=len(points), pair_union_upper=best,
                         triple_capacity=3, uncovered_lower=len(points)-best-3,
                         residual_pair_frontier=frontier))
    return dict(N=2*math.prod(primes), method="coupled axis lines; exact integer enumeration",
                canonical=canonical, uncovered_lower=min(r['uncovered_lower'] for r in rows), cases=rows)


if __name__ == '__main__':
    print(json.dumps(coupled_bound(), indent=2))
