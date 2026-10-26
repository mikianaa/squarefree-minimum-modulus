"""Exact 210-core geometry up to permutations of non-anchor CRT labels.

Prime and twice-prime classes are anchored. Enumerate the three odd pair lines,
then their even doubles. Unused labels are introduced in first-use order, a
lossless quotient by permutations fixing labels 0 and 1. The resulting masks
are residuals BEFORE the final classes 105 and 210.
"""
import argparse
import itertools
import json
import time
from pathlib import Path
from fractions import Fraction

PRIMES = (3,5,7)


def points_for(pattern):
    return [(e,a) for e in (0,1) for a in itertools.product(*(range(1,p) for p in PRIMES))
            if all(a[i] != 1 or pattern[i] != e for i in range(3))]


def pair_residuals(pattern, canonical_even=True):
    points = points_for(pattern)
    whole = (1 << len(points))-1
    fixed = [[i for i in range(3) if i != direction] for direction in range(3)]
    masks = {}
    for direction in range(3):
        for labels in itertools.product(*(range(1,PRIMES[i]) for i in fixed[direction])):
            for e in (0,1):
                masks[direction,e,labels] = sum(1 << j for j,(parity,a) in enumerate(points)
                                                if parity == e and all(a[i]==v for i,v in zip(fixed[direction],labels)))
            masks[direction,None,labels] = masks[direction,0,labels] | masks[direction,1,labels]
    pair_choices = [[(1,1),(1,2),(2,1),(2,2)]+([(2,3)] if p>3 else []) for p in PRIMES]
    residuals = set()
    configurations = 0
    for pairs in itertools.product(*pair_choices):
        covered = 0
        for direction in range(3):
            labels = tuple(pairs[i][[d for d in range(3) if d != i].index(direction)] for i in fixed[direction])
            covered |= masks[direction,None,labels]
        used = tuple(max(pair) for pair in pairs)

        def visit(index,union,maximum):
            nonlocal configurations
            if index == 3:
                residuals.add(whole ^ union)
                configurations += 1
                return
            direction = 2-index  # Moduli 30, 42, 70, in that order.
            coords = fixed[direction]
            for parity in (0,1):
                alphabets = [tuple(v for v in range(1,min(PRIMES[i]-1,maximum[i]+1)+1)
                                   if v != 1 or pattern[i] != parity)
                             if canonical_even else
                             tuple(v for v in range(1,PRIMES[i]) if v != 1 or pattern[i] != parity)
                             for i in coords]
                for labels in itertools.product(*alphabets):
                    updated = list(maximum)
                    for i,v in zip(coords,labels):
                        updated[i] = max(updated[i],v)
                    visit(index+1,union | masks[direction,parity,labels],tuple(updated))
        visit(0,covered,used)
    return points,residuals,configurations


def triple_groups(points):
    groups = {}
    for j,(e,a) in enumerate(points):
        groups[a] = groups.get(a,0) | (1 << j)
    return list(groups.values())


def minimal_triple_residuals(mask,groups):
    """All residuals after nonredundant choices of 105 and then 210.

    A 105 class is a group with one or two points of the anchor box. Empty
    intersections can be moved to a nonempty one without losing coverage.
    A 210 class removes one further point; it can similarly be made nonredundant.
    Every arbitrary pair of triple placements therefore has a result containing
    one of the masks yielded here. Different nonempty groups must all be kept.
    """
    for group in groups:
        intersection = mask & group
        if not intersection:
            continue
        left = mask ^ intersection
        if not left:
            yield 0
            continue
        bits = left
        while bits:
            bit = bits & -bits
            yield left ^ bit
            bits ^= bit


def profile(mask,points):
    counts = [0,0]
    histograms = {(d,e):[0]*d for d in (3,5,7,15,21,35,105) for e in (0,1)}
    for j,(e,a) in enumerate(points):
        if not mask >> j & 1:
            continue
        counts[e] += 1
        # CRT reconstruction in the odd period, using direct small arithmetic.
        n = (70*a[0]+21*a[1]+15*a[2]) % 105
        for d in (3,5,7,15,21,35,105):
            histograms[d,e][n%d] += 1
    caps = [max(histograms[d,e]) for d in (3,5,7,15,21,35,105) for e in (0,1)]
    return tuple(counts+caps)


def optimistic_profiles(points,masks):
    """Keep actual pair-residual intersection capacities, relax the last 3 points.

    The 105 class removes at most one point per parity and 210 removes one more
    in a single parity. Use these optimistic lower counts with the original
    pair-residual capacities, truncated to the lower count. This is valid for
    cover feasibility by the monotonicity of min(r,C)/r, not a realizability claim.
    """
    profiles = set()
    for mask in masks:
        p = profile(mask,points)
        for extra in (0,1):
            r = [max(0,p[e]-1-(extra==e)) for e in (0,1)]
            caps = [min(r[e],p[2+2*j+e]) for j in range(7) for e in (0,1)]
            profiles.add(tuple(r+caps))
    return profiles


def dominates_for_cover(easier,harder):
    for e in (0,1):
        if easier[e] == 0:
            continue
        if harder[e] == 0:
            return False
        if any(easier[2+2*j+e]*harder[e] < harder[2+2*j+e]*easier[e] for j in range(7)):
            return False
    return True


def easiest_profiles(profiles):
    """Keep a dominating set for cover feasibility, using exact rational order.

    A profile with greater c_d(e)/r_e for every d,e is easier to cover. Therefore
    infeasibility for it excludes any dominated profile. A zero-count parity
    imposes no constraint. The retained profiles are still only relaxations.
    """
    def score(p):
        return sum((sum(Fraction(p[2+2*j+e],p[e]) for j in range(7)) if p[e] else Fraction(8)) for e in (0,1))
    kept = []
    for p in sorted(profiles,key=lambda p:(score(p),p),reverse=True):
        if not any(dominates_for_cover(q,p) for q in kept):
            kept.append(p)
    return kept


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='011')
    ap.add_argument('--full-even',action='store_true')
    ap.add_argument('--output',required=True)
    ap.add_argument('--profiles',help='also save optimistic capacity profiles')
    a = ap.parse_args()
    start = time.monotonic()
    pattern = tuple(map(int,a.pattern))
    assert len(pattern)==3 and pattern[0]==0 and all(v in (0,1) for v in pattern)
    points,masks,total = pair_residuals(pattern,not a.full_even)
    result = dict(pattern=pattern,canonical_even=not a.full_even,points=points,
                  configurations=total,distinct_pair_residuals=len(masks),
                  pair_residual_masks=sorted(masks),minimum_pair_residual=min(m.bit_count() for m in masks),
                  seconds=round(time.monotonic()-start,3))
    Path(a.output).write_text(json.dumps(result)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('points','pair_residual_masks')}),flush=True)
    if a.profiles:
        profiles = optimistic_profiles(points,masks)
        easiest = easiest_profiles(profiles)
        result = dict(pattern=pattern,canonical_even=not a.full_even,
                      distinct_profiles=len(profiles),profiles=sorted(profiles),
                      easiest_profiles=easiest,easiest_count=len(easiest),
                      seconds=round(time.monotonic()-start,3),
                      scope='optimistic necessary-capacity profiles; not realizable families')
        Path(a.profiles).write_text(json.dumps(result)+'\n')
        print(json.dumps({k:v for k,v in result.items() if k not in ('profiles','easiest_profiles')}),flush=True)
