"""Sharp size-only bound for the sum of maximal cylinder intersections.

For a product grid with side lengths p_i and |R|=r, the maximum of
sum_{nonempty S} max_a |R intersect cylinder(S,a)| is W(r)-r, where W(r)
is the sum of the r largest weights 2^(number of zero coordinates).
Coordinate compression proves the upper bound; a zero-count/lex initial
segment attains it. This optimizes over arbitrary sets, not cover residuals.
"""
import json
import math
from pathlib import Path


def zero_counts(dimensions):
    counts = [1]
    for p in dimensions:
        assert isinstance(p,int) and p>=1
        following = [0]*(len(counts)+1)
        for z,c in enumerate(counts):
            following[z] += (p-1)*c
            following[z+1] += c
        counts = following
    assert sum(counts)==math.prod(dimensions)
    return counts


def top_sum(counts,r):
    assert 0<=r<=sum(counts)
    total = 0
    remaining = r
    for z in range(len(counts)-1,-1,-1):
        take = min(remaining,counts[z])
        total += take*(1<<z)
        remaining -= take
    assert remaining==0
    return total


def maximum_capacity(dimensions,r):
    return top_sum(zero_counts(dimensions),r)-r


def extension_bound(dimensions,lower,q):
    """Minimum of q*r-W(r) over all integer r>=lower, truncated at zero."""
    counts = zero_counts(dimensions)
    corner = sum(c for z,c in enumerate(counts) if (1<<z)>q)
    minimizing_size = max(lower,corner)
    return dict(uncovered_lower=max(0,q*minimizing_size-top_sum(counts,minimizing_size)),
                minimizing_residual_size=minimizing_size,input_lower=lower,new_prime=q)


def exclusion_threshold(dimensions,q):
    """First positive size r where q*r-W(r)>0; larger sizes also exclude."""
    counts = zero_counts(dimensions)
    size = mass = 0
    for z in range(len(counts)-1,-1,-1):
        weight = 1<<z
        end = size+counts[z]
        end_mass = mass+weight*counts[z]
        if q*end>end_mass:
            assert q>weight
            return max(size,(mass-weight*size)//(q-weight))+1
        size,mass = end,end_mass
    return size+1


def main():
    examples = []
    for primes,lower,q,kind in [([2,3,5,7],22,11,'integer_lower'),
                           ([2,3,5,7,11],90,13,'integer_lower'),
                           ([2,3,5,7,11],106,13,'integer_lower'),
                           ([2,3,5,7,11],121,13,'solver_certified_lower'),
                           ([2,3,5,7,11,13],792,17,'witness_upper'),
                           ([2,3,5,7,11,13,17],8000,19,'hypothetical'),
                           ([2,3,5,7,11,13,17,19],79979,23,'witness_upper')]:
        result = dict(primes=primes,period=math.prod(primes),
                      sharp_capacity=maximum_capacity(primes,lower),
                      exclusion_threshold=exclusion_threshold(primes,q),
                      **extension_bound(primes,lower,q),
                      input_kind=kind,
                      scope=('valid propagated lower bound' if kind.endswith('_lower') else
                             'conditional only; this input is not an established residual lower bound'))
        examples.append(result)
        print(json.dumps(result),flush=True)
    Path('data/extremal-slice-capacity-examples.json').write_text(json.dumps(examples,indent=2)+'\n')


if __name__=='__main__':
    main()
