"""Exact limitation of one quadratic, factored-moment tail certificate.

This does NOT exclude first moments, stop-loss bounds, exact future anchors,
or improved estimates retaining compatibility of class intersections.
"""
from fractions import Fraction
import json
import math
from pathlib import Path


def primes_to(limit):
    sieve = bytearray(b'\1')*(limit+1)
    sieve[:2] = b'\0\0'
    for p in range(2,math.isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes(len(sieve[p*p::p]))
    return [p for p in range(2,limit+1) if sieve[p]]


def tail_lower(primes,scale=10**40):
    """Lower bound on the optimal finite quadratic loss per initial moment.

    Bellman: F_p(K) = K + (1+sqrt(1+12*p*K))^2/(4*p^2).
    Integer square roots and divisions round down; F_p is increasing in K.
    """
    value = 0
    for p in reversed(primes):
        root = math.isqrt(scale*scale+12*p*value*scale)
        value += (scale+root)**2//(4*p*p*scale)
    return Fraction(value,scale)


def certificate(prefix,cutoff):
    primes = primes_to(cutoff)
    minimum_moment = math.prod(Fraction(p+3,p) for p in primes if p<=prefix)
    lower = tail_lower([p for p in primes if p>prefix])
    product = minimum_moment*lower
    assert product>1
    pair = lambda x:[x.numerator,x.denominator]
    return dict(prefix_last_prime=prefix,finite_tail_last_prime=cutoff,
                minimum_initial_moment=pair(minimum_moment),tail_factor_lower=pair(lower),
                loss_bound_lower=pair(product),loss_bound_lower_decimal=float(product),
                scope='barrier for the factored quadratic catalogue estimate only')


if __name__=='__main__':
    result = [certificate(p,q) for p,q in ((7,109),(11,241),(13,709),(17,3943))]
    Path('data/quadratic-tail-barrier.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in result:
        print(json.dumps({k:v for k,v in row.items() if k not in ('minimum_initial_moment','tail_factor_lower','loss_bound_lower')}))
