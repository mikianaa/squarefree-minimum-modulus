"""Independent finite checks of the sharp compression bound."""
from collections import Counter
import itertools
import math
import unittest
from extremal_slice_capacity import maximum_capacity,extension_bound,exclusion_threshold


def actual_capacity(residual,period):
    return sum(max(Counter(n%d for n in residual).values(),default=0)
               for d in range(2,period+1) if period%d==0)


class CompressionTests(unittest.TestCase):
    def test_all_subsets_of_six(self):
        for r in range(7):
            observed = max(actual_capacity(R,6) for R in itertools.combinations(range(6),r))
            self.assertEqual(observed,maximum_capacity([2,3],r))

    def test_all_small_subsets_of_thirty(self):
        for r in range(5):
            observed = max(actual_capacity(R,30) for R in itertools.combinations(range(30),r))
            self.assertEqual(observed,maximum_capacity([2,3,5],r))

    def test_attaining_construction_at_every_size_210(self):
        primes = [2,3,5,7]
        order = sorted(range(210),key=lambda n:(-sum(n%p==0 for p in primes),tuple(n%p for p in primes)))
        for r in range(211):
            self.assertEqual(actual_capacity(order[:r],210),maximum_capacity(primes,r))

    def test_convex_extension_minimum_and_threshold(self):
        for primes,q in [([2,3],5),([2,3,5],7),([2,3,5,7],11)]:
            period = math.prod(primes)
            values = [(q-1)*r-maximum_capacity(primes,r) for r in range(period+1)]
            self.assertEqual(exclusion_threshold(primes,q),next((r for r,v in enumerate(values) if v>0),period+1))
            for lower in range(period+1):
                self.assertEqual(extension_bound(primes,lower,q)['uncovered_lower'],max(0,min(values[lower:])))


if __name__=='__main__':
    unittest.main()
