"""Independent finite checks of the majorization and plane-deficit formula."""
import itertools
import unittest
import numpy as np
from plane_capacity_bound import balanced_lower,plane_deficits,ALLOCATIONS


def packing_vectors(weights,bins):
    states = {(0,)*bins}
    for w in weights:
        following = set()
        for state in states:
            for i in range(bins):
                v = list(state)
                v[i] += w
                following.add(tuple(sorted(v,reverse=True)))
        states = following
    return states


class PlaneTests(unittest.TestCase):
    def test_majorization_against_all_small_packings(self):
        for length in range(1,7):
            for weights in itertools.combinations_with_replacement(range(1,5),length):
                for bins in range(1,5):
                    lower = balanced_lower(list(weights),bins)
                    self.assertEqual(sum(lower),sum(weights))
                    self.assertEqual(lower,sorted(lower,reverse=True))
                    for actual in packing_vectors(weights,bins):
                        self.assertTrue(all(sum(actual[:k])>=sum(lower[:k]) for k in range(1,bins+1)))

    def test_tensor_values_by_direct_cell_enumeration(self):
        profile = [10,12,10,7,4,9,3,4,4,5,3,2,1,2,1,1]
        arrays,dims,caps = plane_deficits(profile,(11,13,17),(1,1,1))
        indices = [(0,0,0),(127,127,127),(0,127,43),(67,25,99)]
        # Include the minimizing allocation, independently of the tensor formula.
        indices.append(tuple(map(int,np.unravel_index(int((arrays[0]+arrays[1]).argmin()),arrays[0].shape))))
        for index in indices:
            for e in (0,1):
                vectors = []
                for D,k in zip(dims[e],index):
                    bits = ALLOCATIONS[k]
                    weights = [c[e] for c in caps]+[c[e] for b,c in zip(bits,caps) if b==e]
                    vectors.append(balanced_lower(weights,D))
                expected = sum(max(0,profile[e]-sum(load)) for load in itertools.product(*vectors))
                self.assertEqual(int(arrays[e][index]),expected)
        self.assertEqual(int((arrays[0]+arrays[1]).min()),6236)


if __name__=='__main__':
    unittest.main()
