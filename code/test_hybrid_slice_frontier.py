import itertools
import math
import unittest
from hybrid_slice_frontier import individual_caps,pair_outputs
from extremal_slice_capacity import maximum_capacity
from audit_hybrid_slice_bound import dual_bound
from parity_capacity import pareto


class HybridTests(unittest.TestCase):
    def test_all_assignments_on_small_boxes(self):
        dims = ((1,3,2),(2,2,2))
        for counts in itertools.product(range(7),range(9)):
            caps = individual_caps(dims,counts)
            F = [maximum_capacity(d,r) for d,r in zip(dims,counts)]
            for b in (0,1):
                outputs = []
                for allocation in itertools.product((0,1),repeat=7):
                    even = [min(F[e],sum(c[e] for c,a in zip(caps,allocation) if a==e)) for e in (0,1)]
                    outputs.append(tuple(max(0,(6-(b==e))*counts[e]-F[e]-even[e]) for e in (0,1)))
                expected = pareto(outputs)
                actual = pair_outputs(dims,counts,7,b)
                self.assertEqual(actual,[list(p) for p in expected])
                self.assertLessEqual(dual_bound(dims,counts,7,b)['integer_lower'],min(map(sum,expected)))

    def test_monotonicity_of_clipped_layer_bounds(self):
        dims = (2,3,3)
        caps = [math.prod(dims[i] for i in range(3) if not s>>i&1) for s in range(1,8)]
        for allocation in itertools.product((0,1),repeat=7):
            for D in range(1,9):
                values = []
                for r in range(19):
                    F = maximum_capacity(dims,r)
                    even = min(F,sum(min(r,c) for c,a in zip(caps,allocation) if a))
                    values.append(max(0,D*r-F-even))
                self.assertEqual(values,sorted(values))


if __name__=='__main__':
    unittest.main()
