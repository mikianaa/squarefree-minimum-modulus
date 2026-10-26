from collections import Counter
import itertools
import math
import unittest
from core_hit_concentration import enumerate_capacities
from audit_fibre_concentration import gray_capacities


class FibreTests(unittest.TestCase):
    def test_subset_capacity_against_direct_sets(self):
        residual = [8,13,19,59,67,68]
        moduli = [3,5,6,15,30,105,210]
        expected = [[-1]*5 for _ in range(3)]
        for mask in range(1<<len(residual)):
            points = [n for j,n in enumerate(residual) if mask>>j&1]
            counts = [sum(n%2==e for n in points) for e in (0,1)]
            value = sum(max(Counter(n%d for n in points).values(),default=0) for d in moduli)
            expected[counts[0]][counts[1]] = max(expected[counts[0]][counts[1]],value)
        self.assertEqual(enumerate_capacities(residual,moduli,chunk_bits=3)['two_layer_capacities'],expected)
        self.assertEqual(gray_capacities(residual,moduli),expected)

    def test_hit_excess_bounds_for_all_small_projection_choices(self):
        residual = [1,2,4,5]
        moduli = [2,3,6]
        data = enumerate_capacities(residual,moduli)
        H = data['two_layer_capacities']
        for labels in itertools.product(*(range(d) for d in moduli)):
            hits = [sum(n%d==a for d,a in zip(moduli,labels)) for n in residual]
            for t0,t1 in itertools.product(range(4),repeat=2):
                actual = sum(max(0,h-(t0 if n%2==0 else t1)) for n,h in zip(residual,hits))
                upper = max(value-t0*a-t1*b for a,row in enumerate(H) for b,value in enumerate(row))
                self.assertLessEqual(actual,upper)

    def test_plane_product_against_actual_coordinate_coverage(self):
        residual = [8,13,19]
        high = [11,13]
        shapes = [3,5,6]
        # Colliding high-coordinate labels show why counting incidences is a
        # relaxation; actual distinct deleted labels cannot exceed incidences.
        for assignments in itertools.product((2,3),repeat=6):
            classes = [(i,d,1,assignments[3*i+j]) for i in range(2) for j,d in enumerate(shapes)]
            bound = actual = 0
            for n in residual:
                alphabets = [range(1+(n%2==1),q) for q in high]
                counts = [sum(n%d==a for i,d,a,label in classes if i==axis) for axis in range(2)]
                bound += math.prod(max(0,len(A)-c) for A,c in zip(alphabets,counts))
                actual += sum(not any(n%d==a and point[i]==label for i,d,a,label in classes)
                              for point in itertools.product(*alphabets))
            self.assertLessEqual(bound,actual)


if __name__=='__main__':
    unittest.main()
