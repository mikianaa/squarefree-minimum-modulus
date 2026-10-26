import itertools
import json
from pathlib import Path
import unittest
import numpy as np
from audit_fibre_concentration import gray_capacities
from portable_tail_certificates import CoreCapacityUpper,PortableDuals,concentration_bounds


class PortableTailTests(unittest.TestCase):
    def test_joint_upper_dominates_every_small_actual_subset(self):
        for pattern in ((0,0,0),(0,0,1),(0,1,0),(0,1,1)):
            model = CoreCapacityUpper(pattern)
            indices = list(range(0,len(model.points),7))
            mask = sum(1<<j for j in indices)
            R = [model.numbers[j] for j in indices]
            actual = np.array(gray_capacities(R,[d for d in range(3,211) if 210%d==0]))
            self.assertTrue(np.all(model.upper(mask)>=actual))

    def test_concentration_transform_against_direct_maximum(self):
        H = np.array([[0,2,3],[3,4,5],[4,5,6]])
        for shapes in (14,16):
            actual = concentration_bounds(H,shapes)
            for t0,t1 in itertools.product(range(17),repeat=2):
                expected = max(int(v)+(a+b+max(a,b) if shapes==16 else 0)-t0*a-t1*b for a,row in enumerate(H) for b,v in enumerate(row))
                self.assertEqual(actual[t0,t1],expected)

    def test_portable_evaluation_reproduces_audited_fixed_core_bounds(self):
        data = json.loads(Path('data/conditional-tail-core-pattern-certificates.json').read_text())
        H = np.array(json.loads(Path('data/core-hit-concentration-7-15.json').read_text())['two_layer_capacities'])
        actual = PortableDuals(data).bounds(H).tolist()
        self.assertEqual(actual,[r['certificate']['integer_lower'] for r in data['cases']])

    def test_optimistic_deletions_are_contained_in_every_actual_triple_result(self):
        model = CoreCapacityUpper((0,1,1))
        R = model.numbers[::3]
        counts = [sum(n%2==e for n in R) for e in (0,1)]
        targets = [[max(0,counts[e]-1-(e==extra)) for e in (0,1)] for extra in (0,1)]
        for a in range(105):
            for b in range(210):
                remaining = [sum(n%2==e and n%105!=a and n!=b for n in R) for e in (0,1)]
                self.assertTrue(any(all(t[e]<=remaining[e] for e in (0,1)) for t in targets))


if __name__=='__main__':
    unittest.main()
