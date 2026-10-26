import copy
import itertools
import math
import unittest
import numpy as np
from search_joint_core_projections import JointProjections, search, polish
from audit_joint_core_projections import audit


class JointProjectionTests(unittest.TestCase):
    def test_multilinear_update_matches_direct_formula(self):
        source = dict(core_residual=[1, 8, 13, 59], high_primes=[11, 13, 17], high_pattern=[0, 1, 0])
        model = JointProjections(source)
        rng = np.random.default_rng(32)
        for trial in range(30):
            counts = rng.integers(0, 15, size=(7, 4))
            direct = []
            for x, point in enumerate(source['core_residual']):
                free = [max(0, q-1-(b == point % 2)-int(counts[i, x])) for i, (q, b) in
                        enumerate(zip(source['high_primes'], source['high_pattern']))]
                direct.append(math.prod(free)-sum(int(counts[h, x])*math.prod(free[i] for i in range(3) if i not in s)
                    for h, s in enumerate(model.groups[3:], 3)))
            np.testing.assert_array_equal(model.evaluate(counts), direct)
            for h in range(7):
                changed = counts.copy(); changed[h] += 1
                np.testing.assert_array_equal(model.evaluate(changed)-model.evaluate(counts), model.increment(counts, h))
                table = model.group_costs(counts, h)
                self.assertTrue(np.all(table[:-1] >= table[1:]))

    def test_independent_audit_and_corruption_rejection(self):
        source = dict(core_residual=[1], high_primes=[11], high_pattern=[0])
        result = search(source, seconds=0, initial=[[1]*14])
        self.assertTrue(audit(result)['simultaneous_projection_relaxation_feasible'])
        self.assertFalse(audit(result)['covering_constructed'])
        for field in ('residue', 'F', 'objective'):
            bad = copy.deepcopy(result)
            if field == 'residue': bad['group_assignments'][0][0] = 3
            if field == 'F': bad['local_F'][0] += 1
            if field == 'objective': bad['positive_volume'] += 1
            with self.assertRaises(AssertionError): audit(bad)

    def test_whole_group_search_has_independently_checked_incumbent(self):
        source = dict(core_residual=[1], high_primes=[11], high_pattern=[0])
        result = search(source, seconds=0, initial=[[0]*14])
        self.assertGreater(audit(result)['positive_volume'], 0)
        improved = polish(result, seconds=3, seconds_per_group=1)
        self.assertEqual(improved['positive_volume'], 0)
        self.assertTrue(audit(improved)['simultaneous_projection_relaxation_feasible'])


if __name__ == '__main__': unittest.main()
