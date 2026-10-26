"""Checks of the candidate generator; the mathematical certificate has separate verifiers."""
import unittest
import numpy as np

from optimize_stoploss_breakpoints import candidate_probabilities, stage_costs, compress_rounded_schedule
from audit_verify import D_for_stage


class BreakpointTests(unittest.TestCase):
    def test_cutoff_six_adds_modulus_seven(self):
        r=np.array([1/7])
        old=stage_costs(7,r,4,8)[0]
        new=stage_costs(7,r,4,6)[0]
        np.testing.assert_allclose(new-old,np.full(4,1/7),rtol=0,atol=1e-15)

    def test_one_coordinate_minimum_dominates_dense_grid(self):
        # A nontrivial distribution and continuation; checks threshold handling
        # at the exceptional stages 5,7 and the ordinary stages thereafter.
        w = np.array([.2, .35, .3, .1, .04, .01])
        v = np.array([0., .01, .03, .1, .5, 2.])
        for p in (5, 7, 11, 31, 101):
            options = candidate_probabilities(p)
            grid = np.linspace(1/p, 1, 10001)
            def evaluate(rs):
                return (stage_costs(p, rs, len(w)) + (1-rs[:, None])*v + rs[:, None]*np.r_[v[1:], 0.]) @ w
            self.assertLessEqual(float(evaluate(options).min()), float(evaluate(grid).min()) + 1e-12)

    def test_export_and_lookup_at_all_boundaries(self):
        source = [0., 0., .2000001, .2000002, .6, .5, .333333]
        edges, params = compress_rounded_schedule(source)
        values = [round(x*1000000) for x in params]
        expected = [0, 0, 200000, 200000, 500000, 500000, 333333]
        self.assertEqual([D_for_stage(k, edges, values) for k in range(1, 8)], expected)
        with self.assertRaises(ValueError):
            D_for_stage(8, edges, values)


if __name__ == '__main__':
    unittest.main()
