from contextlib import redirect_stdout
import io
import unittest
from audit_joint_core_projections import audit
from search_joint_core_projections import JointProjections, search
from optimize_joint_projections import solve


class CoupledProjectionTests(unittest.TestCase):
    def initial(self, primes):
        source = dict(core_residual=[1, 2], high_primes=primes, high_pattern=[0]*len(primes))
        problem = JointProjections(source)
        return search(source, seconds=0, initial=[[0]*len(s) for s in problem.shapes])

    def run_model(self, *args, **kwargs):
        with redirect_stdout(io.StringIO()):
            result = solve(*args, seconds=5, **kwargs)
        self.assertTrue(audit(result)['verified'])
        self.assertFalse(audit(result)['covering_constructed'])
        return result

    def test_one_axis_matches_direct_incidence_bound(self):
        # These two points disagree modulo every resource. 14 resources
        # cannot meet their total demand 10+9=19; allocating 9+5 attains 5.
        initial = self.initial([11])
        for tables in (True, False):
            result = self.run_model(initial, tables=tables)
            self.assertEqual(result['optimization']['status'], 'OPTIMAL')
            self.assertEqual(result['positive_volume'], 5)

    def test_all_tails_fixed_planes_matches_total_incidence(self):
        # The zero plane residues hit neither point. The tail modulus 1
        # hits both, and the other 15 tail resources hit at most one each.
        # Thus the exact objective is 10*12+9*11-(2+15)=202.
        initial = self.initial([11, 13])
        for tables in (True, False):
            result = self.run_model(initial, axes=[2], tables=tables)
            self.assertEqual(result['optimization']['status'], 'OPTIMAL')
            self.assertEqual(result['positive_volume'], 202)
            self.assertEqual(result['group_assignments'][:2], initial['group_assignments'][:2])

    def test_coupled_planes_and_full_polynomial_agree(self):
        initial = self.initial([11, 13])
        for axes, tables in (([0, 1], True), ([0, 1], False), (None, False)):
            result = self.run_model(initial, axes=axes, tables=tables)
            self.assertEqual(result['positive_volume'], 0)
            self.assertTrue(audit(result)['simultaneous_projection_relaxation_feasible'])


if __name__ == '__main__': unittest.main()
