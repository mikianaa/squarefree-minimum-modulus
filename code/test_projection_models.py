"""Small independent controls for the new projection and SAT audit models."""
import itertools
import json
from pathlib import Path
import unittest
from pysat.solvers import Solver

from core_projection import solve as projected_solve
from layered_core_projection import solve as layered_solve
from sat_uncovered_audit import build as sat_build, check_rup


class ProjectionTests(unittest.TestCase):
    def test_capacity_projection_recovers_known_finite_obstruction(self):
        result = projected_solve([11,13],[0,1,1,1,1],[8,14],seconds=10,workers=1)
        self.assertEqual(result['status'],'INFEASIBLE')

    def test_zero_residual_is_a_feasible_relaxation(self):
        # A deliberately unattainable empty core is a positive control for the
        # relaxation implementation, not a covering-system construction.
        result = projected_solve([11],[0,1,1,1],[0,0],seconds=5,workers=1)
        self.assertIn(result['status'],('OPTIMAL','FEASIBLE'))
        result = projected_solve([11],[0,1,1,1],None,seconds=5,workers=1,minimize=True,frontier=[(0,0)])
        self.assertEqual(result['status'],'OPTIMAL')
        self.assertEqual(result['relaxation_uncovered'],0)
        result = layered_solve([11],[0,1,1,1],[(0,0)],[],seconds=5,workers=1)
        self.assertIn(result['status'],('OPTIMAL','FEASIBLE'))

    def test_new_optimal_witnesses_by_direct_integer_coverage(self):
        root = Path(__file__).resolve().parents[1]
        entries = json.loads((root/'data/min-uncovered-2310-five-patterns-optimal.json').read_text())
        for entry in entries:
            family = [(v['modulus'],v['residue']) for v in entry['system']]
            self.assertEqual(len(set(d for d,a in family)),len(family))
            self.assertTrue(all(d>=3 and 2310%d==0 for d,a in family))
            uncovered = sum(not any(n%d==a for d,a in family) for n in range(2310))
            self.assertEqual(uncovered,entry['uncovered'])
            self.assertEqual(entry['status'],'OPTIMAL')
            # This last equality checks the saved report, not its lower-bound
            # proof. The coverage count independently checks the upper witness.
            self.assertEqual(entry['lower_bound'],uncovered)


class SatAuditTests(unittest.TestCase):
    def test_sat_encoding_at_known_thresholds(self):
        for pattern,minimum in (([0,0],8),([0,1],7)):
            for symmetry in (False,True):
                for limit in (minimum-1,minimum):
                    clauses,_ = sat_build([3,5],pattern,limit,symmetry)
                    with Solver(name='minisat22',bootstrap_with=clauses) as solver:
                        self.assertEqual(solver.solve(),limit>=minimum)

    def test_rup_checker_accepts_refutation_and_rejects_invalid_step(self):
        clauses = [[1,2],[-1,2],[1,-2],[-1,-2]]
        self.assertTrue(check_rup(clauses,['2 0','0'])['verified'])
        self.assertFalse(check_rup([[1,2]],['1 0','0'])['verified'])


if __name__ == '__main__':
    unittest.main()
