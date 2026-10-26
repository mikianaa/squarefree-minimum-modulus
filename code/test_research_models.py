"""Checks of mathematical reductions, including an unquotiented enumeration."""
import itertools
import json
import pathlib
import sys
import unittest
from unittest.mock import patch

from three_axis_bound import coupled_bound


class ThreeAxisTests(unittest.TestCase):
    def test_full_enumeration_agrees_with_symmetry_quotient(self):
        reduced = coupled_bound()
        full = coupled_bound(canonical=False)
        self.assertEqual(reduced['cases'], full['cases'])
        self.assertEqual(reduced['uncovered_lower'], 22)

    def test_existing_witness_attains_bound(self):
        root = pathlib.Path(__file__).resolve().parents[1]
        witness = json.loads((root/'data/min-uncovered-lcm210-witness22.json').read_text())
        family = [(e['modulus'],e['residue']) for e in witness['system']]
        self.assertEqual(len({d for d,a in family}), len(family))
        self.assertTrue(all(d >= 3 and 210 % d == 0 for d,a in family))
        missed = [x for x in range(210) if not any(x % d == a for d,a in family)]
        self.assertEqual(len(missed),22)

    def test_slice_extension_lower_bound(self):
        ds = [d for d in range(2,211) if 210 % d == 0]
        lower = min(10*r-sum(min(r,210//d) for d in ds) for r in range(22,211))
        self.assertEqual(lower,26)

    def test_2310_witness_respects_new_frontier(self):
        root = pathlib.Path(__file__).resolve().parents[1]
        witness = json.loads((root/'data/min-uncovered-lcm2310-feasible121.json').read_text())
        family = {e['modulus']:e['residue'] for e in witness['system']}
        pattern = [family[2*p] % 2 for p in (3,5,7,11)]
        counts = [sum(not any(n%d == a for d,a in family.items()) for n in range(e,2310,2)) for e in (0,1)]
        if pattern[0]:
            pattern = [b^1 for b in pattern]
            counts.reverse()
        rows = json.loads((root/'data/parity-capacity-30030-certificate.json').read_text())
        row = next(r for r in rows if r['pattern'] == pattern)
        self.assertTrue(any(a <= counts[0] and b <= counts[1] for a,b in row['frontier']))
        self.assertEqual(sum(counts),121)

    def test_convex_minimum_matches_integer_scan(self):
        from parity_capacity import minimum_deficit
        for caps in itertools.combinations_with_replacement(range(1,5),4):
            for coefficient in range(1,6):
                for lower in range(7):
                    expected = min(max(0,coefficient*r-sum(min(r,c) for c in caps))
                                   for r in range(lower,max(lower,max(caps))+2))
                    self.assertEqual(minimum_deficit(lower,coefficient,caps),expected)


class SolverTests(unittest.TestCase):
    def test_small_model_with_and_without_label_symmetry(self):
        from parity_cover import solve
        for symmetry in (False,True):
            for pattern in ((0,0),(0,1)):
                r = solve(30, pattern, seconds=10, workers=1, minimize=True, symmetry=symmetry)
                self.assertEqual(r['status'],'OPTIMAL')
                self.assertEqual(r['uncovered'], 8 if pattern == (0,0) else 7)

    def test_unknown_slice_assignment_cannot_be_discarded(self):
        import few_progression_residual as module
        # The first-stage problem is feasible with these parameters. Force its
        # second-stage solver to time out, and ensure the loop stops immediately.
        argv = ['few_progression_residual.py','30','7','2','1','--time','10','--workers','1','--iters','2']
        with patch.object(sys, 'argv', argv), patch.object(module,'assign_slices',return_value=(None,'UNKNOWN')) as assignment:
            module.main()
        self.assertEqual(assignment.call_count,1)


if __name__ == '__main__':
    unittest.main()
