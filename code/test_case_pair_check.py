from fractions import Fraction
import unittest
from audit_case_pair_check import score


class CasePairCheckTests(unittest.TestCase):
    def test_fractional_pair_violation_with_unit_resource(self):
        scenario=dict(fixed=[],count_distributions=[[[1,2],[0,1],[1,2]],[[0,1],[1,1],[0,1]]])
        value,capacity=score([1,2],[1,2,3],scenario,0,1,3,2)
        self.assertEqual(value,Fraction(3,2));self.assertEqual(capacity,4)
        with self.assertRaises(AssertionError):score([1,2],[1,2,3],scenario,0,1,2,2)

    def test_fixed_residues_tighten_the_pair_capacity(self):
        scenario=dict(fixed=[[2,0],[3,0]],count_distributions=[[[1,1]],[[1,1]]])
        value,capacity=score([1,2],[1,2,3],scenario,0,1,2,2,True)
        self.assertEqual(capacity,3);self.assertEqual(value,1)
        with self.assertRaises(AssertionError):score([1,2],[1,2,3],scenario,0,1,2,2,False)


if __name__=='__main__':unittest.main()
