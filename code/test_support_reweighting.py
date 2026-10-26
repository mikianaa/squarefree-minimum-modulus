import copy
from contextlib import redirect_stdout
import io
import unittest
import test_weighted_coherent_dual as fixtures
from positive_state_dual import strengthen
from reweight_nonpositive_support import optimize
from audit_support_reweighting import audit


class SupportReweightingTests(unittest.TestCase):
    def fixture(self,small=False):
        data=fixtures.WeightedCoherentDualTests().fixture()
        if small:
            data['input']['high_primes']=[3]
            c=data['result']['pricing_history'][0]
            c.update(point_multipliers=[2,1],capacity_numerator=-11,pool_numerator=-11)
        with redirect_stdout(io.StringIO()):data['result']['positive_state_pricing']=strengthen(data,support_only=True)
        data['result']['support_reweighting']=optimize(data,seconds=5)
        return data

    def test_recovers_positive_known_bound(self):
        report=audit(self.fixture())
        self.assertTrue(report['verified'])
        self.assertEqual(report['pool_numerator'],5*report['denominator'])
        self.assertFalse(report['fixed_point_reweighting_barrier_verified'])

    def test_exact_mixture_certifies_reweighting_barrier(self):
        data=self.fixture(small=True);report=audit(data)
        self.assertTrue(report['fixed_point_reweighting_barrier_verified'])
        self.assertGreater(report['minimum_score_excess_numerator'],0)
        self.assertEqual(report['pool_numerator'],0)
        invalid=copy.deepcopy(data)
        invalid['result']['support_reweighting']['barrier_witness']['point_score_excess_numerators'][0]+=1
        with self.assertRaises(AssertionError):audit(invalid)


if __name__=='__main__':unittest.main()
