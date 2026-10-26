import copy
import unittest
from audit_weighted_coherent_dual import audit
from event_score_majorant import upper_bound


class WeightedCoherentDualTests(unittest.TestCase):
    def fixture(self):
        source = dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        shape = [d for d in range(3,211) if 210%d==0]
        events = [dict(point=x,threshold=t,weight=1) for x in range(2) for t in range(1,15)]
        cap = upper_bound([1,2],shape,[(r['point'],r['threshold']) for r in events],[1]*len(events))
        certificate = dict(point_weight_denominator=1,event_scale=1,denominator=1,
                           point_weights=[1,1],event_multipliers=[[[0]+[-1]*14 for x in range(2)]],
                           point_multipliers=[10,9],capacity_numerator=5,pool_numerator=5,
                           pool_upper_bounds=[14],pool_sizes=[2],
                           group_pricing=[dict(events=events,capacity_certificate=cap)])
        return dict(input=source,result=dict(group_columns=[[[x%d for d in shape] for x in (1,2)]],
                                             pricing_history=[certificate]))

    def test_recovers_known_integer_obstruction_five(self):
        report = audit(self.fixture())
        self.assertTrue(report['verified'])
        self.assertEqual(report['unrestricted_capacity_numerator'],5)
        self.assertEqual(report['pool_numerator'],5)
        self.assertTrue(report['excludes_unrestricted_core_pattern'])

    def test_point_weight_change_without_repricing_is_rejected(self):
        data = self.fixture();data['result']['pricing_history'][0]['point_weights'][0] = 0
        with self.assertRaises(AssertionError):audit(data)

    def test_inadequate_capacity_majorant_is_rejected(self):
        data = self.fixture()
        data['result']['pricing_history'][0]['group_pricing'][0]['capacity_certificate']['slopes'][0] = [1,2]
        with self.assertRaises(AssertionError):audit(data)


if __name__=='__main__':unittest.main()
