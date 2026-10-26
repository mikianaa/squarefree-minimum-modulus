import copy
from contextlib import redirect_stdout
import io
import unittest
from fractional_group_scenarios import maximize
from case_relaxation_barrier import generate
from audit_case_relaxation_barrier import audit


class CaseRelaxationBarrierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=dict(core_residual=[1,2],high_primes=[3],high_pattern=[0])
        shapes=[d for d in range(3,211) if 210%d==0]
        events=[[[0]+[-1]*len(shapes) for p in source['core_residual']]]
        data=dict(input=source,result=dict(pricing_history=[dict(event_multipliers=events)]))
        seed=dict(input=source,result=dict(support_seed=[dict(point_index=0,counts=[2]),dict(point_index=1,counts=[1])]))
        rows=[maximize([1,2],shapes,[list(range(15))]*2,coef,[[3,a]]) for a in (1,2) for coef in ([1,0],[0,1])]
        scenarios=dict(input=source,baseline_event_multipliers=events,group_scenarios=[dict(axis=0,scenarios=rows)])
        trees=dict(input=source,group_bounds=[dict(axis=0,certificate=dict(nodes=[
            dict(fixed=[],resource=3,children=[1,2]),dict(fixed=[[3,1]]),dict(fixed=[[3,2]])]))])
        with redirect_stdout(io.StringIO()):cls.certificate=generate(data,seed,scenarios,trees,5)

    def test_exact_dominating_mixtures_are_independently_verified(self):
        report=audit(self.certificate)
        self.assertTrue(report['all_nonnegative_group_point_coefficients'])
        self.assertGreater(report['minimum_normalized_surplus'][0],0)
        self.assertTrue(report['integer_rounding_cannot_restore_positive_margin'])
        self.assertFalse(report['all_nonnegative_event_weights_barrier_verified'])
        self.assertGreater(report['deficient_thresholds'],0)
        self.assertFalse(report['actual_assignment_barrier_verified'])

    def test_rejects_noncovering_local_state_and_bad_probability(self):
        bad=copy.deepcopy(self.certificate);bad['point_mixtures'][0]['states'][0]['counts']=[0]
        with self.assertRaises(AssertionError):audit(bad)
        bad=copy.deepcopy(self.certificate);bad['group_mixtures'][0]['scenarios'][0]['count_distributions'][0][0]=[2,1]
        with self.assertRaises(AssertionError):audit(bad)

    def test_rejects_missing_case_and_false_surplus(self):
        bad=copy.deepcopy(self.certificate);bad['case_trees'][0]['nodes'][0]['children'].pop()
        with self.assertRaises(AssertionError):audit(bad)
        bad=copy.deepcopy(self.certificate);bad['coordinate_surpluses'][0][0]=[10**9,1]
        with self.assertRaises(AssertionError):audit(bad)

    def test_threshold_deficit_can_defeat_only_the_stated_mixture(self):
        deficit=audit(self.certificate)['threshold_deficits'][0]
        new=[[[0]*len(row) for row in axis] for axis in self.certificate['baseline_event_multipliers']]
        new[deficit['axis']][deficit['point_index']][deficit['threshold']]=-1
        report=audit(self.certificate,new)
        self.assertTrue(report['comparison_defeats_this_witness'])
        old=audit(self.certificate,self.certificate['baseline_event_multipliers'])
        self.assertLess(old['comparison_local_minus_group'][0],0)


if __name__=='__main__':unittest.main()
