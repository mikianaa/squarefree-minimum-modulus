import copy
from contextlib import redirect_stdout
import io
import unittest
from unittest.mock import patch
from audit_support_failure import audit
from expand_support_catalogue import expand
from price_curve_basis_dual import price
from probe_curve_basis_cuts import probe
import test_group_support_weights as fixtures


class SupportFailureTests(unittest.TestCase):
    def fixture(self):
        return dict(input=dict(core_residual=[1,2],high_primes=[11],high_pattern=[0]),denominator=1,
                    event_multipliers=[[[0]*14+[-1],[0]*14+[-1]]],
                    point_states=[dict(point_index=0,counts=[10]),dict(point_index=1,counts=[9])],
                    group_assignments=[[1]*14],claimed_score_difference=-1)

    def test_explicit_assignments_and_local_states_certify_failure(self):
        report=audit(self.fixture())
        self.assertTrue(report['verified']);self.assertTrue(report['these_weights_cannot_exclude_unrestricted_core'])
        self.assertFalse(report['covering_constructed'])
        for change in ('counts','assignment','difference'):
            bad=self.fixture()
            if change=='counts':bad['point_states'][0]['counts']=[0]
            if change=='assignment':bad['group_assignments'][0][0]=3
            if change=='difference':bad['claimed_score_difference']=0
            with self.assertRaises(AssertionError):audit(bad)

    def test_positive_difference_does_not_certify_failure(self):
        data=self.fixture();data['event_multipliers']=[[[0]+[-1]*14 for x in range(2)]]
        data['claimed_score_difference']=5
        self.assertFalse(audit(data)['these_weights_cannot_exclude_unrestricted_core'])

    def test_catalogue_expansion_preserves_old_columns_and_remaps_case_indices(self):
        data=fixtures.GroupSupportWeightTests().fixture([11],[0])
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=r['point_index'],counts=r['attaining_state'])
                    for r in data['result']['integer_support_refinement']['points']]))
        scenarios=dict(input=data['input'],group_scenarios=[dict(axis=0,scenarios=[])])
        with redirect_stdout(io.StringIO()):
            numerical=probe(data,seed,scenarios,[],seconds=5)
            certificate=price(data,seed,scenarios,[],numerical,10**6)
        certificate['result']['threshold_support_optimization']['group_constraint_rows']=[[0,2]]
        candidate=[1]*13+[2]
        with patch('expand_support_catalogue.maximize',return_value=dict(status='FEASIBLE',candidate_assignments=[dict(residues=candidate)])):
            with redirect_stdout(io.StringIO()):source,continued,witness,report=expand(data,certificate,.1)
        self.assertEqual(source['result']['group_columns'][0][:2],data['result']['group_columns'][0])
        self.assertEqual(source['result']['group_columns'][0][2],candidate)
        self.assertEqual(continued['result']['threshold_support_optimization']['group_constraint_rows'],[[0,3]])
        self.assertEqual(continued['result']['integer_support_refinement'],certificate['result']['integer_support_refinement'])
        self.assertEqual(source['result']['pricing_history'][0]['event_multipliers'],data['result']['pricing_history'][0]['event_multipliers'])
        self.assertTrue(audit(witness)['verified'])


if __name__=='__main__':unittest.main()
