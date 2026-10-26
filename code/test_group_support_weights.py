from contextlib import redirect_stdout
import io
import unittest
from optimize_group_support_weights import optimize
from audit_refined_support_slopes import audit
from retarget_support_certificate import retarget
from branched_event_bound import bound
import copy


class GroupSupportWeightTests(unittest.TestCase):
    def fixture(self,high,pattern):
        source=dict(core_residual=[1,2],high_primes=high,high_pattern=pattern)
        n=len(high);G=2**n-1
        shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
        base=dict(event_scale=1,event_multipliers=[[[0]+[-1]*len(shape) for x in source['core_residual']] for shape in shapes])
        states=[dict(point_index=x,attaining_state=[min(14,q-1-(b==p%2)) for q,b in zip(high,pattern)]+[16]*(G-n))
                for x,p in enumerate(source['core_residual'])]
        return dict(input=source,result=dict(pricing_history=[base],
                    group_columns=[[[p%d for d in shape] for p in source['core_residual']] for shape in shapes],
                    integer_support_refinement=dict(points=states)))

    def test_recovers_known_obstruction_with_auditable_new_weights(self):
        with redirect_stdout(io.StringIO()):result=optimize(self.fixture([11],[0]),seconds=5)
        report=audit(result)
        self.assertTrue(report['exact_integer_support_minima'])
        self.assertEqual(report['pool_numerator'],5*report['denominator'])
        self.assertTrue(report['excludes_unrestricted_core_pattern'])

    def test_three_group_weights_emit_valid_support_certificate(self):
        with redirect_stdout(io.StringIO()):result=optimize(self.fixture([3,5],[0,1]),seconds=5,min_weight=.5)
        self.assertIsNotNone(result)
        report=audit(result)
        self.assertTrue(report['verified'])
        self.assertEqual(report['pool_numerator'],result['result']['group_support_optimization']['exact_pool_numerator'])
        self.assertEqual(len(result['result']['group_support_optimization']['weights']),3)

    def test_expanded_pool_preserves_independent_support_proof(self):
        with redirect_stdout(io.StringIO()):result=optimize(self.fixture([11],[0]),seconds=5)
        pool=copy.deepcopy(result)
        pool['result']['group_columns'][0].append([1]*7+[2]*7)
        expanded=retarget(result,pool,True);report=audit(expanded)
        self.assertEqual(report['pool_numerator'],5*report['denominator'])
        self.assertEqual(expanded['result']['retargeted_bounds']['configurations'],3)
        shapes=[d for d in range(3,211) if 210%d==0]
        events=result['result']['pricing_history'][0]['group_pricing'][0]['events']
        tree=bound([1,2],shapes,events,1)
        branches=dict(input=result['input'],complete=True,group_bounds=[dict(axis=0,certificate=tree)])
        branched=retarget(result,pool,False,branches)
        self.assertEqual(branched['result']['pricing_history'][0]['group_pricing'][0]['capacity_certificate']['method'],'branch_affine')
        self.assertTrue(audit(branched)['verified'])

    def test_fractional_scenarios_do_not_enter_actual_pool_certificate(self):
        data=self.fixture([11],[0])
        scenarios=dict(input=data['input'],baseline_event_multipliers=data['result']['pricing_history'][0]['event_multipliers'],
                       group_scenarios=[dict(axis=0,scenarios=[dict(payoffs=[14.,14.])])])
        with redirect_stdout(io.StringIO()):result=optimize(data,seconds=5,min_weight=.5,scenarios=scenarios)
        report=audit(result)
        self.assertTrue(report['verified'])
        self.assertEqual(result['result']['group_columns'],data['result']['group_columns'])
        self.assertEqual(result['result']['group_support_optimization']['scenario_count'],1)
        self.assertLess(result['result']['group_support_optimization']['numerical_scenario_margin'],report['pool_numerator']/report['denominator'])


if __name__=='__main__':unittest.main()
