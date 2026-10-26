import copy
from contextlib import redirect_stdout
import io
import unittest
from probe_curve_basis_barrier import probe,projection_representatives
from probe_curve_basis_cuts import probe as cut_probe
from certify_curve_basis_barrier import certify
from audit_curve_basis_barrier import audit
from escape_curve_basis_barrier import escape
from price_curve_basis_dual import price
from audit_refined_support_slopes import audit as audit_price
import test_group_support_weights as fixtures


class CurveBasisBarrierTests(unittest.TestCase):
    def test_projection_deduplication_never_rounds_nearby_values(self):
        self.assertEqual(projection_representatives([[1,2],[1,2],[1,2+1e-12],[0,2]]).tolist(),[0,2,3])

    def test_identical_projected_columns_preserve_lp_and_clp_agrees(self):
        data=fixtures.GroupSupportWeightTests().fixture([3,5],[0,1])
        for h,axis in enumerate(data['result']['pricing_history'][0]['event_multipliers']):
            for x,row in enumerate(axis):axis[x]=[0,-len(row)+1]+[0]*(len(row)-2)
        seed=dict(input=data['input'],result=dict(support_seed=[]))
        for row in data['result']['integer_support_refinement']['points']:
            for tail in (15,16):
                seed['result']['support_seed'].append(dict(point_index=row['point_index'],counts=row['attaining_state'][:-1]+[tail]))
        for columns in data['result']['group_columns']:columns.append(list(columns[0]))
        scenarios=dict(input=data['input'],group_scenarios=[dict(axis=h,scenarios=[]) for h in range(3)])
        with redirect_stdout(io.StringIO()):
            original=probe(data,seed,scenarios,[],group_pool='all',seconds=5)
            reduced=probe(data,seed,scenarios,[],group_pool='all',seconds=5,deduplicate=True)
            alternate=probe(data,seed,scenarios,[],group_pool='all',seconds=5,deduplicate=True,solver_name='CLP')
        for result in (original,reduced,alternate):self.assertEqual(result['status'],'OPTIMAL')
        self.assertAlmostEqual(original['numerical_margin'],reduced['numerical_margin'])
        self.assertAlmostEqual(original['numerical_margin'],alternate['numerical_margin'])
        self.assertLess(reduced['local_state_count'],original['local_state_count'])
        self.assertLess(reduced['group_state_count'],original['group_state_count'])

    def fixture(self,high):
        data=fixtures.GroupSupportWeightTests().fixture([high],[0])
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        scenarios=dict(input=data['input'],group_scenarios=[dict(axis=0,scenarios=[])])
        return data,seed,scenarios

    def test_positive_numerical_mixture_is_not_automatically_a_proof(self):
        data,seed,scenarios=self.fixture(3)
        basis=copy.deepcopy(data)
        basis['result']['pricing_history'][0]['event_multipliers']=[[[0]*14+[-14] for x in range(2)]]
        with redirect_stdout(io.StringIO()):result=probe(data,seed,scenarios,[basis],seconds=5,group_pool='all')
        self.assertEqual(result['status'],'OPTIMAL')
        self.assertGreater(result['numerical_margin'],0)
        self.assertEqual(result['coordinates'],4)
        self.assertFalse(result['verified'])

    def test_known_separating_curves_prevent_dominating_mixture(self):
        data,seed,scenarios=self.fixture(11)
        events=[[[0]*15 for x in range(2)]]
        events[0][0][10]=-14;events[0][1][9]=-14
        data['result']['pricing_history'][0]['event_multipliers']=events
        with redirect_stdout(io.StringIO()):result=probe(data,seed,scenarios,[],seconds=5,group_pool='all',dual_simplex=False)
        self.assertEqual(result['status'],'OPTIMAL')
        self.assertAlmostEqual(result['numerical_margin'],-.5)
        self.assertAlmostEqual(result['coordinate_dual_weight_sum'],1)
        self.assertTrue(all(row['weight']>0 for row in result['coordinate_dual_weights']))
        self.assertFalse(result['verified'])

    def test_constraint_generation_closes_omitted_group_rows_and_recovers_mixtures(self):
        data,seed,scenarios=self.fixture(11)
        events=[[[0]*15 for x in range(2)]];events[0][0][10]=-14;events[0][1][9]=-14
        data['result']['pricing_history'][0]['event_multipliers']=events
        with redirect_stdout(io.StringIO()):
            result=cut_probe(data,seed,scenarios,[],seconds=5,group_pool='all',batch=1)
        self.assertEqual(result['status'],'OPTIMAL');self.assertTrue(result['finite_pool_closed'])
        self.assertAlmostEqual(result['numerical_margin'],-.5)
        self.assertGreater(result['history'][0]['group_violation'],0)
        for row in result['point_mixtures']+result['group_mixtures']:
            self.assertAlmostEqual(sum(s['weight'] for s in row['states']),1)
        with redirect_stdout(io.StringIO()):
            resumed=cut_probe(data,seed,scenarios,[],seconds=5,group_pool='all',warm_start=result)
        self.assertAlmostEqual(result['numerical_margin'],resumed['numerical_margin'])
        self.assertEqual(len(resumed['history']),1)
        self.assertAlmostEqual(resumed['finite_pool_numerical_duality_gap'],0)
        suggestion=copy.deepcopy(result);suggestion['status']='FEASIBLE'
        with redirect_stdout(io.StringIO()):priced=price(data,seed,scenarios,[],suggestion,10**9)
        self.assertTrue(audit_price(priced)['excludes_unrestricted_core_pattern'])

    def test_generated_positive_mixture_passes_independent_rational_audit(self):
        data,seed,scenarios=self.fixture(3)
        with redirect_stdout(io.StringIO()):result=cut_probe(data,seed,scenarios,[],seconds=5)
        self.assertTrue(result['finite_pool_closed']);self.assertGreater(result['numerical_margin'],0)
        exact=certify(data,scenarios,[],result)
        self.assertTrue(audit(exact)['curve_cone_barrier_verified'])

    def test_dual_direction_produces_independently_audited_integer_certificate(self):
        data,seed,scenarios=self.fixture(11)
        events=[[[0]*15 for x in range(2)]];events[0][0][10]=-14;events[0][1][9]=-14
        data['result']['pricing_history'][0]['event_multipliers']=events
        with redirect_stdout(io.StringIO()):
            numerical=probe(data,seed,scenarios,[],seconds=5,group_pool='all')
            result=price(data,seed,scenarios,[],numerical,10**9)
        report=audit_price(result)
        self.assertTrue(report['verified'])
        self.assertTrue(report['excludes_unrestricted_core_pattern'])
        self.assertGreater(result['result']['curve_basis_dual_pricing']['normalized_pool_gap'],0)

    def exact_fixture(self):
        data,seed,scenarios=self.fixture(3);basis=copy.deepcopy(data)
        basis['result']['pricing_history'][0]['event_multipliers']=[[[0]*14+[-14] for x in range(2)]]
        with redirect_stdout(io.StringIO()):numerical=probe(data,seed,scenarios,[basis],seconds=5,group_pool='all')
        return certify(data,scenarios,[basis],numerical)

    def test_exact_multiple_curve_barrier_and_threshold_escape(self):
        certificate=self.exact_fixture();report=audit(certificate)
        self.assertTrue(report['curve_cone_barrier_verified'])
        self.assertEqual(report['basis_count'],2)
        self.assertTrue(report['actual_assignment_only_barrier'])
        self.assertGreater(report['minimum_normalized_surplus'][0],0)
        self.assertFalse(report['all_allowed_nonnegative_threshold_weights_barrier_verified'])
        new=escape(certificate,report)
        checked=audit(certificate,new['result']['pricing_history'][0]['event_multipliers'])
        self.assertTrue(checked['comparison_defeats_this_witness'])
        self.assertEqual(checked['comparison_local_minus_group'],new['escape_direction']['expected_local_minus_group'])

    def test_integer_case_repair_and_pair_violation_rejection(self):
        data,seed,scenarios=self.fixture(3)
        shape=[d for d in range(3,211) if 210%d==0]
        row=[0.]*15;row[0]=.125;row[8]=.875
        scenario=dict(fixed=[],count_distributions=[row,list(row)],resource_distributions=[
            dict(resource=d,residues=[1,2],probabilities=[.5,.5]) for d in shape])
        scenarios['group_scenarios'][0]['scenarios']=[scenario]
        numerical=dict(input=data['input'],numerical_margin=.1,
            point_mixtures=[dict(point_index=s['point_index'],states=[dict(counts=s['counts'],weight=1.)]) for s in seed['result']['support_seed']],
            group_mixtures=[dict(axis=0,states=[dict(kind='scenario',index=0,weight=1.)])])
        certificate=certify(data,scenarios,[],numerical)
        report=audit(certificate)
        self.assertGreater(report['conditional_pair_inequalities_checked'],0)
        self.assertEqual(report['fractional_group_states'],1)
        bad=copy.deepcopy(certificate);D=bad['probability_denominator']
        row=[0]*15;row[0]=D//8;row[8]=7*D//8
        # These marginals retain mean 7, exactly their resource supply, but
        # P(count_x>=8)+P(count_y>=8)=7/4 violates the pair capacity 14.
        bad['group_mixtures'][0]['states'][0]['count_numerators']=[row,list(row)]
        with self.assertRaises(AssertionError):audit(bad)

    def test_rejects_noncovering_state_bad_assignment_and_false_surplus(self):
        certificate=self.exact_fixture()
        bad=copy.deepcopy(certificate);bad['point_mixtures'][0]['states'][0]['counts']=[0]
        with self.assertRaises(AssertionError):audit(bad)
        bad=copy.deepcopy(certificate);bad['group_mixtures'][0]['states'][0]['assignment'][0]=3
        with self.assertRaises(AssertionError):audit(bad)
        bad=copy.deepcopy(certificate);bad['coordinate_surpluses'][0]['surplus']=[999,1]
        with self.assertRaises(AssertionError):audit(bad)


if __name__=='__main__':unittest.main()
