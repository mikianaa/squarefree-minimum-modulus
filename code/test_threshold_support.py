from contextlib import redirect_stdout
import io
import copy
from fractions import Fraction
import unittest
import numpy as np
import test_group_support_weights as fixtures
from optimize_threshold_support import optimize,scenario_tails,valid_state_pool
from audit_refined_support_slopes import audit


class ThresholdSupportTests(unittest.TestCase):
    def test_compact_master_retains_omitted_seed_states_and_audits(self):
        data=fixtures.GroupSupportWeightTests().fixture([3,5],[0,1])
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        seed['result']['support_seed'] += [dict(point_index=0,counts=[0,0,k]) for k in range(6,17)]
        with redirect_stdout(io.StringIO()):result=optimize(data,seed,seconds=5,iterations=1,compact_master=True,dual_simplex=True)
        saved={(row['point_index'],tuple(row['counts'])) for row in result['result']['support_seed']}
        self.assertTrue(all((row['point_index'],tuple(row['counts'])) in saved for row in seed['result']['support_seed']))
        history=result['result']['threshold_support_optimization']['history'][-1]
        self.assertLess(history['active_point_constraints'],history['point_states'])
        self.assertTrue(audit(result)['verified'])

    def test_cached_slopes_are_recomputed_and_saved_for_independent_audit(self):
        data=fixtures.GroupSupportWeightTests().fixture([11],[0])
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        with redirect_stdout(io.StringIO()):first=optimize(data,seed,seconds=5,iterations=1)
        # Cache entries are arbitrary nonnegative slope suggestions; the old
        # proof is not reused. Force a new value to exercise its saved audit.
        first['result']['integer_support_refinement']['points'][0]['refined_planes'].append(dict(slope=[999,7]))
        with redirect_stdout(io.StringIO()):resumed=optimize(data,first,seconds=5,iterations=1,cached_slopes=12,dual_simplex=True)
        slopes=resumed['result']['positive_state_pricing']['point_slope_numerators'][0]
        old_mass=sum(-sum(axis[0]) for axis in first['result']['pricing_history'][-1]['event_multipliers'])
        new_mass=sum(-sum(axis[0]) for axis in resumed['result']['pricing_history'][-1]['event_multipliers'])
        self.assertGreater(new_mass,0)
        value=Fraction(999*new_mass,7*old_mass)
        self.assertTrue(value.numerator//value.denominator in slopes or -(-value.numerator//value.denominator) in slopes)
        self.assertTrue(audit(resumed)['verified'])

    def test_warm_group_rows_require_identical_matrices(self):
        data=fixtures.GroupSupportWeightTests().fixture([11],[0])
        data['result']['group_columns'][0]=[[1]*k+[2]*(14-k) for k in range(15)]
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        with redirect_stdout(io.StringIO()):
            first=optimize(data,seed,seconds=5)
            resumed=optimize(data,first,seconds=5,diverse=1,share_parity=True,cached_slopes=12)
            changed=copy.deepcopy(data)
            changed['result']['group_columns'][0].reverse()
            cold=optimize(changed,first,seconds=5,diverse=1)
        first_rows=first['result']['threshold_support_optimization']['group_constraint_rows']
        self.assertGreater(sum(map(len,first_rows)),0)
        self.assertEqual(resumed['result']['threshold_support_optimization']['reused_group_constraints'],sum(map(len,first_rows)))
        self.assertEqual(cold['result']['threshold_support_optimization']['reused_group_constraints'],0)
        reports=[audit(result) for result in (first,resumed,cold)]
        self.assertTrue(all(report['verified'] for report in reports))
        self.assertEqual([report['pool_numerator'] for report in reports],[reports[0]['pool_numerator']]*3)

    def test_new_curve_can_improve_audited_obstruction(self):
        data=fixtures.GroupSupportWeightTests().fixture([11],[0])
        # For these two points each resource hits exactly one point, so all
        # fifteen count splits exhaust the occupied-residue score profiles.
        data['result']['group_columns'][0]=[[1]*k+[2]*(14-k) for k in range(15)]
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        with redirect_stdout(io.StringIO()):result=optimize(data,seed,seconds=5)
        report=audit(result)
        self.assertTrue(report['verified'])
        self.assertTrue(report['excludes_unrestricted_core_pattern'])
        self.assertGreater(report['pool_numerator']/report['denominator'],5)
        self.assertEqual(result['result']['group_columns'],data['result']['group_columns'])

    def test_scenario_tail_probabilities_include_fixed_and_unit_hits(self):
        row=dict(fixed=[[3,1]],count_distributions=[[0.25,0.75],[0.5,0.5]])
        tails=scenario_tails([1,2],[1,2,3],row)
        np.testing.assert_allclose(tails,[[1,1,.75],[1,.5,0]])

    def test_scenario_search_preserves_exact_actual_pool(self):
        data=fixtures.GroupSupportWeightTests().fixture([11],[0])
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        # A deliberately generous numerical search row must not become a proof
        # pool entry or an asserted unrestricted upper bound.
        scenario=dict(fixed=[],count_distributions=[[0]*14+[1],[0]*14+[1]])
        scenarios=dict(input=data['input'],group_scenarios=[dict(axis=0,scenarios=[scenario])])
        with redirect_stdout(io.StringIO()):result=optimize(data,seed,scenarios,seconds=5)
        self.assertTrue(audit(result)['verified'])
        self.assertEqual(result['result']['group_columns'],data['result']['group_columns'])
        self.assertLessEqual(result['result']['threshold_support_optimization']['numerical_scenario_margin'],1e-6)

    def test_known_pool_keeps_only_valid_nonpositive_states(self):
        source=dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        rows=[dict(point_index=x,counts=[k]) for x in range(2) for k in (-1,0,9,10,11,14,15)]
        valid=valid_state_pool(source,rows)
        self.assertEqual(valid[0].tolist(),[[10]])
        self.assertEqual(valid[1].tolist(),[[9]])

    def test_parity_sharing_preserves_domains_and_does_not_cross_parities(self):
        source=dict(core_residual=[1,2,3,4],high_primes=[11],high_pattern=[0])
        rows=[dict(point_index=0,counts=[10]),dict(point_index=1,counts=[9])]
        plain=valid_state_pool(source,rows)
        shared=valid_state_pool(source,rows,True)
        self.assertEqual([array.tolist() for array in plain],[[[10]],[[9]],[],[]])
        self.assertEqual([array.tolist() for array in shared],[[[10]],[[9]],[[10]],[[9]]])
        self.assertIs(shared[0],shared[2])
        self.assertIs(shared[1],shared[3])
        # Sharing also works when one parity has no core points.
        source['core_residual']=[1,3]
        shared=valid_state_pool(source,[rows[0]],True)
        self.assertEqual([array.tolist() for array in shared],[[[10]],[[10]]])

    def test_two_curve_basis_recovers_audited_threshold_improvement(self):
        data=fixtures.GroupSupportWeightTests().fixture([11],[0])
        data['result']['group_columns'][0]=[[1]*k+[2]*(14-k) for k in range(15)]
        seed=dict(input=data['input'],result=dict(support_seed=[dict(point_index=row['point_index'],counts=row['attaining_state'])
            for row in data['result']['integer_support_refinement']['points']]))
        new=copy.deepcopy(data)
        new['result']['pricing_history'][0]['event_multipliers']=[[[0]*15 for x in range(2)]]
        for x,t in enumerate((10,9)):new['result']['pricing_history'][0]['event_multipliers'][0][x][t]=-14
        with redirect_stdout(io.StringIO()):result=optimize(data,seed,seconds=5,basis=[new])
        report=audit(result)
        self.assertTrue(report['excludes_unrestricted_core_pattern'])
        self.assertEqual(report['pool_numerator']/report['denominator'],14)
        self.assertEqual(result['result']['threshold_support_optimization']['curve_basis_size'],2)


if __name__=='__main__':unittest.main()
