import copy
import itertools
import unittest
from audit_refined_support_slopes import audit
from refine_support_slopes import refine_point
from support_tail_knapsack import minimum
from multitail_lp import LocalStates
import test_support_reweighting as fixtures
from reweight_nonpositive_support import optimize
from audit_support_reweighting import audit as audit_reweighting


class RefinedSupportTests(unittest.TestCase):
    def test_adaptive_refinement_matches_full_plane_enumeration(self):
        model=LocalStates([3,5],[0,1]);potentials=[[0]+[-k*k for k in range(1,15)],
              [0]+[-3*k for k in range(1,15)],[0]+[-2*k*k for k in range(1,17)]]
        result=refine_point(model,0,potentials,[0,1,2,4],4)
        candidates=[]
        for a,b in itertools.product(range(2),range(5)):
            row=minimum((1-a)*(4-b),[1],[[-v for v in potentials[2]]],-potentials[0][a]-potentials[1][b])
            if row:candidates.append(row['value'])
        value=min(candidates)
        self.assertEqual(result['lower_numerator'],value.numerator//value.denominator)

    def test_complete_certificate_and_modified_bound(self):
        data=fixtures.SupportReweightingTests().fixture(small=True)
        base=data['result']['pricing_history'][0];old=data['result']['positive_state_pricing']
        model=LocalStates(data['input']['high_primes'],data['input']['high_pattern']);rows=[]
        for x,point in enumerate(data['input']['core_residual']):
            potentials=[[base['point_weight_denominator']*sum(g[x][:t+1]) for t in range(len(g[x]))] for g in base['event_multipliers']]
            row=refine_point(model,point%2,potentials,old['point_slope_numerators'][x],old['denominator'])
            row['point_index']=x;rows.append(row)
        data['result']['refined_support_slopes']=dict(points=rows,complete=True)
        self.assertTrue(audit(data)['verified'])
        data['result']['support_reweighting']=optimize(data,seconds=5)
        self.assertTrue(audit_reweighting(data)['fixed_point_reweighting_barrier_verified'])
        invalid=copy.deepcopy(data)
        invalid['result']['refined_support_slopes']['points'][0]['lower_numerator']+=1
        with self.assertRaises(AssertionError):audit(invalid)


if __name__=='__main__':unittest.main()
