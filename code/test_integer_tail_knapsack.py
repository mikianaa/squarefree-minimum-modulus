import itertools
import random
import unittest
from integer_tail_knapsack import minimum
from refine_integer_support import refine_point,diverse_support_states
from multitail_lp import LocalStates
from audit_refined_support_slopes import integer_minimum,audit
import test_support_reweighting as fixtures
import copy


class IntegerTailTests(unittest.TestCase):
    def test_diverse_states_are_integer_feasible_and_costed_after_rounding(self):
        model=LocalStates([3,5,7],[0,1,0])
        potentials=[[-(h+1)*k*k for k in range(15 if h<3 else 17)] for h in range(7)]
        rows=diverse_support_states(model,0,potentials,[0,1,2,4],4,2)
        self.assertGreater(len(rows),32)
        for row in rows:
            state=row['counts']
            self.assertTrue(all(type(k) is int for k in state))
            self.assertLessEqual(model.cost(0,state),0)
            self.assertEqual(row['penalty'],-sum(potentials[h][k] for h,k in enumerate(state)))
        plain=refine_point(model,0,potentials,[0,1,2,4],4)
        enriched=refine_point(model,0,potentials,[0,1,2,4],4,2)
        self.assertEqual(plain,{k:v for k,v in enriched.items() if k!='diverse_support'})

    def test_diverse_states_handle_no_tails_and_infeasible_slices(self):
        model=LocalStates([11],[0]);potentials=[[-k for k in range(15)]]
        for parity,expected in ((0,9),(1,10)):
            rows=diverse_support_states(model,parity,potentials,[0,1],1,3)
            self.assertEqual(rows,[dict(counts=[expected],penalty=expected)])
        self.assertEqual(diverse_support_states(model,0,potentials,[0,1],1,0),[])
        impossible=LocalStates([23],[0])
        self.assertEqual(diverse_support_states(impossible,0,potentials,[0,1],1,1),[])

    def test_requirement_dp_matches_exhaustive_integer_choices(self):
        rng=random.Random(919)
        for trial in range(60):
            volumes=[rng.randrange(5) for h in range(3)]
            penalties=[[0]+[rng.randrange(20) for k in range(4)] for h in volumes]
            target=rng.randrange(4*sum(volumes)+2)
            values=[sum(row[k] for row,k in zip(penalties,state)) for state in itertools.product(range(5),repeat=3)
                    if sum(v*k for v,k in zip(volumes,state))>=target]
            result=minimum(target,volumes,penalties)
            checked=integer_minimum(target,volumes,penalties)
            if not values:
                self.assertIsNone(result);self.assertIsNone(checked)
            else:
                self.assertEqual(result['value'],min(values))
                self.assertEqual(checked,min(values))
                self.assertGreaterEqual(sum(v*k for v,k in zip(volumes,result['counts'])),target)

    def test_full_support_pricing_matches_integer_enumeration(self):
        model=LocalStates([3,5],[0,1])
        potentials=[[-k*k for k in range(15)],[-3*k for k in range(15)],[-2*k*k for k in range(17)]]
        result=refine_point(model,0,potentials,[0,1,2,4],4)
        values=[-sum(potentials[h][k] for h,k in enumerate(state)) for state in itertools.product(range(2),range(5),range(17))
                if model.cost(0,state)<=0]
        self.assertEqual(result['lower_numerator'],min(values))
        self.assertLessEqual(model.cost(0,result['attaining_state']),0)

    def test_integer_certificate_audit_rejects_changed_attaining_state(self):
        data=fixtures.SupportReweightingTests().fixture(small=True)
        base=data['result']['pricing_history'][0];old=data['result']['positive_state_pricing']
        model=LocalStates(data['input']['high_primes'],data['input']['high_pattern']);rows=[]
        for x,point in enumerate(data['input']['core_residual']):
            potentials=[[base['point_weight_denominator']*sum(g[x][:t+1]) for t in range(len(g[x]))] for g in base['event_multipliers']]
            row=refine_point(model,point%2,potentials,old['point_slope_numerators'][x],old['denominator'])
            row['point_index']=x;rows.append(row)
        data['result']['integer_support_refinement']=dict(points=rows,complete=True)
        self.assertTrue(audit(data)['exact_integer_support_minima'])
        invalid=copy.deepcopy(data);invalid['result']['integer_support_refinement']['points'][0]['attaining_state'][0]=0
        with self.assertRaises(AssertionError):audit(invalid)


if __name__=='__main__':unittest.main()
