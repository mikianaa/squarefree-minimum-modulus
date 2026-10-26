from contextlib import redirect_stdout
import copy
import io
import itertools
import random
import unittest
from multitail_lp import LocalStates
from positive_state_dual import lower_bound,strengthen
from audit_positive_state_dual import minimum,audit
import test_weighted_coherent_dual as fixtures


class PositiveStateDualTests(unittest.TestCase):
    def test_envelope_is_below_exhaustive_positive_cost_and_above_raw_cost(self):
        model = LocalStates([3,5],[0,1])
        rng = random.Random(371)
        strict = False
        for parity in (0,1):
            D = [q-1-(b==parity) for q,b in zip(model.high,model.pattern)]
            for trial in range(8):
                potentials = []
                for limit in (14,14,16):
                    row = [0]
                    for k in range(limit):row.append(row[-1]-rng.randrange(8))
                    potentials.append(row)
                coefficients,den = [0,2,4,6,8],8
                bound = lower_bound(model,parity,potentials,coefficients,den)
                independent,count = minimum(D,potentials,coefficients,den)
                self.assertEqual(bound,independent)
                costs,raw = [],[]
                for state in itertools.product(range(D[0]+1),range(D[1]+1),range(17)):
                    cost = model.cost(parity,state)
                    penalty = -sum(potentials[h][v] for h,v in enumerate(state))
                    costs.append(den*max(0,cost)+penalty)
                    raw.append(den*cost+penalty)
                self.assertLessEqual(bound,min(costs))
                self.assertGreaterEqual(bound,min(raw))
                strict |= bound>min(raw)
        self.assertTrue(strict)

    def test_certificate_audit_and_tampering(self):
        data = fixtures.WeightedCoherentDualTests().fixture()
        with redirect_stdout(io.StringIO()):certificate = strengthen(data,steps=4)
        data['result']['positive_state_pricing'] = certificate
        report = audit(data)
        self.assertTrue(report['verified'])
        self.assertEqual(report['pool_numerator'],5)
        invalid = copy.deepcopy(data)
        invalid['result']['positive_state_pricing']['point_multipliers'][0] += 1
        with self.assertRaises(AssertionError):audit(invalid)
        invalid = copy.deepcopy(data)
        invalid['result']['positive_state_pricing']['point_slope_numerators'][0].append(2)
        with self.assertRaises(AssertionError):audit(invalid)

    def test_nonpositive_support_allows_large_slopes_but_is_not_a_cost_bound(self):
        model = LocalStates([3,5],[0,1]);D = [1,4]
        potentials = [[-9*k for k in range(15)],[-5*k for k in range(15)],[-13*k for k in range(17)]]
        coefficients,den = [0,1,2,4,8,16],1
        bound = lower_bound(model,0,potentials,coefficients,den,support_only=True)
        self.assertEqual(bound,minimum(D,potentials,coefficients,den,support_only=True)[0])
        admissible = [-sum(potentials[h][v] for h,v in enumerate(s))
                      for s in itertools.product(range(2),range(5),range(17)) if model.cost(0,s)<=0]
        self.assertLessEqual(bound,min(admissible))
        with self.assertRaises(AssertionError):lower_bound(model,0,potentials,coefficients,den)


if __name__=='__main__':unittest.main()
