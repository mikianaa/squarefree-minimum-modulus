import unittest
from cached_support_slopes import coefficients,safe_maximum
from multitail_lp import LocalStates
from refine_integer_support import refine_point
from positive_state_dual import lower_bound


class CachedSlopeTests(unittest.TestCase):
    def test_rescaling_rounding_budget_and_overflow_filter(self):
        planes=[dict(slope=[3,2]),dict(slope=[999999999999999999999,1])]
        self.assertEqual(coefficients([0,1],planes,3,5,12,10),[0,1,2,3])
        many=[dict(slope=[k,1]) for k in range(1,100)]
        result=coefficients([0],many,1,1,12,100)
        self.assertEqual(len(result),13)
        self.assertEqual(coefficients([0,1],planes,3,5,0),[0,1])
        self.assertEqual(coefficients([0],planes[:1],3,5,12,reward_ratio=(2,1)),[0,3])
        self.assertEqual(coefficients([0],planes[:1],3,5,12,reward_ratio=(0,1)),[0])

    def test_more_slopes_preserve_exact_minimum_and_strengthen_envelope(self):
        model=LocalStates([3,5],[0,1])
        potentials=[[-k*k for k in range(15)],[-3*k for k in range(15)],[-2*k*k for k in range(17)]]
        plain=refine_point(model,0,potentials,[0,1],4)
        slopes=coefficients([0,1],plain['refined_planes'],4,4,12,safe_maximum(model,0,potentials))
        enriched=refine_point(model,0,potentials,slopes,4)
        self.assertEqual(enriched['lower_numerator'],plain['lower_numerator'])
        self.assertGreaterEqual(enriched['initial_envelope_minimum'],plain['initial_envelope_minimum'])
        self.assertLessEqual(len(enriched['refined_planes']),len(plain['refined_planes']))
        maximum=safe_maximum(model,0,potentials)
        lower_bound(model,0,potentials,[maximum],4,True)
        with self.assertRaises(AssertionError):lower_bound(model,0,potentials,[maximum+1],4,True)


if __name__=='__main__':unittest.main()
