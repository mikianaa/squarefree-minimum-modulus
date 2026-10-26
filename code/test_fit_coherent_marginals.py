from contextlib import redirect_stdout
import io
import unittest

from fit_coherent_marginals import fit
from two_configuration_mixtures import evaluate,structure


class FitCoherentMarginalsTests(unittest.TestCase):
    def test_known_mixture_is_recovered_from_its_configuration_pool(self):
        source = dict(core_residual=[1,2,4],high_primes=[11,13],high_pattern=[0,1])
        groups,shapes = structure(source)
        A = [[1%d for d in shape] for shape in shapes]
        B = [[2%d for d in shape] for shape in shapes]
        exact = evaluate(source,A,B,save_witness=True)['coherent_primal_witness']
        pool = dict(group_columns=[[a,b] for a,b in zip(A,B)])
        with redirect_stdout(io.StringIO()):
            result = fit(source,exact,pool)
        self.assertTrue(result['audit']['verified'])
        self.assertTrue(all(r['domination']>1-1e-8 for r in result['group_reports']))
        self.assertAlmostEqual(result['audit']['objective_numerator']/exact['denominator'],
                               exact['objective_numerator']/exact['denominator'],places=5)

    def test_zero_demand_and_single_group_are_well_defined(self):
        source = dict(core_residual=[1],high_primes=[11],high_pattern=[0])
        groups,shapes = structure(source)
        witness = dict(denominator=10,support=[dict(point_index=0,counts=[0],weight=10)])
        pool = dict(group_columns=[[[0 for d in shapes[0]]]])
        with redirect_stdout(io.StringIO()):
            result = fit(source,witness,pool)
        self.assertTrue(result['audit']['verified'])
        self.assertEqual(result['group_reports'][0]['domination'],1)


if __name__=='__main__':
    unittest.main()
