from fractions import Fraction
import itertools
import unittest
from ortools.linear_solver import pywraplp

from two_configuration_mixtures import point_coupling,evaluate,structure


class TwoConfigurationTests(unittest.TestCase):
    def test_reduced_transport_matches_full_count_distribution_lp(self):
        source = dict(core_residual=[1],high_primes=[5,7,11],high_pattern=[1,0,1])
        a,b = [1,3,0,1,2,1,0],[3,1,2,2,1,2,2]
        probability = Fraction(1,3)
        reduced = point_coupling(source,1,a,b,probability)
        solver = pywraplp.Solver.CreateSolver('GLOP')
        states = list(itertools.product(*(range(max(v,w)+1) for v,w in zip(a,b))))
        variables = [solver.NumVar(0,1,str(j)) for j in range(len(states))]
        solver.Add(sum(variables)==1)
        for h,(v,w) in enumerate(zip(a,b)):
            for t in range(1,max(v,w)+1):
                supply = float((1-probability)*(v>=t)+probability*(w>=t))
                solver.Add(sum(p for row,p in zip(states,variables) if row[h]>=t)<=supply)
        def explicit_cost(row):
            f0,f1,f2 = 3-row[0],6-row[1],9-row[2]
            return f0*f1*f2-row[3]*f2-row[4]*f1-row[5]*f0-row[6]
        solver.Minimize(sum(explicit_cost(row)*p for row,p in zip(states,variables)))
        self.assertEqual(solver.Solve(),pywraplp.Solver.OPTIMAL)
        self.assertAlmostEqual(reduced['objective'],solver.Objective().Value(),places=7)
        self.assertAlmostEqual(sum(p for row,p in reduced['support']),1,places=10)
        for h,(v,w) in enumerate(zip(a,b)):
            for t in range(1,max(v,w)+1):
                supply = float((1-probability)*(v>=t)+probability*(w>=t))
                used = sum(p for row,p in reduced['support'] if row[h]>=t)
                self.assertLessEqual(used,supply+1e-8)

    def test_explicit_assignment_mixture_passes_independent_integer_audit(self):
        source = dict(core_residual=[1,2,4],high_primes=[11,13],high_pattern=[0,1])
        groups,shapes = structure(source)
        A = [[1%d for d in shape] for shape in shapes]
        B = [[2%d for d in shape] for shape in shapes]
        result = evaluate(source,A,B,Fraction(1,2),save_witness=True)
        self.assertTrue(result['audit']['verified'])
        self.assertEqual(result['audit']['groups_checked'],3)
        self.assertEqual(result['audit']['resource_assignments_checked'],6)
        self.assertEqual(result['audit']['threshold_links_checked'],132)
        exact = result['coherent_primal_witness']
        self.assertAlmostEqual(result['objective'],exact['objective_numerator']/exact['denominator'],places=5)


if __name__=='__main__':
    unittest.main()
