from fractions import Fraction
import itertools
import random
import unittest
from ortools.linear_solver import pywraplp
from support_tail_knapsack import minimum


class SupportTailKnapsackTests(unittest.TestCase):
    def test_matches_independent_fractional_lp_and_bounds_integer_choices(self):
        rng=random.Random(734)
        for iteration in range(25):
            volumes=[rng.randrange(4) for h in range(3)]
            penalties=[]
            for h in volumes:
                row=[0]
                for k in range(4):row.append(row[-1]+rng.randrange(8))
                penalties.append(row)
            target=rng.randrange(4*sum(volumes)+2)
            result=minimum(target,volumes,penalties)
            admissible=[sum(row[k] for row,k in zip(penalties,state)) for state in itertools.product(range(5),repeat=3)
                        if sum(v*k for v,k in zip(volumes,state))>=target]
            if not admissible:
                self.assertIsNone(result);continue
            self.assertLessEqual(result['value'],min(admissible))
            lp=pywraplp.Solver.CreateSolver('GLOP');objective=lp.Objective();objective.SetMinimization()
            credit=lp.Constraint(target,lp.infinity())
            for h,(v,row) in enumerate(zip(volumes,penalties)):
                choices=[lp.NumVar(0,lp.infinity(),f'x{h}_{k}') for k in range(5)]
                lp.Add(sum(choices)==1)
                for k,var in enumerate(choices):credit.SetCoefficient(var,v*k);objective.SetCoefficient(var,row[k])
            self.assertEqual(lp.Solve(),pywraplp.Solver.OPTIMAL)
            self.assertAlmostEqual(float(result['value']),objective.Value(),places=8)

    def test_fractional_solution_is_not_reported_as_an_integer_state(self):
        result=minimum(1,[2],[[0,3]])
        self.assertEqual(result['value'],Fraction(3,2))
        self.assertEqual(result['counts'],[Fraction(1,2)])
        self.assertEqual(result['slope'],Fraction(3,2))


if __name__=='__main__':unittest.main()
