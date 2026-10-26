from collections import defaultdict
import itertools
import random
import unittest
import numpy as np
from ortools.linear_solver import pywraplp
from coherent_transport_repair import clip_distribution,drop_cost_bounds
from decomposed_coherent_groups import PointOracle
from multitail_lp import LocalStates


class CoherentTransportRepairTests(unittest.TestCase):
    def test_diverse_pricing_preserves_every_coordinate_slice_minimum(self):
        model = LocalStates([3,5],[0,1])
        rng = random.Random(19)
        for e in (0,1):
            penalties = [np.array([rng.randrange(-30,1) for k in range(limit+1)]) for limit in (14,14,16)]
            rows = model.price(e,penalties,3,batch=1,diverse=True)
            D = [q-1-(b==e) for q,b in zip(model.high,model.pattern)]
            all_states = list(itertools.product(range(D[0]+1),range(D[1]+1),range(17)))
            def value(s):return 3*model.cost(e,s)-sum(int(p[v]) for p,v in zip(penalties,s))
            self.assertEqual(rows[0][0],min(map(value,all_states)))
            for h in (0,1):
                for v in range(D[h]+1):
                    self.assertEqual(min(c for c,s in rows if s[h]==v),min(value(s) for s in all_states if s[h]==v))

    def test_one_hit_bounds_hold_for_all_small_two_prime_states(self):
        model = LocalStates([5,7],[1,0])
        for e in (0,1):
            D = [q-1-(b==e) for q,b in zip(model.high,model.pattern)]
            bounds = drop_cost_bounds(model.high,model.pattern,e)
            for state in itertools.product(range(D[0]+1),range(D[1]+1),range(17)):
                for h,count in enumerate(state):
                    if count:
                        lower = list(state);lower[h] -= 1
                        self.assertLessEqual(model.cost(e,lower)-model.cost(e,state),bounds[h])

    def test_quantile_repair_preserves_mass_and_achieves_minimal_hit_removal(self):
        rng = random.Random(32)
        for repeat in range(100):
            rows = [(tuple(rng.randrange(5) for h in range(3)),rng.randrange(1,10)) for j in range(8)]
            mass = sum(w for s,w in rows)
            supply = [[0]+sorted([rng.randrange(mass+1) for t in range(4)],reverse=True) for h in range(3)]
            fixed,removed = clip_distribution(rows,supply,mass)
            self.assertEqual(sum(w for s,w in fixed),mass)
            for h in range(3):
                expected = sum(max(0,sum(w for s,w in rows if s[h]>=t)-supply[h][t]) for t in range(1,5))
                self.assertEqual(removed[h],expected)
                for t in range(1,5):
                    self.assertLessEqual(sum(w for s,w in fixed if s[h]>=t),supply[h][t])
                self.assertEqual(sum(w*s[h] for s,w in rows)-sum(w*s[h] for s,w in fixed),expected)

    def test_capped_oracle_agrees_with_unpenalized_exhaustive_lp(self):
        source = dict(core_residual=[1],high_primes=[3,5],high_pattern=[0,1])
        model = LocalStates(source['high_primes'],source['high_pattern'])
        states = list(itertools.product(range(3),range(4),range(17)))
        supplies = [[0]+[.65,.2]+[0]*12,[0]+[.8,.6,.1]+[0]*11,
                    [0]+[.9,.5,.25,.1]+[0]*12]
        solver = pywraplp.Solver.CreateSolver('GLOP')
        variables = [solver.NumVar(0,1,str(j)) for j in range(len(states))]
        solver.Add(sum(variables)==1)
        for h,line in enumerate(supplies):
            for t,p in enumerate(line[1:],1):
                solver.Add(sum(v for s,v in zip(states,variables) if s[h]>=t)<=p)
        solver.Minimize(sum(v*model.cost(1,s) for s,v in zip(states,variables)))
        self.assertEqual(solver.Solve(),pywraplp.Solver.OPTIMAL)
        oracle = PointOracle(source,model,0,[])
        answer = oracle.solve(supplies,passes=100)
        self.assertAlmostEqual(answer['primal'],solver.Objective().Value(),places=6)
        self.assertAlmostEqual(answer['lower'],answer['primal'],places=5)
        for h,line in enumerate(answer['multipliers']):
            self.assertGreaterEqual(min(line),-oracle.bounds[h]*answer['denominator'])


if __name__=='__main__':unittest.main()
