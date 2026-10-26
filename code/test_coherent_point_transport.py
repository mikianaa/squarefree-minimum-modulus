import itertools
import unittest
from ortools.linear_solver import pywraplp
from coherent_point_transport import solve_point
from multitail_lp import LocalStates
import numpy as np
from fractions import Fraction
from two_configuration_mixtures import point_coupling


class CoherentPointTransportTests(unittest.TestCase):
    def test_plane_supply_above_available_coordinate_count_has_zero_extra_penalty(self):
        source = dict(high_primes=[3,5],high_pattern=[0,1])
        supplies = [[0]+[1]*14,[0]+[1]*14,[0]+[1]*16]
        result = solve_point(source,1,supplies)
        self.assertEqual(result['objective'],-16)
        self.assertEqual(result['event_multipliers'][0][3:],[0]*12)
        self.assertEqual(result['event_multipliers'][1][4:],[0]*11)
        self.check_lifted_dual(source,1,supplies,result)

    def test_equal_volume_aggregation_matches_unaggregated_three_plane_transport(self):
        source = dict(high_primes=[5,7,11],high_pattern=[1,0,1])
        a,b = [1,3,0,1,2,1,0],[3,1,2,2,1,2,2]
        probability = Fraction(1,3)
        supplies = [[0]+[float((1-probability)*(v>=t)+probability*(w>=t)) for t in range(1,(14 if h<3 else 16)+1)]
                    for h,(v,w) in enumerate(zip(a,b))]
        result = solve_point(source,1,supplies)
        other = point_coupling(source,1,a,b,probability)
        self.assertAlmostEqual(result['objective'],other['objective'],places=7)
        self.assertLess(result['variables'],result['plane_states']*5)
        self.check_lifted_dual(source,1,supplies,result)

    def check_lifted_dual(self,source,point,supply,result):
        model = LocalStates(source['high_primes'],source['high_pattern'])
        multipliers = result['event_multipliers']
        common = [np.cumsum(np.array(line,dtype=np.int64)) for line in multipliers]
        constant = model.price(point%2,common,result['denominator'])[0][0]
        bound = (constant+sum(v*p for row,prob in zip(multipliers,supply) for v,p in zip(row,prob)))/result['denominator']
        self.assertAlmostEqual(bound,result['objective'],places=4)

    def test_general_marginals_match_full_state_lp(self):
        source = dict(high_primes=[3,5],high_pattern=[0,1])
        supplies = [[0,.65,.2]+[0]*12,[0,.8,.6,.1]+[0]*11,
                    [0,.9,.5,.25,.1]+[0]*12]
        result = solve_point(source,1,supplies)
        states = list(itertools.product(range(3),range(4),range(17)))
        lp = pywraplp.Solver.CreateSolver('GLOP')
        variables = [lp.NumVar(0,1,str(j)) for j in range(len(states))]
        lp.Add(sum(variables)==1)
        for h,line in enumerate(supplies):
            for t,p in enumerate(line[1:],1):
                lp.Add(sum(v for state,v in zip(states,variables) if state[h]>=t)<=p)
        lp.Minimize(sum(((2-s[0])*(3-s[1])-s[2])*v for s,v in zip(states,variables)))
        self.assertEqual(lp.Solve(),pywraplp.Solver.OPTIMAL)
        self.assertAlmostEqual(result['objective'],lp.Objective().Value(),places=7)
        self.assertAlmostEqual(sum(p for s,p in result['support']),1,places=9)
        for h,line in enumerate(supplies):
            for t,p in enumerate(line[1:],1):
                self.assertLessEqual(sum(w for s,w in result['support'] if s[h]>=t),p+1e-9)
        self.check_lifted_dual(source,1,supplies,result)

    def test_gaps_between_supplied_counts_need_no_intermediate_plane_states(self):
        source = dict(high_primes=[7,11],high_pattern=[0,1])
        supplies = [[0]+[.7]*5+[0]*9,[0]+[.8]*3+[.4]*4+[0]*7,
                    [0]+[1]*2+[.6]*4+[.1]*5+[0]*5]
        result = solve_point(source,1,supplies)
        self.assertEqual(result['plane_states'],6)
        # Independently enumerate every intermediate plane and tail count.
        states = list(itertools.product(range(7),range(10),range(17)))
        lp = pywraplp.Solver.CreateSolver('GLOP')
        variables = [lp.NumVar(0,1,str(j)) for j in range(len(states))]
        lp.Add(sum(variables)==1)
        for h,line in enumerate(supplies):
            for t,p in enumerate(line[1:],1):
                lp.Add(sum(v for state,v in zip(states,variables) if state[h]>=t)<=p)
        lp.Minimize(sum(((6-s[0])*(9-s[1])-s[2])*v for s,v in zip(states,variables)))
        self.assertEqual(lp.Solve(),pywraplp.Solver.OPTIMAL)
        self.assertAlmostEqual(result['objective'],lp.Objective().Value(),places=7)
        self.check_lifted_dual(source,1,supplies,result)


if __name__=='__main__':unittest.main()
