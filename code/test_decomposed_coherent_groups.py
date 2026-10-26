from contextlib import redirect_stdout
from fractions import Fraction
import io
import itertools
import unittest
from ortools.linear_solver import pywraplp

from decomposed_coherent_groups import solve,group_prices
from multitail_lp import LocalStates
from optimize_two_configuration_mixtures import optimize
from two_configuration_mixtures import point_coupling,structure


class DecomposedCoherentTests(unittest.TestCase):
    def test_point_dual_is_a_supporting_function_for_independent_probabilities(self):
        source = dict(core_residual=[1],high_primes=[5,7,11],high_pattern=[1,0,1])
        a,b = [1,3,0,1,2,1,0],[3,1,2,2,1,2,2]
        base = point_coupling(source,1,a,b,[.3]*7)
        for theta in ([0]*7,[1]*7,[.8,.1,.6,.2,.9,.1,.4],[.1,.9,.1,.8,.2,.7,.3]):
            value = point_coupling(source,1,a,b,theta)['objective']
            supporting = base['supporting_intercept']+sum(g*p for g,p in zip(base['gradient'],theta))
            self.assertLessEqual(supporting,value+1e-7)

    def test_two_assignment_master_matches_full_joint_lp(self):
        source = dict(core_residual=[1,2],high_primes=[3,5],high_pattern=[0,1])
        groups,shapes = structure(source)
        A = [[1%d for d in shape] for shape in shapes]
        B = [[2%d for d in shape] for shape in shapes]
        result = optimize(source,A,B,seconds=10)
        solver = pywraplp.Solver.CreateSolver('GLOP')
        theta = [solver.NumVar(0,1,str(h)) for h in groups]
        objective = solver.Objective();objective.SetMinimization()
        for x in source['core_residual']:
            da,db = 2-(x%2==0),4-(x%2==1)
            states = list(itertools.product(range(da+1),range(db+1),range(17)))
            variables = [solver.NumVar(0,1,str((x,j))) for j in range(len(states))]
            solver.Add(sum(variables)==1)
            for row,var in zip(states,variables):
                objective.SetCoefficient(var,(da-row[0])*(db-row[1])-row[2])
            for h,shape in enumerate(shapes):
                ca = sum(x%d==a for d,a in zip(shape,A[h]))
                cb = sum(x%d==b for d,b in zip(shape,B[h]))
                for t in range(1,len(shape)+1):
                    solver.Add(sum(var for row,var in zip(states,variables) if row[h]>=t)
                               <=int(ca>=t)+(int(cb>=t)-int(ca>=t))*theta[h])
        self.assertEqual(solver.Solve(),pywraplp.Solver.OPTIMAL)
        self.assertAlmostEqual(result['objective'],objective.Value(),places=6)
        self.assertTrue(result['audit']['verified'])

    def test_finite_pool_master_recovers_known_optimum_five(self):
        source = dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        groups,shapes = structure(source)
        configurations = [[x%d for d in shapes[0]] for x in (1,2)]
        configurations.append([1 if j<7 else 2 for j,d in enumerate(shapes[0])])
        pool = dict(group_columns=[configurations])
        with redirect_stdout(io.StringIO()):
            result = solve(source,pool,seconds=10)
        self.assertTrue(result['audit']['verified'])
        w = result['coherent_primal_witness']
        self.assertAlmostEqual(w['objective_numerator']/w['denominator'],5,places=6)
        self.assertAlmostEqual(result['history'][-1]['pool_lower'],5,places=6)

    def test_group_pricing_adds_missing_balanced_assignments(self):
        source = dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        groups,shapes = structure(source)
        configurations = [[x%d for d in shapes[0]] for x in (1,2)]
        with redirect_stdout(io.StringIO()):
            result = solve(source,dict(group_columns=[configurations]),seconds=10,
                           pricing_seconds=.2,pricing_every=1)
        witness = result['coherent_primal_witness']
        self.assertAlmostEqual(witness['objective_numerator']/witness['denominator'],5,places=6)
        self.assertGreater(len(result['group_columns'][0]),2)
        self.assertTrue(result['audit']['verified'])
        self.assertAlmostEqual(result['pricing_history'][-1]['lower'],5,places=6)

    def test_pointwise_objective_rejects_cancellation_between_fibres(self):
        source = dict(core_residual=[1,2],high_primes=[3,5],high_pattern=[0,1])
        groups,shapes = structure(source)
        pool = dict(group_columns=[[[1%d for d in shape]] for shape in shapes])
        with redirect_stdout(io.StringIO()):
            result = solve(source,pool,seconds=10,passes=5,pointwise=True)
        self.assertEqual(result['objective_kind'],'positive_point_means')
        self.assertAlmostEqual(result['objective_value'],3,places=6)
        self.assertAlmostEqual(result['history'][-1]['pool_lower'],3,places=6)
        self.assertTrue(result['audit']['negative_cost'])
        self.assertFalse(result['audit']['all_point_means_nonpositive'])

    def test_level_stabilization_preserves_a_zero_small_optimum(self):
        source = dict(core_residual=[1,2],high_primes=[3,5],high_pattern=[0,1])
        groups,shapes = structure(source)
        pool = dict(group_columns=[[[x%d for d in shape] for x in (1,2)] for shape in shapes])
        with redirect_stdout(io.StringIO()):
            result = solve(source,pool,seconds=10,passes=5,pointwise=True,level_fraction=.5,compact_states=4096)
        self.assertTrue(result['audit']['verified'])
        self.assertLess(result['objective_value'],.001)
        self.assertAlmostEqual(result['history'][-1]['pool_lower'],0,places=6)
        self.assertTrue(any(row['level_stabilized'] for row in result['history']))
        self.assertTrue(all(row['compact_points_solved']==2 for row in result['history']))

    def test_weighted_global_pricing_matches_independent_state_enumeration(self):
        source = dict(core_residual=[1,2],high_primes=[3,5],high_pattern=[0,1])
        groups,shapes = structure(source)
        model = LocalStates(source['high_primes'],source['high_pattern'])
        class Dual:
            def __init__(self,value):self.value=value
            def dual_value(self):return self.value
        rows = []
        for x,weight in enumerate((.25,.75)):
            multipliers = [[0]+[-(h+t+x+1)*1000000 for t in range(1,len(shape)+1)] for h,shape in enumerate(shapes)]
            rows.append((x,dict(constant=0,multipliers=multipliers,denominator=1000000),Dual(weight)))
        result = group_prices(source,model,shapes,rows,[-16,-16],.5,pointwise=True)
        self.assertEqual(result['point_weights'],[250000,750000])
        for x,v in enumerate(source['core_residual']):
            D = [q-1-(b==v%2) for q,b in zip(source['high_primes'],source['high_pattern'])]
            values = []
            for state in itertools.product(range(D[0]+1),range(D[1]+1),range(17)):
                cost = (D[0]-state[0])*(D[1]-state[1])-state[2]
                penalty = -sum(sum(result['event_multipliers'][h][x][1:count+1]) for h,count in enumerate(state))
                values.append(cost*1000*result['point_weights'][x]+1000000*penalty)
            self.assertEqual(result['point_multipliers'][x],min(values))
        for h,row in enumerate(result['group_pricing']):
            counts = [(k,14-k) for k in range(15)] if h<2 else [(k+1,16-k) for k in range(16)]
            exact = max(sum(e['weight'] for e in row['events'] if pair[e['point']]>=e['threshold']) for pair in counts)
            self.assertEqual(row['upper_bound'],exact)


if __name__=='__main__':
    unittest.main()
