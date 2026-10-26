import itertools
import math
import random
import unittest
import json
from pathlib import Path
import numpy as np
from conditional_tail_lp import price_point,build
from audit_fibre_concentration import gray_capacities


class ConditionalTailTests(unittest.TestCase):
    def test_abstract_count_model_has_the_same_local_state_costs(self):
        actual = json.loads(Path('data/core-fibre-eight-plane.json').read_text())
        H = json.loads(Path('data/core-hit-concentration-7-15.json').read_text())
        abstract = dict(abstract_counts=[7,15],high_primes=actual['high_primes'],high_pattern=actual['high_pattern'],choices=[])
        first = build(actual,H)
        second = build(abstract,dict(two_layer_capacities=H['two_layer_capacities']))
        self.assertEqual(first[1:3],second[1:3])
        self.assertEqual(second[3],[])
        for index in (4,5,6):
            for a,b in zip(first[index],second[index]):
                np.testing.assert_array_equal(a,b)
        for a,b in zip(first[7],second[7]):
            for c,d in zip(a,b):
                np.testing.assert_array_equal(c,d)

    def test_extra_unit_shapes_capacity_identity(self):
        R = [1,2,4,5,8]
        shapes = [3,5,6,10]
        old = gray_capacities(R,shapes)
        new = gray_capacities(R,[1,2]+shapes)
        for a,row in enumerate(old):
            for b,value in enumerate(row):
                self.assertEqual(new[a][b],value+a+b+max(a,b))

    def test_pricing_separability_against_all_local_states(self):
        grid = np.array(list(itertools.product((0,1),repeat=4)),dtype=np.int64)
        free = 1-grid
        cost = np.prod(free,axis=1)
        volumes = [free[:,2]*free[:,3],free[:,0]]
        links = [-2,-1,-3,-2,-1,-2]
        common = [[np.array([0,-i-1]) for i in range(4)]+
                  [np.array([-max(0,n-3)*2 for n in range(17)]),
                   np.array([-max(0,n-4)*3 for n in range(17)])] for e in (0,1)]
        actual,state,_ = price_point(0,links,common,[grid,grid],[cost,cost],[volumes,volumes],3)
        expected = min(3*(int(cost[j])-h0*int(volumes[0][j])-h1*int(volumes[1][j]))-
                       sum(links[i]*n+int(common[0][i][n]) for i,n in enumerate(tuple(ns)+(h0,h1)))
                       for j,ns in enumerate(grid) for h0,h1 in itertools.product(range(17),repeat=2))
        self.assertEqual(actual,expected)
        self.assertEqual(len(state),6)

    def test_conditional_capacity_with_colliding_plane_labels(self):
        rng = random.Random(81025)
        R = [1,2,4,5]
        Q = [3,5,7]
        groups = [S for k in (2,3) for S in itertools.combinations(range(3),k)]
        for trial in range(80):
            planes = [(i,d,rng.randrange(d),rng.randrange(Q[i])) for i in range(3) for d in (2,3,6)]
            tails = [(S,d,rng.randrange(d),tuple(rng.randrange(Q[i]) for i in S)) for S in groups for d in (1,2,3,6)]
            lower = actual = augmented = 0
            for x in R:
                alphabets = [set(range(1+(x%2==1),q)) for q in Q]
                deleted = [{label for i,d,a,label in planes if i==axis and x%d==a and label in alphabets[axis]} for axis in range(3)]
                hits = [sum(i==axis and x%d==a for i,d,a,label in planes) for axis in range(3)]
                augmented_deleted = [set(s) for s in deleted]
                for axis in range(3):
                    wanted = min(len(alphabets[axis]),hits[axis])
                    for label in sorted(alphabets[axis]-augmented_deleted[axis]):
                        if len(augmented_deleted[axis])==wanted:
                            break
                        augmented_deleted[axis].add(label)
                free = [len(A)-len(B) for A,B in zip(alphabets,augmented_deleted)]
                lower += math.prod(free)
                lower -= sum(math.prod(free[i] for i in range(3) if i not in S) for S,d,a,labels in tails if x%d==a)
                for point in itertools.product(*alphabets):
                    tail_hit = any(x%d==a and all(point[i]==label for i,label in zip(S,labels)) for S,d,a,labels in tails)
                    actual += not tail_hit and not any(point[i] in deleted[i] for i in range(3))
                    augmented += not tail_hit and not any(point[i] in augmented_deleted[i] for i in range(3))
            self.assertLessEqual(lower,augmented)
            self.assertLessEqual(augmented,actual)


if __name__=='__main__':
    unittest.main()
