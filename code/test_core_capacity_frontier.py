import itertools
import json
from pathlib import Path
import unittest
from core_capacity_frontier import intersect_frontiers
from weighted_residual_bound import solve_weighted


class CoreCapacityTests(unittest.TestCase):
    def test_frontier_intersection_as_upward_sets(self):
        first = [[0,6],[2,3],[5,0]]
        second = [[1,5],[3,2],[6,0]]
        result = intersect_frontiers(first,second)
        for a,b in itertools.product(range(9),repeat=2):
            expected = all(any(x<=a and y<=b for x,y in f) for f in (first,second))
            self.assertEqual(any(x<=a and y<=b for x,y in result),expected)

    def test_saved_actual_families_respect_new_frontiers(self):
        root = Path(__file__).resolve().parents[1]
        data = json.loads((root/'data/core-capacity-frontiers.json').read_text())
        families = json.loads((root/'data/min-uncovered-2310-five-patterns-optimal.json').read_text())
        families.append(json.loads((root/'data/min-uncovered-30030-witness.json').read_text()))
        for source in families:
            period = source['N']
            rows = data['at2310' if period==2310 else 'at30030']
            family = [(r['modulus'],r['residue']) for r in source['system']]
            counts = [sum(not any(n%d==a for d,a in family) for n in range(e,period,2)) for e in (0,1)]
            for row in rows:
                if row['pattern']==source['pattern']:
                    self.assertTrue(any(a<=counts[0] and b<=counts[1] for a,b in row['frontier']))

    def test_weighted_solver_against_all_small_congruence_placements(self):
        for pattern in ((0,0),(0,1)):
            anchors = [(3,0),(5,0),(6,4),(10,6 if pattern[1]==0 else 1)]
            possible = []
            for odd,even in itertools.product(range(15),range(30)):
                family = anchors+[(15,odd),(30,even)]
                possible.append([sum(not any(n%d==a for d,a in family) for n in range(e,30,2)) for e in (0,1)])
            for weights in ((1,2),(2,1)):
                expected = min(sum(c*w for c,w in zip(pair,weights)) for pair in possible)
                result = solve_weighted(30,pattern,weights,seconds=10,workers=1)
                self.assertEqual(result['status'],'OPTIMAL')
                self.assertEqual(result['weighted_lower'],expected)


if __name__=='__main__':
    unittest.main()
