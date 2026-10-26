"""Independent residue checks for the core profiles and heuristic witnesses."""
import itertools
import json
import math
from pathlib import Path
import unittest

from core_residual_geometry import (points_for,profile,dominates_for_cover,
    easiest_profiles,optimistic_profiles,minimal_triple_residuals,triple_groups)
from extend_residual_witness import extend
from improve_cover_witness import improve
from block_profile_capacity import weighted_obstruction


ROOT = Path(__file__).resolve().parents[1]


class GeometryTests(unittest.TestCase):
    def test_intersection_capacities_using_ordinary_residues(self):
        for tail in itertools.product((0,1),repeat=2):
            points = points_for((0,)+tail)
            # A nonrectangular subset, separate from any generated core family.
            selected = [j for j,(e,a) in enumerate(points) if (sum(a)+e)%3 != 0]
            mask = sum(1<<j for j in selected)
            residues = [[],[]]
            for j in selected:
                e,a = points[j]
                n = next(n for n in range(e,210,2) if tuple(n%p for p in (3,5,7))==a)
                residues[e].append(n)
            expected = tuple(map(len,residues))+tuple(
                max(sum(n%d==a for n in residues[e]) for a in range(d))
                for d in (3,5,7,15,21,35,105) for e in (0,1))
            self.assertEqual(profile(mask,points),expected)

    def test_optimistic_ratios_dominate_actual_triple_removals(self):
        points = points_for((0,1,1))
        mask = sum(1<<j for j,(e,a) in enumerate(points) if sum(a)%3 != 0)
        optimistic = optimistic_profiles(points,[mask])
        actuals = set(minimal_triple_residuals(mask,triple_groups(points)))
        self.assertTrue(actuals)
        for left in actuals:
            self.assertTrue(any(dominates_for_cover(p,profile(left,points)) for p in optimistic))

    def test_zero_layer_dominance_and_triple_edge_case(self):
        zero = (0,0)+(0,)*14
        one = (1,1)+(1,)*14
        self.assertTrue(dominates_for_cover(zero,one))
        self.assertFalse(dominates_for_cover(one,zero))
        self.assertEqual(easiest_profiles([zero,one]),[zero])
        self.assertEqual(list(minimal_triple_residuals(1,[1])),[0])

    def test_saved_profile_audit_and_all_anchor_cases(self):
        data = json.loads((ROOT/'data/core-geometry-profile-audit.json').read_text())
        self.assertEqual([r['easiest_count'] for r in data],[60,66,115,126])
        for row in data:
            self.assertTrue(row['profile_sets_equal'] and row['every_profile_dominated'])
            self.assertEqual(easiest_profiles(map(tuple,row['easiest_profiles'])),list(map(tuple,row['easiest_profiles'])))

    def test_weighted_bound_against_direct_modulus_enumeration(self):
        rows = json.loads((ROOT/'data/core-geometry-profile-audit.json').read_text())
        core_ds = (3,5,7,15,21,35,105)
        for row in rows:
            p = row['easiest_profiles'][0]
            result = weighted_obstruction(p,(11,13),(0,1))
            def gap(weights):
                demand = sum(weights[e]*p[e]*(10-(e==0))*(12-(e==1)) for e in (0,1))
                capacity = 0
                for d in range(3,30031):
                    if 30030%d or 210%d==0 or d in (11,22,13,26):
                        continue
                    core = math.gcd(d,210)
                    if core%2==0:
                        core //= 2
                    pair = [p[e] if core==1 else p[2+2*core_ds.index(core)+e] for e in (0,1)]
                    amounts = [weights[e]*pair[e]*(1 if d%11==0 else 10-(e==0))
                               *(1 if d%13==0 else 12-(e==1)) for e in (0,1)]
                    capacity += max(amounts) if d%2==0 else sum(amounts)
                return demand-capacity
            self.assertEqual(result['gap'],gap(result['weights']))
            for weights in itertools.product(range(4),repeat=2):
                if sum(weights):
                    self.assertGreaterEqual(result['gap']*sum(weights),gap(weights)*sum(result['weights']))


class WitnessTests(unittest.TestCase):
    def test_saved_30030_witness(self):
        source = json.loads((ROOT/'data/min-uncovered-30030-witness.json').read_text())
        family = [(v['modulus'],v['residue']) for v in source['system']]
        self.assertEqual(len({d for d,a in family}),62)
        self.assertTrue(all(d>=3 and 30030%d==0 for d,a in family))
        missed = sum(not any(n%d==a for d,a in family) for n in range(30030))
        self.assertEqual(missed,792)
        self.assertEqual(missed,source['uncovered'])

    def test_extension_and_improvement_full_period(self):
        # Small control with known nonempty residual. Both helpers also execute
        # a complete, ordinary-integer coverage check on the emitted family.
        family = [(3,0),(5,0),(6,4),(10,1),(15,2),(30,7)]
        out = extend(family,30,7,1,seconds=0,seed=7)
        out['pattern'] = [0,1,1]
        self.assertTrue(out['verified_full_period'])
        improved = improve(out,seconds=0,seed=7)
        self.assertLessEqual(improved['uncovered'],out['uncovered'])
        self.assertTrue(improved['verified_full_period'])
        for v in improved['system']:
            self.assertTrue(all(v['modulus']%(p*p) for p in (2,3,5,7)))


if __name__ == '__main__':
    unittest.main()
