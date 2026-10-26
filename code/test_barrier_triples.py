import copy
import itertools
import random
import unittest
import numpy as np
from separate_barrier_triples import strongest, triple_frontier
from audit_barrier_triples import audit


class BarrierTripleTests(unittest.TestCase):
    def test_separation_matches_exhaustive_thresholds_including_empty_pair(self):
        rng = random.Random(917)
        points, shapes = [0, 1, 2], [2, 3, 4]
        for fixed in ([], [[2, 0]], [[3, 1]]):
            remaining = [d for d in shapes if d not in dict(fixed)]
            reachable = []
            for a in itertools.product(*(range(d) for d in remaining)):
                selected = dict(fixed) | dict(zip(remaining, a))
                reachable.append([sum(p % d == selected[d] for d in shapes) for p in points])
            for trial in range(10):
                rows = [[0] * (len(remaining) + 1) for p in points]
                for row in rows:
                    for j in range(8):
                        row[rng.randrange(len(row))] += 1
                best = 0
                for ts in itertools.product(range(1, 4), repeat=3):
                    if any(all(k >= t for k, t in zip(state, ts)) for state in reachable):
                        continue
                    score = sum(sum(v for k, v in enumerate(row)
                                    if k + sum(p % d == a for d, a in fixed) >= t)
                                for p, row, t in zip(points, rows, ts))
                    best = max(best, score - 16)
                cut = strongest(points, shapes, fixed, rows, 8, True)
                self.assertEqual(cut['violation_numerator'] if cut else 0, best)

    def test_frontier_matches_ordinary_residue_enumeration(self):
        for points in ([0, 1, 2], [0, 2, 3], [1, 4, 7]):
            for shapes in ([2, 3], [2, 3, 4], [3, 5, 6]):
                common = separate = 0
                pairs = [0, 0, 0]
                for d in shapes:
                    r = [p % d for p in points]
                    if len(set(r)) == 1:
                        common += 1
                    elif len(set(r)) == 3:
                        separate += 1
                    else:
                        pairs[next(j for j in range(3) if r.count(r[j]) == 1)] += 1
                N = len(shapes)
                actual = np.full((N + 1, N + 1), -1)
                for assignment in itertools.product(*(range(d) for d in shapes)):
                    k = [sum(p % d == a for d, a in zip(shapes, assignment)) for p in points]
                    actual[:k[0] + 1, :k[1] + 1] = np.maximum(actual[:k[0] + 1, :k[1] + 1], k[2])
                np.testing.assert_array_equal(triple_frontier(N, common, separate, tuple(pairs)), actual)

    def fixture(self):
        points = [0, 1, 2]
        shapes = [d for d in range(3, 211) if 210 % d == 0]
        fixed = [[d, 0] for d in shapes if d not in (5, 7, 35)]
        counts = [[0, 4, 0, 0], [0, 4, 0, 0], [3, 0, 1, 0]]
        state = dict(kind='fractional', fixed=fixed, count_numerators=counts)
        barrier = dict(input=dict(core_residual=points, high_primes=[11]),
                       probability_denominator=4, group_mixtures=[dict(axis=0, states=[state])])
        return points, shapes, fixed, counts, barrier

    def test_three_events_can_fail_without_pair_conflict(self):
        points, shapes, fixed, counts, barrier = self.fixture()
        # Remaining resources separate all three points. Their marginal means
        # (1,1,1/2) fit into capacity 3, but events (1,1,2) cannot coexist.
        for exact in (False, True):
            cut = strongest(points, shapes, fixed, counts, 4, exact)
            self.assertIsNotNone(cut)
            report = dict(denominator=4, cuts=[dict(axis=0, state_index=0, **cut)])
            checked = audit(barrier, report)
            self.assertTrue(checked['verified'])
            self.assertEqual(checked['largest_violation_numerator'], 1)
            self.assertFalse(checked['covering_obstruction'])
            bad = copy.deepcopy(report)
            bad['cuts'][0]['violation_numerator'] += 1
            with self.assertRaises(AssertionError):
                audit(barrier, bad)

    def test_actual_assignment_mixture_is_not_rejected(self):
        points, shapes, fixed, counts, barrier = self.fixture()
        remaining = [d for d in shapes if d not in dict(fixed)]
        hist = [[0] * 4 for p in points]
        assignments = list(itertools.product(*(range(d) for d in remaining)))
        for assignment in assignments:
            for x, p in enumerate(points):
                hist[x][sum(p % d == a for d, a in zip(remaining, assignment))] += 1
        self.assertIsNone(strongest(points, shapes, fixed, hist, len(assignments), True))


if __name__ == '__main__':
    unittest.main()
