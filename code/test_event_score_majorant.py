import itertools
import random
import unittest
from event_score_majorant import upper_bound,optimized_upper_bound
from audit_weighted_coherent_dual import capacity
import copy


class EventScoreMajorantTests(unittest.TestCase):
    def test_bound_dominates_every_actual_assignment_with_and_without_unit_resource(self):
        rng = random.Random(9)
        points = [1,2,4]
        for shapes in ([2,3,6],[1,2,3,6]):
            options = [sorted({x%d for x in points}) for d in shapes]
            counts = [[sum(x%d==a for d,a in zip(shapes,assignment)) for x in points]
                      for assignment in itertools.product(*options)]
            events = [(x,t) for x in range(len(points)) for t in range(1,len(shapes)+1)]
            for repeat in range(30):
                weights = [rng.randrange(10) for e in events]
                exact = max(sum(w for (x,t),w in zip(events,weights) if row[x]>=t) for row in counts)
                self.assertGreaterEqual(upper_bound(points,shapes,events,weights)['bound'],exact)
                improved=optimized_upper_bound(points,shapes,events,weights)
                self.assertGreaterEqual(improved['bound'],exact)
                self.assertLessEqual(improved['bound'],upper_bound(points,shapes,events,weights)['bound'])
                listed=[dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,weights)]
                self.assertEqual(capacity(points,shapes,listed,improved),improved['bound'])

    def test_free_intercept_can_improve_bound_and_is_independently_checked(self):
        points=[1,2,4];shapes=[2,3,6];events=[(x,t) for x in range(3) for t in range(1,4)]
        rng=random.Random(9);found=None
        for repeat in range(100):
            weights=[rng.randrange(10) for e in events]
            certificate=optimized_upper_bound(points,shapes,events,weights)
            if 'point_intercepts' in certificate:found=(weights,certificate);break
        self.assertIsNotNone(found)
        weights,certificate=found
        self.assertLess(certificate['bound'],upper_bound(points,shapes,events,weights)['bound'])
        bad=copy.deepcopy(certificate);bad['point_intercepts'][0]=[-1,1]
        with self.assertRaises(AssertionError):capacity(points,shapes,[dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,weights)],bad)


if __name__=='__main__':unittest.main()
