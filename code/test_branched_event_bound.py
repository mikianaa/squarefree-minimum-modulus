import copy
import itertools
import random
import unittest
from branched_event_bound import bound
from audit_branched_event_bound import audit_bound


class BranchedEventTests(unittest.TestCase):
    def test_full_split_matches_enumeration_and_partial_split_dominates(self):
        points=[0,1,4];rng=random.Random(23)
        for shapes in ([2,3],[1,2,3]):
            for repeat in range(6):
                events=[dict(point=x,threshold=t,weight=rng.randrange(10)) for x in range(3) for t in range(1,len(shapes)+1)]
                exact=max(sum(e['weight'] for e in events if sum(points[e['point']]%d==a for d,a in zip(shapes,residues))>=e['threshold'])
                          for residues in itertools.product(*(range(d) for d in shapes)))
                partial=bound(points,shapes,events,1)
                self.assertGreaterEqual(audit_bound(points,shapes,events,partial),exact)
                full=bound(points,shapes,events,len(shapes))
                self.assertEqual(audit_bound(points,shapes,events,full),exact)

    def test_missing_case_and_false_bound_are_rejected(self):
        points=[0,1,4];shapes=[2,3];events=[dict(point=0,threshold=2,weight=9)]
        certificate=bound(points,shapes,events,2)
        bad=copy.deepcopy(certificate);bad['nodes'][0]['children'].pop()
        with self.assertRaises(AssertionError):audit_bound(points,shapes,events,bad)
        bad=copy.deepcopy(certificate);bad['bound']-=1
        with self.assertRaises(AssertionError):audit_bound(points,shapes,events,bad)


if __name__=='__main__':unittest.main()
