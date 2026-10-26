import copy
import itertools
import random
import unittest
from pair_event_majorant import upper_bound
from event_score_majorant import optimized_upper_bound
from audit_weighted_coherent_dual import capacity
from branched_event_bound import bound,shifted


class PairEventMajorantTests(unittest.TestCase):
    def test_pair_multiplier_proves_a_strict_integer_improvement(self):
        points=[1,2];shapes=[1,2,3];events=[(0,3),(1,2)];weights=[2,2]
        old=optimized_upper_bound(points,shapes,events,weights)
        new=upper_bound(points,shapes,events,weights)
        listed=[dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,weights)]
        self.assertEqual(old['bound'],3);self.assertEqual(new['bound'],2)
        self.assertEqual(capacity(points,shapes,listed,new),2)
        for field,value in (('thresholds',[2,2]),('multiplier',[-2,1])):
            bad=copy.deepcopy(new);bad['pair_cuts'][0][field]=value
            with self.assertRaises(AssertionError):capacity(points,shapes,listed,bad)
        bad=copy.deepcopy(new);bad['point_intercepts'][0]=[-1,1]
        with self.assertRaises(AssertionError):capacity(points,shapes,listed,bad)

    def test_random_scores_dominate_exhaustive_assignments_and_shifted_cases(self):
        rng=random.Random(71);points=[1,2,4]
        for shapes in ([2,3,6],[1,2,3,6]):
            for fixed in ([],[[2,1]],[[d,1%d] for d in shapes]):
                for repeat in range(8):
                    events=[dict(point=x,threshold=t,weight=rng.randrange(8)) for x in range(len(points)) for t in range(1,len(shapes)+1)]
                    remaining,reduced,constant=shifted(points,shapes,events,fixed)
                    pairs=[(e['point'],e['threshold']) for e in reduced];weights=[e['weight'] for e in reduced]
                    cert=upper_bound(points,remaining,pairs,weights)
                    self.assertEqual(capacity(points,remaining,reduced,cert),cert['bound'])
                    self.assertLessEqual(cert['bound'],optimized_upper_bound(points,remaining,pairs,weights)['bound'])
                    exact=max(sum(e['weight'] for e in events if sum(points[e['point']]%d==a for d,a in fixed)
                                  +sum(points[e['point']]%d==a for d,a in zip(remaining,choice))>=e['threshold'])
                              for choice in itertools.product(*(range(d) for d in remaining)))
                    self.assertGreaterEqual(constant+cert['bound'],exact)

    def test_complete_pair_case_tree_is_verified_by_existing_tree_auditor(self):
        points=[1,2];shapes=[1,2,3];events=[dict(point=0,threshold=3,weight=2),dict(point=1,threshold=2,weight=2)]
        cert=bound(points,shapes,events,depth=1,pairs=True)
        self.assertEqual(capacity(points,shapes,events,cert),2)
        self.assertEqual(cert['nodes'][0]['affine']['method'],'pair_affine')

    def test_selected_case_wrapper_rejects_changed_constant(self):
        from certify_pair_case_probe import certify
        from audit_pair_case_probe import audit
        source=dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        events=[dict(point=0,threshold=1,weight=2),dict(point=0,threshold=10,weight=2),dict(point=1,threshold=9,weight=2)]
        data=dict(input=source,result=dict(pricing_history=[dict(group_pricing=[dict(events=events)])]))
        probe=dict(input=source,comparisons=[dict(axis=0,fixed=[[3,1]])])
        proof=certify(data,probe)
        self.assertTrue(audit(proof)['verified'])
        self.assertFalse(audit(proof)['complete_group_coverage'])
        proof['cases'][0]['constant']+=1
        with self.assertRaises(AssertionError):audit(proof)


if __name__=='__main__':unittest.main()
