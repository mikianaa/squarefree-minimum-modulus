import copy
import itertools
import random
import unittest
from clique_case_constraints import maximum_clique, strongest, screen
from audit_case_cliques import audit
from clique_event_majorant import upper_bound
from pair_event_majorant import upper_bound as pair_bound
from audit_weighted_coherent_dual import capacity
from fractional_group_scenarios import maximize
from branched_event_bound import bound, shifted
from repair_clique_scenarios import repair_pool


class CliqueCaseTests(unittest.TestCase):
    def test_weighted_clique_matches_exhaustive_subsets(self):
        rng=random.Random(381)
        for n in range(1,8):
            for repeat in range(10):
                caps=[[2]*n for j in range(n)]
                for x,y in itertools.combinations(range(n),2):
                    caps[x][y]=caps[y][x]=rng.randrange(1,3)
                weights=[rng.randrange(8) for j in range(n)]
                expected=max(sum(weights[j] for j in subset)
                             for k in range(n+1) for subset in itertools.combinations(range(n),k)
                             if all(caps[x][y]<2 for x,y in itertools.combinations(subset,2)))
                score,subset=maximum_clique([(j,1) for j in range(n)],weights,caps)
                self.assertEqual(score,expected)
                self.assertEqual(sum(weights[j] for j in subset),score)

    def test_tiny_nonnegative_tails_do_not_become_negative(self):
        score,events=strongest([0,1],[2,3],[[.999999999999,1e-12,1e-30],[1,0,0]])
        self.assertGreaterEqual(score,0)

    def test_clique_strictly_strengthens_pair_lp_and_has_exact_majorant(self):
        points=[0,1,2];shapes=[2,3,6];events=[(x,3) for x in range(3)];weights=[6]*3
        old=pair_bound(points,shapes,events,weights)
        new=upper_bound(points,shapes,events,weights)
        self.assertGreater(old['bound'],6)
        self.assertEqual(new['bound'],6)
        self.assertEqual(new['method'],'clique_affine')
        listed=[dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,weights)]
        self.assertEqual(capacity(points,shapes,listed,new),6)
        bad=copy.deepcopy(new);bad['clique_cuts'][0]['events'][0]['threshold']=1
        with self.assertRaises(AssertionError):capacity(points,shapes,listed,bad)
        bad=copy.deepcopy(new);bad['clique_cuts'][0]['multiplier']=[-1,1]
        with self.assertRaises(AssertionError):capacity(points,shapes,listed,bad)

    def test_random_shifted_bounds_cover_all_actual_assignments(self):
        rng=random.Random(52);points=[0,1,2,4]
        for shapes in ([2,3,6],[1,2,3,6]):
            for fixed in ([],[[2,1]],[[d,0] for d in shapes]):
                for repeat in range(4):
                    events=[dict(point=x,threshold=t,weight=rng.randrange(5)) for x in range(len(points)) for t in range(1,len(shapes)+1)]
                    remaining,reduced,constant=shifted(points,shapes,events,fixed)
                    certificate=upper_bound(points,remaining,[(e['point'],e['threshold']) for e in reduced],[e['weight'] for e in reduced])
                    self.assertEqual(capacity(points,remaining,reduced,certificate),certificate['bound'])
                    actual=max(sum(e['weight'] for e in events if sum(points[e['point']]%d==a for d,a in fixed)
                                   +sum(points[e['point']]%d==a for d,a in zip(remaining,choice))>=e['threshold'])
                               for choice in itertools.product(*(range(d) for d in remaining)))
                    self.assertGreaterEqual(constant+certificate['bound'],actual)
        cert=bound(points,[2,3,6],events=[dict(point=x,threshold=3,weight=6) for x in range(3)],depth=1,cliques=True)
        self.assertEqual(capacity(points,[2,3,6],[dict(point=x,threshold=3,weight=6) for x in range(3)],cert),cert['bound'])

    def test_exact_case_cut_rejects_false_probability_and_compatible_pair(self):
        shape=[d for d in range(3,211) if 210%d==0]
        state=dict(kind='fractional',fixed=[],count_numerators=[[1]+[0]*13+[1] for x in range(3)])
        data=dict(input=dict(core_residual=[0,1,2],high_primes=[11]),probability_denominator=2,
                  group_mixtures=[dict(axis=0,states=[state])])
        report=screen(data);checked=audit(data,report)
        self.assertTrue(checked['verified']);self.assertEqual(checked['conditional_clique_cuts_checked'],1)
        self.assertFalse(checked['covering_obstruction'])
        bad=copy.deepcopy(report);bad['cuts'][0]['probability_sum_numerator']+=1
        with self.assertRaises(AssertionError):audit(data,bad)
        bad=copy.deepcopy(report)
        for e in bad['cuts'][0]['events']:e['threshold']=1
        with self.assertRaises(AssertionError):audit(data,bad)

    def test_repair_preserves_means_and_prior_pair_feasibility(self):
        from repair_pair_scenarios import maximum_pair_score
        points=[0,1,2];shape=[d for d in range(3,211) if 210%d==0]
        rewards=[[0]*14+[6] for p in points]
        scenario=maximize(points,shape,rewards,pair_compatibility=True)
        data=dict(input=dict(core_residual=points,high_primes=[11]),baseline_event_multipliers=[[[0]*14+[-6] for p in points]],
                  group_scenarios=[dict(axis=0,scenarios=[scenario])])
        repaired=repair_pool(data)['group_scenarios'][0]['scenarios'][0]
        self.assertLessEqual(strongest(points,shape,repaired['count_distributions'])[0],1+1e-9)
        self.assertLessEqual(maximum_pair_score(points,shape,[],repaired['count_distributions']),1+1e-9)
        for old,new in zip(scenario['count_distributions'],repaired['count_distributions']):
            self.assertAlmostEqual(sum(new),1)
            self.assertLessEqual(sum(k*v for k,v in enumerate(new)),sum(k*v for k,v in enumerate(old))+1e-9)
        self.assertEqual(scenario['resource_distributions'],repaired['resource_distributions'])


if __name__=='__main__':unittest.main()
