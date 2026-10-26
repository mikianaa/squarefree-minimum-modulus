import itertools
import random
import unittest
from fractional_group_scenarios import maximize
from branched_event_bound import shifted
from event_score_majorant import optimized_upper_bound


class FractionalGroupScenarioTests(unittest.TestCase):
    def test_pair_compatibility_removes_an_unrealizable_fractional_optimum(self):
        points=[1,2];shapes=[1,2,3];rewards=[[0,0,0,2],[0,0,2,2]]
        loose=maximize(points,shapes,rewards)
        tightened=maximize(points,shapes,rewards,pair_compatibility=True)
        self.assertAlmostEqual(loose['objective'],3)
        self.assertAlmostEqual(tightened['objective'],2)
        self.assertGreater(tightened['pair_constraints_added'],0)
        for fixed in ([],[[2,1]],[[1,0]],[[1,0],[2,1],[3,1]]):
            answer=maximize(points,shapes,rewards,fixed=fixed,pair_compatibility=True)
            remaining=[d for d in shapes if d not in dict(fixed)]
            exact=max(sum(row[sum(p%d==a for d,a in fixed)+sum(p%d==a for d,a in zip(remaining,choice))]
                          for p,row in zip(points,rewards)) for choice in itertools.product(*(range(d) for d in remaining)))
            self.assertGreaterEqual(answer['objective']+1e-8,exact)

    def test_matches_affine_dual_and_dominates_exhaustive_assignments(self):
        rng=random.Random(821)
        points=[1,2,4]
        for shapes in ([2,3,6],[1,2,3,6]):
            for fixed in ([],[[2,1]],[[d,1%d] for d in shapes]):
                for repeat in range(6):
                    rewards=[]
                    for x in points:
                        row=[0]
                        for d in shapes:row.append(row[-1]+rng.randrange(10))
                        rewards.append(row)
                    coefficients=[rng.randrange(1,5) for x in points]
                    result=maximize(points,shapes,rewards,coefficients,fixed)
                    events=[dict(point=x,threshold=k,weight=c*(row[k]-row[k-1]))
                            for x,(row,c) in enumerate(zip(rewards,coefficients)) for k in range(1,len(row))]
                    remaining,reduced,constant=shifted(points,shapes,events,fixed)
                    dual=optimized_upper_bound(points,remaining,[(e['point'],e['threshold']) for e in reduced],
                                               [e['weight'] for e in reduced])
                    # The exact dual rounds its final integer bound down.
                    self.assertLessEqual(abs(result['objective']-(constant+dual['bound'])),1.00001)
                    options=[range(d) for d in remaining]
                    exact=max(sum(c*row[sum(p%d==a for d,a in fixed)+sum(p%d==a for d,a in zip(remaining,choice))]
                                  for p,row,c in zip(points,rewards,coefficients)) for choice in itertools.product(*options))
                    self.assertGreaterEqual(result['objective']+1e-7,exact)
                    supplies={r['resource']:dict(zip(r['residues'],r['probabilities'])) for r in result['resource_distributions']}
                    for row in supplies.values():self.assertAlmostEqual(sum(row.values()),1)
                    for x,(p,row) in enumerate(zip(points,result['count_distributions'])):
                        self.assertAlmostEqual(sum(row),1)
                        self.assertGreaterEqual(min(row),-1e-8)
                        self.assertLessEqual(sum(j*v for j,v in enumerate(row)),sum(r.get(p%d,0) for d,r in supplies.items())+1e-7)
                        payoff=sum(v*rewards[x][result['fixed_hits'][x]+result['count_base']+j] for j,v in enumerate(row))
                        self.assertAlmostEqual(result['payoffs'][x],payoff)


if __name__=='__main__':unittest.main()
