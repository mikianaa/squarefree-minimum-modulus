import itertools
import unittest
from fractional_group_scenarios import maximize
from repair_pair_scenarios import repair,maximum_pair_score


class PairScenarioRepairTests(unittest.TestCase):
    def test_repairs_fractional_violation_without_increasing_mean_demands(self):
        points=[1,2];shapes=[1,2,3];rewards=[[0,0,0,2],[0,0,2,2]]
        old=maximize(points,shapes,rewards);new,violation=repair(points,shapes,rewards,old)
        self.assertGreater(violation,1)
        self.assertLessEqual(maximum_pair_score(points,shapes,[],new['count_distributions']),1+1e-10)
        for a,b in zip(old['count_distributions'],new['count_distributions']):
            self.assertAlmostEqual(sum(b),1)
            self.assertLessEqual(sum(k*p for k,p in enumerate(b)),sum(k*p for k,p in enumerate(a)))
        self.assertLessEqual(sum(new['payoffs']),sum(old['payoffs']))

    def test_every_ordinary_assignment_remains_unmodified(self):
        points=[1,2,4];shapes=[1,2,3];rewards=[list(range(4)) for p in points]
        for fixed in ([],[[2,1]],[[1,0],[2,1],[3,1]]):
            remaining=[d for d in shapes if d not in dict(fixed)];base=int(1 in remaining)
            for choice in itertools.product(*(range(d) for d in remaining)):
                rows=[]
                for p in points:
                    hits=sum(p%d==a for d,a in zip(remaining,choice));row=[0.]*(len(remaining)-base+1);row[hits-base]=1.;rows.append(row)
                result,old=repair(points,shapes,rewards,dict(fixed=fixed,count_distributions=rows))
                self.assertEqual(result['pair_repair_factor'],1)
                self.assertEqual(result['count_distributions'],rows)


if __name__=='__main__':unittest.main()
