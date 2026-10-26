"""Exact small-instance checks of the lower-certificate construction."""
import contextlib
from fractions import Fraction as F
import io
from itertools import product
import unittest

from certify_stoploss_half_barrier import initial_cases,certify,SCALE
from audit_stoploss_half_barrier import enumerate_prefixes


def law(rs):
    w=[F(1)]
    for r in rs:
        nxt=[F(0)]*(len(w)+1)
        for n,p in enumerate(w):nxt[n]+=p*(1-r);nxt[n+1]+=p*r
        w=nxt
    return w


def advance(w,p,d):
    charge=sum(v*max(F(0),(F(2**n,p)-d)/(1-d)) for n,v in enumerate(w))
    r=1/(p*(1-d));nxt=[F(0)]*(len(w)+1)
    for n,v in enumerate(w):nxt[n]+=v*(1-r);nxt[n+1]+=v*r
    return charge,nxt


class BarrierChecks(unittest.TestCase):
    def test_independent_prefix_enumeration(self):
        a=initial_cases();b=enumerate_prefixes()
        self.assertEqual(len(a),60)
        for x,y in zip(a,b):
            self.assertEqual(x['deltas'],y['deltas'])
            self.assertEqual(x['prefix_loss'],y['prefix_loss'])
            self.assertEqual(x['seed_numerators'],y['numerators'])
            self.assertEqual(x['seed_denominator'],y['denominator'])

    def test_majorization_convex_comparison_exactly(self):
        ps=[13,17,19]
        for fractions in product((F(0),F(1,2),F(1)),repeat=3):
            rs=[(1+t)/p for p,t in zip(ps,fractions)]
            excess=sum(t/p for t,p in zip(fractions,ps))
            greedy=[]
            for p in ps:
                add=min(excess,F(1,p));greedy.append(F(1,p)+add);excess-=add
            self.assertEqual(excess,0)
            actual=law(rs);lower=law(greedy)
            for threshold in range(9):
                a=sum(w*max(0,2**n-threshold) for n,w in enumerate(actual))
                b=sum(w*max(0,2**n-threshold) for n,w in enumerate(lower))
                self.assertLessEqual(b,a)

    def test_all_small_horizon_cases_against_exact_optima(self):
        with contextlib.redirect_stdout(io.StringIO()):certificate=certify(end=17,H=100,states=10)
        actions13=[F(0),F(1,13),F(2,13),F(4,13),F(1,2)]
        actions17=[F(0),F(1,17),F(2,17),F(4,17),F(8,17),F(1,2)]
        for case in certificate['cases']:
            w=[F(x,case['seed_denominator']) for x in case['seed_numerators']]
            minimum=None
            for d,e in product(actions13,actions17):
                a,nextw=advance(w,13,d);b,_=advance(nextw,17,e)
                value=F(case['prefix_loss'])+a+b
                minimum=value if minimum is None else min(minimum,value)
            self.assertLessEqual(F(case['total_lower'],SCALE),minimum)


if __name__=='__main__':unittest.main()
