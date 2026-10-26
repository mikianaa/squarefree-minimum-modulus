"""Exact rational checks for the new full-range reductions."""
from fractions import Fraction as F
from itertools import product
import unittest

from certify_fullrange_action_pruning import initial_atoms,future_marginal_bounds
from certify_fullrange_stoploss_barrier import mixture_grid,MIX,prefix_cases
from audit_fullrange_stoploss_barrier import initial_cases as audit_cases,padded_mixtures,Q


def distribution(rs):
    w=[F(1)]
    for r in rs:
        out=[F(0)]*(len(w)+1)
        for n,v in enumerate(w):out[n]+=v*(1-r);out[n+1]+=v*r
        w=out
    return w


def fixed_future_cost(n,primes,atoms):
    w=[F(0)]*n+[F(1)];total=F(0)
    for p,a in zip(primes,atoms):
        total+=sum(v*F(max(0,2**k-a),p-a) for k,v in enumerate(w))
        out=[F(0)]*(len(w)+1)
        for k,v in enumerate(w):out[k]+=v*F(p-a-1,p-a);out[k+1]+=v/F(p-a)
        w=out
    return total


class FullRangeChecks(unittest.TestCase):
    def test_large_delta_after_last_atom_is_dominated(self):
        for p in (5,7,11,17,31):
            a=initial_atoms(p)[-1];d=F(a,p)
            for e in (d+(1-d)/3,d+2*(1-d)/3):
                for n in range(8):
                    s=F(2**n-int(p==5),p)
                    self.assertLessEqual(max(0,(s-d)/(1-d)),max(0,(s-e)/(1-e)))
                self.assertLessEqual(F(1,p)/(1-d),min(F(1),F(1,p)/(1-e)))

    def test_marginal_future_lower_against_every_small_schedule(self):
        ps=[5,7,11,13,17];allowed=[initial_atoms(p) for p in ps]
        S=10**6;table=future_marginal_bounds(ps,allowed,8,S)
        for atoms in product(allowed[3],allowed[4]):
            for n in range(6):
                true=fixed_future_cost(n+1,ps[3:],atoms)-fixed_future_cost(n,ps[3:],atoms)
                self.assertLessEqual(F(table[3][n],S),true)

    def test_sorted_unequal_caps_majorize_and_lower_convex_cost(self):
        ps=[13,17,19];lows=[F(1,p) for p in ps]
        original_upper=[F(1,5),F(1),F(1,3)];upper=sorted(original_upper,reverse=True)
        for choices in product((F(0),F(1,2),F(1)),repeat=3):
            actual=[l+(u-l)*t for l,u,t in zip(lows,original_upper,choices)]
            excess=sum(actual)-sum(lows);greedy=[]
            for l,u in zip(lows,upper):
                add=min(excess,u-l);greedy.append(l+add);excess-=add
            self.assertEqual(excess,0)
            for j in range(1,4):self.assertGreaterEqual(sum(greedy[:j]),sum(sorted(actual,reverse=True)[:j]))
            before=distribution(actual);after=distribution(greedy)
            for t in (0,1,2,3,4,5,7):
                self.assertLessEqual(sum(w*max(0,2**n-t) for n,w in enumerate(after)),
                                     sum(w*max(0,2**n-t) for n,w in enumerate(before)))

    def test_two_grid_enclosures_against_exact_greedy_law(self):
        ps=[13,17,19];ds=[1,3,5];H=100
        caps=[F(0)]
        for p,d in zip(ps,ds):caps.append(caps[-1]+F(1,d)-F(1,p))
        width=int(caps[-1]*H)+1
        ts,left,right=mixture_grid(ps,ds,width,H)
        # The independent audit uses its own fixed H=5000, tested separately.
        audit_width=int(caps[-1]*5000)+1
        ats,af=padded_mixtures(ps,ds,audit_width)
        for j in range(width):
            x=F(j,H);t=int(ts[j]);f=(x-caps[t])/(F(1,ds[t])-F(1,ps[t])) if t<3 else F(0)
            self.assertLessEqual(F(int(left[j]),MIX),1-f)
            self.assertLessEqual(F(int(right[j]),MIX),f)
            idx=j*50;t2=int(ats[idx]);fraction=F(int(af[idx]),Q)
            approx=[F(1,ds[i]) if i<t2 else F(1,ps[i])+(F(1,ds[i])-F(1,ps[i]))*fraction if i==t2 else F(1,ps[i]) for i in range(3)]
            true=[F(1,ds[i]) if i<t else F(1,ps[i])+(F(1,ds[i])-F(1,ps[i]))*f if i==t else F(1,ps[i]) for i in range(3)]
            self.assertTrue(all(a<=b for a,b in zip(approx,true)))

    def test_prefix_enumerations_agree(self):
        choices={5:[0,1],7:[0,1,2,4],11:[0,1,2,4]}
        a=prefix_cases(choices);b=audit_cases(choices)
        self.assertEqual(len(a),32)
        for x,y in zip(a,b):
            self.assertEqual(x['deltas'],y['deltas'])
            self.assertEqual(x['prefix_loss'],y['prefix_loss'])
            self.assertEqual(x['seed_numerators'],y['numerators'])
            self.assertEqual(x['seed_denominator'],y['denominator'])


if __name__=='__main__':unittest.main()
