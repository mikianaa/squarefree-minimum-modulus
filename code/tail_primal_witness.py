"""Repair an aggregate LP point distribution to exact rational feasibility.

A feasible nonpositive-cost distribution is a witness only for the relaxation,
not a congruence covering. Shrink rounded positive masses to satisfy every
nonnegative upper constraint, then fill each parity's remaining mass with
the zero-hit state.
"""
import math
from fractions import Fraction


def exact_witness(columns,parities,masses,tail_rows,high,pattern,rank_rows=()):
    dimensions = len(high)
    assert len(pattern)==dimensions
    group_count = (1<<dimensions)-1
    scale = 10**9
    support = [(parities[x],ns,max(0,math.floor(v.solution_value()*scale))) for x,ns,v in columns if v.solution_value()>0]
    support = [(e,ns,w) for e,ns,w in support if w]
    ratio = Fraction(1)
    for e,total in enumerate(masses):
        used = sum(w for parity,ns,w in support if parity==e)
        if used:
            ratio = min(ratio,Fraction(total*scale,used))
    for h,t0,t1,bound,_ in tail_rows:
        used = sum(w*max(0,ns[h]-(t0 if e==0 else t1)) for e,ns,w in support)
        if used:
            ratio = min(ratio,Fraction(bound*scale,used))
    for h,t0,t1,w0,w1,bound,_ in rank_rows:
        used = sum(w*(w0 if e==0 else w1)*int(ns[h]>=(t0 if e==0 else t1)) for e,ns,w in support)
        if used:
            ratio = min(ratio,Fraction(bound*scale,used))
    support = [(e,ns,w*ratio.numerator//ratio.denominator) for e,ns,w in support]
    support = [(e,ns,w) for e,ns,w in support if w]
    for e,total in enumerate(masses):
        left = total*scale-sum(w for parity,ns,w in support if parity==e)
        assert left>=0
        if left:
            support.append((e,(0,)*group_count,left))
    from itertools import combinations
    groups = [S for k in range(2,dimensions+1) for S in combinations(range(dimensions),k)]
    objective = 0
    for e,ns,w in support:
        f = [q-1-(b==e)-n for q,b,n in zip(high,pattern,ns[:dimensions])]
        cost = math.prod(f)-sum(ns[h]*math.prod(f[i] for i in range(dimensions) if i not in S) for h,S in enumerate(groups,dimensions))
        objective += w*cost
    return dict(denominator=scale,objective_numerator=objective,
                support=[dict(parity=e,counts=list(ns),weight=w) for e,ns,w in support],
                scope='exact rational distribution for the aggregate relaxation only; not a covering system')
