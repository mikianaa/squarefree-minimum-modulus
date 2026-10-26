"""Rational affine group bounds strengthened by pair-threshold incompatibility."""
from fractions import Fraction
from fractional_group_scenarios import maximize
from event_score_majorant import optimized_upper_bound


def pair(v):return [v.numerator,v.denominator]


def upper_bound(points,shapes,events,weights):
    assert len(events)==len(weights) and all(type(w) is int and w>=0 for w in weights)
    N=len(shapes);base=int(1 in shapes)
    for x,t in events:assert type(x) is int and 0<=x<len(points) and type(t) is int and 1<=t<=N
    rewards=[[sum(w for (x,t),w in zip(events,weights) if x==j and t<=k) for k in range(N+1)] for j in range(len(points))]
    primal=maximize(points,shapes,rewards,pair_compatibility=True,return_dual=True)
    dual=primal['numerical_dual'];quantum=10**6
    slopes=[Fraction(max(0,round(v*quantum)),quantum) for v in dual['slopes']]
    cuts=[];adjustments=[[Fraction(0) for k in range(N+1)] for p in points]
    for row in dual['pair_cuts']:
        mu=Fraction(max(0,round(row['multiplier']*quantum)),quantum)
        if not mu:continue
        cuts.append(dict(points=row['points'],thresholds=row['thresholds'],multiplier=pair(mu)))
        for x,t in zip(row['points'],row['thresholds']):
            for k in range(t,N+1):adjustments[x][k]+=mu
    offsets=[max(Fraction(rewards[x][k])-(k-base)*slopes[x]-adjustments[x][k] for k in range(base,N+1))
             for x in range(len(points))]
    capacities=[]
    for d in shapes:
        if d!=1:capacities.append(max(sum(a for p,a in zip(points,slopes) if p%d==r) for r in {p%d for p in points}))
    total=sum(offsets)+sum(capacities)+sum(Fraction(*row['multiplier']) for row in cuts)
    bound=min(sum(row[-1] for row in rewards),total.numerator//total.denominator)
    old=optimized_upper_bound(points,shapes,events,weights)
    if bound>=old['bound']:return old
    return dict(method='pair_affine',base=base,bound=bound,slopes=list(map(pair,slopes)),
                point_intercepts=list(map(pair,offsets)),intercept_fraction=pair(sum(offsets)),
                resource_bounds=list(map(pair,capacities)),pair_cuts=cuts)
