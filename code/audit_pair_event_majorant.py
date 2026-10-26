"""Independent exact verifier for affine bounds with pair-event multipliers."""
from fractions import Fraction


def number(v):
    assert len(v)==2 and all(type(k) is int for k in v) and v[1]>0
    return Fraction(*v)


def audit_bound(points,shapes,events,certificate):
    assert points and len(set(shapes))==len(shapes) and all(type(d) is int and d>0 for d in shapes)
    assert certificate['method']=='pair_affine'
    N=len(shapes);base=int(1 in shapes);assert certificate['base']==base
    rewards=[[0]*(N+1) for p in points]
    for e in events:
        x,t,w=e['point'],e['threshold'],e['weight']
        assert type(x) is int and 0<=x<len(points) and type(t) is int and 1<=t<=N and type(w) is int and w>=0
        for k in range(t,N+1):rewards[x][k]+=w
    slopes=list(map(number,certificate['slopes']));offsets=list(map(number,certificate['point_intercepts']))
    assert len(slopes)==len(offsets)==len(points) and min(slopes)>=0
    corrections=[[Fraction(0) for k in range(N+1)] for p in points];cut_total=Fraction(0)
    for row in certificate['pair_cuts']:
        assert len(row['points'])==len(row['thresholds'])==2
        x,y=row['points'];t,u=row['thresholds']
        assert all(type(v) is int for v in (x,y,t,u)) and 0<=x<y<len(points) and 1<=t<=N and 1<=u<=N
        # Each remaining resource contributes at most one hit unless its
        # residue classes can identify the two points, when it contributes two.
        capacity=sum(2 if (points[x]-points[y])%d==0 else 1 for d in shapes)
        assert t+u>capacity
        mu=number(row['multiplier']);assert mu>=0;cut_total+=mu
        for j,threshold in ((x,t),(y,u)):
            for k in range(threshold,N+1):corrections[j][k]+=mu
    for x in range(len(points)):
        for k in range(base,N+1):assert rewards[x][k]<=offsets[x]+slopes[x]*(k-base)+corrections[x][k]
    assert sum(offsets)==number(certificate['intercept_fraction'])
    capacities=[]
    for d in shapes:
        if d!=1:capacities.append(max(sum(slopes[x] for x,p in enumerate(points) if p%d==a) for a in range(d)))
    assert capacities==list(map(number,certificate['resource_bounds']))
    total=sum(offsets)+sum(capacities)+cut_total
    bound=min(sum(row[-1] for row in rewards),total.numerator//total.denominator)
    assert type(certificate['bound']) is int and bound==certificate['bound']
    return bound
