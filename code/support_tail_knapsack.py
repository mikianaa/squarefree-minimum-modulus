"""Exact convexified tail support pricing for a fixed plane tuple.

The fractional solution is not an integer count state. Its rational slope
certifies a lower bound on the event penalty of every F<=0 count state.
"""
from fractions import Fraction


def minimum(target,volumes,penalties,plane_penalty=0):
    assert type(target) is int and target>=0
    assert len(volumes)==len(penalties)
    assert all(type(v) is int and v>=0 for v in volumes)
    assert all(row and row[0]==0 and all(type(v) is int and v>=0 for v in row)
               and all(a<=b for a,b in zip(row,row[1:])) for row in penalties)
    if target>sum(v*(len(row)-1) for v,row in zip(volumes,penalties)):return None
    segments=[]
    for h,(volume,row) in enumerate(zip(volumes,penalties)):
        if not volume:continue
        hull=[]
        for k,q in enumerate(row):
            while len(hull)>=2:
                a,b=hull[-2:]
                if (row[b]-row[a])*(k-b)<(q-row[b])*(b-a):break
                hull.pop()
            hull.append(k)
        for a,b in zip(hull,hull[1:]):
            credit=volume*(b-a);cost=row[b]-row[a]
            segments.append((Fraction(cost,credit),h,a,b,credit,cost))
    remaining=target;value=Fraction(plane_penalty);slope=Fraction(0)
    counts=[Fraction(0) for row in penalties]
    for density,h,a,b,credit,cost in sorted(segments):
        if not remaining:break
        used=min(remaining,credit)
        counts[h]+=Fraction(used,credit)*(b-a)
        value+=Fraction(used,credit)*cost
        remaining-=used;slope=density
    assert remaining==0
    dual=plane_penalty+slope*target+sum(min(q-slope*v*k for k,q in enumerate(row)) for v,row in zip(volumes,penalties))
    assert dual==value
    return dict(value=value,slope=slope,counts=counts)
