"""Solver-free rational upper bound on a weighted threshold-event score."""
from collections import defaultdict
from fractions import Fraction


def upper_bound(points,shapes,events,weights):
    assert len(events)==len(weights) and all(type(w) is int and w>=0 for w in weights)
    limit,base = len(shapes),int(1 in shapes)
    rewards = [[0]*(limit+1) for x in points]
    for (x,t),w in zip(events,weights):
        assert 0<=x<len(points) and 1<=t<=limit
        for k in range(t,limit+1):rewards[x][k] += w
    intercept = sum(row[base] for row in rewards)
    slopes = [max([Fraction(0)]+[Fraction(row[k]-row[base],k-base) for k in range(base+1,limit+1)]) for row in rewards]
    capacities = []
    for d in shapes:
        if d==1:continue
        totals = defaultdict(Fraction)
        for x,w in zip(points,slopes):totals[x%d] += w
        capacities.append(max(totals.values(),default=Fraction(0)))
    rational = intercept+sum(capacities)
    pointwise = sum(row[-1] for row in rewards)
    return dict(bound=min(pointwise,rational.numerator//rational.denominator),
                base=base,intercept=intercept,pointwise_bound=pointwise,
                slopes=[[v.numerator,v.denominator] for v in slopes],
                resource_bounds=[[v.numerator,v.denominator] for v in capacities],
                method='integer linear majorant of each threshold reward followed by resource capacities')


def optimized_upper_bound(points,shapes,events,weights,seconds=2):
    """Choose free affine intercepts; round slopes and rebuild all bounds exactly."""
    from ortools.linear_solver import pywraplp
    assert len(events)==len(weights) and all(type(w) is int and w>=0 for w in weights)
    limit,base=len(shapes),int(1 in shapes)
    rewards=[[sum(w for (x,t),w in zip(events,weights) if x==j and t<=k) for k in range(limit+1)] for j in range(len(points))]
    normalizer=max(1,max(row[-1] for row in rewards))
    solver=pywraplp.Solver.CreateSolver('GLOP');solver.SetTimeLimit(int(seconds*1000))
    slopes=[solver.NumVar(0,solver.infinity(),f'a{x}') for x in points]
    intercepts=[solver.NumVar(0,solver.infinity(),f'b{x}') for x in points]
    objective=solver.Objective();objective.SetMinimization()
    for x,row in enumerate(rewards):
        objective.SetCoefficient(intercepts[x],1)
        for k in range(base,limit+1):solver.Add(intercepts[x]+(k-base)*slopes[x]>=row[k]/normalizer)
    for d in shapes:
        if d==1:continue
        capacity=solver.NumVar(0,solver.infinity(),f'c{d}');objective.SetCoefficient(capacity,1)
        for residue in {x%d for x in points}:solver.Add(capacity>=sum(a for x,a in zip(points,slopes) if x%d==residue))
    if solver.Solve()!=pywraplp.Solver.OPTIMAL:return upper_bound(points,shapes,events,weights)
    denominator=10**6
    rational_slopes=[Fraction(max(0,round(var.solution_value()*normalizer*denominator)),denominator) for var in slopes]
    offsets=[max(Fraction(row[k])-(k-base)*slope for k in range(base,limit+1)) for row,slope in zip(rewards,rational_slopes)]
    capacities=[]
    for d in shapes:
        if d==1:continue
        capacities.append(max(sum(a for x,a in zip(points,rational_slopes) if x%d==residue) for residue in {x%d for x in points}))
    total=sum(offsets)+sum(capacities);bound=min(sum(row[-1] for row in rewards),total.numerator//total.denominator)
    old=upper_bound(points,shapes,events,weights)
    if bound>=old['bound']:return old
    return dict(bound=bound,base=base,point_intercepts=[[v.numerator,v.denominator] for v in offsets],
                intercept_fraction=[sum(offsets).numerator,sum(offsets).denominator],
                slopes=[[v.numerator,v.denominator] for v in rational_slopes],
                resource_bounds=[[v.numerator,v.denominator] for v in capacities],
                method='optimized affine threshold majorants with exact reconstructed intercepts')
