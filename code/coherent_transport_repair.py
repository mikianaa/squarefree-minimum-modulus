"""Integer quantile repair and exact one-hit cost bounds for group mixtures.

Repair changes point counts downwards, preserving each point's mass and
leaving the chosen actual group assignments intact. No solver is used here.
"""
from collections import defaultdict
from fractions import Fraction
import itertools
import math


def drop_cost_bounds(high, pattern, parity):
    D = [q-1-(b==parity) for q,b in zip(high,pattern)]
    groups = [S for k in range(1,len(D)+1)
              for S in itertools.combinations(range(len(D)),k)]
    return [math.prod(d for i,d in enumerate(D) if i not in S) for S in groups]


def clip_distribution(support, supply, mass):
    """Return integer (state,weight) pairs dominated by survival supplies.

    Each supply[h][t] is the available mass with count >=t; entry zero is
    unused. Monotone quantile coupling and taking min(old,available) removes
    exactly sum_t max(demand[t]-supply[t],0) hits on each coordinate.
    """
    rows = {tuple(state): 0 for state,w in support}
    for state,w in support:
        assert type(w) is int and w>=0
        rows[tuple(state)] += w
    assert sum(rows.values())==mass
    removed = []
    for h,line in enumerate(supply):
        assert mass>=line[1]>=0
        assert all(a>=b>=0 for a,b in zip(line[1:],line[2:]))
        demand = [sum(w for state,w in rows.items() if state[h]>=t)
                  for t in range(1,len(line))]
        deficit = sum(max(d-s,0) for d,s in zip(demand,line[1:]))
        removed.append(deficit)
        if not deficit:
            continue
        survival = [mass]+list(line[1:])+[0]
        atoms = [(t,survival[t]-survival[t+1])
                 for t in range(len(line)-1,-1,-1)
                 if survival[t]>survival[t+1]]
        j,left = 0,atoms[0][1]
        repaired = defaultdict(int)
        actual_drop = 0
        for state,w in sorted(rows.items(),key=lambda pair:pair[0][h],reverse=True):
            while w:
                used = min(w,left)
                count = min(state[h],atoms[j][0])
                changed = list(state);changed[h] = count
                repaired[tuple(changed)] += used
                actual_drop += used*(state[h]-count)
                w -= used;left -= used
                if not left:
                    j += 1
                    if j<len(atoms):left = atoms[j][1]
        assert actual_drop==deficit
        rows = dict(repaired)
    return [(state,w) for state,w in rows.items() if w],removed


def rational_witness(source,model,columns,configurations,shapes):
    scale,R = 10**9,source['core_residual']
    supply = [[[0]*(len(shape)+1) for shape in shapes] for x in R]
    groups = []
    for h,configs in enumerate(configurations):
        values = [max(Fraction(0),Fraction(str(var.solution_value())))
                  for residues,counts,var in configs]
        total = sum(values);assert total>0
        weights = [int(value*scale//total) for value in values]
        weights[0] += scale-sum(weights)
        assignments = []
        for (residues,counts,var),w in zip(configs,weights):
            if not w:continue
            assignments.append(dict(residues=list(residues),weight=w))
            for x,v in enumerate(counts):
                for t in range(1,v+1):supply[x][h][t] += w
        groups.append(dict(axis=h,assignments=assignments))
    point_rows = [[] for x in R]
    for x,state,var in columns:
        value = max(Fraction(0),Fraction(str(var.solution_value())))
        if value:point_rows[x].append((state,value))
    support,repair = [],[]
    for x,rows in enumerate(point_rows):
        total = sum(value for state,value in rows);assert total>0
        integer_rows = [(state,int(value*scale//total)) for state,value in rows]
        missing = scale-sum(w for state,w in integer_rows)
        integer_rows.append(((0,)*len(shapes),missing))
        fixed,removed = clip_distribution(integer_rows,supply[x],scale)
        before = sum(w*model.cost(R[x]%2,state) for state,w in integer_rows)
        after = sum(w*model.cost(R[x]%2,state) for state,w in fixed)
        bounds = drop_cost_bounds(source['high_primes'],source['high_pattern'],R[x]%2)
        penalty = sum(b*d for b,d in zip(bounds,removed))
        assert after<=before+penalty
        repair.append(dict(point_index=x,removed_hits=removed,cost_before=before,
                           cost_after=after,cost_increase_bound=penalty))
        support.extend(dict(point_index=x,counts=list(state),weight=w) for state,w in fixed)
    objective = sum(row['weight']*model.cost(R[row['point_index']]%2,row['counts']) for row in support)
    return dict(denominator=scale,objective_numerator=objective,groups=groups,
                support=support,transport_repair=repair,
                positive_point_means_numerator=sum(max(0,row['cost_after']) for row in repair))
