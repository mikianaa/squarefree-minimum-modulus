"""Retain shared odd-plane labels in one high coordinate.

Other high coordinates retain the proved balanced-load lower vectors. The
line/point classes have only their total capacities. This is a necessary
condition for one core profile and one triple of even-plane allocations.
"""
import itertools
import math
from ortools.sat.python import cp_model
from plane_capacity_bound import ALLOCATIONS,balanced_lower


def solve(profile,high,pattern,allocations,axis=0,seconds=10):
    caps = [profile[2+2*j:4+2*j] for j in range(7)]
    dims = [[p-1-(pattern[i]==e) for i,p in enumerate(high)] for e in (0,1)]
    tables = []
    for e in (0,1):
        vectors = []
        for i in range(3):
            if i==axis:
                continue
            bits = ALLOCATIONS[allocations[i]]
            weights = [c[e] for c in caps]+[c[e] for bit,c in zip(bits,caps) if bit==e]
            vectors.append(balanced_lower(weights,dims[e][i]))
        tables.append([sum(max(0,profile[e]-load-x-y) for x in vectors[0] for y in vectors[1]) for load in range(profile[e]+1)])
    mdl = cp_model.CpModel()
    labels = []
    bins = {(e,a):[] for e in (0,1) for a in range(1,high[axis]) if a!=1 or pattern[axis]!=e}
    even_bits = ALLOCATIONS[allocations[axis]]
    choices = []
    for even in (False,True):
        for j,cap in enumerate(caps):
            parity = even_bits[j] if even else None
            allowed = [a for a in range(1,high[axis]) if not even or a!=1 or pattern[axis]!=parity]
            xs = [mdl.NewBoolVar(f'x{even}_{j}_{a}') for a in allowed]
            mdl.AddExactlyOne(xs)
            label = mdl.NewIntVar(1,high[axis]-1,f'label{even}_{j}')
            mdl.Add(label==sum(a*x for a,x in zip(allowed,xs)))
            labels.append(label)
            choices.append((even,j,allowed,xs))
            for a,x in zip(allowed,xs):
                for e in (0,1):
                    if (parity is None or parity==e) and (e,a) in bins:
                        bins[e,a].append(cap[e]*x)
    previous = 1
    for j,label in enumerate(labels):
        mdl.Add(label<=previous+1)
        current = mdl.NewIntVar(1,high[axis]-1,f'max{j}')
        mdl.AddMaxEquality(current,[previous,label])
        previous = current
    deficits = [[],[]]
    for (e,a),terms in bins.items():
        u = mdl.NewIntVar(0,tables[e][0],f'deficit{e}_{a}')
        load = sum(terms)
        # The convex table is the maximum of these supporting affine pieces.
        for k in range(profile[e]):
            slope = tables[e][k+1]-tables[e][k]
            mdl.Add(u>=tables[e][k]+slope*(load-k))
        deficits[e].append(u)
    items = []
    for k in (2,3):
        for subset in itertools.combinations(range(3),k):
            volume = [math.prod(dims[e][i] for i in range(3) if i not in subset) for e in (0,1)]
            items.extend(tuple(c[e]*volume[e] for e in (0,1)) for c in [profile[:2]]+caps)
    odd = [sum(c[e] for c in items) for e in (0,1)]
    moves = [mdl.NewBoolVar(f'move{j}') for j in range(len(items))]
    mdl.Add(sum(deficits[0])<=odd[0]+sum(c[0]*x for c,x in zip(items,moves)))
    mdl.Add(sum(deficits[1])<=odd[1]+sum(c[1]*(1-x) for c,x in zip(items,moves)))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = 1
    status = solver.Solve(mdl)
    out = dict(status=solver.StatusName(status),axis=axis,allocations=list(allocations),seconds=solver.WallTime(),
               scope='one shared-label coordinate plus relaxed capacities; INFEASIBLE remains solver-dependent')
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        chosen = [(even,j,next(a for a,x in zip(allowed,xs) if solver.Value(x))) for even,j,allowed,xs in choices]
        actual = [0,0]
        for e,a in bins:
            load = sum(caps[j][e] for even,j,label in chosen if label==a and (not even or even_bits[j]==e))
            actual[e] += tables[e][min(load,profile[e])]
        assigned = [solver.Value(x) for x in moves]
        capacity = [odd[e]+sum(c[e] for c,x in zip(items,assigned) if x==(1-e)) for e in (0,1)]
        assert all(actual[e]<=capacity[e] for e in (0,1))
        out.update(chosen_planes=chosen,tail_moves=assigned,deficits=actual,capacities=capacity,verified=True)
    return out
