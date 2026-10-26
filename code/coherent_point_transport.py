"""Exact finite-support reduction of a fixed-mixture point problem.

Mathematical reduction is exact; GLOP outputs are numerical candidates and
must still be rounded/repaired and independently audited before use.
"""
import itertools
import math
from ortools.linear_solver import pywraplp
from two_configuration_mixtures import structure,cost
from coherent_transport_repair import drop_cost_bounds


def atom_counts(survival,limit):
    tail = [1.0]+[min(1,max(0,float(v))) for v in survival[1:limit+1]]+[0.0]
    assert all(a+1e-9>=b for a,b in zip(tail,tail[1:]))
    return [t for t in range(limit+1) if tail[t]-tail[t+1]>1e-13]


def lift_potential(endpoints,duals,limit,bound,normalization,scale=10**6):
    """Interpolate cumulative event penalties, then extend with slope -B.

    The one-coordinate cost is affine between endpoints. Beyond the last
    supplied count its downward slope is no steeper than -B.
    """
    values = {0:0.0}
    total = 0.0
    for t in sorted(set(endpoints)-{0}):
        total += min(0,duals.get(t,0.0)*normalization)
        values[t] = total
    line = [0]*(limit+1)
    left = 0
    for right in sorted(set(endpoints)-{0}):
        slope = (values[right]-values[left])/(right-left)
        for t in range(left+1,right+1):line[t] = min(0,round(slope*scale))
        left = right
    for t in range(left+1,limit+1):line[t] = -bound*scale
    return line


def solve_point(source,point,supply,max_states=4096,max_variables=100000):
    groups,shapes = structure(source)
    n = len(source['high_primes'])
    D = [q-1-(b==point%2) for q,b in zip(source['high_primes'],source['high_pattern'])]
    planes = [sorted({0,*atom_counts(supply[h],min(14,D[h]))}) for h in range(n)]
    size = math.prod(map(len,planes))
    if size>max_states:return None
    tails = [atom_counts(supply[h],len(shapes[h])) for h in range(n,len(groups))]
    assert all(tails) if tails else True
    states = list(itertools.product(*planes))
    free_rows = [[d-v for d,v in zip(D,state)] for state in states]
    volumes = [[math.prod(free[i] for i in range(n) if i not in S) for S in groups[n:]] for free in free_rows]
    volume_groups = [dict() for atoms in tails]
    for j,row in enumerate(volumes):
        for h,v in enumerate(row):volume_groups[h].setdefault(v,[]).append(j)
    needed = size+sum(sum(v!=0 for v in table)*(len(atoms)-1) for table,atoms in zip(volume_groups,tails))
    if needed>max_variables:return None
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetTimeLimit(3000)
    objective = solver.Objective();objective.SetMinimization()
    normalization = max(1,math.prod(D))
    mass = [solver.NumVar(0,solver.infinity(),f'm{j}') for j in range(size)]
    mass_constraint = solver.Add(sum(mass)==1)
    plane_links = [dict() for h in range(n)]
    for h,values in enumerate(planes):
        for t in values:
            if t:
                plane_links[h][t] = solver.Add(sum(p for s,p in zip(states,mass) if s[h]>=t)<=supply[h][t])
    marginal = [[solver.Constraint(-solver.infinity(),supply[h][t]) for t in atoms[1:]]
                for h,atoms in enumerate(tails,n)]
    for j,p in enumerate(mass):
        baseline = math.prod(free_rows[j])-sum(atoms[0]*v for atoms,v in zip(tails,volumes[j]))
        objective.SetCoefficient(p,baseline/normalization)
    transport = [dict() for atoms in tails]
    for tail,(atoms,table,constraints) in enumerate(zip(tails,volume_groups,marginal)):
        for v,indices in table.items():
            if not v:continue
            previous = sum(mass[j] for j in indices)
            row = []
            for k,t in enumerate(atoms[1:]):
                y = solver.NumVar(0,solver.infinity(),f'y{tail}_{v}_{t}')
                solver.Add(y<=previous)
                constraints[k].SetCoefficient(y,1)
                objective.SetCoefficient(y,-(t-atoms[k])*v/normalization)
                row.append((t-atoms[k],y))
                previous = y
            transport[tail][v] = row
    status = solver.Solve()
    if status!=pywraplp.Solver.OPTIMAL:return None
    support = []
    group_masses = [{v:sum(max(0,mass[j].solution_value()) for j in indices) for v,indices in table.items()}
                    for table in volume_groups]
    for j,(state,p) in enumerate(zip(states,mass)):
        weight = max(0,p.solution_value())
        if weight<1e-13:continue
        heights = []
        for tail,table in enumerate(transport):
            v = volumes[j][tail]
            if not v:continue
            fraction = weight/group_masses[tail][v]
            heights.extend((n+tail,delta,min(weight,max(0,var.solution_value()*fraction))) for delta,var in table[v])
        endpoints = sorted({0.0,weight,*[height for h,delta,height in heights]})
        for left,right in zip(endpoints,endpoints[1:]):
            if right-left<1e-13:continue
            counts = list(state)+[atoms[0] for atoms in tails]
            for h,delta,height in heights:
                if height>left:counts[h] += delta
            support.append((counts,right-left))
    value = sum(p*cost(source,point,state) for state,p in support)
    assert abs(value-objective.Value()*normalization)<1e-5*max(1,abs(value))
    bounds = drop_cost_bounds(source['high_primes'],source['high_pattern'],point%2)
    multipliers = []
    for h,endpoints in enumerate(planes):
        duals = {t:c.dual_value() for t,c in plane_links[h].items()}
        limit = min(14,D[h])
        # Counts above D_h have already been clipped and cannot occur in a
        # local state. Their constraints have zero demand, so their duals
        # must remain zero even when the original supply has more hits.
        line = lift_potential(endpoints,duals,limit,bounds[h],normalization)
        multipliers.append(line+[0]*(14-limit))
    for h,(endpoints,constraints) in enumerate(zip(tails,marginal),n):
        duals = {t:c.dual_value() for t,c in zip(endpoints[1:],constraints)}
        multipliers.append(lift_potential([0]+endpoints,duals,16,bounds[h],normalization))
    return dict(objective=value,support=support,plane_states=size,
                event_multipliers=multipliers,denominator=10**6,
                unrounded_mass_dual=mass_constraint.dual_value()*normalization,
                variables=solver.NumVariables(),constraints=solver.NumConstraints())
