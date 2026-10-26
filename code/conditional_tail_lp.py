"""Keep tail capacities conditional on all plane projection hit counts.

For fixed core R and four additional primes, augment the plane exclusions
independently in each fibre until exactly min(D_i,n_i) labels are removed.
Any cover remains a cover after this augmentation. For a higher-codimension
shape S the capacity at x is then at most product_{i not in S}(D_i-n_i).
Thus a cover requires F = sum_x [product_i(D_i-n_i) -
sum_{|S|>=2} h_S(x) product_{i not in S}(D_i-n_i)] <= 0.
Here h_S counts core projections from all 16 divisors of 210. This LP links
all 15 hit-count marginals to actual core projection choices, and imposes
their exact core-subset concentration inequalities. It is only a relaxation
for one fixed core. Rounded duals are priced using integer arithmetic.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from core_fibre_plane_model import prepare


def build(incumbent, concentration):
    high, pattern = incumbent['high_primes'], incumbent['high_pattern']
    abstract = 'abstract_counts' in incumbent
    if abstract:
        core = {d:None for d in range(3,211) if 210%d==0}
        R = [0]*incumbent['abstract_counts'][0]+[1]*incumbent['abstract_counts'][1]
        assert np.asarray(concentration['two_layer_capacities']).shape==tuple(n+1 for n in incumbent['abstract_counts'])
    else:
        core, R, dims, plane_resources, _ = prepare(dict(system=incumbent['core_system']), high, pattern)
        assert R==concentration['core_residual']
    assert len(high)==4
    subsets = [list(S) for k in range(2,5) for S in itertools.combinations(range(4),k)]
    groups = [[i] for i in range(4)] + subsets
    moduli = [sorted(core)]*4 + [[1,2]+sorted(core)]*len(subsets)
    resources = []
    for h, shapes in enumerate(moduli if not abstract else []):
        for d in shapes:
            resources.append(dict(axis=h, core_modulus=d,
                options=[(a, [j for j,x in enumerate(R) if x%d==a]) for a in sorted({x%d for x in R})]))
    H = np.asarray(concentration['two_layer_capacities'],dtype=np.int64)
    m0,m1 = np.indices(H.shape)
    tail_H = H+m0+m1+np.maximum(m0,m1)
    grids, costs, volumes = [], [], []
    for e in (0,1):
        D = np.array([q-1-(b==e) for q,b in zip(high,pattern)],dtype=np.int64)
        grid = np.array(list(itertools.product(*(range(min(14,int(d))+1) for d in D))),dtype=np.int64)
        free = D-grid
        grids.append(grid)
        costs.append(np.prod(free,axis=1))
        volumes.append([np.prod(free[:,[i for i in range(4) if i not in S]],axis=1) for S in subsets])
    return R,groups,moduli,resources,[H,tail_H],grids,costs,volumes


def price_point(e, link, common, grids, costs, volumes, scale):
    """Return exact minimum, an attaining state, and point-local pricing size."""
    grid = grids[e]
    values = costs[e]*scale - sum(grid[:,i]*link[i]+common[e][i][grid[:,i]] for i in range(4))
    choices = []
    for h, volume in enumerate(volumes[e],4):
        v = np.arange(int(volume.max())+1,dtype=np.int64)
        hits = np.arange(17,dtype=np.int64)
        table = -scale*v[:,None]*hits[None,:]-link[h]*hits[None,:]-common[e][h][None,:]
        args = table.argmin(axis=1)
        values += table[np.arange(len(v)),args][volume]
        choices.append(args[volume])
    j = int(values.argmin())
    return int(values[j]), tuple(map(int,grid[j]))+tuple(int(a[j]) for a in choices), len(grid)


def solve(incumbent, concentration, seconds=60, aggregate=False, stop_after_positive=False,target_lower=1,verbose=True,stop_on_negative_primal=False):
    abstract = 'abstract_counts' in incumbent
    assert not abstract or aggregate
    R,groups,moduli,resources,Hs,grids,costs,volumes = build(incumbent,concentration)
    original_R = R
    if aggregate:
        R = [0,1]
        resources = []
    masses = [sum(n%2==e for n in original_R) for e in (0,1)] if aggregate else [1]*len(R)
    start = time.monotonic()
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    objective = solver.Objective()
    objective.SetMinimization()
    mass = [solver.Constraint(size,size) for size in masses]
    links = [[solver.Constraint(-solver.infinity(),0) for group in groups] for x in R]
    for j,r in enumerate(resources):
        c = solver.Constraint(1,1)
        for a,points in r['options']:
            v = solver.NumVar(0,solver.infinity(),f'w{j}_{a}')
            c.SetCoefficient(v,1)
            for x in points:
                links[x][r['axis']].SetCoefficient(v,-1)
    tail_rows = []
    for h,shapes in enumerate(moduli):
        H = Hs[int(h>=4)]
        m0,m1 = np.indices(H.shape)
        for t0,t1 in itertools.product(range(len(shapes)+1),repeat=2):
            if t0==t1==len(shapes):
                continue
            bound = int(np.max(H-t0*m0-t1*m1))
            tail_rows.append((h,t0,t1,bound,solver.Constraint(-solver.infinity(),bound)))
    seen = [set() for x in R]
    columns = []
    # Index nonzero excess coefficients once, avoiding all rows per column.
    row_lookup = [[[[] for n in range(len(shapes)+1)] for shapes in moduli] for e in (0,1)]
    for h,t0,t1,b,c in tail_rows:
        for e,t in enumerate((t0,t1)):
            for n in range(t+1,len(moduli[h])+1):
                row_lookup[e][h][n].append((c,n-t))
    def add_column(x, ns):
        if ns in seen[x]:
            return False
        seen[x].add(ns)
        v = solver.NumVar(0,solver.infinity(),f'z{x}_{len(seen[x])}')
        columns.append((x,ns,v))
        mass[x].SetCoefficient(v,1)
        for h,n in enumerate(ns):
            if n:
                if not aggregate:
                    links[x][h].SetCoefficient(v,n)
                for c,coef in row_lookup[R[x]%2][h][n]:
                    c.SetCoefficient(v,coef)
        free = [q-1-(b==R[x]%2)-n for q,b,n in zip(incumbent['high_primes'],incumbent['high_pattern'],ns[:4])]
        value = math.prod(free)-sum(ns[h]*math.prod(free[i] for i in range(4) if i not in S) for h,S in enumerate(groups[4:],4))
        objective.SetCoefficient(v,value)
        return True
    for x in range(len(R)):
        add_column(x,(0,)*len(groups))
    # Seed a feasible integer projection assignment to avoid a long degenerate
    # phase in which all point columns have zero tail coverage.
    seed = [[0]*len(groups) for x in original_R]
    for x,n in enumerate(original_R):
        for i,(q,b) in enumerate(zip(incumbent['high_primes'],incumbent['high_pattern'])):
            seed[x][i] = min(q-1-(b==n%2),sum(n%r['core_modulus']==r['core_residue'] for r in incumbent['choices'] if r['high_prime']==q))
    # In aggregate mode use the original resources only to seed columns; the
    # relaxation deliberately does not enforce these projection identities.
    seed_resources = build(incumbent,concentration)[3] if aggregate else resources
    for r in seed_resources:
        h = r['axis']
        if h<4:
            continue
        weights = [math.prod(q-1-(b==n%2)-seed[x][i] for i,(q,b) in enumerate(zip(incumbent['high_primes'],incumbent['high_pattern'])) if i not in groups[h]) for x,n in enumerate(original_R)]
        _,points = max(r['options'],key=lambda option:sum(weights[x] for x in option[1]))
        for x in points:
            seed[x][h] += 1
    for x,ns in enumerate(seed):
        add_column(original_R[x]%2 if aggregate else x,tuple(ns))
    if abstract:
        for e in (0,1):
            for n,h in itertools.product((0,2,4,6,8,10,12,14),(0,2,4,6,8,16)):
                plane = tuple(min(n,q-1-(b==e)) for q,b in zip(incumbent['high_primes'],incumbent['high_pattern']))
                add_column(e,plane+(h,)*(len(groups)-4))
    best = None
    iterations = 0
    converged = False
    primal_witness = None
    scale = 10**6
    while time.monotonic()-start<seconds:
        solver.SetTimeLimit(max(1,int(1000*(seconds-(time.monotonic()-start)))))
        status = solver.Solve()
        if status!=pywraplp.Solver.OPTIMAL:
            break
        iterations += 1
        if stop_on_negative_primal and aggregate and objective.Value()<0:
            from tail_primal_witness import exact_witness
            candidate = exact_witness(columns,[n%2 for n in R],masses,tail_rows,incumbent['high_primes'],incumbent['high_pattern'])
            if candidate['objective_numerator']<0:
                primal_witness = candidate
                break
        link_integer = np.zeros((len(R),len(groups)),dtype=np.int64) if aggregate else np.minimum(0,np.rint(np.array([[c.dual_value() for c in row] for row in links])*scale).astype(np.int64))
        tail_integer = np.minimum(0,np.rint(np.array([row[4].dual_value() for row in tail_rows])*scale).astype(np.int64))
        common = [[np.zeros(len(shapes)+1,dtype=np.int64) for shapes in moduli] for e in (0,1)]
        for coeff,(h,t0,t1,b,c) in zip(tail_integer,tail_rows):
            if coeff:
                for e,t in enumerate((t0,t1)):
                    common[e][h] += coeff*np.maximum(0,np.arange(len(moduli[h])+1)-t)
        points,prices = [],[]
        violation = 0
        for x,n in enumerate(R):
            p,state,_ = price_point(n%2,link_integer[x],common,grids,costs,volumes,scale)
            points.append(p)
            gap = mass[x].dual_value()-p/scale
            violation = max(violation,gap)
            if gap>0.001:
                prices.append((x,state))
        resource_integer = [min(sum(int(link_integer[x,r['axis']]) for x in pts) for a,pts in r['options']) for r in resources]
        numerator = sum(v*m for v,m in zip(points,masses))+sum(resource_integer)+sum(int(v)*row[3] for v,row in zip(tail_integer,tail_rows))
        certificate = dict(numerator=numerator,denominator=scale,integer_lower=-(-numerator//scale),
                           point_multipliers=points,resource_multipliers=resource_integer,
                           link_multipliers=link_integer.tolist(),
                           tail_multipliers=[dict(axis=row[0],thresholds=list(row[1:3]),bound=row[3],multiplier=int(v))
                                             for v,row in zip(tail_integer,tail_rows) if v])
        if best is None or numerator>best['numerator']:
            best = certificate
        if stop_after_positive and best['integer_lower']>=target_lower:
            break
        if verbose and (iterations%5==0 or not prices):
            print(json.dumps(dict(iteration=iterations,columns=sum(map(len,seen)),master=objective.Value(),
                                  lower=certificate['integer_lower'],best=best['integer_lower'],violation=violation,
                                  seconds=round(time.monotonic()-start,3))),flush=True)
        if not prices:
            converged = True
            break
        if not sum(add_column(x,state) for x,state in prices):
            break
    return dict(core_residual=None if abstract else original_R,high_primes=incumbent['high_primes'],high_pattern=incumbent['high_pattern'],
                groups=groups,certificate=best,iterations=iterations,converged_to_tolerance=converged,
                aggregate=aggregate,abstract_core=abstract,point_masses=masses,
                primal_relaxation_witness=primal_witness,
                master_status=int(status),seconds=round(time.monotonic()-start,3),
                scope=('aggregate conditional-tail bound from supplied concentration upper bounds; no core realizability claim' if abstract else
                       'conditional-tail necessary condition for one specified core, not a global obstruction'))


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--output',required=True)
    ap.add_argument('--aggregate',action='store_true')
    a = ap.parse_args()
    read = lambda name: json.loads(Path('data',name).read_text())
    result = solve(read('core-fibre-eight-plane.json'),read('core-hit-concentration-7-15.json'),a.time,a.aggregate)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('certificate','core_residual','groups')}),flush=True)
