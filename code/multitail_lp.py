"""Aggregate conditional-tail relaxation with up to five higher primes.

Keeps the 210-core, 14 plane projection shapes, and 16 tail shapes. All
nonempty subsets of the supplied higher primes are included. This generator
does not independently verify its certificates or the supplied capacity H.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from tail_primal_witness import exact_witness


class LocalStates:
    def __init__(self,high,pattern):
        self.n = len(high)
        assert 1<=self.n<=5 and len(pattern)==self.n
        self.high,self.pattern = list(high),list(pattern)
        self.groups = [tuple(S) for k in range(1,self.n+1) for S in itertools.combinations(range(self.n),k)]
        self.grids,self.costs,self.volumes = [],[],[]
        for e in (0,1):
            D = [q-1-(b==e) for q,b in zip(high,pattern)]
            shape = tuple(min(14,d)+1 for d in D)
            grid = np.indices(shape,dtype=np.int8).reshape(self.n,-1).T.copy()
            free = np.array(D,dtype=np.int64)-grid
            self.grids.append(grid)
            self.costs.append(np.prod(free,axis=1))
            self.volumes.append([np.prod(free[:,[i for i in range(self.n) if i not in S]],axis=1)
                                 for S in self.groups[self.n:]])

    def cost(self,e,state):
        f = [q-1-(b==e)-state[i] for i,(q,b) in enumerate(zip(self.high,self.pattern))]
        return math.prod(f)-sum(state[h]*math.prod(f[i] for i in range(self.n) if i not in S)
                                for h,S in enumerate(self.groups[self.n:],self.n))

    def price(self,e,common,scale,batch=1,diverse=False):
        grid = self.grids[e]
        # Establish a conservative no-overflow bound before integer arrays.
        worst = abs(scale)*int(self.costs[e].max())+sum(max(map(abs,map(int,c))) for c in common)
        worst += sum(abs(scale)*16*int(v.max()) for v in self.volumes[e])
        assert worst<2**63
        values = scale*self.costs[e]-sum(common[i][grid[:,i]] for i in range(self.n))
        choices = []
        for h,volume in enumerate(self.volumes[e],self.n):
            v = np.arange(int(volume.max())+1,dtype=np.int64)
            table = -scale*v[:,None]*np.arange(17,dtype=np.int64)[None,:]-common[h][None,:]
            arg = table.argmin(axis=1)
            values += table[np.arange(len(v)),arg][volume]
            choices.append(arg)
        size = min(batch,len(values))
        indices = np.argpartition(values,size-1)[:size]
        if diverse:
            # Adjacent globally cheapest grid points can have almost the same
            # counts. Include a cheapest point on every coordinate slice;
            # the globally exact minimum is retained without approximation.
            indices = set(map(int,indices))
            shape = tuple(int(grid[:,i].max())+1 for i in range(self.n))
            tensor = values.reshape(shape)
            for axis,length in enumerate(shape):
                for v in range(length):
                    selection = [slice(None)]*self.n;selection[axis] = v
                    section = tensor[tuple(selection)]
                    coordinates = list(np.unravel_index(int(np.argmin(section)),section.shape))
                    coordinates.insert(axis,v)
                    indices.add(int(np.ravel_multi_index(tuple(coordinates),shape)))
        indices = sorted(map(int,indices),key=lambda i:(int(values[i]),i))
        rows = []
        for i in indices:
            state = tuple(map(int,grid[i]))+tuple(int(arg[vol[i]]) for arg,vol in zip(choices,self.volumes[e]))
            rows.append((int(values[i]),state))
        return rows


def solve(source,seconds=120,batch=12,target_lower=1,rank_cuts=False,seed=None):
    start = time.monotonic()
    high,pattern,counts = source['high_primes'],source['high_pattern'],source['counts']
    model = LocalStates(high,pattern)
    n,groups = model.n,model.groups
    H = np.array(source['H_upper'],dtype=np.int64)
    assert H.shape==tuple(v+1 for v in counts)
    a,b = np.indices(H.shape)
    tables = [H,H+a+b+np.maximum(a,b)]
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    objective = solver.Objective()
    objective.SetMinimization()
    mass = [solver.Constraint(v,v) for v in counts]
    rows = []
    limits = [14 if h<n else 16 for h in range(len(groups))]
    lookup = [[[[] for v in range(limit+1)] for limit in limits] for e in (0,1)]
    for h,limit in enumerate(limits):
        table = tables[int(h>=n)]
        for t0,t1 in itertools.product(range(limit+1),repeat=2):
            if t0==t1==limit:
                continue
            bound = int(np.max(table-t0*a-t1*b))
            c = solver.Constraint(-solver.infinity(),bound)
            rows.append((h,t0,t1,bound,c))
            for e,t in enumerate((t0,t1)):
                for v in range(t+1,limit+1):
                    lookup[e][h][v].append((c,v-t))
    columns = []
    seen = [set(),set()]
    def add(e,state):
        if state in seen[e]:
            return 0
        seen[e].add(state)
        variable = solver.NumVar(0,solver.infinity(),f'z{e}_{len(seen[e])}')
        mass[e].SetCoefficient(variable,1)
        for h,v in enumerate(state):
            for c,coefficient in lookup[e][h][v]:
                c.SetCoefficient(variable,coefficient)
        objective.SetCoefficient(variable,model.cost(e,state))
        columns.append((e,state,variable))
        return 1
    for e in (0,1):
        add(e,(0,)*len(groups))
        for plane,tail in itertools.product(range(0,15,2),(0,2,4,6,8,16)):
            state = tuple(min(plane,q-1-(p==e)) for q,p in zip(high,pattern))+(tail,)*(len(groups)-n)
            add(e,state)
    if seed:
        for row in seed['support']:
            add(row['parity'],tuple(row['counts']))
    rank_rows,rank_seen = [],set()
    rank_catalogue = []
    if rank_cuts:
        from rank_capacity import catalogue,violations
        rank_catalogue = catalogue(H,n)
    def add_rank(row):
        if row in rank_seen:
            return
        rank_seen.add(row)
        h,t0,t1,w0,w1,bound = row
        c = solver.Constraint(-solver.infinity(),bound)
        rank_rows.append((*row,c))
        for e,t,w in ((0,t0,w0),(1,t1,w1)):
            if w:
                for v in range(t,limits[h]+1):
                    lookup[e][h][v].append((c,w))
        for e,state,variable in columns:
            coefficient = (w0 if e==0 else w1)*int(state[h]>=(t0 if e==0 else t1))
            if coefficient:
                c.SetCoefficient(variable,coefficient)
    scale = 10**6
    best,witness = None,None
    iterations,status = 0,None
    while time.monotonic()-start<seconds:
        solver.SetTimeLimit(max(1,int(1000*(seconds-(time.monotonic()-start)))))
        status = solver.Solve()
        if status!=pywraplp.Solver.OPTIMAL:
            break
        iterations += 1
        if objective.Value()<0:
            candidate = exact_witness(columns,[0,1],counts,rows,high,pattern,rank_rows)
            if candidate['objective_numerator']<0:
                violated = violations(candidate,n,rank_catalogue) if rank_cuts else []
                if violated:
                    for gap,row in violated:
                        add_rank(row)
                    print(json.dumps(dict(stage='rank_cuts',added=len(violated),total=len(rank_rows),seconds=round(time.monotonic()-start,3))),flush=True)
                    continue
                witness = candidate
                break
        multipliers = np.minimum(0,np.rint(np.array([r[4].dual_value() for r in rows])*scale).astype(np.int64))
        common = [[np.zeros(limit+1,dtype=np.int64) for limit in limits] for e in (0,1)]
        for y,(h,t0,t1,bound,c) in zip(multipliers,rows):
            if y:
                for e,t in enumerate((t0,t1)):
                    common[e][h] += y*np.maximum(0,np.arange(limits[h]+1)-t)
        rank_multipliers = [min(0,round(row[-1].dual_value()*scale)) for row in rank_rows]
        for y,(h,t0,t1,w0,w1,bound,c) in zip(rank_multipliers,rank_rows):
            for e,t,w in ((0,t0,w0),(1,t1,w1)):
                common[e][h] += y*w*(np.arange(limits[h]+1)>=t)
        prices = [model.price(e,common[e],scale,batch) for e in (0,1)]
        points = [r[0][0] for r in prices]
        numerator = sum(p*v for p,v in zip(points,counts))+sum(int(y)*r[3] for y,r in zip(multipliers,rows))
        numerator += sum(y*r[5] for y,r in zip(rank_multipliers,rank_rows))
        certificate = dict(numerator=numerator,denominator=scale,integer_lower=-(-numerator//scale),point_multipliers=points,
                           tail_multipliers=[dict(axis=r[0],thresholds=list(r[1:3]),bound=r[3],multiplier=int(y)) for y,r in zip(multipliers,rows) if y],
                           rank_multipliers=[dict(axis=r[0],thresholds=list(r[1:3]),weights=list(r[3:5]),bound=r[5],multiplier=y) for y,r in zip(rank_multipliers,rank_rows) if y])
        if best is None or numerator>best['numerator']:
            best = certificate
        if iterations%5==0 or best['integer_lower']>=target_lower:
            print(json.dumps(dict(iteration=iterations,columns=len(columns),master=objective.Value(),
                                  lower=best['integer_lower'],seconds=round(time.monotonic()-start,3))),flush=True)
        if best['integer_lower']>=target_lower:
            break
        additions = 0
        mass_duals = [c.dual_value() for c in mass]
        for e in (0,1):
            for value,state in prices[e]:
                if mass_duals[e]-value/scale>0.001:
                    additions += add(e,state)
        if not additions:
            break
    return dict(certificate=best,primal_relaxation_witness=witness,iterations=iterations,
                master_status=status,columns=len(columns),seconds=round(time.monotonic()-start,3),
                groups=groups,plane_grid_sizes=[len(g) for g in model.grids],rank_cuts_enabled=rank_cuts,
                rank_constraints=[list(row[:-1]) for row in rank_rows],
                column_states=[dict(parity=e,counts=list(state)) for e,state,var in columns],
                scope='aggregate conditional-tail relaxation only; supplied capacity validity requires a separate audit')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=120)
    ap.add_argument('--batch',type=int,default=12)
    ap.add_argument('--target-lower',type=int,default=1)
    ap.add_argument('--rank-cuts',action='store_true')
    ap.add_argument('--seed-source')
    args = ap.parse_args()
    source = json.loads(Path(args.input).read_text())
    seed = None
    if args.seed_source:
        saved = json.loads(Path(args.seed_source).read_text())['result']
        seed = saved.get('primal_relaxation_witness') or dict(support=saved['column_states'])
    result = solve(source,args.time,args.batch,args.target_lower,args.rank_cuts,seed)
    Path(args.output).write_text(json.dumps(dict(input=source,result=result),separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('certificate','groups','primal_relaxation_witness','column_states','rank_constraints')}|
                     dict(lower=result['certificate']['integer_lower'] if result['certificate'] else None,negative_primal=bool(result['primal_relaxation_witness']))),flush=True)
