"""Capacity relaxation retaining all coordinates beyond a fixed 210 core.

The normalized core family leaves r_0,r_1 points in its parity layers. For each
core divisor d, a class can meet at most min(r_e, product D_p(e) over p not in d)
of these points. Keep the actual placements in the additional prime coordinates
and require their capacities to cover r_e in every surviving high-coordinate
cell. Base residue choices are relaxed, so FEASIBLE does not construct a cover.

Coverability of this relaxation is monotone downwards in (r_0,r_1): each ratio
min(r,C)/r decreases with r. Therefore it suffices to test a componentwise lower
frontier proved for the core. UNSAT for one fixed pair/pattern excludes that case
only, never the entire problem.
"""
import argparse
import itertools
import json
import math
import time
from ortools.sat.python import cp_model


def build(high_primes, pattern, residual, symmetry=True, minimize=False, frontier=None, core_capacities=None):
    core = (3,5,7)
    assert len(pattern) == 3+len(high_primes) and pattern[0] == 0
    assert all(b in (0,1) for b in pattern)
    assert list(high_primes) == sorted(set(high_primes))
    assert all(p>7 and all(p%d for d in range(2,math.isqrt(p)+1)) for p in high_primes)
    mdl = cp_model.CpModel()
    core_dims = [[p-1-(pattern[i] == e) for i,p in enumerate(core)] for e in (0,1)]
    if frontier is not None:
        assert minimize
        residual = [math.prod(d) for d in core_dims]
        targets = [mdl.NewIntVar(0,residual[e],f'core_residual{e}') for e in (0,1)]
        mdl.AddAllowedAssignments(targets,frontier)
    else:
        assert residual is not None and len(residual) == 2
        assert all(0<=residual[e]<=math.prod(core_dims[e]) for e in (0,1))
        targets = residual
    alphabets = [[tuple(v for v in range(1,p) if v != 1 or pattern[i+3] != e)
                  for i,p in enumerate(high_primes)] for e in (0,1)]
    cells = [(e,a) for e in (0,1) if residual[e] for a in itertools.product(*alphabets[e])]
    constraints = [[] for _ in cells]
    ds = [d for d in range(1,211) if 210 % d == 0]
    resources = []
    for k in range(1,len(high_primes)+1):
        for subset in itertools.combinations(range(len(high_primes)),k):
            for d in ds:
                if k == 1 and d in (1,2):
                    continue  # These are the already normalized high anchors.
                resources.append((d*math.prod(high_primes[i] for i in subset),d,subset))
    resources.sort()
    labels_by_prime = [[] for _ in high_primes]
    placements_all = []
    for modulus,d,subset in resources:
        cap = [min(residual[e], math.prod(core_dims[e][i] for i,p in enumerate(core) if d % p)) for e in (0,1)]
        if core_capacities is not None:
            assert frontier is None, 'fixed-family capacities cannot be used for arbitrary core families'
            supplied = core_capacities[d]
            assert all(0<=supplied[e]<=cap[e] for e in (0,1))
            cap = supplied
        if not any(cap):
            continue
        even = d % 2 == 0
        placements = []
        lookup = {}
        for parity in ((0,1) if even else (None,)):
            alphabet = [alphabets[parity][i] if even else range(1,high_primes[i]) for i in subset]
            for labels in itertools.product(*alphabet):
                if even and not cap[parity]:
                    continue
                if not even and not any(cap[e] and all(a in alphabets[e][i] for i,a in zip(subset,labels)) for e in (0,1)):
                    continue
                v = mdl.NewBoolVar(f'x{modulus}_{parity}_{labels}')
                placements.append((parity,labels,v))
                lookup[(parity,labels)] = v
        mdl.AddExactlyOne(v for _,_,v in placements)
        for j,i in enumerate(subset):
            label = mdl.NewIntVar(1,high_primes[i]-1,f'label{modulus}_{high_primes[i]}')
            mdl.Add(label == sum(labels[j]*v for _,labels,v in placements))
            labels_by_prime[i].append(label)
        for j,(e,a) in enumerate(cells):
            if cap[e]:
                v = lookup[(e if even else None,tuple(a[i] for i in subset))]
                constraints[j].append(cap[e]*v)
        placements_all.append((modulus,d,subset,cap,placements))
    deficits = []
    for j,((e,a), terms) in enumerate(zip(cells,constraints)):
        if minimize:
            u = mdl.NewIntVar(0,residual[e],f'uncovered{j}')
            deficits.append(u)
            mdl.Add(sum(terms)+u >= targets[e])
        else:
            mdl.Add(sum(terms) >= residual[e])
    if minimize:
        mdl.Minimize(sum(deficits))
    if symmetry:
        for i,labels in enumerate(labels_by_prime):
            prev = 1
            for j,label in enumerate(labels):
                mdl.Add(label <= prev+1)
                current = mdl.NewIntVar(1,high_primes[i]-1,f'max{i}_{j}')
                mdl.AddMaxEquality(current,[prev,label])
                prev = current
    return mdl,cells,placements_all,targets


def solve(high_primes,pattern,residual,seconds=60,workers=8,minimize=False,frontier=None,core_capacities=None):
    start = time.monotonic()
    mdl,cells,resources,targets = build(high_primes,pattern,residual,minimize=minimize,frontier=frontier,core_capacities=core_capacities)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.linearization_level = 2
    built = time.monotonic()
    status = solver.Solve(mdl)
    out = dict(N=210*math.prod(high_primes),pattern=pattern,core_residual=None if frontier is not None else residual,
               status=solver.StatusName(status),cells=len(cells),resources=len(resources),
               build_seconds=round(built-start,3),seconds=round(time.monotonic()-start,3),
               minimize=minimize,variable_core_frontier=frontier is not None,
               fixed_family_capacities=core_capacities is not None,
               scope='capacity relaxation at one fixed anchor pattern')
    if minimize:
        out['uncovered_lower_bound'] = math.ceil(solver.BestObjectiveBound()-1e-7)
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        chosen = []
        for modulus,d,subset,cap,placements in resources:
            parity,labels,_ = next(t for t in placements if solver.Value(t[2]))
            chosen.append(dict(modulus=modulus,core_divisor=d,high_subset=list(subset),
                               parity=parity,labels=list(labels),capacity=cap))
        # Independently check every aggregate cell, without reading constraints.
        selected = [solver.Value(v) for v in targets] if frontier is not None else residual
        deficit = 0
        for e,a in cells:
            capacity = sum(v['capacity'][e] for v in chosen
                           if (v['parity'] is None or v['parity'] == e)
                           and all(a[i] == b for i,b in zip(v['high_subset'],v['labels'])))
            deficit += max(0,selected[e]-capacity)
            if not minimize:
                assert capacity >= selected[e]
        if minimize:
            assert deficit <= round(solver.ObjectiveValue())
            out['relaxation_uncovered'] = deficit
            out['selected_core_residual'] = selected
        out['placements'] = chosen
        out['verdict'] = 'relaxation feasible; not a covering system'
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--high',type=int,nargs='+',required=True)
    ap.add_argument('--pattern',required=True)
    ap.add_argument('--residual',type=int,nargs=2,help='required unless --core-frontier is used')
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--output',required=True)
    ap.add_argument('--minimize',action='store_true')
    ap.add_argument('--core-frontier',action='store_true',help='optimize over the certified 210 residual frontier')
    a = ap.parse_args()
    pattern = list(map(int,a.pattern))
    frontier = None
    if a.core_frontier:
        from pathlib import Path
        cert = json.loads((Path(__file__).resolve().parents[1]/'data/three-axis-210-certificate.json').read_text())
        frontier = next(row['residual_pair_frontier'] for row in cert['cases'] if row['pattern'] == pattern[:3])
    elif a.residual is None:
        ap.error('--residual is required without --core-frontier')
    result = solve(a.high,pattern,a.residual,a.time,a.workers,a.minimize,frontier)
    print(json.dumps({k:v for k,v in result.items() if k != 'placements'}),flush=True)
    with open(a.output,'w') as f:
        json.dump(result,f,indent=2)
