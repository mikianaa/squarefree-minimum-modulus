"""Numerical joint core-210 research: actual uncovered sets and tail envelopes.

The LP uses heuristic separation for thresholds 2..7. Its priced values
are NOT certified upper bounds. Thresholds >=8 are exact: (16-j)*max(q).
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation

DIVS=[d for d in range(1,211) if 210%d==0]
ALLOWED=[d for d in DIVS if d>6]


def make_prefix(output,seconds=30):
    model=cp_model.CpModel()
    z={(d,r):model.NewBoolVar(f'z{d}_{r}') for d in ALLOWED for r in range(d)}
    for d in ALLOWED:model.AddExactlyOne(z[d,r] for r in range(d))
    model.Add(z[7,0]==1);model.Add(z[10,0]==1)
    covered=[model.NewBoolVar(f'b{x}') for x in range(210)]
    for x,b in enumerate(covered):model.Add(b<=sum(z[d,x%d] for d in ALLOWED))
    model.Maximize(sum(covered))
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=seconds;solver.parameters.num_workers=4
    status=solver.Solve(model)
    assert status in [cp_model.OPTIMAL,cp_model.FEASIBLE]
    family={d:next(r for r in range(d) if solver.Value(z[d,r])) for d in ALLOWED}
    support=[x for x in range(210) if all(x%d!=r for d,r in family.items())]
    result=dict(family=family,support=support,uncovered=len(support),solver_status=solver.StatusName(status),
                solver_covered_bound=solver.BestObjectiveBound(),scope='Actual finite prefix; not a covering')
    Path(output).write_text(json.dumps(result,indent=2));print(result,flush=True)
    return result


def tail_basis(cutoff=100000):
    cache=Path(f'tmp/core210-tail-basis-{cutoff}.json')
    if cache.exists():return np.array(json.loads(cache.read_text())['G'])
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 11<=p<=cutoff]
    schedule=json.loads(Path('data/core30-case4-numerical-3m.json').read_text())
    byprime=dict(zip(schedule['primes'],schedule['deltas']))
    count=45;dist=np.zeros(count);dist[0]=1;G=np.zeros(17)
    for p in ps:
        r=1/(p*(1-byprime[p]))
        G+=np.maximum(0,1+r*(np.arange(17)[:,None]*2.**np.arange(count)-p))@dist
        dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
    k,_=quadratic_continuation(allps,cutoff)
    G+=k*(dist@4.**np.arange(count))*np.arange(17)**2
    cache.write_text(json.dumps(dict(G=G.tolist(),scope='numerical only',cutoff=cutoff)))
    return G


def price(points,weights,threshold,rng,random_starts=30,sweeps=6,ceiling=16):
    # All aligned catalogues, plus random ones. Batch coordinate ascent.
    residues=np.array([[x%d for d in DIVS] for x in range(210)]+[
        [int(rng.integers(d)) for d in DIVS] for _ in range(random_starts)])
    counts=sum((points[None,:]%d==residues[:,i,None]).astype(np.int8) for i,d in enumerate(DIVS))
    for sweep in range(sweeps):
        changes=0
        for i in rng.permutation(len(DIVS)):
            d=DIVS[i];rem=points%d
            counts-=(rem[None,:]==residues[:,i,None]).astype(np.int8)
            bins=(rem[None,:]+np.arange(len(residues))[:,None]*d).ravel()
            gain=weights[None,:]*(np.maximum(np.minimum(counts+1,ceiling)-threshold,0)-np.maximum(np.minimum(counts,ceiling)-threshold,0))
            scores=np.bincount(bins,weights=gain.ravel(),minlength=len(residues)*d).reshape(len(residues),d)
            best=np.argmax(scores,axis=1)
            changes+=np.count_nonzero(best!=residues[:,i]);residues[:,i]=best
            counts+=(rem[None,:]==best[:,None]).astype(np.int8)
        if not changes:break
    rows=np.maximum(np.minimum(counts,ceiling)-threshold,0)
    values=rows@weights
    indices=np.argsort(values)[-20:]
    return float(values.max()),rows[indices],residues[indices]


def optimize(prefix_file,output,iterations=70,cutoff=100000,basis_file=None,warm_file=None,cap11=False):
    start=time.time();prefix=json.loads(Path(prefix_file).read_text())
    pts=np.array(prefix['support']);G=tail_basis(cutoff) if basis_file is None else np.array(json.loads(Path(basis_file).read_text())['G'])
    fullG=G.copy();a11=None
    if cap11:
        schedule=json.loads(Path(basis_file or 'data/core30-case4-numerical-3m.json').read_text())
        delta=dict(zip(schedule['primes'],schedule['deltas']))[11]
        a11=round(11*delta);assert abs(11*delta-a11)<1e-7
        G=G-np.maximum(0,(np.arange(17)-a11)/(11-a11))
    c=np.r_[0,G[2]-G[1],np.diff(G[1:],2)]
    assert min(c)>-1e-10
    model=pywraplp.Solver.CreateSolver('GLOP')
    model.SetSolverSpecificParametersAsString('use_dual_simplex: true\nprimal_feasibility_tolerance: 1e-8\ndual_feasibility_tolerance: 1e-8')
    q=[model.NumVar(0,1,f'q{x}') for x in pts];model.Add(sum(q)==1)
    top=model.NumVar(0,1,'maxq')
    for v in q:model.Add(top>=v)
    caps=[]
    for d in DIVS[1:]:
        cap=model.NumVar(0,1,f'cap{d}');caps.append(cap)
        for r in range(d):
            inds=np.flatnonzero(pts%d==r)
            if len(inds):model.Add(cap>=sum(q[int(i)] for i in inds))
    hs={j:model.NumVar(0,16,f'h{j}') for j in range(2,8)}
    high=float(sum(c[j]*(16-j) for j in range(8,16)))
    capped=model.NumVar(0,11,'capped11') if cap11 else None
    model.Minimize(float(G[1])+float(c[1])*sum(caps)+sum(float(c[j])*hs[j] for j in hs)+high*top+(capped/(11-a11) if cap11 else 0))
    seen=set();pool=[]
    def add(j,row,residue,ceiling=16):
        key=(j,ceiling,tuple(int(v) for v in row))
        if key not in seen:
            seen.add(key);model.Add((capped if ceiling==11 else hs[j])>=sum(int(v)*z for v,z in zip(row,q)))
            pool.append(dict(threshold=j,ceiling=ceiling,residues=[int(v) for v in residue]))
    for x in range(210):
        residues=[x%d for d in DIVS]
        counts=sum((pts%d==x%d).astype(np.int8) for d in DIVS)
        for j in hs:add(j,np.maximum(counts-j,0),residues)
        if cap11:add(a11,np.maximum(np.minimum(counts,11)-a11,0),residues,11)
    if warm_file:
        warm=json.loads(Path(warm_file).read_text());assert warm['support']==prefix['support']
        for item in warm['pool']:
            residues=item['residues'];counts=sum((pts%d==r).astype(np.int8) for d,r in zip(DIVS,residues))
            if item.get('ceiling',16)==16 and item['threshold'] in hs:
                add(item['threshold'],np.maximum(counts-item['threshold'],0),residues)
            if cap11:add(a11,np.maximum(np.minimum(counts,11)-a11,0),residues,11)
    rng=np.random.default_rng(71810);history=[]
    for it in range(iterations):
        status=model.Solve()
        if status!=pywraplp.Solver.OPTIMAL:
            model.SetSolverSpecificParametersAsString('use_dual_simplex: false\nuse_preprocessing: false\nprimal_feasibility_tolerance: 1e-7\ndual_feasibility_tolerance: 1e-7')
            status=model.Solve()
        if status!=pywraplp.Solver.OPTIMAL:raise RuntimeError(f'LP status {status}; not an infeasibility conclusion')
        weights=np.array([v.solution_value() for v in q]);lower=model.Objective().Value()
        H=np.zeros(17);H[1]=sum(max(np.bincount(pts%d,weights=weights,minlength=d)) for d in DIVS[1:])
        H[0]=H[1]+1
        for j in range(8,17):H[j]=(16-j)*max(weights)
        for j in hs:
            value,rows,residues=price(pts,weights,j,rng)
            H[j]=value
            for row,residue in zip(rows,residues):add(j,row,residue)
        capvalue=0
        if cap11:
            capvalue,rows,residues=price(pts,weights,a11,rng,ceiling=11)
            for row,residue in zip(rows,residues):add(a11,row,residue,11)
        priced=float(G[1]+c@H[:16]+(capvalue/(11-a11) if cap11 else 0));history.append(dict(iteration=it,restricted_lp=lower,heuristic_priced=priced,rows=len(seen)))
        if it%5==0 or priced-lower<1e-8:print(it,lower,priced,'rows',len(seen),'seconds',round(time.time()-start,1),flush=True)
        if priced-lower<1e-8:break
    result=dict(prefix_file=prefix_file,family=prefix['family'],support=pts.tolist(),weights=weights.tolist(),
                envelope=H.tolist(),G=fullG.tolist(),coefficients=c.tolist(),restricted_lp=lower,
                heuristic_priced=priced,history=history,pool=pool,p_exact=cutoff,basis_file=basis_file,seconds=time.time()-start,
                cap11=cap11,cap11_atom=a11,cap11_stoploss=capvalue,
                scope='Heuristic pricing. Neither priced value nor convergence certifies an upper bound.')
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print('FINAL',lower,priced,'support',len(pts),flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prefix','optimize']);ap.add_argument('output')
    ap.add_argument('--prefix',default='data/core210-dense-prefix.json');ap.add_argument('--iterations',type=int,default=70)
    ap.add_argument('--p-exact',type=int,default=100000);ap.add_argument('--seconds',type=int,default=30)
    ap.add_argument('--basis')
    ap.add_argument('--warm');ap.add_argument('--cap11',action='store_true')
    a=ap.parse_args()
    if a.mode=='prefix':make_prefix(a.output,a.seconds)
    else:optimize(a.prefix,a.output,a.iterations,a.p_exact,a.basis,a.warm,a.cap11)
