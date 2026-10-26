"""Actual stage-11 extension and a survivor-only core-2310 comparison.

All families and residual sets are actual. Envelope pricing and initial
measure optimization are heuristic, so scores are not proof certificates.
"""
import argparse,json,time,math
from pathlib import Path
import numpy as np
from ortools.linear_solver import pywraplp
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation


def extension(output,restarts=100):
    old=json.loads(Path('data/core210-hard-prefix.json').read_text());base={int(d):r for d,r in old['family'].items()}
    core=np.array(old['support']);points=np.array([x for x in range(2310) if x%210 in set(core)])
    divs=[d for d in range(1,211) if 210%d==0];rng=np.random.default_rng(1113210)
    best=None
    for trial in range(restarts):
        family={11*d:(0 if d==1 else int(rng.integers(11*d))) for d in divs}
        counts=sum((points%d==r).astype(int) for d,r in family.items())
        target=np.ones(len(points)) if trial%3==0 else rng.uniform(.8,1.2,len(points))
        for sweep in range(20):
            changed=False
            for d in rng.permutation([11*m for m in divs if m>1]):
                prev=family[d];counts-=(points%d==prev)
                score=np.bincount(points%d,weights=target*(counts==0),minlength=d)
                r=int(np.argmax(score));family[d]=r;counts+=(points%d==r);changed|=prev!=r
            if not changed:break
        remaining=points[counts==0].tolist()
        if best is None or len(remaining)<best['uncovered']:
            best=dict(family={**base,**{int(d):int(r) for d,r in family.items()}},extension={int(d):int(r) for d,r in family.items()},
                      support=remaining,uncovered=len(remaining),base_support=core.tolist(),period=2310,
                      scope='An actual stage-11 extension, not a minimum-residual theorem or a covering')
    Path(output).write_text(json.dumps(best,indent=2));print('extension residual',best['uncovered'],flush=True)
    return best


def tail_basis(cutoff=100000):
    cache=Path(f'tmp/core2310-tail-basis-{cutoff}.json')
    if cache.exists():return np.array(json.loads(cache.read_text())['G'])
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text());new=json.loads(Path('data/core210-hard2-continuous-hull-probe.json').read_text())
    ds=dict(zip(old['primes'],old['deltas']));ds.update(zip(new['primes'],new['deltas']))
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 13<=p<=cutoff];N=45
    dist=np.zeros(N);dist[0]=1;G=np.zeros(33);loads=np.arange(33)[:,None]*2.**np.arange(N)
    for p in ps:
        r=1/(p*(1-ds[p]));G+=np.maximum(0,1+r*(loads-p))@dist;dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
    k,_=quadratic_continuation(allps,cutoff);G+=k*(dist@4.**np.arange(N))*np.arange(33)**2
    cache.write_text(json.dumps(dict(G=G.tolist(),scope='numerical continuation')));return G


def price(points,weights,t,rng,starts=96,sweeps=4):
    period=2310;divs=[d for d in range(1,period+1) if period%d==0]
    leaders=points[np.argsort(weights)[-starts//2:]]
    seeds=np.r_[leaders,rng.choice(period,starts-len(leaders),replace=False)]
    residues=np.array([[x%d for d in divs] for x in seeds])
    residues=np.r_[residues,np.array([[int(rng.integers(d)) for d in divs] for _ in range(starts//4)])]
    counts=sum((points[None,:]%d==residues[:,i,None]).astype(np.int8) for i,d in enumerate(divs))
    for sweep in range(sweeps):
        changes=0
        for i in rng.permutation(len(divs)):
            d=divs[i];rem=points%d;counts-=(rem[None,:]==residues[:,i,None]).astype(np.int8)
            gain=weights[None,:]*(np.maximum(counts+1-t,0)-np.maximum(counts-t,0))
            bins=(rem[None,:]+np.arange(len(residues))[:,None]*d).ravel()
            scores=np.bincount(bins,weights=gain.ravel(),minlength=len(residues)*d).reshape(len(residues),d)
            best=np.argmax(scores,axis=1);changes+=np.count_nonzero(best!=residues[:,i]);residues[:,i]=best
            counts+=(rem[None,:]==best[:,None]).astype(np.int8)
        if not changes:break
    rows=np.maximum(counts-t,0);values=rows@weights;pick=np.argsort(values)[-8:]
    return float(values.max()),rows[pick],residues[pick]


def optimize(prefix_file,output,iterations=35,cutoff=100000):
    start=time.time();prefix=json.loads(Path(prefix_file).read_text());points=np.array(prefix['support'])
    base=sorted(int(x) for x in set(points%210));index={x:i for i,x in enumerate(base)};fi=np.array([index[x%210] for x in points]);sizes=np.bincount(fi)
    divs=[d for d in range(1,2311) if 2310%d==0];G=tail_basis(cutoff);c=np.r_[0,G[2]-G[1],np.diff(G[1:],2)]
    feature=lambda row:np.bincount(fi,weights=row,minlength=len(base))/sizes
    model=pywraplp.Solver.CreateSolver('CLP');q=[model.NumVar(0,1,f'q{x}') for x in base];model.Add(sum(q)==1)
    top=model.NumVar(0,1,'maxpoint')
    for v,n in zip(q,sizes):model.Add(top>=v/int(n))
    caps=[]
    for d in divs[1:]:
        cap=model.NumVar(0,1,f'cap{d}');caps.append(cap)
        for r in sorted(set(points%d)):
            row=feature(points%d==r);model.Add(cap>=sum(float(a)*v for a,v in zip(row,q) if a))
    hs={j:model.NumVar(0,32,f'h{j}') for j in range(2,16)}
    high=float(sum(c[j]*(32-j) for j in range(16,32)))
    model.Minimize(float(G[1])+float(c[1])*sum(caps)+sum(float(c[j])*hs[j] for j in hs)+high*top)
    seen=set();pool=[]
    def add(j,row,residues):
        f=feature(row);key=(j,tuple(np.round(f,12)))
        if key not in seen:
            seen.add(key);model.Add(hs[j]>=sum(float(a)*v for a,v in zip(f,q) if a));pool.append(dict(threshold=j,residues=[int(r) for r in residues]))
    for x in points[::max(1,len(points)//128)]:
        residues=[x%d for d in divs];counts=sum((points%d==r).astype(int) for d,r in zip(divs,residues))
        for j in hs:add(j,np.maximum(counts-j,0),residues)
    rng=np.random.default_rng(231013);history=[]
    def snapshot(fibre,weights,H,lower,priced):
        return dict(prefix_file=prefix_file,period=2310,support=prefix['support'],family=prefix['family'],base=base,
                    fibre_weights=fibre.tolist(),weights=weights.tolist(),envelope=H.tolist(),G=G.tolist(),
                    restricted_lp=lower,heuristic_priced=priced,history=history,pool=pool,p_exact=cutoff,seconds=time.time()-start,
                    scope='Heuristic core-2310 pricing. Optimizes fibre-constant survivor weights for one actual stage-11 family only; not a universal or certified bound.')
    for it in range(iterations):
        if model.Solve()!=pywraplp.Solver.OPTIMAL:raise RuntimeError('Numerical LP did not solve; no infeasibility inference')
        lower=model.Objective().Value();fibre=np.array([v.solution_value() for v in q]);weights=fibre[fi]/sizes[fi]
        H=np.zeros(33);H[1]=sum(max(np.bincount(points%d,weights=weights,minlength=d)) for d in divs[1:]);H[0]=H[1]+1
        for j in range(16,33):H[j]=(32-j)*max(weights)
        for j in hs:
            val,rows,residues=price(points,weights,j,rng);H[j]=val
            for row,residue in zip(rows,residues):add(j,row,residue)
        priced=float(G[1]+c@H[:32]);history.append(dict(iteration=it,restricted_lp=lower,heuristic_priced=priced,rows=len(seen)))
        print(it,lower,priced,'rows',len(seen),'seconds',round(time.time()-start,1),flush=True)
        if it%5==0:Path(output).write_text(json.dumps(snapshot(fibre,weights,H,lower,priced),separators=(',',':')))
        if priced-lower<1e-7:break
    result=snapshot(fibre,weights,H,lower,priced)
    Path(output).write_text(json.dumps(result,separators=(',',':')));return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['extension','optimize']);ap.add_argument('output');ap.add_argument('--prefix',default='data/core2310-overlap-prefix.json');ap.add_argument('--iterations',type=int,default=35);ap.add_argument('--p-exact',type=int,default=100000)
    a=ap.parse_args()
    if a.mode=='extension':extension(a.output)
    else:optimize(a.prefix,a.output,a.iterations,a.p_exact)
