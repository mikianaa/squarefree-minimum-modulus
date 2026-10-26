"""Numerical exact stop-loss hull of a finite catalogue pool.

The hull has noninteger atoms. These introduce control breakpoints absent
from the older integer-envelope optimization. The pool is not exhaustive.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation
from core210_joint_probe import DIVS


def hull(data):
    pts=np.array(data['support']);q=np.maximum(np.array(data['weights']),0);q/=sum(q)
    catalogues=[a['residues'] for a in data['pool']]+[[x%d for d in DIVS] for x in range(210)]
    catalogues.append([int(np.argmax(np.bincount(pts%d,weights=q,minlength=d))) for d in DIVS])
    hist=[]
    for residues in catalogues:
        counts=sum((pts%d==r).astype(int) for d,r in zip(DIVS,residues))
        hist.append(np.bincount(counts,weights=q,minlength=17))
    hist=np.array(hist);segments=[]
    for j in range(16):
        slopes=-hist[:,j+1:].sum(axis=1);intercepts=hist[:,j+1:]@np.arange(j+1,17)
        unique={}
        for m,b in zip(slopes,intercepts):
            key=round(float(m),13)
            if key not in unique or b>unique[key][1]:unique[key]=(float(m),float(b))
        lines=[];starts=[]
        for m,b in sorted(unique.values()):
            x=-float('inf')
            while lines:
                oldm,oldb=lines[-1]
                if abs(m-oldm)<1e-13:
                    if b<=oldb:break
                    lines.pop();starts.pop();continue
                x=(oldb-b)/(m-oldm)
                if x<=starts[-1]+1e-13:lines.pop();starts.pop()
                else:break
            if lines and abs(m-lines[-1][0])<1e-13:continue
            if not lines:x=-float('inf')
            lines.append((m,b));starts.append(x)
        for i,(m,b) in enumerate(lines):
            lo=max(float(j),starts[i]);hi=min(float(j+1),starts[i+1] if i+1<len(starts) else float('inf'))
            if hi>lo+1e-12:segments.append((lo,hi,m,b))
    atoms=[];weights=[];previous=-1.
    for lo,hi,m,b in segments:
        jump=m-previous
        assert jump>-1e-9
        if jump>1e-10:atoms.append(lo);weights.append(jump)
        previous=m
    if -previous>1e-10:atoms.append(16.);weights.append(-previous)
    atoms=np.array(atoms);weights=np.array(weights);weights/=sum(weights)
    assert abs(weights@atoms-max(hist@np.arange(17)))<1e-9
    for t in np.linspace(0,16,129):
        actual=max(hist@np.maximum(np.arange(17)-t,0));claimed=weights@np.maximum(atoms-t,0)
        assert abs(actual-claimed)<1e-9,(t,actual,claimed)
    return atoms,weights


def terminal(atoms,control=10000,states=45):
    grid=np.array([1,1.25,1.5,1.75,2,2.5,3,3.5,4,5,6,7,8,10,12,14,16.])
    cache=Path(f'tmp/core210-real-fixed-{control}-{states}.npz')
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text());mapping=dict(zip(old['primes'],old['deltas']))
    ps=primes_to(16000000)[:1000000];controlled=[p for p in ps if 11<=p<=control]
    if cache.exists():d=np.load(cache);V=d['values'];quad=float(d['quadratic'])
    else:
        k,_=quadratic_continuation(ps,3000000);V=np.zeros((len(grid),states));quad=k
        counts=grid[:,None]*2.**np.arange(states)
        for p in reversed([p for p in ps if control<p<=3000000]):
            a=round(p*mapping[p]);assert abs(a-p*mapping[p])<1e-7
            if a:assert a//(a&-a)<=7
            r=1/(p-a);quad*=1+3*r
            V=np.maximum(0,1+r*(counts-p))+(1-r)*V+r*np.c_[V[:,1:],np.zeros(len(grid))]
        np.savez(cache,values=V,quadratic=quad)
    interpolated=np.array([np.interp(atoms,grid,V[:,j]) for j in range(states)]).T
    interpolated+=quad*atoms[:,None]**2*4.**np.arange(states)
    return controlled,mapping,interpolated


def run(measure,initial,output,sweeps=5):
    start=time.time();data=json.loads(Path(measure).read_text());atoms,weights=hull(data)
    ps,mapping,basis=terminal(atoms);previous=json.loads(Path(initial).read_text());mapping.update(zip(previous['primes'],previous['deltas']))
    rs=np.array([1/(p*(1-mapping[p])) for p in ps]);states=basis.shape[1];end=weights@basis
    counts=atoms[:,None]*2.**np.arange(states)
    def cost(p,options):
        ans=np.zeros((len(options),states))
        for c,w in zip(counts,weights):ans+=w*np.maximum(0,1+options[:,None]*(c-p))
        return ans
    print('hull atoms',list(zip(atoms,weights)),flush=True)
    for sweep in range(sweeps):
        values=np.zeros((len(ps)+1,states));values[-1]=end
        for i in range(len(ps)-1,-1,-1):
            v=values[i+1];r=rs[i];values[i]=cost(ps[i],np.array([r]))[0]+(1-r)*v+r*np.r_[v[1:],0]
        dist=np.zeros(states);dist[0]=1;total=0
        for i,p in enumerate(ps):
            knots=(atoms[:,None]*2.**np.arange(int(math.log2(p))+1)).ravel();knots=knots[(knots>=0)&(knots<=p-1)]
            options=1/(p-np.unique(np.r_[0,knots,p-1]));v=values[i+1];charges=cost(p,options)
            objective=charges+(1-options[:,None])*v+options[:,None]*np.r_[v[1:],0]
            idx=np.argmin(objective@dist);r=options[idx];rs[i]=r
            total+=charges[idx]@dist;dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
        total+=dist@end
        print('continuous hull sweep',sweep,'start',values[0,0],'end',total,flush=True)
    out=dict(measure_file=measure,atoms=atoms.tolist(),weights=weights.tolist(),primes=ps,deltas=(1-1/(np.array(ps)*rs)).tolist(),
             numerical_total=float(total),seconds=time.time()-start,
             scope='Finite catalogue hull, numerical noninteger-breakpoint search; no global pricing or covering certificate')
    Path(output).write_text(json.dumps(out,separators=(',',':')));return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('measure');ap.add_argument('initial');ap.add_argument('output');ap.add_argument('--sweeps',type=int,default=5)
    a=ap.parse_args();run(a.measure,a.initial,a.output,a.sweeps)
