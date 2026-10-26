"""Fast numerical core/tail alternation with controls through a chosen prime.

The later continuation is cached as 16 separate core-load value functions.
Neither floating arithmetic nor heuristic core pricing is a certificate.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation
from optimize_core210_tail import envelope,cost
from core210_joint_probe import optimize as optimize_measure
from core210_joint_probe import DIVS


def setup(control_cutoff=10000,exact_cutoff=3000000,states=45):
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text());oldmap=dict(zip(old['primes'],old['deltas']))
    cache=Path(f'tmp/core210-fixed-basis-{control_cutoff}-{exact_cutoff}-{states}.npy')
    allps=primes_to(16000000)[:1000000]
    ps=[p for p in allps if 11<=p<=control_cutoff]
    if cache.exists():terminal=np.load(cache)
    else:
        k,_=quadratic_continuation(allps,exact_cutoff)
        counts=np.arange(1,17)[:,None]*2.**np.arange(states)
        terminal=k*counts**2
        later=[p for p in allps if control_cutoff<p<=exact_cutoff]
        for p in reversed(later):
            r=1/(p*(1-oldmap[p]))
            terminal=np.maximum(0,1+r*(counts-p))+(1-r)*terminal+r*np.c_[terminal[:,1:],np.zeros(16)]
        np.save(cache,terminal)
    rs=np.array([1/(p*(1-oldmap[p])) for p in ps])
    return ps,rs,terminal


def refine(measure,output,control_cutoff=10000,sweeps=6,initial_basis=None):
    start=time.time();data=json.loads(Path(measure).read_text());law=envelope(data)
    ps,rs,basis=setup(control_cutoff);states=basis.shape[1];terminal=law@basis
    if initial_basis:
        prev=json.loads(Path(initial_basis).read_text());prevmap=dict(zip(prev['primes'],prev['deltas']))
        rs=np.array([1/(p*(1-prevmap[p])) for p in ps])
    active=np.flatnonzero(law>1e-10)+1
    capped_H=None
    if data.get('cap11'):
        pts=np.array(data['support']);w=np.maximum(0,np.array(data['weights']));w/=sum(w);capped_H=np.zeros(12)
        for item in data['pool']:
            counts=sum((pts%d==r).astype(int) for d,r in zip(DIVS,item['residues']))
            capped_H=np.maximum(capped_H,np.maximum(np.minimum(counts[:,None],11)-np.arange(12),0).T@w)
    for sweep in range(sweeps):
        V=np.zeros((len(ps)+1,states));V[-1]=terminal
        for i in range(len(ps)-1,-1,-1):
            r=rs[i];v=V[i+1];V[i]=cost(ps[i],np.array([r]),states,law)[0]+(1-r)*v+r*np.r_[v[1:],0]
        dist=np.zeros(states);dist[0]=1;total=0
        for i,p in enumerate(ps):
            atoms=np.unique((active[:,None]*2.**np.arange(int(math.log2(p))+1)).ravel());atoms=atoms[atoms<p]
            opts=1/(p-np.r_[0,atoms]);v=V[i+1];charges=cost(p,opts,states,law)
            if p==11 and capped_H is not None:
                opts=1/(11-np.arange(11));charges=cost(p,opts,states,law)
                charges[:,0]=capped_H[:11]/(11-np.arange(11))
            values=charges+(1-opts[:,None])*v+opts[:,None]*np.r_[v[1:],0]
            best=np.argmin(values@dist);r=opts[best];rs[i]=r
            total+=charges[best]@dist;dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
        total+=dist@terminal
    G=np.zeros(17);dist=np.zeros(states);dist[0]=1
    for p,r in zip(ps,rs):
        G+=np.maximum(0,1+r*(np.arange(17)[:,None]*2.**np.arange(states)-p))@dist
        dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
    G[1:]+=basis@dist
    out=dict(measure_file=measure,G=G.tolist(),law=law.tolist(),primes=ps,deltas=(1-1/(np.array(ps)*rs)).tolist(),
             control_cutoff=control_cutoff,p_exact=3000000,numerical_total=float(total),seconds=time.time()-start,
             scope='Floating, heuristic core envelope; fixed controls after control_cutoff from core30 case4')
    Path(output).write_text(json.dumps(out,separators=(',',':')));print('REFINED',float(total),'seconds',time.time()-start,flush=True)
    return out


def alternate(prefix,measure,stem,rounds=4,cap11=False):
    current=measure;previous=None;rows=[]
    for i in range(rounds):
        out=f'{stem}-tail-{i}.json';next_measure=f'{stem}-measure-{i}.json'
        r=refine(current,out,initial_basis=previous)
        m=optimize_measure(prefix,next_measure,iterations=70,cutoff=3000000,basis_file=out,warm_file=current,cap11=cap11)
        rows.append(dict(round=i,tail_value=r['numerical_total'],restricted_lp=m['restricted_lp'],priced=m['heuristic_priced']))
        Path(f'{stem}-summary.json').write_text(json.dumps(rows,indent=2))
        current=next_measure;previous=out
    return rows

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('prefix');ap.add_argument('measure');ap.add_argument('stem');ap.add_argument('--rounds',type=int,default=4);ap.add_argument('--cap11',action='store_true')
    a=ap.parse_args();alternate(a.prefix,a.measure,a.stem,a.rounds,a.cap11)
