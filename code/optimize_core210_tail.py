"""Numerical continuation search for a priced core-210 envelope.

Input catalogue pricing is heuristic. Outputs are not proof certificates.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation
from core210_joint_probe import DIVS


def envelope(data):
    pts=np.array(data['support']);w=np.maximum(0,np.array(data['weights']));w/=sum(w)
    catalogues=[a['residues'] for a in data['pool']]
    catalogues.extend([[x%d for d in DIVS] for x in range(210)])
    catalogues.append([int(np.argmax(np.bincount(pts%d,weights=w,minlength=d))) for d in DIVS])
    H=np.zeros(17)
    for residues in catalogues:
        counts=sum((pts%d==r).astype(int) for d,r in zip(DIVS,residues))
        H=np.maximum(H,np.maximum(counts[:,None]-np.arange(17),0).T@w)
    law=np.r_[np.diff(H,2),H[-2]]
    assert min(law)>-1e-8
    law=np.maximum(law,0);law/=sum(law)
    return law


def cost(p,rs,states,law):
    out=np.zeros((len(rs),states));twos=2.**np.arange(states)
    for i,w in enumerate(law,1):
        if w>1e-12:out+=w*np.maximum(0,1+rs[:,None]*(i*twos-p))
    return out


def optimize(measure,output,cutoff=100000,sweeps=4):
    start=time.time();data=json.loads(Path(measure).read_text());law=envelope(data)
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 11<=p<=cutoff];states=40
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text());oldmap=dict(zip(old['primes'],old['deltas']))
    rs=np.array([1/(p*(1-oldmap[p])) for p in ps]);k,far=quadratic_continuation(allps,cutoff)
    terminal=k*(law@np.arange(1,17)**2)*4.**np.arange(states)
    active=np.flatnonzero(law>1e-10)+1
    for sweep in range(sweeps):
        V=np.zeros((len(ps)+1,states));V[-1]=terminal
        for i in range(len(ps)-1,-1,-1):
            v=V[i+1];r=rs[i];V[i]=cost(ps[i],np.array([r]),states,law)[0]+(1-r)*v+r*np.r_[v[1:],0]
        dist=np.zeros(states);dist[0]=1;total=0
        print('sweep',sweep,'start',V[0,0],flush=True)
        for i,p in enumerate(ps):
            atoms=np.unique((active[:,None]*2.**np.arange(int(math.log2(p))+1)).ravel());atoms=atoms[atoms<p]
            opts=1/(p-np.r_[0,atoms]);v=V[i+1];charges=cost(p,opts,states,law)
            q=charges+(1-opts[:,None])*v+opts[:,None]*np.r_[v[1:],0]
            best=np.argmin(q@dist);r=opts[best];rs[i]=r
            total+=charges[best]@dist;dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
        total+=dist@terminal
    G=np.zeros(17);dist=np.zeros(states);dist[0]=1
    for p,r in zip(ps,rs):
        G+=np.maximum(0,1+r*(np.arange(17)[:,None]*2.**np.arange(states)-p))@dist
        dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
    G+=k*(dist@4.**np.arange(states))*np.arange(17)**2
    result=dict(measure_file=measure,law=law.tolist(),G=G.tolist(),primes=ps,deltas=(1-1/(np.array(ps)*rs)).tolist(),
                p_exact=cutoff,numerical_total=float(total),seconds=time.time()-start,
                scope='Numerical candidate based on heuristic core envelope; not a certificate')
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print('FINAL',total,'seconds',time.time()-start,flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('measure');ap.add_argument('output');ap.add_argument('--p-exact',type=int,default=100000);ap.add_argument('--sweeps',type=int,default=4)
    a=ap.parse_args();optimize(a.measure,a.output,a.p_exact,a.sweeps)
