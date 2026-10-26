"""Numerical conditional-comparison diagnostic retaining stage-11 fibre sizes.

The inequality is valid, but pricing used to evaluate it here is heuristic.
"""
import json,sys
from pathlib import Path
import numpy as np
from core210_joint_probe import DIVS,price


def mixture_price(points,w0,w1,t,rng,sweeps=6):
    seeds=np.arange(210);res=np.array([[x%d for d in DIVS] for x in seeds])
    res=np.r_[res,np.array([[int(rng.integers(d)) for d in DIVS] for _ in range(40)])]
    counts=sum((points[None,:]%d==res[:,i,None]).astype(np.int8) for i,d in enumerate(DIVS))
    for sweep in range(sweeps):
        changed=0
        for i in rng.permutation(16):
            d=DIVS[i];r=points%d;counts-=(r[None,:]==res[:,i,None]).astype(np.int8)
            gains=w0[None,:]*(np.maximum(counts+1-t,0)-np.maximum(counts-t,0))
            gains+=w1[None,:]*(np.maximum(2*counts+2-t,0)-np.maximum(2*counts-t,0))
            bins=(r[None,:]+np.arange(len(res))[:,None]*d).ravel()
            scores=np.bincount(bins,weights=gains.ravel(),minlength=len(res)*d).reshape(len(res),d)
            new=np.argmax(scores,axis=1);changed+=np.count_nonzero(new!=res[:,i]);res[:,i]=new
            counts+=(r[None,:]==new[:,None]).astype(np.int8)
        if not changed:break
    vals=np.maximum(counts-t,0)@w0+np.maximum(2*counts-t,0)@w1
    return float(max(vals))


def run(measure,output):
    data=json.loads(Path(measure).read_text());points=np.array(data['base']);q=np.array(data['fibre_weights'])
    support=np.array(data['support']);sizes=np.array([sum(support%210==x) for x in points]);r=1/sizes
    H=np.array(data['envelope']);new=H.copy();rng=np.random.default_rng(132310)
    for t in range(2,16):
        first=mixture_price(points,q*(1-r),q*r/2,t,rng)
        second=price(points,q*r,t/2,rng)[0]
        new[t]=first+second
    G=np.array(data['G']);c=np.r_[0,G[2]-G[1],np.diff(G[1:],2)]
    out=dict(heuristic_conditional_value=float(G[1]+c@new[:32]),direct_priced=data['heuristic_priced'],
             conditional_envelope=new.tolist(),scope='Diagnostic, not a certified upper bound; the conditional Jensen inequality itself is valid')
    Path(output).write_text(json.dumps(out,indent=2));print(out,flush=True)

if __name__=='__main__':run(*sys.argv[1:])
