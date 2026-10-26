"""Integer lower certificate: cutoff 6 fails for all schedules 0<=delta<=1/2.

This does NOT exclude larger distortions, a different comparison distribution,
or bounds retaining cross-stage compatibility. The proof is a finite-horizon
majorization relaxation, with all state and arithmetic rounding DOWN.
No optimization solver or floating-point arithmetic is used in the bound.
NumPy int64 operations are used only where explicit overflow bounds apply.
"""
import argparse
import bisect
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from verify_squarefree import primes_to

SCALE = 10**9
MIX_SCALE = 10**9


def initial_cases():
    options = [[F(0), F(1,5), F(1,2)],
               [F(0), F(1,7), F(2,7), F(1,2)],
               [F(0), F(1,11), F(2,11), F(4,11), F(1,2)]]
    cases = []
    for ds in itertools.product(*options):
        w = [F(1,3), F(1,2), F(1,6)]
        loss = F(0)
        for p,d in zip((5,7,11),ds):
            shift = int(p == 5)
            loss += sum(prob*max(F(0), (F(2**n-shift,p)-d)/(1-d)) for n,prob in enumerate(w))
            r = 1/(p*(1-d))
            new = [F(0)]*(len(w)+1)
            for n,prob in enumerate(w):
                new[n] += prob*(1-r)
                new[n+1] += prob*r
            w = new
        den = math.lcm(*(v.denominator for v in w))
        nums = [int(v*den) for v in w]
        cases.append(dict(deltas=[str(d) for d in ds], seed_numerators=nums,
                          seed_denominator=den, prefix_loss=str(loss),
                          prefix_lower=loss.numerator*SCALE//loss.denominator))
    return cases


def bernoulli_lower(w,p,numerator):
    out = w*(p-numerator)
    out[1:] += w[:-1]*numerator
    return out//p


def grid_metadata(primes,H):
    caps = [F(0)]
    for p in primes: caps.append(caps[-1]+F(1,p))
    sizes = [int(x*H)+1 for x in caps]
    indices=[]; flo=[]; fhi=[]
    for j in range(sizes[-1]):
        x=F(j,H)
        t=bisect.bisect_right(caps,x)-1
        f=(x-caps[t])*primes[t] if t<len(primes) else F(0)
        assert 0<=f<=1
        lower=f.numerator*MIX_SCALE//f.denominator
        upper=-((-f.numerator*MIX_SCALE)//f.denominator)
        indices.append(t);flo.append(lower);fhi.append(upper)
    return sizes,np.array(indices),np.array(flo,dtype=np.int64),np.array(fhi,dtype=np.int64)


def certify(end=5003,H=4000,states=20,output=None):
    started=time.time()
    assert end<=10000 and H<=10000 and states<=25
    primes=[p for p in primes_to(end) if p>=13]
    assert primes
    cases=initial_cases()
    seeds=np.array([c['seed_numerators'] for c in cases],dtype=np.int64)
    denominators=np.array([c['seed_denominator'] for c in cases],dtype=np.int64)
    assert max(denominators)<100000
    sizes,indices,flo,fhi=grid_metadata(primes,H)
    upper_prefix=[np.r_[np.int64(SCALE),np.zeros(states-1,dtype=np.int64)]]
    moment=F(1)
    for p in primes:
        upper_prefix.append(bernoulli_lower(upper_prefix[-1],p,2))
        moment*=F(p+2,p)
    # Bounds used to guarantee int64 safety in every convolution and cost.
    assert moment<64 and SCALE*SCALE<2**63
    assert 4*max(H,end)*64*SCALE*32<2**63
    assert 500*SCALE*max(denominators)<2**63
    powers=np.array([1<<n for n in range(states)],dtype=np.int64)
    value=np.zeros((len(cases),sizes[-1]),dtype=np.int64)
    checkpoints=[]
    for k in range(len(primes)-1,-1,-1):
        p=primes[k];m=sizes[k]
        suffix=[None]*(k+1)
        suffix[k]=upper_prefix[0]
        for t in range(k-1,-1,-1):
            suffix[t]=bernoulli_lower(suffix[t+1],primes[t],1)
        extrema=np.array([np.convolve(upper_prefix[t],suffix[t])[:states]//SCALE for t in range(k+1)])
        assert np.all(extrema>=0) and np.all(extrema.sum(axis=1)<=SCALE)
        ts=indices[:m]
        left=extrema[ts]
        right=extrema[np.minimum(ts+1,k)]
        # Both interpolation weights are rounded down; their sum is <=1.
        w=((MIX_SCALE-fhi[:m,None])*left+flo[:m,None]*right)//MIX_SCALE
        tailP=np.cumsum(w[:,::-1],axis=1)[:,::-1]
        tailA=np.cumsum((w*powers)[:,::-1],axis=1)[:,::-1]
        assert int(tailP.max())<=SCALE and int(tailA.max())<=64*SCALE

        def cost(nu_num,nu_den):
            conditional=[]
            for a in range(seeds.shape[1]):
                n=0
                while n<states and (1<<(n+a))*nu_num<=p*(nu_num-nu_den): n+=1
                if n==states:
                    conditional.append(np.zeros(m,dtype=np.int64))
                else:
                    numerator=nu_num*(1<<a)*tailA[:,n]-(nu_num-nu_den)*p*tailP[:,n]
                    assert np.all(numerator>=0)
                    conditional.append(numerator//(nu_den*p))
            conditional=np.array(conditional)
            assert int(conditional.max())<=500*SCALE
            return (seeds@conditional)//denominators[:,None]

        knot_bins={}
        a=1
        while 2*a<=p:
            q=H*a//(p*(p-a))
            knot_bins.setdefault(q,[]).append(a)
            a*=2
        out=np.full((len(cases),m),10000*SCALE,dtype=np.int64)
        low=cost(1,1)
        for q in range(H//p+1):
            high=cost(H+p*(q+1),H) if p*(q+1)<=H else cost(2,1)
            here=np.minimum(low,high)
            for a in knot_bins.get(q,[]): here=np.minimum(here,cost(p,p-a))
            out=np.minimum(out,here+value[:,q:q+m])
            low=high
        value=out
        if k%50==0:
            row=dict(k=k,p=p,states=m,smallest_continuation=int(value[:,0].min()),seconds=round(time.time()-started,2))
            checkpoints.append(row)
            print(json.dumps(row),flush=True)
    for i,case in enumerate(cases):
        case['continuation_lower']=int(value[i,0])
        case['total_lower']=case['prefix_lower']+case['continuation_lower']
    minimum=min(c['total_lower'] for c in cases)
    result=dict(status='CERTIFIED_HALF_RANGE_BARRIER' if minimum>SCALE else 'INCONCLUSIVE',
                scope='cutoff 6; original first/stop/quadratic comparison; all prime schedules 0<=delta<=1/2',
                excluded_extensions=['delta>1/2','stronger or conditional comparison bounds'],
                end_prime=primes[-1],grid_denominator=H,count_terms=states,scale=SCALE,mix_scale=MIX_SCALE,
                minimum_lower_integer=minimum,minimum_lower=minimum/SCALE,
                cases=cases,checkpoints=checkpoints,seconds=round(time.time()-started,2))
    if output:Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','checkpoints')},indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    parser.add_argument('--end',type=int,default=5003)
    parser.add_argument('--grid',type=int,default=4000)
    args=parser.parse_args()
    certify(args.end,args.grid,output=args.output)
