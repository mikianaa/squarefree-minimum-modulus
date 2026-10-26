"""Separate integer recomputation of the half-range stop-loss barrier.

Uses a different grid (1/5000), probability scale (10^8), count truncation
(18 unknown successes), explicit small-prime subset enumeration, and full
convolution with the initial distribution before evaluating costs. Does not
import the primary certificate generator or the covering verifiers.
"""
import argparse
from bisect import bisect_right
from fractions import Fraction
from itertools import product
import json
from math import isqrt,lcm
from pathlib import Path
import time
import numpy as np

Q=10**8
H=5000
TERMS=18


def prime_list(n):
    flags=[True]*(n+1)
    flags[0]=flags[1]=False
    for d in range(2,isqrt(n)+1):
        if flags[d]:
            for multiple in range(d*d,n+1,d): flags[multiple]=False
    return [i for i,yes in enumerate(flags) if yes]


def small_configurations(primes,deltas):
    out=[]
    for bits in product((0,1),repeat=len(primes)):
        probability=Fraction(1)
        divisors=[1]
        for p,d,b in zip(primes,deltas,bits):
            r=1/(p*(1-d))
            probability*=r if b else 1-r
            if b:divisors += [p*a for a in divisors]
        out.append((sum(bits),divisors,probability))
    return out


def enumerate_prefixes():
    primes=[2,3,5,7,11]
    option_lists=[]
    for k in (2,3,4):
        p=primes[k]
        values={Fraction(0),Fraction(1,2)}
        for _,divisors,_ in small_configurations(primes[:k],[Fraction(0)]*k):
            s=Fraction(sum(d*p>6 for d in divisors),p)
            if 0<=s<=Fraction(1,2):values.add(s)
        option_lists.append(sorted(values))
    results=[]
    for ds in product(*option_lists):
        deltas=[Fraction(0),Fraction(0)]+list(ds)
        cost=Fraction(0)
        for k in (2,3,4):
            for _,divisors,w in small_configurations(primes[:k],deltas[:k]):
                s=Fraction(sum(d*primes[k]>6 for d in divisors),primes[k])
                cost+=w*max(Fraction(0),(s-deltas[k])/(1-deltas[k]))
        probabilities=[Fraction(0)]*6
        for n,_,w in small_configurations(primes,deltas):probabilities[n]+=w
        denominator=lcm(*(w.denominator for w in probabilities))
        results.append(dict(deltas=[str(d) for d in ds],prefix_loss=str(cost),
                            prefix_lower=int(cost*Q),denominator=denominator,
                            numerators=[int(w*denominator) for w in probabilities]))
    assert len(results)==60
    return results


def step_polynomial(v,p,a):
    u=np.zeros(len(v),dtype=np.int64)
    for n in range(len(v)):
        u[n]=((p-a)*int(v[n])+(a*int(v[n-1]) if n else 0))//p
    return u


def run(output):
    start=time.time()
    cases=enumerate_prefixes()
    seeds=np.array([x['numerators'] for x in cases],dtype=np.int64)
    den=np.array([x['denominator'] for x in cases],dtype=np.int64)
    primes=[p for p in prime_list(5003) if p>=13]
    caps=[Fraction(0)]
    for p in primes:caps.append(caps[-1]+Fraction(1,p))
    lengths=[int(H*c)+1 for c in caps]
    tindex=[];mixleft=[];mixright=[]
    for j in range(lengths[-1]):
        x=Fraction(j,H);t=bisect_right(caps,x)-1
        f=(x-caps[t])*primes[t] if t<len(primes) else Fraction(0)
        tindex.append(t)
        mixleft.append(int((1-f)*Q));mixright.append(int(f*Q))
    tindex=np.array(tindex)
    mixleft=np.array(mixleft,dtype=np.int64);mixright=np.array(mixright,dtype=np.int64)
    up=[np.array([Q]+[0]*(TERMS-1),dtype=np.int64)]
    for p in primes:up.append(step_polynomial(up[-1],p,2))
    current=np.zeros((60,lengths[-1]),dtype=np.int64)
    full_terms=TERMS+5
    powers=np.array([2**n for n in range(full_terms)],dtype=np.int64)
    # All arrays below are nonnegative probability subdistributions.
    assert Q*Q<2**63 and int(den.max())*Q<2**63
    for k in range(len(primes)-1,-1,-1):
        p=primes[k];width=lengths[k]
        down=[None]*(k+1);down[k]=up[0]
        for i in range(k-1,-1,-1):down[i]=step_polynomial(down[i+1],primes[i],1)
        laws=[]
        for i in range(k+1):
            laws.append(np.convolve(down[i],up[i])[:TERMS]//Q)
        laws=np.array(laws)
        t=tindex[:width]
        unknown=(laws[t]*mixleft[:width,None]+laws[np.minimum(t+1,k)]*mixright[:width,None])//Q
        assert np.all(unknown.sum(axis=1)<=Q)
        combined=np.zeros((60,width,full_terms),dtype=np.int64)
        for a in range(6):
            combined[:,:,a:a+TERMS]+=seeds[:,a,None,None]*unknown[None,:,:]
        combined//=den[:,None,None]
        assert np.all(combined.sum(axis=2)<=Q)
        tail_probability=np.flip(np.cumsum(np.flip(combined,axis=2),axis=2),axis=2)
        tail_power=np.flip(np.cumsum(np.flip(combined*powers,axis=2),axis=2),axis=2)
        # The full mean E[2^N] is <64*32; this also bounds all cost intermediates.
        assert int(tail_power.max())<2048*Q
        assert 4*max(H,5003)*2048*Q<2**63

        def evaluate(a,b):
            n=0
            while n<full_terms and a*(1<<n)<=p*(a-b):n+=1
            if n==full_terms:return np.zeros((60,width),dtype=np.int64)
            numerator=a*tail_power[:,:,n]-(a-b)*p*tail_probability[:,:,n]
            assert np.all(numerator>=0)
            return numerator//(b*p)

        atoms=[];atom=1
        while 2*atom<=p:
            atoms.append(atom);atom*=2
        next_value=np.full((60,width),10000*Q,dtype=np.int64)
        for jump in range(H//p+1):
            low_num,low_den=H+p*jump,H
            high_num,high_den=(H+p*(jump+1),H) if p*(jump+1)<=H else (2,1)
            local=np.minimum(evaluate(low_num,low_den),evaluate(high_num,high_den))
            for atom in atoms:
                # Include exactly the hinge breakpoints lying in this closed control interval.
                if low_num*(p-atom)<=p*low_den and p*high_den<=high_num*(p-atom):
                    local=np.minimum(local,evaluate(p,p-atom))
            next_value=np.minimum(next_value,local+current[:,jump:jump+width])
        current=next_value
        if k%100==0:print('k',k,'p',p,'seconds',round(time.time()-start,1),flush=True)
    for i,row in enumerate(cases):
        row['continuation_lower']=int(current[i,0])
        row['total_lower']=row['prefix_lower']+row['continuation_lower']
    minimum=min(row['total_lower'] for row in cases)
    result=dict(status='VERIFIED_HALF_RANGE_BARRIER' if minimum>Q else 'INCONCLUSIVE',
                scale=Q,grid_denominator=H,unknown_count_terms=TERMS,end_prime=5003,
                minimum_lower_integer=minimum,minimum_lower=minimum/Q,cases=cases,
                scope='original comparison, cutoff 6, all schedules with 0<=delta<=1/2',
                seconds=round(time.time()-start,2))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    run(parser.parse_args().output)
