"""Separate lower-barrier recomputation with upward-padded capacity thresholds.

The padding creates a stochastically smaller comparison law. This differs
from the primary program's two-sided interpolation coefficient enclosures.
Imports only helpers from the earlier independent HALF-range audit.
"""
import argparse,hashlib,itertools,json,math,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from audit_stoploss_half_barrier import prime_list,small_configurations,step_polynomial

Q=10**8
H=5000
R=10**12
TERMS=18


def initial_cases(actions):
    primes=[2,3,5,7,11];rows=[]
    for atoms in itertools.product(*(actions[p] for p in (5,7,11))):
        ds=[F(0),F(0)]+[F(a,p) for a,p in zip(atoms,(5,7,11))]
        cost=F(0)
        for k in range(2,5):
            for _,divisors,prob in small_configurations(primes[:k],ds[:k]):
                load=F(sum(d*primes[k]>6 for d in divisors),primes[k])
                cost+=prob*max(F(0),(load-ds[k])/(1-ds[k]))
        law=[F(0)]*6
        for n,_,prob in small_configurations(primes,ds):law[n]+=prob
        den=math.lcm(*(x.denominator for x in law))
        rows.append(dict(atoms=list(atoms),deltas=[str(x) for x in ds[2:]],prefix_loss=str(cost),
                         prefix_lower=int(cost*Q),denominator=den,numerators=[int(x*den) for x in law]))
    return rows


def padded_mixtures(primes,denominators,width):
    thresholds=[0]
    for p,d in zip(primes,denominators):
        assert d<p
        thresholds.append(thresholds[-1]+(R*(p-d)+d*p-1)//(d*p))
    assert thresholds[-1]<4*R
    x=np.arange(width,dtype=np.int64)*(R//H)
    indices=np.searchsorted(np.array(thresholds,dtype=np.int64),x,side='right')-1
    factors=[]
    for z,t in zip(x,indices):
        t=int(t)
        if t==len(primes):factors.append(0);continue
        p,d=primes[t],denominators[t]
        f=min(Q,(int(z)-thresholds[t])*d*p*Q//(R*(p-d)))
        assert 0<=f<=Q
        factors.append(f)
    return indices,np.array(factors,dtype=np.int64)


def audit(pruning_path,pruning_audit_path,output):
    started=time.time();raw=Path(pruning_path).read_bytes();source=json.loads(raw)
    checked=json.loads(Path(pruning_audit_path).read_text())
    source_hash=hashlib.sha256(raw).hexdigest()
    assert checked['status']=='VERIFIED_FULL_RANGE_ACTION_DOMINANCE' and checked['source_sha256']==source_hash
    assert checked['eliminations_checked']==len(source['eliminations'])
    assert [p for p in prime_list(source['horizon_prime']) if p>=5]==source['primes']
    actions=dict(zip(source['primes'],source['remaining_atoms']))
    rows=initial_cases(actions)
    seeds=np.array([x['numerators'] for x in rows],dtype=np.int64)
    seed_den=np.array([x['denominator'] for x in rows],dtype=np.int64)
    primes=[p for p in prime_list(5003) if p>=13]
    upper_den=[p-actions[p][-1] for p in primes]
    sizes=[1];capacity=F(0);moment=F(1)
    for p,d in zip(primes,upper_den):
        capacity+=F(1,d)-F(1,p);sizes.append(int(H*capacity)+1);moment*=F(d+1,d)
    assert moment<64 and R%H==0
    maxnum=max(max(primes),max((H*p+d-1)//d for p,d in zip(primes,upper_den)))
    assert Q*Q<2**63 and maxnum*(2048+5003)*Q<2**63
    assert int(seed_den.max())*Q<2**63
    base=np.array([Q]+[0]*(TERMS-1),dtype=np.int64)
    total_states=TERMS+5;powers=np.array([1<<n for n in range(total_states)],dtype=np.int64)
    values=np.zeros((len(rows),sizes[-1]),dtype=np.int64)
    for k in range(len(primes)-1,-1,-1):
        p=primes[k];last_den=upper_den[k];width=sizes[k]
        ordered=sorted(upper_den[:k])
        upper=[base]
        for d in ordered:upper.append(step_polynomial(upper[-1],d,1))
        lower=[None]*(k+1);lower[k]=base
        for i in range(k-1,-1,-1):lower[i]=step_polynomial(lower[i+1],primes[i],1)
        laws=np.array([np.convolve(lower[t],upper[t])[:TERMS]//Q for t in range(k+1)])
        index,f=padded_mixtures(primes[:k],ordered,width)
        unknown=((Q-f[:,None])*laws[index]+f[:,None]*laws[np.minimum(index+1,k)])//Q
        assert np.all(unknown.sum(axis=1)<=Q)
        complete=np.zeros((len(rows),width,total_states),dtype=np.int64)
        for a in range(6):complete[:,:,a:a+TERMS]+=seeds[:,a,None,None]*unknown[None,:,:]
        complete//=seed_den[:,None,None]
        assert np.all(complete.sum(axis=2)<=Q)
        tailP=np.flip(np.cumsum(np.flip(complete,axis=2),axis=2),axis=2)
        tailA=np.flip(np.cumsum(np.flip(complete*powers,axis=2),axis=2),axis=2)
        assert int(tailA.max())<2048*Q
        def charge(a,b):
            n=0
            while n<total_states and a*(1<<n)<=p*(a-b):n+=1
            if n==total_states:return np.zeros((len(rows),width),dtype=np.int64)
            numerator=a*tailA[:,:,n]-(a-b)*p*tailP[:,:,n]
            assert np.all(numerator>=0)
            return numerator//(b*p)
        atoms=[];a=1
        while a<=p-last_den:atoms.append(a);a*=2
        assert atoms==actions[p][1:]
        new=np.full((len(rows),width),10000*Q,dtype=np.int64)
        for jump in range(H*(p-last_den)//(p*last_den)+1):
            low_a,low_b=H+p*jump,H
            high_a,high_b=(H+p*(jump+1),H) if last_den*(H+p*(jump+1))<=H*p else (p,last_den)
            now=np.minimum(charge(low_a,low_b),charge(high_a,high_b))
            for a in atoms:
                if low_a*(p-a)<=p*low_b and p*high_b<=high_a*(p-a):now=np.minimum(now,charge(p,p-a))
            new=np.minimum(new,now+values[:,jump:jump+width])
        values=new
        if k%100==0:print('k',k,'p',p,'seconds',round(time.time()-started,1),flush=True)
    for i,row in enumerate(rows):row['continuation_lower']=int(values[i,0]);row['total_lower']=row['prefix_lower']+int(values[i,0])
    bound=min(x['total_lower'] for x in rows)
    result=dict(status='VERIFIED_FULL_RANGE_BARRIER' if bound>Q else 'INCONCLUSIVE',
                scale=Q,grid_denominator=H,count_terms=TERMS,capacity_padding_scale=R,
                original_objective_horizon=source['horizon_prime'],relaxed_prefix_last_prime=5003,
                minimum_lower_integer=bound,minimum_lower=bound/Q,cases=rows,
                pruning_file=str(pruning_path),pruning_sha256=source_hash,pruning_audit_file=str(pruning_audit_path),
                seconds=round(time.time()-started,2))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('pruning');p.add_argument('pruning_audit');p.add_argument('output')
    args=p.parse_args();audit(args.pruning,args.pruning_audit,args.output)
