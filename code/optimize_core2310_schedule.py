"""Numerical schedule search for an exact core envelope's power-of-two majorant.

This script exports a candidate, not a certificate. Its scope is one fixed
prefix. Do not pass the schedule to an unconditional cutoff-6 verifier.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation,compress_rounded_schedule


def run(envelope,output,sweeps=3,cutoff=3000000):
    start_time=time.time();core=json.loads(Path(envelope).read_text())
    allps=primes_to(16000000)[:1000000];ps=np.array([p for p in allps if 13<=p<=cutoff])
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text())
    new=json.loads(Path('data/core210-hard2-continuous-hull-probe.json').read_text())
    ds=dict(zip(old['primes'],old['deltas']));ds.update(zip(new['primes'],new['deltas']))
    r=1/(ps*(1-np.array([ds[int(p)] for p in ps])))
    K,later=quadratic_continuation(allps,cutoff);N=50;A=2.**np.arange(N);terminal=K*A*A
    seed=np.zeros(N);nums=core['power_seed_numerators'];seed[:len(nums)]=np.array(nums)/core['power_seed_denominator']
    history=[]
    for sweep in range(sweeps):
        V=np.zeros((len(ps)+1,N));V[-1]=terminal
        for i in range(len(ps)-1,-1,-1):
            v=V[i+1];V[i]=np.maximum(0,1+r[i]*(A-ps[i]))+(1-r[i])*v+r[i]*np.r_[v[1:],0.]
        initial=float(seed@V[0]);dist=seed.copy();cost=0.
        for i,p in enumerate(ps):
            # All hinges in the exported admissible interval, including endpoints.
            deltas=np.unique(np.r_[0.,A[A<=p/2]/p,.5]);options=1/(p*(1-deltas));v=V[i+1]
            charges=np.maximum(0,1+options[:,None]*(A-p))
            values=charges+(1-options[:,None])*v+options[:,None]*np.r_[v[1:],0.]
            pick=np.argmin(values@dist);r[i]=options[pick];cost+=float(charges[pick]@dist)
            dist=(1-r[i])*dist+r[i]*np.r_[0.,dist[:-1]]
        total=cost+float(dist@terminal);history.append(dict(sweep=sweep,initial=initial,total=total))
        print(history[-1], 'seconds',round(time.time()-start_time,1),flush=True)
    raw=[0.]*5+(1-1/(ps*r)).tolist()+later;edges,params=compress_rounded_schedule(raw)
    result=dict(comparison='fixed-core2310',envelope_file=envelope,seed_numerators=nums,seed_denominator=core['power_seed_denominator'],
                edges=edges,params=params,p_exact=cutoff,n_stages=1000000,history=history,
                numerical_total_before_rounding=total,scope='Numerical candidate for one fixed core, not an upper-bound-6 certificate')
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    # Reprice every possible initial integer load using the rounded schedule.
    used=np.array(params)[np.searchsorted(edges,np.arange(6,len(ps)+6))]
    dist=np.zeros(N);dist[0]=1;loads=np.arange(33)[:,None]*A;G=np.zeros(33)
    for p,d in zip(ps,used):
        prob=1/(p*(1-d));G+=np.maximum(0,1+prob*(loads-p))@dist
        dist=(1-prob)*dist+prob*np.r_[0.,dist[:-1]]
    G+=K*float(dist@(A*A))*np.arange(33)**2
    result['G']=G.tolist();result['numerical_rounded_exact_prefix']=float(seed[:6]@G[[1,2,4,8,16,32]])
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print('rounded-prefix total',result['numerical_rounded_exact_prefix'],flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('envelope');ap.add_argument('output');ap.add_argument('--sweeps',type=int,default=3)
    a=ap.parse_args();run(a.envelope,a.output,a.sweeps)
