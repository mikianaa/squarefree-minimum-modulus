"""Numerical survivor-mass-aware union estimate.

For survivor mass m and stage loss Y in [0,1],
 E[1_U Y] <= inf_t (m*t + E(Y-t)+).
The positive-part term admits the same increasing-convex comparison.
The finite-pool core law used here is numerical, not globally certified.
"""
import json,time,sys
from pathlib import Path
import numpy as np
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation


def run(output):
    start=time.time();core=json.loads(Path('data/core210-hard2-continuous-hull-probe.json').read_text())
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text())
    schedule=dict(zip(old['primes'],old['deltas']));schedule.update(zip(core['primes'],core['deltas']))
    f=np.array(core['atoms']);w=np.array(core['weights']);N=50
    loads=(f[:,None]*2.**np.arange(N)).ravel();order=np.argsort(-loads);loads=loads[order]
    atom_indices=np.repeat(np.arange(len(f)),N)[order];count_indices=np.tile(np.arange(N),len(f))[order]
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 11<=p<=3000000]
    law=np.zeros(N);law[0]=1.;m=1.;oldsum=0.;credit=0.;rows=[];first_failure=None
    for p in ps:
        delta=schedule[p];r=1/(p*(1-delta));weights=w[atom_indices]*law[count_indices]
        losses=np.maximum(0,1+r*(loads-p));ordinary=float(weights@losses);oldsum+=ordinary
        if m>0:
            cumulative=np.cumsum(weights);before=np.r_[0,cumulative[:-1]]
            allocated=np.maximum(0,np.minimum(weights,m-before))
            newloss=float(allocated@losses)
            # Y<=1 separately bounds the mass removed by m.
            newloss=min(m,newloss);gain=ordinary-newloss;credit+=gain;m-=newloss
            if m<=1e-18 and first_failure is None:first_failure=p
            if gain>1e-10 or p in [11,13,17,101,1009,10007]:
                rows.append(dict(p=p,ordinary=ordinary,union_charge=newloss,credit=gain,survivor=m))
        law=(1-r)*law+r*np.r_[0,law[:-1]]
    k,_=quadratic_continuation(allps,3000000)
    tail=k*(law@4.**np.arange(N))*(w@f**2)
    result=dict(ordinary_total=oldsum+tail,survivor_before_tail=m,quadratic_remaining=tail,
                survivor_after_tail=m-tail,first_failure_prime=first_failure,credit_before_failure=credit,
                rows=rows,seconds=time.time()-start,
                scope='Numerical probe using an incompletely priced core law; no covering certificate')
    Path(output).write_text(json.dumps(result,indent=2));print({k:v for k,v in result.items() if k!='rows'},flush=True)
    return result

if __name__=='__main__':run(sys.argv[1])
