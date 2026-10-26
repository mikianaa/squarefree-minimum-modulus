"""Exact dominance eliminations for the unrestricted stop-loss comparison.

All lower bounds use signed integer floor division. A removed largest atom
can be replaced by the preceding one without increasing the FULL finite
objective through horizon_prime, for any other currently allowed controls.
This file alone is not a minimum-modulus barrier.
"""
import argparse,json,time
from pathlib import Path
from verify_squarefree import primes_to

SCALE=10**12


def initial_atoms(p):
    if p==5:return [0,1,3]
    out=[0];a=1
    while a<p:out.append(a);a*=2
    return out


def future_marginal_bounds(primes,atoms,cap,scale):
    values=[[0]*(cap+1) for _ in range(len(primes)+1)]
    for k in range(len(primes)-1,-1,-1):
        p=primes[k];future=values[k+1]
        shift=int(p==5)
        for n in range(cap+1):
            A=(1<<n)-shift;B=(1<<(n+1))-shift
            following=future[min(n+1,cap)]
            values[k][n]=min(((max(0,B-a)-max(0,A-a))*scale+(p-a-1)*future[n]+following)//(p-a) for a in atoms[k])
    return values


def dominance_lower(primes,atoms,k,future,cap,scale):
    p=primes[k];a=atoms[k][-1];shift=int(p==5)
    values=[]
    for n in range(cap+1):
        A=(1<<n)-shift
        values.append((A-p)*scale*(A>=a)+future[n])
    for j in range(k-1,-1,-1):
        low_den=primes[j];high_den=primes[j]-atoms[j][-1]
        values=[min(((low_den-1)*values[n]+values[min(n+1,cap)])//low_den,
                    ((high_den-1)*values[n]+values[min(n+1,cap)])//high_den) for n in range(cap+1)]
    return (2*values[0]+3*values[1]+values[2])//6


def certify(output,horizon=10007,cap=24,max_passes=4):
    started=time.time();primes=[p for p in primes_to(horizon) if p>=5]
    atoms=[initial_atoms(p) for p in primes]
    assert (1<<cap)>2*primes[-1]
    evidence=[];passes=[]
    for sweep in range(max_passes):
        future=future_marginal_bounds(primes,atoms,cap,SCALE)
        changed=0
        for k,p in enumerate(primes):
            if 2*atoms[k][-1]<=p:continue
            bound=dominance_lower(primes,atoms,k,future[k+1],cap,SCALE)
            if bound>=0:
                old=atoms[k].pop()
                evidence.append(dict(pass_index=sweep,p=p,removed_atom=old,replacement_atom=atoms[k][-1],lower_integer=bound))
                changed+=1
        row=dict(pass_index=sweep,removed=changed,seconds=round(time.time()-started,2))
        passes.append(row);print(json.dumps(row),flush=True)
        if not changed:break
    result=dict(status='CERTIFIED_ACTION_DOMINANCE',horizon_prime=primes[-1],state_cap=cap,scale=SCALE,
                primes=primes,remaining_atoms=atoms,eliminations=evidence,passes=passes,
                scope='Preserves the minimum of the full finite pure-stop-loss objective for cutoff 6; not by itself a barrier.',
                seconds=round(time.time()-started,2))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print('certified eliminations',len(evidence),'seconds',result['seconds'],flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output');parser.add_argument('--horizon',type=int,default=10007)
    args=parser.parse_args();certify(args.output,args.horizon)
