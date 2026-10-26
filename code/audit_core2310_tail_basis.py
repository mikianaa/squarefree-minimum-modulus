"""Separate tail-basis audit using direct sums and different precision/helpers."""
import bisect,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
from audit_verify import nth_prime_list,log_enclosure,exp_upper


def audit(schedule,output):
    started=time.time();raw=Path(schedule).read_bytes();data=json.loads(raw);unit=10**50;den=10**6;limit=100
    primes=nth_prime_list(data['n_stages']);cut=data['p_exact'];bands=data['edges'];deltas=[round(v*den) for v in data['params']]
    assert len(bands)==len(deltas) and all(a<b for a,b in zip(bands,bands[1:])) and bands[-1]>=len(primes)
    assert all(0<=d<=den//2 for d in deltas) and primes[-1]>=cut>=13
    up=lambda a,b:(a+b-1)//b
    loads=[1,2,4,8,16,32];values=[0]*6;law=[unit]+[0]*limit;third=fourth=unit
    for stage,p in enumerate(primes,1):
        if p<13:continue
        delta=deltas[bisect.bisect_left(bands,stage)];b=p*(den-delta)
        if p<=cut:
            charges=[0]*6
            for j,prob in enumerate(law):
                for k,f in enumerate(loads):charges[k]+=prob*max(0,(f<<j)*den-p*delta)
            for k,f in enumerate(loads):values[k]+=up(charges[k],b)+up(f*third*2**limit*den,3**limit*b)
            for j in range(limit,0,-1):law[j]=up(law[j]*(b-den)+law[j-1]*den,b)
            law[0]=up(law[0]*(b-den),b)
        else:
            assert delta>0
            for k,f in enumerate(loads):values[k]+=up(f*f*fourth*den**2,4*delta*(den-delta)*p*p)
        third=up(third*(b+2*den),b);fourth=up(fourth*(b+3*den),b)
        if stage%200000==0:print(stage,'basis audit seconds',round(time.time()-started,1),flush=True)
    lo,hi=log_enclosure(primes[-1]);assert lo>3 and len(primes)>=61
    polynomial=sum(Fraction(c)*hi**j for j,c in enumerate([120,120,60,20,5,1]))
    end=2*(1+Fraction(3,2)/lo)*Fraction(fourth,unit)*exp_upper(6/lo**2)*polynomial/(lo**6*primes[-1])
    for k,f in enumerate(loads):values[k]+=up(end.numerator*unit*f*f,end.denominator)
    result=dict(status='VERIFIED',loads=loads,upper_integers=[str(v) for v in values],scale=str(unit),upper_floats=[v/unit for v in values],
                schedule_file=schedule,schedule_sha256=hashlib.sha256(raw).hexdigest(),seconds=time.time()-started,
                scope='Independent convex later-prime coefficients; core and initial free classes require separate verification')
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit(*sys.argv[1:])
