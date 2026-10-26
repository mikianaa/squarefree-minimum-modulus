"""Separate continuation audit: keep core load outside the tail count law.

Uses independent prime/log/tail helpers, a different scale and truncation.
The core inequality is independently checked by audit_core2310_power_envelope.
"""
import bisect,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
from audit_verify import nth_prime_list,log_enclosure,exp_upper

UNIT=10**50;DEN=10**6;LIMIT=100


def audit(schedule,output):
    start=time.time();raw=Path(schedule).read_bytes();data=json.loads(raw)
    assert data['comparison']=='fixed-core2310'
    core_raw=Path(data['envelope_file']).read_bytes();core=json.loads(core_raw)
    nums=data['seed_numerators'];q=data['seed_denominator'];fs=[1,2,4,8,16,32]
    assert nums==core['power_seed_numerators'] and q==core['power_seed_denominator'] and sum(nums)==q
    assert len(nums)==6 and all(w>=0 for w in nums)
    m1=sum(f*w for f,w in zip(fs,nums));m2=sum(f*f*w for f,w in zip(fs,nums))
    up=lambda a,b:(a+b-1)//b
    bands=data['edges'];ds=[round(x*DEN) for x in data['params']]
    assert len(bands)==len(ds) and all(a<b for a,b in zip(bands,bands[1:])) and all(0<=d<=DEN//2 for d in ds)
    ps=nth_prime_list(data['n_stages']);cut=data['p_exact'];assert ps[-1]>=cut>=13 and bands[-1]>=len(ps)
    law=[UNIT]+[0]*LIMIT;moment2=moment3=moment4=UNIT;exact=far=0
    for index,p in enumerate(ps):
        if p<13:continue
        d=ds[bisect.bisect_left(bands,index+1)];b=p*(DEN-d);bounds=[up(m1*moment2,q*p)]
        if p<=cut:
            charge=sum(prob*sum(w*max(0,f*(1<<j)*DEN-p*d) for f,w in zip(fs,nums)) for j,prob in enumerate(law))
            remainder=up(m1*moment3*2**LIMIT*DEN,q*3**LIMIT*b)
            bounds.append(up(charge,q*b)+remainder)
        if d:bounds.append(up(m2*moment4*DEN**2,q*4*d*(DEN-d)*p*p))
        if p<=cut:exact+=min(bounds)
        else:far+=min(bounds)
        if p<=cut:
            for j in range(LIMIT,0,-1):law[j]=up(law[j]*(b-DEN)+law[j-1]*DEN,b)
            law[0]=up(law[0]*(b-DEN),b)
        moment2=up(moment2*(b+DEN),b);moment3=up(moment3*(b+2*DEN),b);moment4=up(moment4*(b+3*DEN),b)
        if (index+1)%200000==0:print(index+1,'audit partial',(exact+far)/UNIT,flush=True)
    lo,hi=log_enclosure(ps[-1]);assert lo>3 and len(ps)>=61
    poly=sum(Fraction(c)*hi**j for j,c in enumerate([120,120,60,20,5,1]))
    bound=2*(1+Fraction(3,2)/lo)*Fraction(m2*moment4,q*UNIT)*exp_upper(6/lo**2)*poly/(lo**6*ps[-1])
    tail=up(bound.numerator*UNIT,bound.denominator);total=exact+far+tail
    result=dict(status='VERIFIED' if total<UNIT else 'NOT_VERIFIED',total_integer=str(total),scale=str(UNIT),total_upper=total/UNIT,
                part_integers=dict(exact=str(exact),far=str(far),tail=str(tail)),
                schedule_sha256=hashlib.sha256(raw).hexdigest(),core_sha256=hashlib.sha256(core_raw).hexdigest(),
                parameters=dict(count_cap=LIMIT,p_exact=cut,n_stages=len(ps),last_prime=ps[-1]),seconds=time.time()-start,
                scope='No covering extension with all moduli >6 of the specified fixed core only; not the global upper bound 6')
    if core.get('initial_free_moduli'):
        result['scope']='Later-prime loss only; free initial classes must be charged separately, and global upper bound 6 is not established'
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit(*sys.argv[1:])
