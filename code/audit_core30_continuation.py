"""Separate integer audit using the unnormalized core measure and a tail-only count.

Initial covered mass is 1/10, with core comparison masses (6,13,7,1)/30.
Unlike the primary verifier, the core count is not convolved into the running
Poisson-binomial law. Remainders and scaling are also different.
"""
import bisect,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
from audit_verify import nth_prime_list,log_enclosure,exp_upper

UNIT=10**50;DEN=10**6;LIMIT=100


def audit(schedule,output):
    start=time.time();raw=Path(schedule).read_bytes();data=json.loads(raw)
    assert data['comparison']=='core30-no10'
    assert data['seed_numerators']==[6,13,7,1] and data['seed_denominator']==27
    ceiling=lambda a,b:(a+b-1)//b
    bands=data['edges'];ds=[round(x*DEN) for x in data['params']]
    assert len(bands)==len(ds) and all(a<b for a,b in zip(bands,bands[1:]))
    assert all(0<=d<=DEN//2 for d in ds)
    ps=nth_prime_list(data['n_stages']);cut=data['p_exact'];assert ps[-1]>=cut and bands[-1]>=len(ps)
    law=[UNIT]+[0]*LIMIT;moment2=UNIT;moment3=UNIT;moment4=UNIT
    mass_weights=[6,13,7,1];fvalues=[1,2,4,8]
    early=UNIT//10;exact=far=0
    for index,p in enumerate(ps):
        if p<7:continue
        d=ds[bisect.bisect_left(bands,index+1)];b=p*(DEN-d)
        bounds=[ceiling(68*moment2,30*p)]
        if p<=cut:
            charge=0
            for j,prob in enumerate(law):
                for f,w in zip(fvalues,mass_weights):
                    charge+=prob*w*max(0,f*(1<<j)*DEN-p*d)
            truncation=ceiling(68*moment3*2**LIMIT*DEN,30*3**LIMIT*b)
            bounds.append(ceiling(charge,30*b)+truncation)
        if d:bounds.append(ceiling(234*moment4*DEN**2,30*4*d*(DEN-d)*p*p))
        loss=min(bounds)
        if p<=cut:exact+=loss
        else:far+=loss
        if p<=cut:
            for j in range(LIMIT,0,-1):law[j]=ceiling(law[j]*(b-DEN)+law[j-1]*DEN,b)
            law[0]=ceiling(law[0]*(b-DEN),b)
        moment2=ceiling(moment2*(b+DEN),b)
        moment3=ceiling(moment3*(b+2*DEN),b)
        moment4=ceiling(moment4*(b+3*DEN),b)
        if (index+1)%100000==0:print(index+1,'audit partial',(early+exact+far)/UNIT,flush=True)
    low,high=log_enclosure(ps[-1]);assert low>3 and len(ps)>=61
    q=sum(Fraction(c)*high**j for j,c in enumerate([120,120,60,20,5,1]))
    bound=2*(1+Fraction(3,2)/low)*Fraction(234*moment4,30*UNIT)*exp_upper(6/low**2)*q/(low**6*ps[-1])
    tail=ceiling(bound.numerator*UNIT,bound.denominator)
    total=early+exact+far+tail
    # Convert the unnormalized total back to the primary normalized scale.
    normalized=Fraction(total-early,UNIT)*Fraction(10,9)
    result=dict(passed=total<UNIT,claim='No cover when all moduli >6 and modulus 10 is absent',
        subprobability_total_integer=str(total),scale=str(UNIT),subprobability_total_upper=total/UNIT,
        normalized_upper=str(normalized),normalized_float=float(normalized),
        part_integers=dict(early=str(early),exact=str(exact),far=str(far),tail=str(tail)),
        parameters=dict(L=DEN,count_cap=LIMIT,p_exact=cut,n_stages=len(ps),last_prime=ps[-1]),
        schedule_sha256=hashlib.sha256(raw).hexdigest(),seconds=time.time()-start)
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit(*sys.argv[1:])
