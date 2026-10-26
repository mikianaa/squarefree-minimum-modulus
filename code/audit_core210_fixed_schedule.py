"""Solver-free rational audit of a fixed-schedule core comparison obstruction.

Recomputes finite-prime coefficients with a different downward-rounded
count recurrence, then checks all dual budgets and all 104 point scores.
"""
import hashlib,json,sys,time
from fractions import Fraction as F
from pathlib import Path


def audit(dual_file,output):
    start=time.time();certificate=json.loads(Path(dual_file).read_text())
    raw=Path(certificate['coefficient_file']).read_bytes();data=json.loads(raw)
    assert hashlib.sha256(raw).hexdigest()==certificate['coefficient_sha256']
    limit=data['p_max'];sieve=bytearray([1])*(limit+1);sieve[:2]=b'\0\0'
    for k in range(2,int(limit**.5)+1):
        if sieve[k]:sieve[k*k::k]=bytes(len(range(k*k,limit+1,k)))
    primes=[p for p in range(11,limit+1) if sieve[p]];assert len(primes)==len(data['atoms'])
    S=10**24;N=50;prob=[S]+[0]*N;const=mean=top=0;hinges={};powers=[1<<j for j in range(N+1)]
    for p,a in zip(primes,data['atoms']):
        assert isinstance(a,int) and 0<=a<p
        den=p-a;cnum=mnum=tnum=0
        for j,two in enumerate(powers):
            if a<=two:
                cnum+=prob[j]*(two-a);mnum+=prob[j]*two
            elif a<8*two:
                key=F(a,two);hinges[key]=hinges.get(key,0)+(prob[j]*two)//den
            elif a<16*two:tnum+=prob[j]*(16*two-a)
        const+=cnum//den;mean+=mnum//den;top+=tnum//den
        # Floor the two incoming flows separately, unlike the generator.
        prob=[prob[j]*(den-1)//den+(prob[j-1]//den if j else 0) for j in range(N+1)]
    family={int(d):r for d,r in certificate['family'].items()}
    divs=[d for d in range(1,211) if 210%d==0]
    assert sorted(family)==[d for d in divs if d>6]
    assert family[7]==0 and family[10]==0
    support=[x for x in range(210) if all(x%d!=r for d,r in family.items())]
    assert support==certificate['support'] and len(support)==104
    D=certificate['dual_denominator'];assert D>0
    cylinder_budget={d:0 for d in divs[1:]};hinge_budget={};point_budget=0
    scores={x:F(const,S) for x in support}
    for term in certificate['terms']:
        w=term['weight'];assert isinstance(w,int) and w>0
        if term['kind']=='point':
            x=term['point'];assert x in scores;point_budget+=w;scores[x]+=F(w,D)
        elif term['kind']=='cylinder':
            d=term['modulus'];r=term['residue'];assert d in cylinder_budget and 0<=r<d
            cylinder_budget[d]+=w
            for x in support:
                if x%d==r:scores[x]+=F(w,D)
        else:
            assert term['kind']=='catalogue';t=F(term['threshold']);assert 1<t<8
            residues=term['residues'];assert len(residues)==16
            assert all(isinstance(r,int) and 0<=r<d for d,r in zip(divs,residues))
            hinge_budget[t]=hinge_budget.get(t,0)+w
            for x in support:
                f=sum(x%d==r for d,r in zip(divs,residues))
                scores[x]+=F(w,D)*max(F(0),f-t)
    assert point_budget*S<=top*D
    assert all(w*S<=mean*D for w in cylinder_budget.values())
    assert all(w*S<=hinges.get(t,0)*D for t,w in hinge_budget.items())
    worst=min(scores,key=scores.get);bound=scores[worst];assert bound>1
    result=dict(passed=True,lower_bound_exact=str(bound),lower_bound_float=float(bound),
                safe_lower_numerator=bound.numerator*10**9//bound.denominator,safe_lower_denominator=10**9,
                worst_point=worst,point_checks=len(support),positive_dual_terms=len(certificate['terms']),
                prime_stages=len(primes),count_cap=N,coefficient_scale=str(S),seconds=time.time()-start,
                scope='For this fixed schedule and the specified stop-loss core-catalogue functional, every initial probability on this actual 104-point core has value >1. Not an all-schedule barrier and not a covering.')
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':audit(*sys.argv[1:])
