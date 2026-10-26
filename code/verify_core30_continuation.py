"""Integer continuation certificate for the independently checked no-10 core.

The normalized initial comparison is 2**K, K has masses (6,13,7,1)/27.
Only primes >=7 are processed. This does NOT prove the unconditional bound 6.
"""
import bisect,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
from verify_squarefree import ceildiv,nth_primes,log_bounds,exp_upper

S=10**60;L=10**6;CAP=120


def verify(schedule,output):
    start=time.time();raw=Path(schedule).read_bytes();d=json.loads(raw)
    assert d['comparison']=='core30-no10' and d['seed_numerators']==[6,13,7,1] and d['seed_denominator']==27
    edges=d['edges'];Ds=[round(float(v)*L) for v in d['params']]
    assert len(edges)==len(Ds) and all(a<b for a,b in zip(edges,edges[1:]))
    assert all(0<=v<=L//2 for v in Ds)
    n=d['n_stages'];cut=d['p_exact'];assert edges[-1]>=n
    ps=nth_primes(n)[:n];assert ps[n-1]>=cut
    seed=[6,13,7,1]
    dist=[ceildiv(v*S,27) for v in seed]+[0]*(CAP-3)
    C=ceildiv(68*S,27);T=ceildiv(234*S,27);E3=5*S
    parts={'exact':0,'far':0};rows=[]
    for k,p in enumerate(ps,1):
        if p<7:continue
        D=Ds[bisect.bisect_left(edges,k)];den=p*(L-D)
        options=[ceildiv(C,p)]
        if p<=cut:
            numerator=sum(prob*max(0,(1<<j)*L-p*D) for j,prob in enumerate(dist))
            stop=ceildiv(numerator,den)+ceildiv(2**CAP*E3*L,3**CAP*den)
            options.append(stop)
        if D:options.append(ceildiv(T*L*L,4*D*(L-D)*p*p))
        value=min(options);parts['exact' if p<=cut else 'far']+=value
        if p<=31:rows.append(dict(p=p,D=D,upper=value/S))
        C=ceildiv(C*(den+L),den);T=ceildiv(T*(den+3*L),den);E3=ceildiv(E3*(den+2*L),den)
        if p<=cut:
            dist=[ceildiv(dist[j]*(den-L)+(dist[j-1]*L if j else 0),den) for j in range(CAP+1)]
        if k%100000==0:print(k,'partial',sum(parts.values())/S,'seconds',round(time.time()-start,1),flush=True)
    lo,hi=log_bounds(ps[-1]);assert lo>3 and n>=61
    polynomial=hi**5+5*hi**4+20*hi**3+60*hi**2+120*hi+120
    tail=2*(1+Fraction(3,2)/lo)*Fraction(T,S)*exp_upper(6/lo**2)/lo**6*polynomial/ps[-1]
    parts['tail']=ceildiv(tail.numerator*S,tail.denominator)
    total=sum(parts.values())
    result=dict(status='CERTIFIED' if total<S else 'NOT_CERTIFIED',
        claim='No distinct squarefree covering system with all moduli >6 and no modulus 10',
        comparison='core30-no10',total_integer=str(total),scale=str(S),total_upper=total/S,
        part_integers={k:str(v) for k,v in parts.items()},parts={k:v/S for k,v in parts.items()},
        parameters=dict(L=L,count_cap=CAP,p_exact=cut,n_stages=n,last_prime=ps[-1]),
        schedule_file=schedule,schedule_sha256=hashlib.sha256(raw).hexdigest(),rows=rows,seconds=time.time()-start)
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    verify(*sys.argv[1:])
