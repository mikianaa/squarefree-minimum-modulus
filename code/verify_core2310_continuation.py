"""Integer infinite continuation for ONE explicitly specified core modulo 2310.

The core comparison must be audited separately. This is not a universal
minimum-modulus theorem: every certificate identifies its fixed prefix.
"""
import bisect,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
from verify_squarefree import ceildiv,nth_primes,log_bounds,exp_upper

S=10**60;L=10**6;CAP=120


def verify(schedule,output):
    start=time.time();raw=Path(schedule).read_bytes();d=json.loads(raw)
    assert d['comparison']=='fixed-core2310'
    core_raw=Path(d['envelope_file']).read_bytes();core=json.loads(core_raw)
    seed=d['seed_numerators'];denom=d['seed_denominator']
    assert seed==core['power_seed_numerators'] and denom==core['power_seed_denominator']
    assert len(seed)==6 and all(isinstance(v,int) and v>=0 for v in seed) and sum(seed)==denom
    edges=d['edges'];Ds=[round(float(v)*L) for v in d['params']]
    assert len(edges)==len(Ds) and all(a<b for a,b in zip(edges,edges[1:])) and all(0<=v<=L//2 for v in Ds)
    n=d['n_stages'];cut=d['p_exact'];ps=nth_primes(n)[:n];assert edges[-1]>=n and ps[-1]>=cut>=13
    dist=[ceildiv(v*S,denom) for v in seed]+[0]*(CAP+1-len(seed))
    C=ceildiv(sum(v*2**j for j,v in enumerate(seed))*S,denom)
    T=ceildiv(sum(v*4**j for j,v in enumerate(seed))*S,denom)
    E3=ceildiv(sum(v*3**j for j,v in enumerate(seed))*S,denom)
    parts={'exact':0,'far':0}
    for index,p in enumerate(ps,1):
        if p<13:continue
        D=Ds[bisect.bisect_left(edges,index)];den=p*(L-D);options=[ceildiv(C,p)]
        if p<=cut:
            numerator=sum(prob*max(0,(1<<j)*L-p*D) for j,prob in enumerate(dist))
            options.append(ceildiv(numerator,den)+ceildiv(2**CAP*E3*L,3**CAP*den))
        if D:options.append(ceildiv(T*L*L,4*D*(L-D)*p*p))
        parts['exact' if p<=cut else 'far']+=min(options)
        C=ceildiv(C*(den+L),den);T=ceildiv(T*(den+3*L),den);E3=ceildiv(E3*(den+2*L),den)
        if p<=cut:dist=[ceildiv(dist[j]*(den-L)+(dist[j-1]*L if j else 0),den) for j in range(CAP+1)]
        if index%200000==0:print(index,'partial',sum(parts.values())/S,flush=True)
    lo,hi=log_bounds(ps[-1]);assert lo>3 and n>=61
    poly=hi**5+5*hi**4+20*hi**3+60*hi**2+120*hi+120
    tail=2*(1+Fraction(3,2)/lo)*Fraction(T,S)*exp_upper(6/lo**2)/lo**6*poly/ps[-1]
    parts['tail']=ceildiv(tail.numerator*S,tail.denominator);total=sum(parts.values())
    result=dict(status='CERTIFIED' if total<S else 'NOT_CERTIFIED',claim='The specified fixed core-2310 family has no distinct squarefree covering extension with all moduli >6',
                total_integer=str(total),scale=str(S),total_upper=total/S,part_integers={k:str(v) for k,v in parts.items()},
                schedule_file=schedule,schedule_sha256=hashlib.sha256(raw).hexdigest(),core_file=d['envelope_file'],core_sha256=hashlib.sha256(core_raw).hexdigest(),
                parameters=dict(p_exact=cut,n_stages=n,last_prime=ps[-1],count_cap=CAP),seconds=time.time()-start,
                scope='Conditional on this fixed family and the audited core comparison; does not prove the global upper bound 6')
    if core.get('initial_free_moduli'):
        result['claim']='Later-prime loss bound for the specified initial measure'
        result['scope']='Free initial classes must be charged separately before any no-cover conclusion; global upper bound 6 is not established'
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    verify(*sys.argv[1:])
