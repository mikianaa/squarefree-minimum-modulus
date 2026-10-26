"""Upward-rounded convex tail coefficients for loads 1,2,4,8,16,32.

Use pure stop-loss through the cutoff and pure quadratic continuation after
it. No per-load minimum of different bounds: that could destroy convexity.
The resulting linear combination works for every audited convex core law.
"""
import bisect,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
from verify_squarefree import ceildiv,nth_primes,log_bounds,exp_upper


def certify(schedule,output):
    t=time.time();raw=Path(schedule).read_bytes();s=json.loads(raw);S=10**60;L=10**6;cap=120
    ps=nth_primes(s['n_stages'])[:s['n_stages']];cut=s['p_exact'];edges=s['edges'];ds=[round(x*L) for x in s['params']]
    assert len(ds)==len(edges) and all(a<b for a,b in zip(edges,edges[1:])) and edges[-1]>=len(ps)
    assert all(0<=d<=L//2 for d in ds) and ps[-1]>=cut>=13
    fs=[1,2,4,8,16,32];exact=[0]*6;far=[0]*6;law=[S]+[0]*cap;m3=m4=S
    for index,p in enumerate(ps,1):
        if p<13:continue
        d=ds[bisect.bisect_left(edges,index)];b=p*(L-d)
        if p<=cut:
            mass=[0]*(cap+2);moment=[0]*(cap+2)
            for j in range(cap,-1,-1):mass[j]=mass[j+1]+law[j];moment[j]=moment[j+1]+(law[j]<<j)
            for k,f in enumerate(fs):
                threshold=min(cap+1,(p*d//(f*L)).bit_length())
                numerator=f*L*moment[threshold]-p*d*mass[threshold];assert numerator>=0
                exact[k]+=ceildiv(numerator,b)+ceildiv(f*2**cap*m3*L,3**cap*b)
            law=[ceildiv(law[j]*(b-L)+(law[j-1]*L if j else 0),b) for j in range(cap+1)]
        else:
            assert d>0
            for k,f in enumerate(fs):far[k]+=ceildiv(f*f*m4*L*L,4*d*(L-d)*p*p)
        m3=ceildiv(m3*(b+2*L),b);m4=ceildiv(m4*(b+3*L),b)
    lo,hi=log_bounds(ps[-1]);assert lo>3 and len(ps)>=61
    polynomial=hi**5+5*hi**4+20*hi**3+60*hi**2+120*hi+120
    tail=2*(1+Fraction(3,2)/lo)*Fraction(m4,S)*exp_upper(6/lo**2)*polynomial/(lo**6*ps[-1])
    tails=[ceildiv(tail.numerator*S*f*f,tail.denominator) for f in fs]
    totals=[a+b+c for a,b,c in zip(exact,far,tails)]
    result=dict(status='CERTIFIED',loads=fs,upper_integers=[str(x) for x in totals],scale=str(S),upper_floats=[x/S for x in totals],
                schedule_file=schedule,schedule_sha256=hashlib.sha256(raw).hexdigest(),parameters=dict(p_exact=cut,n_stages=len(ps),count_cap=cap,last_prime=ps[-1]),
                seconds=time.time()-t,scope='Convex later-prime coefficients only; combine with audited core laws and initial free-class charges')
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    certify(*sys.argv[1:])
