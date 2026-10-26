"""Separate signed-integer verification of every action elimination.

Does not import the pruning generator. Uses a different state cap, scale,
and rounding decomposition; only the sign of each certified bound is used.
"""
import argparse,hashlib,json,math,time
from pathlib import Path

SCALE=10**10
CAP=28


def primes(limit):
    sieve=[True]*(limit+1);sieve[0]=sieve[1]=False
    for d in range(2,math.isqrt(limit)+1):
        if sieve[d]:
            for x in range(d*d,limit+1,d):sieve[x]=False
    return [x for x in range(5,limit+1) if sieve[x]]


def choices(p):
    if p==5:return [0,1,3]
    return [0]+[2**n for n in range(p.bit_length()) if 2**n<p]


def continuation(ps,allowed):
    table=[[0]*(CAP+1) for _ in range(len(ps)+1)]
    for i in range(len(ps)-1,-1,-1):
        p=ps[i];nextrow=table[i+1]
        for n in range(CAP+1):
            count=2**n-int(p==5);count_plus=2**(n+1)-int(p==5)
            possibilities=[]
            for atom in allowed[p]:
                denominator=p-atom
                immediate=(max(0,count_plus-atom)-max(0,count-atom))*SCALE//denominator
                later=((denominator-1)*nextrow[n]+nextrow[min(n+1,CAP)])//denominator
                possibilities.append(immediate+later)
            table[i][n]=min(possibilities)
    return table


def check_elimination(ps,allowed,index,future):
    p=ps[index];atom=allowed[p][-1]
    values=[future[n]+((2**n-int(p==5)-p)*SCALE if 2**n-int(p==5)>=atom else 0) for n in range(CAP+1)]
    for q in reversed(ps[:index]):
        transformed=[]
        for n,v in enumerate(values):
            change=values[min(n+1,CAP)]-v
            denominator=q if change>=0 else q-allowed[q][-1]
            transformed.append(v+change//denominator)
        values=transformed
    return (2*values[0]+3*values[1]+values[2])//6


def audit(source,output):
    started=time.time();raw=Path(source).read_bytes();data=json.loads(raw)
    ps=primes(data['horizon_prime']);assert ps==data['primes'] and 2**CAP>2*ps[-1]
    allowed={p:choices(p) for p in ps};positions={p:i for i,p in enumerate(ps)}
    evidence=data['eliminations'];audited=[]
    previous_pass=-1
    for pass_index in sorted({row['pass_index'] for row in evidence}):
        assert pass_index>previous_pass;previous_pass=pass_index
        table=continuation(ps,allowed)
        last=-1
        for row in (r for r in evidence if r['pass_index']==pass_index):
            p=row['p'];i=positions[p];assert i>last;last=i
            assert allowed[p][-1]==row['removed_atom'] and 2*allowed[p][-1]>p
            assert len(allowed[p])>=2 and allowed[p][-2]==row['replacement_atom']
            lower=check_elimination(ps,allowed,i,table[i+1])
            assert lower>=0,(p,lower)
            allowed[p].pop()
            audited.append(dict(pass_index=pass_index,p=p,lower_integer=lower))
            if len(audited)%100==0:print('verified',len(audited),'seconds',round(time.time()-started,1),flush=True)
    assert [allowed[p] for p in ps]==data['remaining_atoms']
    result=dict(status='VERIFIED_FULL_RANGE_ACTION_DOMINANCE',source=str(source),source_sha256=hashlib.sha256(raw).hexdigest(),
                horizon_prime=ps[-1],state_cap=CAP,scale=SCALE,eliminations_checked=len(audited),
                smallest_verified_lower_integer=min(r['lower_integer'] for r in audited),checks=audited,
                seconds=round(time.time()-started,2))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('source');parser.add_argument('output')
    args=parser.parse_args();audit(args.source,args.output)
