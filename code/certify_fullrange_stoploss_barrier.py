"""Integer barrier after exact full-range action dominance.

The conclusion concerns the ORIGINAL objective through the dominance horizon,
not necessarily its shorter prefix used for the final relaxed lower bound.
All arithmetic is integer/rational, with downward rounding of lower bounds.
"""
import argparse,bisect,hashlib,itertools,json,math,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from verify_squarefree import primes_to
from certify_fullrange_action_pruning import initial_atoms

SCALE=10**9
MIX=10**9
BOUND_SCALE=10**18


def load_actions(path):
    raw=Path(path).read_bytes();data=json.loads(raw)
    primes=[p for p in primes_to(data['horizon_prime']) if p>=5]
    assert primes==data['primes']
    remaining={p:initial_atoms(p) for p in primes}
    for row in data['eliminations']:
        p=row['p']
        assert row['lower_integer']>=0 and remaining[p][-1]==row['removed_atom']
        remaining[p].pop()
        assert remaining[p][-1]==row['replacement_atom']
    assert [remaining[p] for p in primes]==data['remaining_atoms']
    return data,remaining,hashlib.sha256(raw).hexdigest()


def prefix_cases(actions):
    out=[]
    for aa in itertools.product(*(actions[p] for p in (5,7,11))):
        w=[F(1,3),F(1,2),F(1,6)];cost=F(0)
        for p,a in zip((5,7,11),aa):
            cost+=sum(v*F(max(0,(1<<n)-int(p==5)-a),p-a) for n,v in enumerate(w))
            new=[F(0)]*(len(w)+1)
            for n,v in enumerate(w):new[n]+=v*F(p-a-1,p-a);new[n+1]+=v/F(p-a)
            w=new
        den=math.lcm(*(v.denominator for v in w))
        out.append(dict(atoms=list(aa),deltas=[str(F(a,p)) for a,p in zip(aa,(5,7,11))],
                        seed_numerators=[int(v*den) for v in w],seed_denominator=den,
                        prefix_loss=str(cost),prefix_lower=int(cost*SCALE)))
    return out


def bernoulli(w,den):
    v=w*(den-1);v[1:]+=w[:-1]
    return v//den


def mixture_grid(primes,upper_den,width,H):
    """Rigorous mean-budget interpolation; directed rational-sum enclosures."""
    k=len(primes);R=BOUND_SCALE
    assert R%H==0
    lower=[0];upper=[0]
    for p,d in zip(primes,upper_den):
        assert d<p
        lower.append(lower[-1]+R//d-((R+p-1)//p))
        upper.append(upper[-1]+((R+d-1)//d)-R//p)
    assert upper[-1]<4*R<2**63
    xs=np.arange(width,dtype=np.int64)*(R//H)
    left_index=np.searchsorted(np.array(upper,dtype=np.int64),xs,side='right')-1
    right_index=np.searchsorted(np.array(lower,dtype=np.int64),xs,side='right')-1
    exact_caps=None
    indices=[];weights_left=[];weights_right=[]
    for j,(x,t,u) in enumerate(zip(xs,left_index,right_index)):
        x=int(x);t=int(t);u=int(u)
        if t!=u:
            if exact_caps is None:
                exact_caps=[F(0)]
                for p,d in zip(primes,upper_den):exact_caps.append(exact_caps[-1]+F(1,d)-F(1,p))
            t=bisect.bisect_right(exact_caps,F(j,H))-1
            f=(F(j,H)-exact_caps[t])/(F(1,upper_den[t])-F(1,primes[t])) if t<k else F(0)
            lo=int(f*MIX);hi=-((-f.numerator*MIX)//f.denominator)
        elif t==k:
            lo=hi=0
        else:
            p,d=primes[t],upper_den[t]
            denominator=R*(p-d)
            lo=max(0,(x-upper[t])*d*p*MIX//denominator)
            hi=min(MIX,-((-(x-lower[t])*d*p*MIX)//denominator))
        assert 0<=lo<=hi<=MIX
        indices.append(t);weights_left.append(MIX-hi);weights_right.append(lo)
    return np.array(indices),np.array(weights_left,dtype=np.int64),np.array(weights_right,dtype=np.int64)


def certify(pruning_path,output,end=5003,H=4000,states=20):
    started=time.time();pruning,actions,source_hash=load_actions(pruning_path)
    assert end<pruning['horizon_prime']<=20000 and end<=10000 and H<=10000 and states<=25
    primes=[p for p in primes_to(end) if p>=13]
    actual_den=[p-actions[p][-1] for p in primes]
    cases=prefix_cases(actions)
    seeds=np.array([x['seed_numerators'] for x in cases],dtype=np.int64)
    seed_den=np.array([x['seed_denominator'] for x in cases],dtype=np.int64)
    caps=[F(0)];moment=F(1)
    for p,d in zip(primes,actual_den):caps.append(caps[-1]+F(1,d)-F(1,p));moment*=F(d+1,d)
    sizes=[int(x*H)+1 for x in caps]
    assert moment<64 and max(seed_den)<100000
    max_numerator=max(max(primes),max((H*p+d-1)//d for p,d in zip(primes,actual_den)))
    assert SCALE*SCALE<2**63
    assert max_numerator*(64*32+end)*SCALE<2**63
    assert 1000*SCALE*max(seed_den)<2**63
    base=np.array([SCALE]+[0]*(states-1),dtype=np.int64)
    powers=np.array([1<<n for n in range(states)],dtype=np.int64)
    V=np.zeros((len(cases),sizes[-1]),dtype=np.int64)
    checkpoints=[]
    for k in range(len(primes)-1,-1,-1):
        p=primes[k];dmax=actual_den[k];width=sizes[k]
        sorted_den=sorted(actual_den[:k])
        prefix=[base]
        for d in sorted_den:prefix.append(bernoulli(prefix[-1],d))
        suffix=[None]*(k+1);suffix[k]=base
        for i in range(k-1,-1,-1):suffix[i]=bernoulli(suffix[i+1],primes[i])
        endpoints=np.array([np.convolve(prefix[t],suffix[t])[:states]//SCALE for t in range(k+1)])
        assert np.all(endpoints.sum(axis=1)<=SCALE)
        index,left,right=mixture_grid(primes[:k],sorted_den,width,H)
        w=(left[:,None]*endpoints[index]+right[:,None]*endpoints[np.minimum(index+1,k)])//MIX
        P=np.cumsum(w[:,::-1],axis=1)[:,::-1]
        A=np.cumsum((w*powers)[:,::-1],axis=1)[:,::-1]
        assert int(P.max())<=SCALE and int(A.max())<=64*SCALE
        def cost(num,den):
            rows=[]
            for a in range(6):
                n=0
                while n<states and num*(1<<(n+a))<=p*(num-den):n+=1
                if n==states:rows.append(np.zeros(width,dtype=np.int64));continue
                value=num*(1<<a)*A[:,n]-(num-den)*p*P[:,n]
                assert np.all(value>=0)
                rows.append(value//(den*p))
            matrix=np.array(rows)
            assert int(matrix.max())<=1000*SCALE
            return (seeds@matrix)//seed_den[:,None]
        bins={}
        for a in actions[p]:
            if a:bins.setdefault(H*a//(p*(p-a)),[]).append(a)
        result=np.full((len(cases),width),10000*SCALE,dtype=np.int64)
        previous=cost(1,1)
        maxjump=H*(p-dmax)//(p*dmax)
        for jump in range(maxjump+1):
            if dmax*(H+p*(jump+1))<=H*p:
                following=cost(H+p*(jump+1),H)
            else:following=cost(p,dmax)
            here=np.minimum(previous,following)
            for a in bins.get(jump,[]):here=np.minimum(here,cost(p,p-a))
            result=np.minimum(result,here+V[:,jump:jump+width])
            previous=following
        V=result
        if k%100==0:
            row=dict(k=k,p=p,width=width,seconds=round(time.time()-started,2))
            checkpoints.append(row);print(json.dumps(row),flush=True)
    for i,c in enumerate(cases):c['continuation_lower']=int(V[i,0]);c['total_lower']=c['prefix_lower']+int(V[i,0])
    minimum=min(c['total_lower'] for c in cases)
    result=dict(status='CERTIFIED_FULL_RANGE_BARRIER' if minimum>SCALE else 'INCONCLUSIVE',
                scope='cutoff 6; original independent Bernoulli catalogue; every 0<=delta_p<1, using the trivial probability cap at 1',
                original_objective_horizon=pruning['horizon_prime'],relaxed_prefix_last_prime=primes[-1],
                scale=SCALE,mix_scale=MIX,grid_denominator=H,count_terms=states,
                minimum_lower_integer=minimum,minimum_lower=minimum/SCALE,cases=cases,
                pruning_file=str(pruning_path),pruning_sha256=source_hash,checkpoints=checkpoints,
                seconds=round(time.time()-started,2))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','checkpoints')},indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pruning');parser.add_argument('output')
    parser.add_argument('--grid',type=int,default=4000)
    args=parser.parse_args();certify(args.pruning,args.output,H=args.grid)
