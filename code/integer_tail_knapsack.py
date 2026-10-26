"""Exact integer minimum penalty for meeting a tail-volume requirement."""
import numpy as np


def minimum(target,volumes,penalties):
    assert type(target) is int and target>=0 and len(volumes)==len(penalties)
    assert all(type(v) is int and v>=0 for v in volumes)
    assert all(row and row[0]==0 and all(type(v) is int and v>=0 for v in row) for row in penalties)
    if sum(v*(len(row)-1) for v,row in zip(volumes,penalties))<target:return None
    inf=1+sum(max(row) for row in penalties)
    assert 2*inf<2**63 and all(len(row)<256 for row in penalties)
    dp=np.full(target+1,inf,dtype=np.int64);dp[0]=0;choices=[]
    for volume,row in zip(volumes,penalties):
        new=dp.copy();arg=np.zeros(target+1,dtype=np.uint8)
        for k,cost in enumerate(row[1:],1):
            credit=min(target,volume*k)
            better=cost<new[:credit+1]
            new[:credit+1][better]=cost;arg[:credit+1][better]=k
            if credit<target:
                candidate=dp[1:target-credit+1]+cost
                better=candidate<new[credit+1:]
                new[credit+1:][better]=candidate[better];arg[credit+1:][better]=k
        dp=new;choices.append(arg)
    assert int(dp[target])<inf
    remaining=target;counts=[0]*len(volumes)
    for h in range(len(volumes)-1,-1,-1):
        k=int(choices[h][remaining]);counts[h]=k
        remaining=max(0,remaining-volumes[h]*k)
    assert remaining==0 and sum(row[k] for row,k in zip(penalties,counts))==int(dp[target])
    return dict(value=int(dp[target]),counts=counts)
