"""Exact conditional/Jensen upper envelope reducing the 2310-core to 30.

All core cylinder catalogues modulo 30 are enumerated in integer arithmetic.
The resulting envelope bounds every core-2310 catalogue, not just a pool.
The final continuation evaluation printed here remains floating point.
"""
import json,sys,time,math
from fractions import Fraction
from pathlib import Path
import numpy as np
from exact_two_stage import stage7_counts


def reduce_groups(weights,order=(11,7)):
    groups=[{1:4*weights}]
    period=len(weights)
    for prime in order:
        nxt=[];size=period//prime
        for group in groups:
            left={};right={}
            for k,v in group.items():
                cells=v.reshape(prime,size);total=cells.sum(axis=0);cap=cells.max(axis=0)
                assert np.all(cap%2==0)
                left[k]=left.get(k,np.zeros(size,dtype=np.int64))+total-cap
                left[2*k]=left.get(2*k,np.zeros(size,dtype=np.int64))+cap//2
                right[2*k]=right.get(2*k,np.zeros(size,dtype=np.int64))+cap//2
            nxt.extend([left,right])
        groups=nxt;period=size
    assert period==30
    return groups


def run(measure,output,order=(11,7)):
    start=time.time();data=json.loads(Path(measure).read_text());weights=np.zeros(2310,dtype=np.int64)
    for x,q in zip(data['support'],data['weights']):weights[x]=max(0,round(q*10**8))
    den=4*int(weights.sum());assert den>0
    groups=reduce_groups(weights,order);C=1+stage7_counts();maxima=np.zeros((4,33),dtype=np.int64)
    matrices={k:np.array([g.get(k,np.zeros(30,dtype=np.int64)) for g in groups]).T for k in [1,2,4]}
    for lo in range(0,len(C),4096):
        block=C[lo:lo+4096];hist=np.zeros((len(block),4,33),dtype=np.int64)
        for k,W in matrices.items():
            for n in range(1,9):hist[:,:,k*n]+=(block==n).astype(np.int64)@W
        mass=np.zeros((len(block),4),dtype=np.int64);loss=mass.copy()
        for j in range(31,-1,-1):
            mass+=hist[:,:,j+1];loss+=mass;maxima[:,j]=np.maximum(maxima[:,j],loss.max(axis=0))
    H=maxima.sum(axis=0);assert H[0]-H[1]==den
    law=np.r_[np.diff(H,2),H[-2]];assert min(law)>=0 and sum(law)==den
    knots=[0,1,2,4,8,16,32];slopes=[Fraction(int(H[b]-H[a]),b-a) for a,b in zip(knots,knots[1:])]
    powerlaw=[slopes[i+1]-slopes[i] for i in range(5)]+[-slopes[-1]]
    factor=math.lcm(*(v.denominator for v in powerlaw));powernums=[int(v*factor) for v in powerlaw]
    assert min(powernums)>=0 and sum(powernums)==den*factor
    G=np.array(data['G'])
    numerical=float(law@G[1:]/den);power_numerical=float(np.array(powernums)@G[[1,2,4,8,16,32]]/(den*factor))
    result=dict(measure_file=measure,period=2310,elimination_order=list(order),weights_numerators=weights.tolist(),weights_denominator=int(weights.sum()),
                group_scale=4,group_envelope_numerators=maxima.tolist(),envelope_numerators=H.tolist(),envelope_denominator=den,
                law_numerators=law.tolist(),power_seed_numerators=powernums,power_seed_denominator=den*factor,
                numerical_continuation=numerical,power_numerical_continuation=power_numerical,seconds=time.time()-start,
                scope='Exact core upper envelope by conditional comparison and exhaustive integer enumeration; continuation numbers are numerical')
    if data.get('support_mode')=='fixed_subfamily':
        result['initial_free_moduli']=data['free_moduli']
        result['scope']='Exact core comparison; free initial classes must be charged in addition to the later-prime continuation'
    Path(output).write_text(json.dumps(result,indent=2));print({k:v for k,v in result.items() if k not in ['weights_numerators','group_envelope_numerators','envelope_numerators','law_numerators']},flush=True)
    return result

if __name__=='__main__':run(sys.argv[1],sys.argv[2],(7,11) if '--reverse' in sys.argv else (11,7))
