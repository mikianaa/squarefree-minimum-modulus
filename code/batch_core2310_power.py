"""Batch exact-integer power-envelope enumeration using binary64 dot products.

All matrix operands are nonnegative integers. The explicit bound on every
dot product is <2**52, so products and all partial sums are exactly
representable in IEEE binary64. This is exact integer arithmetic, not a
floating error estimate. A separate int64 lookup-table audit is provided.
"""
import argparse,json,math,time
from fractions import Fraction
from pathlib import Path
import numpy as np
from exact_two_stage import stage7_counts
from core2310_conditional_envelope import reduce_groups


def seed(H,den):
    knots=[0,1,2,4,8,16,32];slopes=[Fraction(H[i+1]-H[i],knots[i+1]-knots[i]) for i in range(6)]
    law=[slopes[i+1]-slopes[i] for i in range(5)]+[-slopes[-1]];assert min(law)>=0 and sum(law)==den
    factor=math.lcm(*(v.denominator for v in law));return [int(v*factor) for v in law],den*factor


def build(summary_file,manifest_file):
    started=time.time();rows=json.loads(Path(summary_file).read_text());items=[];columns=[];constants=[]
    for row in rows:
        source=row['measure_file'];data=json.loads(Path(source).read_text());w=np.zeros(2310,dtype=np.int64)
        for x,v in zip(data['support'],data['weights']):w[x]=max(0,round(v*10**8))
        den=int(sum(w));assert den>0;groups=reduce_groups(w);assert sum(int(sum(v)) for g in groups for v in g.values())==4*den
        for group in groups:
            for j in [1,2,4,8,16]:
                coefficient=np.zeros(90,dtype=np.int64);constant=0
                for k,v in group.items():
                    if j<=k:coefficient[:30]+=k*v;constant+=j*int(sum(v))
                    elif j==2*k:coefficient[30:60]+=k*v
                    elif j==4*k:coefficient[60:]+=k*v
                    else:assert j>=8*k
                assert int(sum(coefficient))*8+constant<2**52
                columns.append(coefficient);constants.append(constant)
        items.append(dict(measure_file=source,weights_numerators=w.tolist(),weights_denominator=den,period=2310,
                          group_scale=4,elimination_order=[11,7],initial_free_moduli=data['free_moduli'],
                          power_knots=[0,1,2,4,8,16,32],format='power-knots-only',
                          scope='Exact core comparison; free initial charges and continuation require separate verification'))
    A=np.array(columns,dtype=np.float64).T;constant=np.array(constants,dtype=np.int64);maximum=np.zeros(len(columns),dtype=np.int64);arguments=np.zeros(len(columns),dtype=np.int64)
    C=1+stage7_counts()
    for lo in range(0,len(C),4096):
        B=C[lo:lo+4096].astype(np.int64);features=np.concatenate([B,np.maximum(B-2,0),np.maximum(B-4,0)],axis=1).astype(np.float64)
        products=features@A;indices=np.argmax(products,axis=0);values=products[indices,np.arange(len(columns))];assert np.all(values==np.floor(values))
        values=values.astype(np.int64)-constant;better=values>maximum
        maximum[better]=values[better];arguments[better]=lo+indices[better]
    manifest=[]
    for i,item in enumerate(items):
        maxima=maximum[20*i:20*(i+1)].reshape(4,5);h=[int(v) for v in maxima.sum(axis=0)];H=[h[0]+4*item['weights_denominator']]+h+[0]
        nums,den=seed(H,4*item['weights_denominator']);item.update(power_envelope_numerators=H,group_power_envelope_numerators=maxima.tolist(),
               envelope_denominator=4*item['weights_denominator'],power_seed_numerators=nums,power_seed_denominator=den)
        item['maximizing_catalogues']=[[g,j,int(arguments[20*i+5*g+t])] for g in range(4) for t,j in enumerate([1,2,4,8,16])]
        path=item['measure_file'].replace('-measure.json','-power-envelope.json');assert path!=item['measure_file']
        Path(path).write_text(json.dumps(item,indent=2));manifest.append(dict(core_file=path,measure_file=item['measure_file']))
    Path(manifest_file).write_text(json.dumps(manifest,indent=2));print('BATCH EXACT CORES',len(items),'seconds',time.time()-started,flush=True)


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    ap=argparse.ArgumentParser();ap.add_argument('summary');ap.add_argument('manifest');a=ap.parse_args();build(a.summary,a.manifest)
