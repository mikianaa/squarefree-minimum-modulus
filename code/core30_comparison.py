"""Retain the actual stage-5 core and remove its covered points.

Finite enumeration is integer-exact; schedule search/evaluation is numerical.
For any core cylinder catalogue f (one cylinder for each divisor of 30),
compute max E_Q[(f-j)+], j=0,...,8. Q is the subprobability on points
uncovered after stage 5. The convex polygon through these maxima is the
stop-loss transform of an upper comparison distribution supported on 1,...,8.
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from exact_two_stage import X30, p5_weights, loss5, stage7_counts
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation


def build(output, d5=200000):
    source=json.loads(Path('data/union-accounting-C0-8.json').read_text())['part_a']['all_orbits']
    rows=[]; weights=[]
    for old in source:
        tab=old['table']; w=p5_weights(tab,d5)
        q=[Fraction(0) if tab[X30[x][0]*3+X30[x][1]] >> X30[x][2] & 1 else w[x] for x in range(30)]
        den=math.lcm(*(v.denominator for v in q))
        weights.append([int(v*den) for v in q])
        rows.append(dict(table=tab,denominator=den,loss5=str(loss5(tab,d5)),mass=str(sum(q))))
    W=np.array(weights,dtype=np.int64).T
    counts=1+stage7_counts()
    maxima=np.zeros((32,9),dtype=np.int64); args=np.zeros((32,9),dtype=np.int64)
    for j in range(8):
        for lo in range(0,len(counts),8192):
            vals=np.maximum(counts[lo:lo+8192]-j,0).astype(np.int64)@W
            indices=np.argmax(vals,axis=0); tops=vals[indices,np.arange(32)]
            better=tops>maxima[:,j]; maxima[better,j]=tops[better]; args[better,j]=lo+indices[better]
        print('threshold',j,'done',flush=True)
    for k,row in enumerate(rows):
        h=maxima[k].tolist()
        law=[h[j-1]-2*h[j]+h[j+1] for j in range(1,8)]+[h[7]]
        assert min(law)>=0 and sum(law)==sum(weights[k])
        row.update(stoploss_numerators=h,law_numerators=law,argmax_indices=args[k].tolist())
    result=dict(C0=6,D5=d5,L=1000000,n_configs=810000,rows=rows,
                scope='Integer-exact core envelopes; input orbit completeness requires separate check.')
    Path(output).write_text(json.dumps(result,indent=2))
    return result


def costs(p,rs,states,law):
    counts=2.**np.arange(states)
    out=np.zeros((len(rs),states))
    for f,w in enumerate(law,1):
        if w: out+=w*np.maximum(0,1+rs[:,None]*(f*counts-p))
    return out


def options(p):
    atoms=np.unique((np.arange(1,9)[:,None]*2.**np.arange(int(math.log2(p))+1)).ravel())
    atoms=atoms[atoms<p]
    return 1/(p-np.r_[0,atoms])


def evaluate(envelope_path,cutoff=100000):
    env=json.loads(Path(envelope_path).read_text())
    initial=json.loads(Path('data/schedule6-breakpoints-probe.json').read_text())
    allps=primes_to(16000000)[:1000000]; ps=[p for p in allps if 7<=p<=cutoff]
    ds=np.array(initial['params'])[np.searchsorted(initial['edges'],np.arange(4,len(ps)+4))]
    laws=np.array([np.array(r['law_numerators'])/r['denominator'] for r in env['rows']])
    states=50; dist=np.zeros(states);dist[0]=1
    totals=np.array([float(Fraction(r['loss5'])) for r in env['rows']])
    for p,d in zip(ps,ds):
        r=1/(p*(1-d)); count=2.**np.arange(states)
        charges=np.maximum(0,1+r*(np.arange(1,9)[:,None]*count-p))@dist
        totals+=laws@charges
        dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
    k,_=quadratic_continuation(allps,cutoff)
    totals+=k*(dist@4.**np.arange(states))*(laws@np.arange(1,9)**2)
    print('fixed-schedule totals',sorted(enumerate(totals),key=lambda x:-x[1]),flush=True)
    return totals.tolist()


def optimize(envelope_path,index,output,cutoff=100000,sweeps=4):
    env=json.loads(Path(envelope_path).read_text());row=env['rows'][index]
    law=np.array(row['law_numerators'])/row['denominator']; loss=float(Fraction(row['loss5']))
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 7<=p<=cutoff];states=40
    initial=json.loads(Path('data/schedule6-breakpoints-probe.json').read_text())
    ds=np.array(initial['params'])[np.searchsorted(initial['edges'],np.arange(4,len(ps)+4))]
    rs=1/(np.array(ps)*(1-ds));k,far=quadratic_continuation(allps,cutoff)
    terminal=k*(law@np.arange(1,9)**2)*4.**np.arange(states)
    for sweep in range(sweeps):
        V=np.zeros((len(ps)+1,states));V[-1]=terminal
        for i in range(len(ps)-1,-1,-1):
            v=V[i+1];r=rs[i]
            V[i]=costs(ps[i],np.array([r]),states,law)[0]+(1-r)*v+r*np.r_[v[1:],0]
        dist=np.zeros(states);dist[0]=1;total=loss
        print('case',index,'sweep',sweep,'start',loss+V[0,0],flush=True)
        for i,p in enumerate(ps):
            opts=options(p);v=V[i+1]
            charges=costs(p,opts,states,law)
            objectives=charges+(1-opts[:,None])*v+opts[:,None]*np.r_[v[1:],0]
            j=np.argmin(objectives@dist);r=opts[j];rs[i]=r
            total+=float(charges[j]@dist);dist=(1-r)*dist+r*np.r_[0,dist[:-1]]
        total+=float(dist@terminal)
    result=dict(C0=6,envelope_file=str(envelope_path),case=index,table=row['table'],p_exact=cutoff,
                numerical_total=total,scope='floating candidate, not a certificate',
                primes=ps,deltas=(1-1/(np.array(ps)*rs)).tolist())
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print('FINAL',index,total,flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['build','evaluate','optimize']);ap.add_argument('file')
    ap.add_argument('--d5',type=int,default=200000);ap.add_argument('--case',type=int,default=0)
    ap.add_argument('--p-exact',type=int,default=100000);ap.add_argument('--sweeps',type=int,default=4)
    ap.add_argument('--output',default='tmp/core30-optimized.json');a=ap.parse_args()
    if a.mode=='build':build(a.file,a.d5)
    elif a.mode=='evaluate':evaluate(a.file,a.p_exact)
    else:optimize(a.file,a.case,a.output,a.p_exact,a.sweeps)
