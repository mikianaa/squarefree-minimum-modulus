"""Use actual dyadic stop-loss thresholds, avoiding integer interpolation.

All continuation arithmetic and catalogue separation here are numerical.
No reported objective is an independently certified covering bound.
"""
import argparse,json,time,math
from fractions import Fraction
from pathlib import Path
import numpy as np
from ortools.linear_solver import pywraplp
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation
from core210_joint_probe import DIVS,price


def basis(schedule_file,output):
    start=time.time();new=json.loads(Path(schedule_file).read_text())
    old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text())
    controls=dict(zip(old['primes'],old['deltas']));controls.update(zip(new['primes'],new['deltas']))
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 11<=p<=3000000]
    n=50;dist=np.zeros(n);dist[0]=1;twos=2.**np.arange(n)
    constant=mean=top=0.;coeff={};moment4=1.
    for p in ps:
        a=round(p*controls[p]);assert abs(a-p*controls[p])<1e-7
        den=p-a;r=1/den
        # Contributions with threshold at most 1 are affine in f.
        lo=a.bit_length() if a and a&(a-1) else max(0,a.bit_length()-1)
        meanpart=float(dist[lo:]@twos[lo:])/den
        mean+=meanpart;constant+=meanpart-a*float(dist[lo:].sum())/den
        for j in range(max(0,lo-4),lo):
            t=Fraction(a,1<<j)
            if 8<=t<16:top+=dist[j]*((16<<j)-a)/den
            elif 1<t<8:
                coeff[t]=coeff.get(t,0.)+dist[j]*(1<<j)/den
        dist=(1-r)*dist+r*np.r_[0,dist[:-1]];moment4*=1+3*r
    k,_=quadratic_continuation(allps,3000000);quad=k*moment4
    constant+=quad;mean+=3*quad;top+=72*quad
    for j in range(2,8):coeff[Fraction(j)]=coeff.get(Fraction(j),0.)+2*quad
    out=dict(schedule_file=schedule_file,constant=constant,mean_coefficient=mean,top_coefficient=top,
             thresholds=[str(t) for t in sorted(coeff)],coefficients=[coeff[t] for t in sorted(coeff)],
             seconds=time.time()-start,scope='Numerical finite-prime sum and quadratic continuation')
    Path(output).write_text(json.dumps(out,indent=2));print('basis',len(coeff),'thresholds',out['seconds'],flush=True)
    return out


def optimize(measure_file,basis_file,output,iterations=70):
    start=time.time();old=json.loads(Path(measure_file).read_text());pts=np.array(old['support'])
    b=json.loads(Path(basis_file).read_text());thresholds=[float(Fraction(t)) for t in b['thresholds']];coeff=b['coefficients']
    model=pywraplp.Solver.CreateSolver('CLP');q=[model.NumVar(0,1,f'q{x}') for x in pts];model.Add(sum(q)==1)
    top=model.NumVar(0,1,'maxq')
    for v in q:model.Add(top>=v)
    caps=[]
    for d in DIVS[1:]:
        cap=model.NumVar(0,1,f'cap{d}');caps.append(cap)
        for r in range(d):
            inds=np.flatnonzero(pts%d==r)
            if len(inds):model.Add(cap>=sum(q[int(i)] for i in inds))
    h=[model.NumVar(0,16,f'h{i}') for i in range(len(thresholds))]
    model.Minimize(b['constant']+b['mean_coefficient']*sum(caps)+b['top_coefficient']*top+sum(c*v for c,v in zip(coeff,h)))
    seen=set();pool=[]
    def add(i,residues):
        counts=sum((pts%d==r).astype(int) for d,r in zip(DIVS,residues))
        # Dyadic coefficients are exactly representable in binary.
        row=np.maximum(counts-thresholds[i],0);key=(i,tuple(row.tolist()))
        if key not in seen:
            seen.add(key);model.Add(h[i]>=sum(float(v)*z for v,z in zip(row,q)))
            pool.append(dict(threshold=thresholds[i],residues=[int(v) for v in residues]))
    for x in range(210):
        residues=[x%d for d in DIVS]
        for i in range(len(thresholds)):add(i,residues)
    for item in old['pool']:
        if item.get('ceiling',16)!=16:continue
        for i,t in enumerate(thresholds):
            if abs(t-item['threshold'])<1:add(i,item['residues'])
    rng=np.random.default_rng(71026);history=[]
    for it in range(iterations):
        status=model.Solve()
        if status!=pywraplp.Solver.OPTIMAL:raise RuntimeError(f'LP status {status}, not an impossibility result')
        lower=model.Objective().Value();weights=np.array([v.solution_value() for v in q])
        vals=[]
        for i,t in enumerate(thresholds):
            v,rows,residues=price(pts,weights,t,rng)
            vals.append(v)
            for a in residues:add(i,a)
        mean=sum(max(np.bincount(pts%d,weights=weights,minlength=d)) for d in DIVS[1:])
        priced=b['constant']+b['mean_coefficient']*mean+b['top_coefficient']*max(weights)+np.dot(coeff,vals)
        history.append(dict(iteration=it,restricted_lp=lower,heuristic_priced=priced,rows=len(seen)))
        if it%5==0 or priced-lower<1e-8:print(it,lower,priced,'rows',len(seen),'seconds',time.time()-start,flush=True)
        if priced-lower<1e-8:break
    result=dict(prefix_file=old.get('prefix_file'),family=old['family'],support=old['support'],weights=weights.tolist(),
                thresholds=b['thresholds'],envelope=vals,restricted_lp=lower,heuristic_priced=float(priced),
                pool=pool,history=history,basis_file=basis_file,seconds=time.time()-start,
                scope='Dyadic-threshold numerical probe with heuristic pricing; not a certified upper bound')
    Path(output).write_text(json.dumps(result,separators=(',',':')));print('DYADIC FINAL',lower,priced,flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['basis','optimize']);ap.add_argument('input');ap.add_argument('output');ap.add_argument('--basis');ap.add_argument('--iterations',type=int,default=70)
    a=ap.parse_args()
    if a.mode=='basis':basis(a.input,a.output)
    else:optimize(a.input,a.basis,a.output,a.iterations)
