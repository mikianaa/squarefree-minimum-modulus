"""Exact lower certificate for one fixed stop-loss schedule, all core measures.

The finite continuation coefficients are rounded DOWN. A rational dual of
a finite catalogue LP then lower-bounds the full catalogue functional.
Scope: this schedule and this functional, not every possible schedule.
"""
import argparse,json,math,hashlib
from fractions import Fraction
from pathlib import Path
from verify_squarefree import primes_to


def coefficients(schedule_file,output):
    new=json.loads(Path(schedule_file).read_text());old=json.loads(Path('data/core30-case4-numerical-3m.json').read_text())
    ds=dict(zip(old['primes'],old['deltas']));ds.update(zip(new['primes'],new['deltas']))
    S=10**30;N=60;law=[S]+[0]*N;constant=mean=top=0;hinges={};atoms=[]
    ps=[p for p in primes_to(3000000) if p>=11]
    for index,p in enumerate(ps):
        a=round(p*ds[p]);assert abs(a-p*ds[p])<1e-7 and 0<=a<p
        atoms.append(a);den=p-a
        lo=a.bit_length() if a and a&(a-1) else max(0,a.bit_length()-1)
        m=sum(law[j]*(1<<j) for j in range(lo,N+1))
        mean+=m//den;constant+=(m-a*sum(law[lo:]))//den
        tn=0
        for j in range(max(0,lo-4),lo):
            if j>N:continue
            t=Fraction(a,1<<j)
            if 8<=t<16:tn+=law[j]*((16<<j)-a)
            elif 1<t<8:hinges[t]=hinges.get(t,0)+(law[j]*(1<<j))//den
        top+=tn//den
        for j in range(N,0,-1):law[j]=(law[j]*(den-1)+law[j-1])//den
        law[0]=law[0]*(den-1)//den
    keys=sorted(hinges)
    result=dict(scale=str(S),count_cap=N,p_max=3000000,first_prime=11,atoms=atoms,
                constant=constant/S,mean_coefficient=mean/S,top_coefficient=top/S,
                thresholds=[str(t) for t in keys],coefficients=[hinges[t]/S for t in keys],
                lower_integers=dict(constant=str(constant),mean=str(mean),top=str(top),hinges=[str(hinges[t]) for t in keys]),
                schedule_file=schedule_file,scope='Exact lower coefficients for finite primes; discarded tail is nonnegative')
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print('Exact finite coefficients',len(ps),'primes',len(keys),'thresholds',flush=True)
    return result


def dual(measure_file,coeff_file,output):
    from ortools.linear_solver import pywraplp
    data=json.loads(Path(measure_file).read_text());b=json.loads(Path(coeff_file).read_text())
    pts=data['support'];divs=[d for d in range(1,211) if 210%d==0]
    ts=[Fraction(t) for t in b['thresholds']];index={t:i for i,t in enumerate(ts)}
    model=pywraplp.Solver.CreateSolver('CLP');q=[model.NumVar(0,1,f'q{x}') for x in pts];model.Add(sum(q)==1)
    top=model.NumVar(0,1,'top');records=[]
    for x,v in zip(pts,q):records.append((model.Add(top-v>=0),dict(kind='point',point=x)))
    caps=[]
    for d in divs[1:]:
        cap=model.NumVar(0,1,f'cap{d}');caps.append(cap)
        for r in range(d):
            ids=[i for i,x in enumerate(pts) if x%d==r]
            if ids:records.append((model.Add(cap-sum(q[i] for i in ids)>=0),dict(kind='cylinder',modulus=d,residue=r)))
    hs=[model.NumVar(0,16,f'h{i}') for i in range(len(ts))]
    seen=set()
    for item in data['pool']:
        t=Fraction(str(item['threshold']))
        if t not in index:continue
        i=index[t];residues=item['residues']
        row=tuple(max(Fraction(0),sum(x%d==r for d,r in zip(divs,residues))-t) for x in pts)
        key=(i,row)
        if key in seen:continue
        seen.add(key)
        records.append((model.Add(hs[i]-sum(float(v)*z for v,z in zip(row,q))>=0),
                        dict(kind='catalogue',threshold=str(t),residues=residues)))
    model.Minimize(b['constant']+b['mean_coefficient']*sum(caps)+b['top_coefficient']*top+sum(c*h for c,h in zip(b['coefficients'],hs)))
    assert model.Solve()==pywraplp.Solver.OPTIMAL
    lower=model.Objective().Value();print('restricted numerical LP',lower,flush=True)
    D=10**12;terms=[]
    for constraint,description in records:
        value=max(0,constraint.dual_value());weight=math.floor(value*(1-1e-8)*D)
        if weight:terms.append(dict(**description,weight=weight))
    result=dict(family=data['family'],support=pts,dual_denominator=D,terms=terms,
                coefficient_file=coeff_file,coefficient_sha256=hashlib.sha256(Path(coeff_file).read_bytes()).hexdigest(),
                numerical_restricted_lp=lower,
                scope='Candidate rational dual for the specified finite stop-loss functional and schedule only')
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print('dual terms',len(terms),flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['coefficients','dual']);ap.add_argument('input');ap.add_argument('output');ap.add_argument('--coefficients')
    a=ap.parse_args()
    if a.mode=='coefficients':coefficients(a.input,a.output)
    else:dual(a.input,a.coefficients,a.output)
