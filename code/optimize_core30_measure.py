"""Numerical LP search for an arbitrary initial probability on the uncovered core.

Uses exhaustive separation of 810000 core cylinder catalogues. This is a
candidate search only: neither floating LP status nor output is a certificate.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from ortools.linear_solver import pywraplp
from exact_two_stage import X30,stage7_counts
from verify_squarefree import primes_to
from optimize_stoploss_breakpoints import quadratic_continuation


def continuation(cutoff=100000):
    allps=primes_to(16000000)[:1000000];ps=[p for p in allps if 7<=p<=cutoff]
    schedule=json.loads(Path('data/schedule6-breakpoints-probe.json').read_text())
    ds=np.array(schedule['params'])[np.searchsorted(schedule['edges'],np.arange(4,len(ps)+4))]
    n=45;law=np.zeros(n);law[0]=1;G=np.zeros(9)
    for p,d in zip(ps,ds):
        r=1/(p*(1-d))
        G+=np.maximum(0,1+r*(np.arange(9)[:,None]*2.**np.arange(n)-p))@law
        law=(1-r)*law+r*np.r_[0,law[:-1]]
    k,_=quadratic_continuation(allps,cutoff)
    G+=k*(law@4.**np.arange(n))*np.arange(9)**2
    coefficients=np.r_[0,G[2]-G[1],np.diff(G[1:],2)]
    assert min(coefficients)>-1e-12
    return G,coefficients


def run(index,output,iterations=100,cutoff=100000):
    start=time.time()
    source=json.loads(Path('data/core30-comparison-d5-200000.json').read_text())
    table=source['rows'][index]['table']
    support=[x for x in range(30) if not table[X30[x][0]*3+X30[x][1]] >> X30[x][2]&1]
    counts=1+stage7_counts()[:,support]
    G,c=continuation(cutoff)
    solver=pywraplp.Solver.CreateSolver('GLOP')
    q=[solver.NumVar(0,1,f'q{x}') for x in support]
    hs=[None]+[solver.NumVar(0,8,f'h{j}') for j in range(1,8)]
    solver.Add(sum(q)==1)
    solver.Minimize(G[1]+sum(float(c[j])*hs[j] for j in range(1,8)))
    seen=set();history=[]
    for it in range(iterations):
        status=solver.Solve();assert status==pywraplp.Solver.OPTIMAL
        weights=np.array([v.solution_value() for v in q]);lower=solver.Objective().Value();H=np.zeros(9);H[8]=0
        new=0
        for j in range(1,8):
            values=np.maximum(counts-j,0)@weights
            best=int(np.argmax(values));H[j]=values[best]
            # Add several distinct near-worst rows to reduce degenerate pricing.
            picks=np.argpartition(values,-min(20,len(values)))[-20:]
            for a in np.r_[best,picks]:
                row=tuple(int(v) for v in np.maximum(counts[a]-j,0));key=(j,row)
                if key not in seen:
                    solver.Add(hs[j]>=sum(float(v)*z for v,z in zip(row,q)))
                    seen.add(key);new+=1
        H[0]=H[1]+1
        upper=float(G[1]+c@H[:8])
        history.append(dict(iteration=it,lower=lower,upper=upper,rows=len(seen)))
        if it%5==0 or upper-lower<1e-9:print(index,it,lower,upper,'rows',len(seen),flush=True)
        if upper-lower<1e-9:break
    law=np.r_[np.diff(H,2),H[7]]
    assert min(law)>-1e-10 and abs(sum(law)-1)<1e-10
    result=dict(case=index,table=table,support=support,weights=weights.tolist(),
                envelope=H.tolist(),law=law.tolist(),numerical_total=upper,
                numerical_lp_lower=lower,history=history,p_exact=cutoff,seconds=time.time()-start,
                scope='floating LP and schedule evaluation; not a certificate')
    Path(output).write_text(json.dumps(result,indent=2))
    print('FINAL',index,upper,'seconds',time.time()-start,flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,default=4)
    ap.add_argument('--output',default='tmp/core30-measure-case4.json')
    ap.add_argument('--iterations',type=int,default=100);ap.add_argument('--p-exact',type=int,default=100000)
    a=ap.parse_args();run(a.case,a.output,a.iterations,a.p_exact)
