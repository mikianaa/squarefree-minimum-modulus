"""Small LP for density-30 measures, optimizing the certifiable power seed.

Reduction maxima are linear in the density on each 30-fibre. Only 30
density variables and the five power-knot envelopes are necessary.
Exhaustive floating pricing is candidate generation, not integer proof.
"""
import json,time
from pathlib import Path
import numpy as np
from ortools.linear_solver import pywraplp
from exact_two_stage import stage7_counts
from core2310_conditional_envelope import reduce_groups

COUNTS=None


def optimize(data,output,schedule,iterations=25,quiet=False,price=True,discount=False):
    global COUNTS
    if COUNTS is None:COUNTS=1+stage7_counts()
    C=COUNTS;t=time.time();support=data['support'];free=data['free_moduli'];G=np.array(schedule['G']);knots=[1,2,4,8,16]
    slopes=np.diff(G[[1,2,4,8,16,32]])/np.diff([1,2,4,8,16,32]);coeff=np.r_[slopes[0],np.diff(slopes)];assert min(coeff)>=0
    indicator=np.zeros(2310,dtype=np.int64);indicator[support]=1;groups=reduce_groups(indicator)
    groups=[{k:v/4 for k,v in group.items()} for group in groups]
    sizes=np.bincount(np.array(support)%30,minlength=30);model=pywraplp.Solver.CreateSolver('CLP')
    w=[model.NumVar(0,1 if sizes[z] else 0,f'w{z}') for z in range(30)];model.Add(sum(int(n)*v for n,v in zip(sizes,w))==1)
    expr=lambda row:sum(float(x)*v for x,v in zip(row,w) if x)
    alpha=w
    if discount:
        alpha=[model.NumVar(0,1 if sizes[z] else 0,f'a{z}') for z in range(30)]
        for a,v in zip(alpha,w):model.Add(a<=v)
    alpha_expr=lambda row:sum(float(x)*v for x,v in zip(row,alpha) if x)
    h={(g,j):model.NumVar(0,32,f'h{g}_{j}') for g in range(4) for j in knots};free_rows=[];caps=[]
    base=set(x%210 for x in support);multiplicity=np.array([sum(y%30==z for y in base) for z in range(30)])
    for d in free:
        assert d%11==0
        m=d//11;projection=m//7 if m%7==0 else m
        factor=(multiplicity>0).astype(float) if m%7==0 else multiplicity
        M=np.array([factor*(np.arange(30)%projection==r) for r in range(projection)])
        free_rows.append(M);cap=model.NumVar(0,1,f'cap{d}');caps.append(cap)
        for row in M:model.Add(cap>=alpha_expr(row))
    beta=1-alpha_expr(sizes) if discount else 0
    model.Minimize(float(G[1])+sum(float(coeff[knots.index(j)])*v for (g,j),v in h.items())+sum(caps)+beta)
    seen=set()
    def add(g,j,row):
        if j not in knots or (g,j,row) in seen:return
        seen.add((g,j,row));feature=sum(v*np.maximum(k*C[row].astype(float)-j,0) for k,v in groups[g].items());model.Add(h[g,j]>=expr(feature))
    for g,j in h:add(g,j,0)
    for g,j,row in data.get('pricing_cuts',[]):add(g,j,row)
    history=[];best=None
    for it in range(iterations):
        assert model.Solve()==pywraplp.Solver.OPTIMAL
        density=np.array([v.solution_value() for v in w]);alpha_density=np.array([v.solution_value() for v in alpha]);maxima=np.zeros((4,33));args=np.zeros((4,33),dtype=int)
        if not price:
            H=np.zeros(33)
            for (g,j),v in h.items():H[j]+=v.solution_value()
            value=model.Objective().Value();best=(value,density,H,alpha_density)
            history.append(dict(iteration=it,restricted_lp=value,cuts=len(seen)))
            break
        matrices={k:np.array([g.get(k,np.zeros(30))*density for g in groups]).T for k in [1,2,4]}
        for lo in range(0,len(C),8192):
            B=C[lo:lo+8192];hist=np.zeros((len(B),4,33))
            for k,W in matrices.items():
                for m in range(1,9):hist[:,:,k*m]+=(B==m)@W
            mass=np.zeros((len(B),4));loss=mass.copy()
            for j in range(31,-1,-1):
                mass+=hist[:,:,j+1];loss+=mass;idx=np.argmax(loss,axis=0);val=loss[idx,np.arange(4)];better=val>maxima[:,j]
                maxima[better,j]=val[better];args[better,j]=lo+idx[better]
        H=maxima.sum(axis=0);value=float(G[1]+coeff@H[knots])+sum(float(max(M@alpha_density)) for M in free_rows)+(float(1-sizes@alpha_density) if discount else 0)
        if best is None or value<best[0]:best=(value,density.copy(),H.copy(),alpha_density.copy())
        lower=model.Objective().Value();history.append(dict(iteration=it,lower=lower,value=value,best=best[0],cuts=len(seen)))
        if not quiet:print(history[-1],flush=True)
        old=len(seen)
        for g,j in h:add(g,j,int(args[g,j]))
        if best[0]-lower<1e-7 or len(seen)==old:break
    value,density,H,alpha_density=best;weights=[float(density[x%30]) for x in support];base=data['base']
    result={k:data[k] for k in ['prefix_file','period','family','support','base','support_mode','fixed_family','free_moduli']}
    result.update(weights=weights,fibre_weights=[sum(v for x,v in zip(support,weights) if x%210==b) for b in base],G=G.tolist(),
                  pricing_cuts=[list(row) for row in sorted(seen)],density_period=30,point_weights=True,envelope=H.tolist(),
                  history=history,seconds=time.time()-t,exhaustive_pricing=price,scope='Numerical small-LP candidate; integer envelope and tail audit required')
    if price:result['numerical_power_bound']=value
    else:
        result['numerical_restricted_lp']=value;result['restricted_envelope']=result.pop('envelope')
        result['scope']='Finite-pool LP only: objective is not an upper bound; exhaustive integer verification required'
    if discount:result['union_density']=alpha_density.tolist()
    Path(output).write_text(json.dumps(result,separators=(',',':')));return result
