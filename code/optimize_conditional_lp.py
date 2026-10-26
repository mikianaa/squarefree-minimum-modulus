"""Cutting-plane LP for the exhaustive conditional block comparison.

Every pricing call enumerates all 810000 core-30 catalogues. Numerical LP
output still requires integer core and continuation verification.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from ortools.linear_solver import pywraplp
from exact_two_stage import stage7_counts


def optimize(measure,output,iterations=60,schedule=None,free_moduli=(),point_weights=False,density_period=None,union_discount=False):
    t=time.time();data=json.loads(Path(measure).read_text());base=data['base'];support=data['support'];n=len(support) if point_weights else len(base);nc=len(base) if point_weights else 0;nv=n+nc+90
    G=np.array(json.loads(Path(schedule).read_text())['G'] if schedule else data['G'])
    coeff=np.r_[0,G[2]-G[1],np.diff(G[1:],2)];assert min(coeff)>-1e-9
    first={k:np.zeros((210,nv)) for k in [1,2]}
    model=pywraplp.Solver.CreateSolver('CLP');variables=[model.NumVar(0,1,f'v{i}') for i in range(nv)]
    model.Add(sum(variables[:n])==1)
    union_variables=variables[:n]
    if union_discount:
        assert not point_weights
        union_variables=[model.NumVar(0,1,f'alpha{i}') for i in range(n)]
        for a,v in zip(union_variables,variables):model.Add(a<=v)
    if density_period is not None:
        assert point_weights and 2310%density_period==0
        leaders={}
        for j,y in enumerate(support):
            residue=y%density_period
            if residue in leaders:model.Add(variables[j]==variables[leaders[residue]])
            else:leaders[residue]=j
    expr=lambda a:sum(float(c)*v for c,v in zip(a,variables) if c)
    if point_weights:
        for i,x in enumerate(base):
            indices=[j for j,y in enumerate(support) if y%210==x]
            first[1][x,indices]=1;first[1][x,n+i]=-1;first[2][x,n+i]=.5
            for j in indices:model.Add(variables[n+i]>=variables[j])
            model.Add(variables[n+i]<=sum(variables[j] for j in indices))
    else:
        for i,x in enumerate(base):
            size=sum(y%210==x for y in support);first[1][x,i]=1-1/size;first[2][x,i]=.5/size
    initial=[first,{2:first[2]}];groups=[];capindex=n+nc
    for group in initial:
        left={};right={}
        for k,M in group.items():
            T=M.reshape(7,30,nv).sum(axis=0);cap=np.zeros((30,nv))
            for x in range(30):
                cap[x,capindex]=1
                for y in range(x,210,30):model.Add(variables[capindex]>=expr(M[y]))
                model.Add(variables[capindex]<=expr(T[x]));capindex+=1
            left[k]=left.get(k,np.zeros((30,nv)))+T-cap
            left[2*k]=left.get(2*k,np.zeros((30,nv)))+cap/2
            right[2*k]=right.get(2*k,np.zeros((30,nv)))+cap/2
        groups.extend([left,right])
    assert capindex==nv
    h={(g,j):model.NumVar(0,32,f'h{g}_{j}') for g in range(4) for j in range(1,32)}
    free_caps=[];free_rows=[]
    index={x:i for i,x in enumerate(base)};sizes={x:sum(y%210==x for y in support) for x in base}
    for d in free_moduli:
        assert d>6 and 2310%d==0
        M=np.zeros((d,nv))
        for j,y in enumerate(support):
            if point_weights:M[y%d,j]=1
            else:M[y%d,index[y%210]]+=1/sizes[y%210]
        cap=model.NumVar(0,1,f'free{d}');free_caps.append(cap);free_rows.append(M)
        for row in M:model.Add(cap>=sum(float(c)*v for c,v in zip(row[:n],union_variables) if c))
    model.Minimize(float(G[1])+sum(float(coeff[j])*h[g,j] for g,j in h)+sum(free_caps)+(1-sum(union_variables) if union_discount else 0))
    C=1+stage7_counts();seen=set();history=[];best=None
    def add(g,j,row):
        key=(g,j,int(row))
        if key in seen:return
        seen.add(key);a=sum(np.maximum(k*C[row].astype(float)-j,0)@M for k,M in groups[g].items())
        model.Add(h[g,j]>=expr(a))
    for g,j in h:add(g,j,0)
    for g,j,row in data.get('pricing_cuts',[]):add(g,j,row)
    for it in range(iterations):
        assert model.Solve()==pywraplp.Solver.OPTIMAL
        solution=np.array([v.solution_value() for v in variables]);maxima=np.zeros((4,33));args=np.zeros((4,33),dtype=int)
        matrices={k:np.array([g.get(k,np.zeros((30,nv)))@solution for g in groups]).T for k in [1,2,4]}
        assert min(M.min() for M in matrices.values())>-1e-7
        for lo in range(0,len(C),8192):
            B=C[lo:lo+8192];hist=np.zeros((len(B),4,33))
            for k,W in matrices.items():
                for m in range(1,9):hist[:,:,k*m]+=(B==m)@W
            mass=np.zeros((len(B),4));loss=mass.copy()
            for j in range(31,-1,-1):
                mass+=hist[:,:,j+1];loss+=mass;idx=np.argmax(loss,axis=0);val=loss[idx,np.arange(4)];better=val>maxima[:,j]
                maxima[better,j]=val[better];args[better,j]=lo+idx[better]
        union_solution=np.array([v.solution_value() for v in union_variables])
        H=maxima.sum(axis=0);value=float(G[1]+coeff@H[:32])+sum(float(np.max(M[:,:n]@union_solution)) for M in free_rows)+(1-sum(union_solution) if union_discount else 0);lower=model.Objective().Value()
        if best is None or value<best[0]:best=(value,solution[:n].copy(),H.copy(),union_solution.copy())
        history.append(dict(iteration=it,lower=lower,value=value,best=best[0],cuts=len(seen)))
        print(history[-1],'seconds',round(time.time()-t,1),flush=True)
        old_cut_count=len(seen)
        for g,j in h:add(g,j,args[g,j])
        value,q,H,alpha=best;sizes={x:sum(y%210==x for y in support) for x in base};index={x:i for i,x in enumerate(base)}
        result={k:data[k] for k in ['prefix_file','period','support','family','base']}
        if 'support_mode' in data:result['support_mode']=data['support_mode']
        if 'fixed_family' in data:result['fixed_family']=data['fixed_family']
        weights=q.tolist() if point_weights else [float(q[index[y%210]]/sizes[y%210]) for y in support]
        fibre=[sum(weights[j] for j,y in enumerate(support) if y%210==x) for x in base]
        result.update(fibre_weights=fibre,weights=weights,point_weights=point_weights,density_period=density_period,pricing_cuts=[list(k) for k in sorted(seen)],
                      G=G.tolist(),envelope=H.tolist(),numerical_conditional_bound=value,free_moduli=list(free_moduli),history=history,
                      scope='Numerical LP with exhaustive catalogue pricing; one fixed prefix, requires integer verification')
        if union_discount:result['union_fibre_weights']=alpha.tolist()
        Path(output).write_text(json.dumps(result,separators=(',',':')))
        if best[0]-lower<1e-7 or len(seen)==old_cut_count:break
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('measure');ap.add_argument('output');ap.add_argument('--iterations',type=int,default=60);ap.add_argument('--schedule');ap.add_argument('--free-moduli',type=int,nargs='*',default=[]);ap.add_argument('--point-weights',action='store_true');ap.add_argument('--density-period',type=int);ap.add_argument('--union-discount',action='store_true')
    a=ap.parse_args();optimize(a.measure,a.output,a.iterations,a.schedule,a.free_moduli,a.point_weights,a.density_period,a.union_discount)
