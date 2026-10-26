"""Numerical weight optimization for an exhaustively evaluated block comparison.

Search derivatives pass through conditional maxima and catalogue maxima.
Final candidates must be rationalized and independently verified.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from exact_two_stage import stage7_counts


def forward(q,support,base):
    index={x:i for i,x in enumerate(base)};fi=np.array([index[x%210] for x in support]);sizes=np.bincount(fi)
    fine=np.zeros(2310);fine[np.array(support)]=q[fi]/sizes[fi]
    groups=[{1:fine}];layers=[];period=2310
    for prime in [11,7]:
        out=[];records=[];size=period//prime
        for group in groups:
            left={};right={};record={}
            for k,v in group.items():
                cells=v.reshape(prime,size);cap=cells.max(axis=0);total=cells.sum(axis=0)
                tied=(cells==cap);selector=tied/tied.sum(axis=0)
                left[k]=left.get(k,np.zeros(size))+total-cap
                left[2*k]=left.get(2*k,np.zeros(size))+cap/2
                right[2*k]=right.get(2*k,np.zeros(size))+cap/2
                record[k]=selector
            out.extend([left,right]);records.append(record)
        layers.append((period,records));groups=out;period=size
    return groups,layers,fi,sizes


def evaluate(q,support,base,C,G):
    groups,layers,fi,sizes=forward(q,support,base);c=np.r_[0,G[2]-G[1],np.diff(G[1:],2)]
    matrices={k:np.array([g.get(k,np.zeros(30)) for g in groups]).T for k in [1,2,4]}
    maxima=np.zeros((4,33));args=np.zeros((4,33),dtype=int)
    for lo in range(0,len(C),8192):
        block=C[lo:lo+8192];hist=np.zeros((len(block),4,33))
        for k,W in matrices.items():
            if np.max(W)==0:continue
            for n in range(1,9):hist[:,:,k*n]+=(block==n)@W
        mass=np.zeros((len(block),4));loss=mass.copy()
        for j in range(31,-1,-1):
            mass+=hist[:,:,j+1];loss+=mass;idx=np.argmax(loss,axis=0);v=loss[idx,np.arange(4)];better=v>maxima[:,j]
            maxima[better,j]=v[better];args[better,j]=lo+idx[better]
    H=maxima.sum(axis=0);value=float(G[1]+c@H[:32])
    gradients=[]
    for g in range(4):
        out={}
        for k in groups[g]:
            v=np.full(30,G[1])
            for j in range(1,32):v+=c[j]*np.maximum(k*C[args[g,j]].astype(float)-j,0)
            out[k]=v
        gradients.append(out)
    for period,records in reversed(layers):
        incoming=[]
        for i,record in enumerate(records):
            left,right=gradients[2*i:2*i+2];out={}
            for k,selector in record.items():
                a=left[k];b=left[2*k];cc=right[2*k]
                out[k]=(a[None,:]+selector*((b+cc)/2-a)[None,:]).ravel()
            incoming.append(out)
        gradients=incoming
    gradfine=gradients[0][1][np.array(support)]
    gradient=np.bincount(fi,weights=gradfine,minlength=len(base))/sizes
    return value,gradient,H


def run(measure,output,iterations=35,rate=.01):
    start=time.time();data=json.loads(Path(measure).read_text());support=data['support'];base=data['base'];q=np.array(data['fibre_weights'])
    q=np.maximum(q,1e-6);q/=sum(q);C=1+stage7_counts();G=np.array(data['G']);best=None;history=[]
    for it in range(iterations):
        value,gradient,H=evaluate(q,support,base,C,G)
        if best is None or value<best[0]:best=(value,q.copy(),H.copy())
        history.append(dict(iteration=it,value=value,best=best[0]));print(it,value,'best',best[0],'seconds',round(time.time()-start,1),flush=True)
        eta=rate/(1+it/20)**.6
        q*=np.exp(-eta*(gradient-np.min(gradient)));q/=sum(q)
    value,q,H=best;idx={x:i for i,x in enumerate(base)};sizes={x:sum(y%210==x for y in support) for x in base}
    weights=[float(q[idx[x%210]]/sizes[x%210]) for x in support]
    result=dict(prefix_file=data['prefix_file'],period=2310,support=support,family=data['family'],base=base,
                fibre_weights=q.tolist(),weights=weights,envelope=H.tolist(),G=G.tolist(),numerical_conditional_bound=value,
                history=history,seconds=time.time()-start,scope='Numerical optimization of exhaustive conditional/Jensen comparison; requires exact envelope and continuation verification')
    Path(output).write_text(json.dumps(result,separators=(',',':')));return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('measure');ap.add_argument('output');ap.add_argument('--iterations',type=int,default=35);ap.add_argument('--rate',type=float,default=.01)
    a=ap.parse_args();run(a.measure,a.output,a.iterations,a.rate)
