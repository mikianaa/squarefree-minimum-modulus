"""Search actual core families containing normalized moduli 7 and 10.

All reported residual sets are re-counted exactly. Numerical scores use
heuristic catalogue pricing and must not be interpreted as certificates.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from core210_joint_probe import ALLOWED,DIVS,price,tail_basis,optimize


def search(output,starts=200,keep=8):
    rng=np.random.default_rng(7102026);xs=np.arange(210);free=[d for d in ALLOWED if d not in (7,10)]
    found={};G=tail_basis();c=np.r_[0,G[2]-G[1],np.diff(G[1:],2)]
    for trial in range(starts):
        family={d:(0 if d in (7,10) else int(rng.integers(d))) for d in ALLOWED}
        counts=sum((xs%d==r).astype(int) for d,r in family.items())
        target=np.ones(210) if trial%3==0 else rng.uniform(.65,1.35,210)
        for sweep in range(15):
            changed=False
            for d in rng.permutation(free):
                old=family[d];counts-=(xs%d==old)
                scores=np.bincount(xs%d,weights=target*(counts==0),minlength=d)
                best=int(np.argmax(scores));family[d]=best;counts+=(xs%d==best);changed|=old!=best
            if not changed:break
        support=tuple(int(x) for x in xs[counts==0])
        if support not in found:found[support]={int(d):int(r) for d,r in family.items()}
    print('distinct prefixes',len(found),flush=True)
    ranked=[]
    for i,(support,family) in enumerate(found.items()):
        pts=np.array(support);w=np.ones(len(pts))/len(pts);H=np.zeros(17)
        H[1]=sum(max(np.bincount(pts%d,weights=w,minlength=d)) for d in DIVS[1:]);H[0]=1+H[1]
        for j in range(2,8):H[j]=price(pts,w,j,rng,random_starts=0,sweeps=3)[0]
        for j in range(8,17):H[j]=(16-j)/len(pts)
        score=float(G[1]+c@H[:16]);ranked.append(dict(family=family,support=support,uncovered=len(pts),uniform_priced=score))
    ranked.sort(key=lambda row:-row['uniform_priced'])
    result=dict(scope='Heuristic numerical candidates, no universal bound',generated=len(found),prefixes=ranked)
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    for i,row in enumerate(ranked[:keep]):
        prefix=f'tmp/core210-prefix-{i}.json';measure=f'tmp/core210-measure-{i}.json'
        Path(prefix).write_text(json.dumps(row));print('CANDIDATE',i,'uncovered',row['uncovered'],'uniform',row['uniform_priced'],flush=True)
        opt=optimize(prefix,measure,iterations=60)
        row.update(measure_file=measure,optimized_priced=opt['heuristic_priced'],restricted_lp=opt['restricted_lp'])
        Path(output).write_text(json.dumps(result,separators=(',',':')))
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('output');ap.add_argument('--starts',type=int,default=200);ap.add_argument('--keep',type=int,default=8)
    a=ap.parse_args();search(a.output,a.starts,a.keep)
