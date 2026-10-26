"""Repair numerical case scenarios by mixing count distributions with their minima.

This preserves mean-resource feasibility and enforces pair inequalities to
numerical tolerance. These are exploratory states, not rational certificates.
"""
import argparse
import copy
import itertools
import json
from pathlib import Path
import numpy as np


def maximum_pair_score(points,shapes,fixed,distributions):
    N=len(shapes);remaining=[d for d in shapes if d not in dict(fixed)];base=int(1 in remaining)
    hits=[sum(p%d==a for d,a in fixed) for p in points]
    assert len(distributions)==len(points)
    hist=np.zeros((len(points),N+1))
    for x,row in enumerate(distributions):
        assert len(row)==len(remaining)-base+1 and min(row)>=-1e-10 and abs(sum(row)-1)<1e-8
        hist[x,hits[x]+base:hits[x]+len(remaining)+1]=row
    tails=np.cumsum(hist[:,::-1],axis=1)[:,::-1]
    pairs=list(itertools.combinations(range(len(points)),2))
    if not pairs or not N:return 0.0
    xs=np.asarray([x for x,y in pairs]);ys=np.asarray([y for x,y in pairs])
    caps=np.asarray([hits[x]+hits[y]+len(remaining)+sum(points[x]%d==points[y]%d for d in remaining) for x,y in pairs])
    ts=np.arange(1,N+1)[None,:];us=caps[:,None]+1-ts;valid=(us>=1)&(us<=N)
    scores=tails[xs[:,None],ts]+tails[ys[:,None],np.clip(us,0,N)]
    return float(np.where(valid,scores,0).max())


def repair(points,shapes,rewards,scenario):
    fixed=scenario['fixed'];assert len(dict(fixed))==len(fixed) and all(d in shapes and 0<=a<d for d,a in fixed)
    old=maximum_pair_score(points,shapes,fixed,scenario['count_distributions'])
    factor=(1-1e-12)/old if old>1+1e-10 else 1.0
    result=copy.deepcopy(scenario)
    if factor<1:
        result['count_distributions']=[[1-factor+factor*row[0]]+[factor*v for v in row[1:]] for row in scenario['count_distributions']]
    new=maximum_pair_score(points,shapes,fixed,result['count_distributions']);assert new<=1+1e-9
    remaining=[d for d in shapes if d not in dict(fixed)];base=int(1 in remaining)
    hits=[sum(p%d==a for d,a in fixed) for p in points]
    result['payoffs']=[sum(v*rewards[x][hits[x]+base+j] for j,v in enumerate(row)) for x,row in enumerate(result['count_distributions'])]
    for key in ('objective','numerical_dual'):result.pop(key,None)
    result.update(pair_compatibility=True,pair_repair_factor=factor,pair_check_maximum=new,
                  scope='Numerically pair-compatible repaired state; not an optimum, rational certificate, or actual assignment')
    return result,old


def repair_pool(data):
    source=data['input'];R=source['core_residual'];n=len(source['high_primes'])
    output=copy.deepcopy(data);changed=0;total=0;largest=0.;smallest=1.
    for group in output['group_scenarios']:
        h=group['axis'];shapes=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        rewards=[[-sum(row[:k+1]) for k in range(len(row))] for row in data['baseline_event_multipliers'][h]]
        for j,s in enumerate(group['scenarios']):
            new,old=repair(R,shapes,rewards,s);group['scenarios'][j]=new
            total+=1;changed+=new['pair_repair_factor']<1;largest=max(largest,old);smallest=min(smallest,new['pair_repair_factor'])
    output['pair_repair']=dict(scenarios=total,repaired=changed,maximum_old_pair_score=largest,minimum_retained_mass=smallest,
                               scope='Numerical compatibility repair, preserving resource means by reducing counts; no proof or optimality claim')
    output['scope']='Pair-compatible numerical scenario pool; excludes the superseded incompatible rows'
    return output


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args();result=repair_pool(json.loads(Path(args.source).read_text()))
    result['pair_repair']['source']=args.source
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps(result['pair_repair']))
