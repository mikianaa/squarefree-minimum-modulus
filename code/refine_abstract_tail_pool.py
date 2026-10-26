"""Find reusable aggregate tail duals for the remaining 001 capacity inputs.

This is certificate generation only. A positive generated bound needs an
independent point-inequality and input-capacity audit before a theorem claim.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from conditional_tail_lp import build,price_point,solve
from portable_tail_certificates import concentration_bounds


def portable(certificate,source):
    common = [[np.zeros(15 if h<4 else 17,dtype=np.int64) for h in range(15)] for e in (0,1)]
    weights = np.zeros((2,17,17),dtype=np.int64)
    for row in certificate['tail_multipliers']:
        h = row['axis']
        t0,t1 = row['thresholds']
        weights[int(h>=4),t0,t1] += row['multiplier']
        for e,t in enumerate((t0,t1)):
            common[e][h] += row['multiplier']*np.maximum(0,np.arange(len(common[e][h]))-t)
    points = []
    for pattern in itertools.product((0,1),repeat=4):
        incumbent = dict(abstract_counts=source['counts'],high_primes=[11,13,17,19],high_pattern=list(pattern),choices=[])
        _,_,_,_,_,grids,costs,volumes = build(incumbent,dict(two_layer_capacities=source['H_upper']))
        points.append([price_point(e,[0]*15,common,grids,costs,volumes,certificate['denominator'])[0] for e in (0,1)])
    return dict(certificate=certificate,point_multipliers_by_pattern=points),weights,np.array(points,dtype=np.int64)


def run(rounds,seconds,resume=False,target_lower=5000):
    start = time.monotonic()
    read = lambda name: json.loads(Path('data',name).read_text())
    inputs = read('anchor001-unresolved-capacity-inputs.json')['cases']
    patterns = list(itertools.product((0,1),repeat=4))
    indices = np.array([patterns.index(tuple(r['high_pattern'])) for r in inputs])
    counts = np.array([r['counts'] for r in inputs],dtype=np.int64)
    Ks = np.array([np.stack([concentration_bounds(np.array(r['H_upper'],dtype=np.int64),s) for s in (14,16)]) for r in inputs])
    best = np.array([r['old_lower'] for r in inputs],dtype=np.int64)
    chosen = [-1]*len(inputs)
    pool = []
    if resume:
        pool = read('anchor001-refined-dual-pool.json')['pool']
        for index,item in enumerate(pool):
            weights = np.zeros((2,17,17),dtype=np.int64)
            for row in item['certificate']['tail_multipliers']:
                weights[int(row['axis']>=4),*row['thresholds']] += row['multiplier']
            points = np.array(item['point_multipliers_by_pattern'],dtype=np.int64)
            numerator = np.sum(points[indices]*counts,axis=1)+np.sum(Ks*weights[None,:,:,:],axis=(1,2,3))
            values = -((-numerator)//item['certificate']['denominator'])
            for j,value in enumerate(values):
                if value>best[j]:
                    best[j],chosen[j] = value,index
        print(json.dumps(dict(stage='resumed',existing_duals=len(pool),remaining=int(np.sum(best<=0)))),flush=True)
    for iteration in range(rounds):
        if np.all(best>0):
            break
        if not pool:
            initial = read('anchor001-worst-capacity-refinement.json')
            certificate,source = initial['result']['certificate'],initial['input']
        else:
            source = inputs[int(np.argmin(best))]
            incumbent = dict(abstract_counts=source['counts'],high_primes=[11,13,17,19],high_pattern=source['high_pattern'],choices=[])
            result = solve(incumbent,dict(two_layer_capacities=source['H_upper']),seconds,aggregate=True,
                           stop_after_positive=True,target_lower=target_lower,verbose=False)
            if result['certificate'] is None:
                break
            certificate = result['certificate']
        item,weights,point_values = portable(certificate,source)
        numerator = np.sum(point_values[indices]*counts,axis=1)+np.sum(Ks*weights[None,:,:,:],axis=(1,2,3))
        values = -((-numerator)//certificate['denominator'])
        for j,value in enumerate(values):
            if value>best[j]:
                best[j] = value
                chosen[j] = len(pool)
        item['source_input'] = source
        pool.append(item)
        report = dict(core_pattern=[0,0,1],patterns=[list(p) for p in patterns],pool=pool,
                      best_lower_bounds=best.tolist(),selected_pool_indices=chosen,
                      independent_audit=False,seconds=round(time.monotonic()-start,3),
                      scope='generated portable dual pool for remaining optimistic capacity inputs; not yet a global branch certificate')
        Path('data/anchor001-refined-dual-pool.json').write_text(json.dumps(report,indent=2)+'\n')
        remaining = int(np.sum(best<=0))
        print(json.dumps(dict(round=len(pool),remaining_inputs=remaining,total_inputs=len(inputs),minimum_lower=int(best.min()),seconds=report['seconds'])),flush=True)
        if not remaining or iteration+1==rounds:
            break
    return dict(pool=pool,best_lower_bounds=best.tolist())


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rounds',type=int,default=8)
    ap.add_argument('--time',type=float,default=30)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--target-lower',type=int,default=5000)
    a = ap.parse_args()
    run(a.rounds,a.time,a.resume,a.target_lower)
