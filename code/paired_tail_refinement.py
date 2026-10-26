"""Generate paired-capacity refinements of the unresolved 010 profiles.

All orbit reductions here are proposals until checked against ordinary
residuals. Stores explicit matching witnesses for independent verification.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from core_residual_geometry import pair_residuals
from core_orbit_reduction import CoreOrbits
from joint_tail_profiles import joint_profile
from paired_slice_capacity import paired_capacity
from portable_tail_certificates import CoreCapacityUpper,concentration_bounds
from compressed_tail_branch import evaluate
from conditional_tail_lp import solve
from refine_abstract_tail_pool import portable


def run(rounds,seconds,resume,pattern=(0,1,0),target_lower=2000):
    start = time.monotonic()
    pattern = list(pattern)
    tag = ''.join(map(str,pattern))
    read = lambda name:json.loads(Path('data',name).read_text())
    path = Path(f'data/paired-tail-refinement-{tag}.json')
    if resume:
        data = json.loads(path.read_text())
        pool,capacities,cases = data['pool'],data['capacities'],data['cases']
    else:
        cat = read(f'joint-tail-profiles-{tag}.json')
        pool = read('critical-010-coupled-pair-screen.json')['pool'][:]
        initial_counts = np.array([r['counts'] for r in cat['cases']])
        initial_K = np.array([np.stack([concentration_bounds(np.array(r['H_upper']),s) for s in (14,16)]) for r in cat['cases']])
        initial_best = np.max([evaluate(p,initial_counts,initial_K) for p in pool],axis=0)
        selected = [i for i,r in enumerate(cat['sources']) if any(np.any(initial_best[j]<=0) for j in r['target_cases'])]
        keys = {(*r['counts_before'],tuple(tuple(tuple(p) for p in shape) for shape in r['slice_pairs'])):i
                for i,r in enumerate(cat['sources']) if i in selected}
        model = CoreCapacityUpper(pattern)
        _,masks,_ = pair_residuals(pattern)
        matching = [m for m in masks if joint_profile(model,m) in keys]
        orbits = CoreOrbits(pattern)
        canonical = sorted(set(orbits.canonical(matching)))
        print(json.dumps(dict(stage='selected',profiles=len(selected),raw_masks=len(matching),orbits=len(canonical))),flush=True)
        capacities,cases = [],[]
        for i,mask in enumerate(canonical):
            R = [n for j,n in enumerate(model.numbers) if mask>>j&1]
            cap = paired_capacity(R,[d for d in range(3,211) if 210%d==0])
            cap['profile_index'] = keys[joint_profile(model,mask)]
            cap['joint_upper'] = model.upper(mask).tolist()
            capacities.append(cap)
            H = np.minimum(np.array(cap['H_upper']),np.array(cap['joint_upper']))
            r0,r1 = cap['layer_sizes']
            for extra in (0,1):
                a,b = max(0,r0-1-(extra==0)),max(0,r1-1-(extra==1))
                cases.append(dict(capacity_index=i,extra_parity=extra,counts=[a,b],H_upper=H[:a+1,:b+1].tolist()))
            if (i+1)%20==0:
                print(json.dumps(dict(stage='capacities',done=i+1,total=len(canonical),seconds=round(time.monotonic()-start,2))),flush=True)
        data = dict(core_pattern=pattern,source_profiles=selected,raw_matching_masks=len(matching),
                    capacities=capacities,cases=cases,pool=pool,attempts=[],independent_audit=False)
    counts = np.array([r['counts'] for r in cases])
    Ks = np.array([np.stack([concentration_bounds(np.array(r['H_upper']),s) for s in (14,16)]) for r in cases])
    values = [evaluate(p,counts,Ks) for p in pool]
    patterns = list(itertools.product((0,1),repeat=4))
    blocked = {(r['case_index'],r['pattern_index']) for r in data['attempts'] if r.get('negative_primal')}
    for iteration in range(rounds+1):
        matrix = np.array(values)
        best = matrix.max(axis=0)
        data.update(pool=pool,best_bounds=best.tolist(),selection=matrix.argmax(axis=0).tolist(),
                    remaining=int(np.sum(best<=0)),minimum=int(best.min()),seconds=round(time.monotonic()-start,3))
        path.write_text(json.dumps(data,separators=(',',':'))+'\n')
        print(json.dumps(dict(stage='bounds',pool=len(pool),remaining=data['remaining'],minimum=data['minimum'],seconds=data['seconds'])),flush=True)
        if iteration==rounds or np.all(best>0):
            break
        eligible = [(int(best[i,j]),i,j) for i,j in zip(*np.where(best<=0)) if (int(i),int(j)) not in blocked]
        if not eligible:
            break
        _,i,j = min(eligible)
        i,j = int(i),int(j)
        row = cases[i]
        source = dict(counts=row['counts'],H_upper=row['H_upper'],high_pattern=list(patterns[j]))
        inc = dict(abstract_counts=row['counts'],high_primes=[11,13,17,19],high_pattern=list(patterns[j]),choices=[])
        result = solve(inc,dict(two_layer_capacities=row['H_upper']),seconds,aggregate=True,
                       stop_after_positive=True,target_lower=target_lower,verbose=False,stop_on_negative_primal=True)
        attempt = dict(case_index=i,pattern_index=j,input=source,result=result,
                       negative_primal=bool(result['primal_relaxation_witness']))
        Path(f'data/paired-tail-{tag}-attempt-{len(data["attempts"])}.json').write_text(json.dumps(attempt,separators=(',',':'))+'\n')
        data['attempts'].append(attempt)
        if attempt['negative_primal']:
            blocked.add((i,j))
        item,_,_ = portable(result['certificate'],source)
        pool.append(item)
        values.append(evaluate(item,counts,Ks))
        print(json.dumps(dict(stage='solved',case=i,pattern=patterns[j],lower=result['certificate']['integer_lower'],
                              negative_primal=attempt['negative_primal'])),flush=True)
    return data


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rounds',type=int,default=8)
    ap.add_argument('--time',type=float,default=45)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--pattern',default='010')
    ap.add_argument('--target-lower',type=int,default=2000)
    a = ap.parse_args()
    run(a.rounds,a.time,a.resume,list(map(int,a.pattern)),a.target_lower)
