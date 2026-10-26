"""Generate a whole-branch certificate using only the core count frontier.

Core-specific subset geometry is deliberately unnecessary here: each frontier
pair uses the sharp compression capacity for arbitrary subsets of the anchor
boxes. Certificates are proposals until audit_compressed_tail_branch passes.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from conditional_tail_lp import solve
from portable_tail_certificates import CoreCapacityUpper,concentration_bounds
from refine_abstract_tail_pool import portable


def evaluate(item,counts,Ks):
    c = item['certificate']
    weights = np.zeros((2,17,17),dtype=np.int64)
    for row in c['tail_multipliers']:
        weights[int(row['axis']>=4),*row['thresholds']] += row['multiplier']
    nums = counts@np.array(item['point_multipliers_by_pattern'],dtype=np.int64).T
    nums += np.sum(Ks*weights[None,:,:,:],axis=(1,2,3))[:,None]
    return -((-nums)//c['denominator'])


def run(pattern,rounds,seconds,resume=False,envelopes=None):
    start = time.monotonic()
    read = lambda name: json.loads(Path('data',name).read_text())
    tag = ''.join(map(str,pattern))
    path = Path(f'data/{"joint" if envelopes else "compressed"}-tail-branch-{tag}.json')
    if envelopes:
        source = json.loads(Path(envelopes).read_text())
        assert source['core_pattern']==pattern
        frontier = [r['counts'] for r in source['cases']]
        tables = [np.array(r['H_upper'],dtype=np.int64) for r in source['cases']]
    else:
        source = next(r for r in read('three-axis-210-certificate.json')['cases'] if r['pattern']==pattern)
        frontier = source['residual_pair_frontier']
        C = CoreCapacityUpper(pattern).compression
        tables = [C[:a+1,:b+1] for a,b in frontier]
    counts = np.array(frontier,dtype=np.int64)
    Ks = np.array([np.stack([concentration_bounds(H,s) for s in (14,16)]) for H in tables])
    high_patterns = [list(p) for p in itertools.product((0,1),repeat=4)]
    if resume:
        pool = json.loads(path.read_text())['pool']
    elif envelopes and Path(f'data/compressed-tail-branch-{tag}.json').exists():
        pool = read(f'compressed-tail-branch-{tag}.json')['pool']
    else:
        pool = read('anchor001-refined-dual-pool.json')['pool'][:]
        prior_path = Path('data/compressed-tail-branch-001.json')
        if pattern!=[0,0,1] and prior_path.exists():
            prior = json.loads(prior_path.read_text())
            if prior.get('complete_positive_bounds'):
                pool.extend(prior['pool'])
        initial = read('anchor001-compression-refinement.json')
        item,_,_ = portable(initial['result']['certificate'],initial['input'])
        pool.append(item)
    values = [evaluate(item,counts,Ks) for item in pool]
    for iteration in range(rounds+1):
        matrix = np.array(values)
        best = matrix.max(axis=0)
        complete = bool(np.all(best>0))
        report = dict(core_pattern=pattern,period=9699690,frontier=frontier,high_patterns=high_patterns,
                      capacity_mode='joint_envelope' if envelopes else 'compression',capacity_source=envelopes,
                      pool=pool,best_bounds=best.tolist(),complete_positive_bounds=complete,
                      independent_audit=False,seconds=round(time.monotonic()-start,3),
                      scope='generated whole-branch compression certificate; requires independent audit')
        if complete:
            # Reduce audit work by selecting a small set of certificates that
            # together cover every count/pattern cell with a positive bound.
            uncovered = np.ones(best.shape,dtype=bool)
            selected = []
            while np.any(uncovered):
                gain = np.sum((matrix>0)&uncovered[None,:,:],axis=(1,2))
                j = int(gain.argmax())
                assert gain[j]>0
                selected.append(j)
                uncovered &= matrix[j]<=0
            restricted = matrix[selected]
            selection = restricted.argmax(axis=0)
            report['pool'] = [pool[j] for j in selected]
            report['selection'] = selection.tolist()
            report['bounds'] = restricted.max(axis=0).tolist()
            report['uniform_lower'] = int(restricted.max(axis=0).min())
        path.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(pattern=pattern,pool=len(pool),remaining=int(np.sum(best<=0)),minimum=int(best.min()),seconds=report['seconds'])),flush=True)
        if complete or iteration==rounds:
            return report
        i,j = map(int,np.unravel_index(best.argmin(),best.shape))
        input_row = dict(counts=frontier[i],high_pattern=high_patterns[j],H_upper=tables[i].tolist())
        incumbent = dict(abstract_counts=frontier[i],high_primes=[11,13,17,19],high_pattern=high_patterns[j],choices=[])
        result = solve(incumbent,dict(two_layer_capacities=tables[i].tolist()),seconds,aggregate=True,
                       stop_after_positive=True,target_lower=5000,verbose=False)
        if result['certificate'] is None:
            return report
        item,_,_ = portable(result['certificate'],input_row)
        item['source_input'] = input_row
        item['solver_summary'] = {k:result[k] for k in ('iterations','converged_to_tolerance','master_status','seconds')}
        pool.append(item)
        values.append(evaluate(item,counts,Ks))


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='001')
    ap.add_argument('--rounds',type=int,default=12)
    ap.add_argument('--time',type=float,default=30)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--envelopes')
    a = ap.parse_args()
    run(list(map(int,a.pattern)),a.rounds,a.time,a.resume,a.envelopes)
