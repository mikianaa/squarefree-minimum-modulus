"""Screen the three selected critical pair residuals using coupled capacities.

This is a local probe, not a completeness assertion about the core catalogue.
"""
import json
import itertools
from pathlib import Path
import numpy as np
from paired_slice_capacity import paired_capacity
from portable_tail_certificates import CoreCapacityUpper, concentration_bounds
from compressed_tail_branch import evaluate
from refine_abstract_tail_pool import portable
from conditional_tail_lp import solve


def run():
    read = lambda name: json.loads(Path('data', name).read_text())
    model = CoreCapacityUpper([0,1,0])
    moduli = [d for d in range(3,211) if 210%d==0]
    pool = read('joint-tail-branch-010.json')['pool'][:]
    prior = read('critical-010-paired-bound.json')
    item,_,_ = portable(prior['result']['certificate'],prior['input'])
    pool.append(item)
    rows = []
    attempted = set()
    patterns = list(itertools.product((0,1),repeat=4))
    for k,mask in enumerate(read('critical-010-profile302.json')['pair_orbit_representatives']):
        R = [n for j,n in enumerate(model.numbers) if mask>>j&1]
        cap = paired_capacity(R,moduli)
        path = Path(f'data/critical-010-pair{k}-coupled-capacity.json')
        path.write_text(json.dumps(cap,separators=(',',':'))+'\n')
        r0,r1 = cap['layer_sizes']
        for extra in (0,1):
            a,b = max(0,r0-1-(extra==0)),max(0,r1-1-(extra==1))
            H = np.array(cap['H_upper'],dtype=np.int64)[:a+1,:b+1]
            Ks = np.stack([concentration_bounds(H,s) for s in (14,16)])[None,:,:,:]
            best = np.max([evaluate(p,np.array([[a,b]]),Ks)[0] for p in pool],axis=0)
            row = dict(pair_index=k,extra_parity=extra,counts=[a,b],capacity_source=str(path),bounds=best.tolist())
            print(json.dumps(row),flush=True)
            for iteration in range(6):
                j = int(best.argmin())
                key = (a,b,H.tobytes(),j)
                if best[j]>0 or key in attempted:
                    break
                attempted.add(key)
                pattern = list(patterns[j])
                input_row = dict(counts=[a,b],H_upper=H.tolist(),high_pattern=pattern)
                inc = dict(abstract_counts=[a,b],high_primes=[11,13,17,19],high_pattern=pattern,choices=[])
                result = solve(inc,dict(two_layer_capacities=H.tolist()),45,aggregate=True,
                               stop_after_positive=True,target_lower=1000,verbose=False,stop_on_negative_primal=True)
                output = dict(input=input_row,result=result)
                boundpath = Path(f'data/critical-010-pair{k}-extra{extra}-high{j}-coupled-bound.json')
                boundpath.write_text(json.dumps(output,separators=(',',':'))+'\n')
                row.setdefault('bound_sources',[]).append(str(boundpath))
                print(json.dumps(dict(pair=k,extra=extra,high=pattern,lower=result['certificate']['integer_lower'],
                                      primal=bool(result['primal_relaxation_witness']))),flush=True)
                item,_,_ = portable(result['certificate'],input_row)
                pool.append(item)
                best = np.maximum(best,evaluate(item,np.array([[a,b]]),Ks)[0])
                if result['primal_relaxation_witness']:
                    break
            row['bounds'] = best.tolist()
            rows.append(row)
    Path('data/critical-010-coupled-pair-screen.json').write_text(json.dumps(dict(cases=rows,pool=pool),indent=2)+'\n')


if __name__=='__main__':
    run()
