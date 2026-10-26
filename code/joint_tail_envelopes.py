"""Count-indexed concentration envelopes over every normalized pair residual.

For each residual before 105 and 210, form the joint-slice/compression upper
bound H. Restrict to each optimistic final count pair and take the pointwise
maximum over all inputs with that same pair. This loses information but gives
a small catalogue covering every pair placement; no inclusion reduction is
used. Generation is not an independent audit.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from core_residual_geometry import pair_residuals
from portable_tail_certificates import CoreCapacityUpper,concentration_bounds
from compressed_tail_branch import evaluate


def run(pattern,full_even=False):
    start = time.monotonic()
    model = CoreCapacityUpper(pattern)
    _,masks,configurations = pair_residuals(pattern,canonical_even=not full_even)
    envelopes = {}
    multiplicities = {}
    for index,mask in enumerate(masks):
        H = model.upper(mask)
        n0,n1 = [s-1 for s in H.shape]
        for extra in (0,1):
            a,b = max(0,n0-1-(extra==0)),max(0,n1-1-(extra==1))
            key = a,b
            restricted = H[:a+1,:b+1]
            if key not in envelopes:
                envelopes[key] = restricted.copy()
                multiplicities[key] = 0
            else:
                np.maximum(envelopes[key],restricted,out=envelopes[key])
            multiplicities[key] += 1
        if index and index%10000==0:
            print(json.dumps(dict(stage='envelope',pattern=pattern,processed=index,total=len(masks),seconds=round(time.monotonic()-start,3))),flush=True)
    rows = [dict(counts=list(key),H_upper=H.tolist(),source_branches=multiplicities[key]) for key,H in sorted(envelopes.items())]
    pool = json.loads(Path('data/compressed-tail-branch-010.json').read_text())['pool']
    counts = np.array([r['counts'] for r in rows],dtype=np.int64)
    Ks = np.array([np.stack([concentration_bounds(np.array(r['H_upper']),s) for s in (14,16)]) for r in rows])
    values = np.array([evaluate(item,counts,Ks) for item in pool])
    best = values.max(axis=0)
    result = dict(core_pattern=pattern,full_even=full_even,pair_configurations=configurations,
                  distinct_pair_residuals=len(masks),cases=rows,best_initial_bounds=best.tolist(),
                  remaining=int(np.sum(best<=0)),seconds=round(time.monotonic()-start,3),
                  scope='generated whole-core count envelopes; independent enumeration audit required')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','best_initial_bounds')}),flush=True)
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='010')
    ap.add_argument('--full-even',action='store_true')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    Path(a.output).write_text(json.dumps(run(list(map(int,a.pattern)),a.full_even),indent=2)+'\n')
