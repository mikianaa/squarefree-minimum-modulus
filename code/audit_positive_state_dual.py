"""Independent integer check of the positive-state slope-envelope bound."""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
from audit_weighted_coherent_dual import audit as audit_original


def minimum(D,potentials,coefficients,denominator,support_only=False,return_envelope=False):
    assert denominator>0 and coefficients and all(type(v) is int and v>=0 and (support_only or v<=denominator) for v in coefficients)
    n = len(D)
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    ranges = [np.arange(max(0,d-14),d+1,dtype=np.int64) for d in D]
    free = [v.reshape(-1) for v in np.meshgrid(*ranges,indexing='ij')]
    base = np.prod(np.stack(free),axis=0)
    risk = max(coefficients)*math.prod(D)+sum(max(map(abs,row)) for row in potentials)
    risk += sum(16*max(coefficients)*math.prod(D[i] for i in range(n) if i not in S) for S in groups[n:])
    assert risk<2**63
    answer = np.full(len(base),-risk,dtype=np.int64)
    for coefficient in coefficients:
        value = coefficient*base
        for h in range(n):value -= np.asarray(potentials[h],dtype=np.int64)[D[h]-free[h]]
        for h,S in enumerate(groups[n:],n):
            outside = [i for i in range(n) if i not in S]
            volumes = np.arange(math.prod(D[i] for i in outside)+1,dtype=np.int64)
            best = np.full(len(volumes),-potentials[h][0],dtype=np.int64)
            for k in range(1,17):best = np.minimum(best,-coefficient*k*volumes-potentials[h][k])
            actual = np.ones(len(base),dtype=np.int64)
            for i in outside:actual *= free[i]
            value += best[actual]
        answer = np.maximum(answer,value)
    return (answer,free,groups) if return_envelope else (int(answer.min()),len(base))


def audit(data):
    original = audit_original(data)
    source = data['input'];c = data['result']['pricing_history'][-1]
    certificate = data['result']['positive_state_pricing']
    den = c['denominator'];assert certificate['denominator']==den
    mode = certificate['mode'];assert mode in ('nonpositive_support','positive_state_cost')
    slopes = certificate['point_slope_numerators']
    assert len(slopes)==len(source['core_residual'])
    minima,checked = [],0
    for x,point in enumerate(source['core_residual']):
        D = [q-1-(b==point%2) for q,b in zip(source['high_primes'],source['high_pattern'])]
        potentials = [[c['point_weight_denominator']*sum(g[x][:t+1]) for t in range(len(g[x]))] for g in c['event_multipliers']]
        value,count = minimum(D,potentials,slopes[x],den,mode=='nonpositive_support')
        minima.append(value);checked += count
    assert minima==certificate['point_multipliers']
    gain = sum(minima)-sum(c['point_multipliers'])
    capacity = original['unrestricted_capacity_numerator']+gain
    pool = original['pool_numerator']+gain if original['pool_numerator'] is not None else None
    assert capacity==certificate['unrestricted_capacity_numerator'] and pool==certificate['pool_numerator']
    return dict(verified=True,mode=mode,denominator=den,unrestricted_capacity_numerator=capacity,pool_numerator=pool,
                gain_over_pointwise_dual_numerator=gain,plane_states_checked=checked,
                configurations_checked=original['configurations_checked'],
                excludes_unrestricted_core_pattern=capacity>0,excludes_stated_pool=pool is not None and pool>0,
                scope=('Necessary event-score excess on nonpositive support' if mode=='nonpositive_support'
                       else 'Expected positive cost of individual states')+'; same fixed core/pattern and finite pool as input')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args = ap.parse_args()
    result = audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
