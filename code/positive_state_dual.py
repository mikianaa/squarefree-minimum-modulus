"""Strengthen a saved dual using the positive cost of each individual state.

For each plane tuple, maximize finitely many affine lower bounds to max(0,F),
after independently minimizing the tail counts. This is a certified lower
bound, not an exact solution of the positive-state pricing problem.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from multitail_lp import LocalStates


def lower_bound(model,parity,potentials,coefficients,denominator,support_only=False,return_envelope=False):
    assert coefficients and all(type(v) is int and v>=0 and (support_only or v<=denominator) for v in coefficients)
    assert denominator>0
    grid = model.grids[parity]
    worst = max(coefficients)*int(model.costs[parity].max())+sum(max(map(abs,row)) for row in potentials)
    worst += sum(max(coefficients)*16*int(v.max()) for v in model.volumes[parity])
    assert worst<2**63
    plane_penalty = -sum(np.asarray(potentials[h],dtype=np.int64)[grid[:,h]] for h in range(model.n))
    envelope = np.full(len(grid),-worst,dtype=np.int64)
    for coefficient in sorted(set(coefficients)):
        values = coefficient*model.costs[parity]+plane_penalty
        for h,volume in enumerate(model.volumes[parity],model.n):
            v = np.arange(int(volume.max())+1,dtype=np.int64)
            choices = -coefficient*v[:,None]*np.arange(17,dtype=np.int64)[None,:]-np.asarray(potentials[h])[None,:]
            values += choices.min(axis=1)[volume]
        envelope = np.maximum(envelope,values)
    return envelope if return_envelope else int(envelope.min())


def strengthen(data,steps=8,support_only=False,slope_max=16):
    assert type(steps) is int and steps>0
    source = data['input'];c = data['result']['pricing_history'][-1]
    model = LocalStates(source['high_primes'],source['high_pattern'])
    denominator,weight_den = c['denominator'],c['point_weight_denominator']
    minima,slopes = [],[]
    for x,point in enumerate(source['core_residual']):
        coefficients = sorted({k*denominator//steps for k in range(steps+1)}|{c['event_scale']*c['point_weights'][x]})
        if support_only:
            assert type(slope_max) is int and slope_max>=1
            coefficients = sorted(set(coefficients)|{denominator*k for k in range(2,slope_max+1)})
        potentials = [[weight_den*sum(group[x][:t+1]) for t in range(len(group[x]))] for group in c['event_multipliers']]
        minima.append(lower_bound(model,point%2,potentials,coefficients,denominator,support_only))
        slopes.append(coefficients)
        assert minima[-1]>=c['point_multipliers'][x]
        print(json.dumps(dict(point_index=x,old_lower=c['point_multipliers'][x]/denominator,
                              strengthened_lower=minima[-1]/denominator)),flush=True)
    return dict(point_multipliers=minima,point_slope_numerators=slopes,denominator=denominator,
                mode='nonpositive_support' if support_only else 'positive_state_cost',
                unrestricted_capacity_numerator=sum(minima)-weight_den*sum(r['capacity_certificate']['bound'] for r in c['group_pricing']),
                pool_numerator=sum(minima)-weight_den*sum(c['pool_upper_bounds']) if c.get('pool_upper_bounds') else None,
                scope=('Necessary event-score excess for distributions supported on F<=0; not a cost bound on arbitrary distributions'
                       if support_only else 'Lower bound on sum of expected positive state costs; slope envelope may be weaker than exact state pricing'))


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--steps',type=int,default=8)
    ap.add_argument('--support-only',action='store_true')
    ap.add_argument('--slope-max',type=int,default=16)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    result = strengthen(data,args.steps,args.support_only,args.slope_max)
    data['result']['positive_state_pricing'] = result
    data['positive_state_provenance'] = dict(source=args.source,steps=args.steps,support_only=args.support_only,slope_max=args.slope_max)
    Path(args.output).write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(json.dumps(result),flush=True)
