"""Refine finite slope grids using exact fractional tail pricing.

Stored rational slopes are proof candidates for independent verification.
Fractional tail choices are never claimed to be covering-system states.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from multitail_lp import LocalStates
from positive_state_dual import lower_bound
from support_tail_knapsack import minimum


def refine_point(model,parity,potentials,coefficients,denominator):
    envelope=lower_bound(model,parity,potentials,coefficients,denominator,True,True)
    feasible=model.costs[parity]<=16*sum(model.volumes[parity])
    indices=np.flatnonzero(feasible)
    indices=indices[np.argsort(envelope[indices],kind='stable')]
    penalties=[[-v for v in row] for row in potentials]
    best=None;checked=[]
    for j in indices:
        j=int(j)
        if best is not None and int(envelope[j])>=best.numerator//best.denominator:break
        plane=list(map(int,model.grids[parity][j]))
        volumes=[int(v[j]) for v in model.volumes[parity]]
        price=minimum(int(model.costs[parity][j]),volumes,penalties[model.n:],
                      sum(penalties[h][k] for h,k in enumerate(plane)))
        assert price is not None
        if best is None or price['value']<best:best=price['value']
        slope=price['slope']
        checked.append(dict(plane=plane,slope=[slope.numerator,slope.denominator]))
    assert best is not None
    bound=best.numerator//best.denominator
    return dict(lower_numerator=bound,minimum_fraction=[best.numerator,best.denominator],
                refined_planes=checked,feasible_plane_count=int(feasible.sum()),
                excluded_plane_count=int((~feasible).sum()),
                scope='Floor of the convexified tail support minimum; rational slopes certify the lower bound')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--point',type=int)
    args=ap.parse_args();data=json.loads(Path(args.source).read_text())
    source=data['input'];base=data['result']['pricing_history'][-1];old=data['result']['positive_state_pricing']
    assert old['mode']=='nonpositive_support'
    model=LocalStates(source['high_primes'],source['high_pattern']);results=[]
    for x,point in enumerate(source['core_residual']):
        if args.point is not None and x!=args.point:continue
        potentials=[[base['point_weight_denominator']*sum(g[x][:t+1]) for t in range(len(g[x]))] for g in base['event_multipliers']]
        # Retain the original grid for provenance and independent checks.
        row=refine_point(model,point%2,potentials,old['point_slope_numerators'][x],old['denominator'])
        row['point_index']=x;assert row['lower_numerator']>=old['point_multipliers'][x]
        results.append(row)
        data['result']['refined_support_slopes']=dict(points=results,complete=len(results)==len(source['core_residual']),
                                                      source=args.source,independently_audited=False)
        Path(args.output).write_text(json.dumps(data,separators=(',',':'))+'\n')
        print(json.dumps(dict(point_index=x,lower=row['lower_numerator']/old['denominator'],
                              gain=(row['lower_numerator']-old['point_multipliers'][x])/old['denominator'],
                              refined_planes=len(row['refined_planes']))),flush=True)
