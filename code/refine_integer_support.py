"""Integer support pricing with coarse-envelope and rational-slope pruning."""
import argparse
import json
from pathlib import Path
import numpy as np
from multitail_lp import LocalStates
from positive_state_dual import lower_bound
from support_tail_knapsack import minimum as fractional_minimum
from integer_tail_knapsack import minimum as integer_minimum


def diverse_support_states(model,parity,potentials,coefficients,denominator,per_slice=1,envelope=None):
    """Construct valid integer states, not lower bounds or optimality claims.

For varied plane tuples, convexified tail counts supply enough volume.
Rounding each count upward preserves that property because volumes are
nonnegative. The resulting penalty is evaluated at the integer counts.
"""
    assert type(per_slice) is int and per_slice>=0
    if not per_slice:return []
    if envelope is None:envelope=lower_bound(model,parity,potentials,coefficients,denominator,True,True)
    feasible=model.costs[parity]<=16*sum(model.volumes[parity])
    indices=np.flatnonzero(feasible);indices=indices[np.argsort(envelope[indices],kind='stable')]
    selected=set(map(int,indices[:32]))
    shape=tuple(min(14,q-1-int(b==parity))+1 for q,b in zip(model.high,model.pattern))
    infinity=np.iinfo(np.int64).max
    tensor=np.where(feasible,envelope,infinity).reshape(shape)
    for axis,length in enumerate(shape):
        for value in range(length):
            selection=[slice(None)]*model.n;selection[axis]=value;section=tensor[tuple(selection)]
            size=min(per_slice,section.size)
            candidates=np.argpartition(section.reshape(-1),size-1)[:size]
            for j in candidates:
                if int(section.reshape(-1)[j])==infinity:continue
                coordinates=list(np.unravel_index(int(j),section.shape));coordinates.insert(axis,value)
                selected.add(int(np.ravel_multi_index(tuple(coordinates),shape)))
    penalties=[[-int(v) for v in row] for row in potentials];result=[]
    for j in sorted(selected):
        plane=list(map(int,model.grids[parity][j]));target=int(model.costs[parity][j])
        volumes=[int(v[j]) for v in model.volumes[parity]]
        plane_cost=sum(penalties[h][k] for h,k in enumerate(plane))
        relaxed=fractional_minimum(target,volumes,penalties[model.n:],plane_cost);assert relaxed is not None
        tail=[-(-v.numerator//v.denominator) for v in relaxed['counts']]
        assert all(0<=k<=16 for k in tail) and sum(v*k for v,k in zip(volumes,tail))>=target
        state=plane+tail;assert model.cost(parity,state)<=0
        result.append(dict(counts=state,penalty=sum(row[k] for row,k in zip(penalties,state))))
    return sorted(result,key=lambda row:(row['penalty'],row['counts']))


def refine_point(model,parity,potentials,coefficients,denominator,diverse=0):
    envelope=lower_bound(model,parity,potentials,coefficients,denominator,True,True)
    feasible=model.costs[parity]<=16*sum(model.volumes[parity])
    indices=np.flatnonzero(feasible);indices=indices[np.argsort(envelope[indices],kind='stable')]
    penalties=[[-v for v in row] for row in potentials]
    best=None;best_state=None;checked=[];dp_calls=0
    for j in indices:
        j=int(j)
        if best is not None and int(envelope[j])>=best:break
        plane=list(map(int,model.grids[parity][j]));target=int(model.costs[parity][j])
        volumes=[int(v[j]) for v in model.volumes[parity]]
        plane_cost=sum(penalties[h][k] for h,k in enumerate(plane))
        fractional=fractional_minimum(target,volumes,penalties[model.n:],plane_cost)
        assert fractional is not None
        slope=fractional['slope'];value=fractional['value']
        row=dict(plane=plane,slope=[slope.numerator,slope.denominator])
        if best is None or -(-value.numerator//value.denominator)<best:
            exact=integer_minimum(target,volumes,penalties[model.n:]);assert exact is not None
            dp_calls+=1;row['integer_cost']=plane_cost+exact['value'];row['tail_counts']=exact['counts']
            if best is None or row['integer_cost']<best:
                best=row['integer_cost'];best_state=plane+exact['counts']
        checked.append(row)
    assert best is not None and model.cost(parity,best_state)<=0
    result=dict(lower_numerator=best,attaining_state=best_state,refined_planes=checked,integer_dp_calls=dp_calls,
                initial_envelope_minimum=int(envelope.min()),
                feasible_plane_count=int(feasible.sum()),excluded_plane_count=int((~feasible).sum()),
                scope='Exact integer event-penalty minimum over all states with F<=0, for these fixed penalties')
    if diverse:result['diverse_support']=diverse_support_states(model,parity,potentials,coefficients,denominator,diverse,envelope)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True);ap.add_argument('--point',type=int)
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());source=data['input']
    base=data['result']['pricing_history'][-1];old=data['result']['positive_state_pricing']
    assert old['mode']=='nonpositive_support'
    model=LocalStates(source['high_primes'],source['high_pattern']);results=[]
    for x,point in enumerate(source['core_residual']):
        if args.point is not None and x!=args.point:continue
        potentials=[[base['point_weight_denominator']*sum(g[x][:t+1]) for t in range(len(g[x]))] for g in base['event_multipliers']]
        row=refine_point(model,point%2,potentials,old['point_slope_numerators'][x],old['denominator'])
        row['point_index']=x;results.append(row)
        data['result']['integer_support_refinement']=dict(points=results,complete=len(results)==len(source['core_residual']),source=args.source)
        Path(args.output).write_text(json.dumps(data,separators=(',',':'))+'\n')
        print(json.dumps(dict(point_index=x,lower=row['lower_numerator']/old['denominator'],
                              refined_planes=len(row['refined_planes']),integer_dp_calls=row['integer_dp_calls'])),flush=True)
