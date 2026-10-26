"""Reprice saved weights with cached slopes and compare exact minima.

This benchmarks pruning; it neither changes weights nor proves a new obstruction.
"""
import argparse
import json
from pathlib import Path
import time
from cached_support_slopes import coefficients,safe_maximum
from multitail_lp import LocalStates
from refine_integer_support import refine_point


def benchmark(data,indices,limit=12):
    source=data['input'];c=data['result']['pricing_history'][-1]
    old=data['result']['positive_state_pricing']
    saved={row['point_index']:row for row in data['result']['integer_support_refinement']['points']}
    model=LocalStates(source['high_primes'],source['high_pattern']);results=[]
    for x in indices:
        parity=source['core_residual'][x]%2
        potentials=[[c['point_weight_denominator']*sum(axis[x][:k+1]) for k in range(len(axis[x]))]
                    for axis in c['event_multipliers']]
        base=old['point_slope_numerators'][x]
        selected=coefficients(base,saved[x]['refined_planes'],old['denominator'],old['denominator'],limit,
                              safe_maximum(model,parity,potentials))
        start=time.monotonic();new=refine_point(model,parity,potentials,selected,old['denominator'])
        assert new['lower_numerator']==saved[x]['lower_numerator']
        row=dict(point_index=x,point=source['core_residual'][x],coefficients=selected,
                 old_exceptions=len(saved[x]['refined_planes']),new_exceptions=len(new['refined_planes']),
                 old_dp_calls=saved[x]['integer_dp_calls'],new_dp_calls=new['integer_dp_calls'],
                 lower_numerator=new['lower_numerator'],seconds=time.monotonic()-start)
        results.append(row);print(json.dumps(row),flush=True)
    return dict(input=source,denominator=old['denominator'],points=results,
                scope='Same-weight exact local repricing agrees with saved minima; pruning benchmark, not a new covering obstruction.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--points',type=int,nargs='+',required=True);ap.add_argument('--limit',type=int,default=12)
    args=ap.parse_args();data=json.loads(Path(args.source).read_text())
    result=benchmark(data,args.points,args.limit);result['source']=args.source
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
