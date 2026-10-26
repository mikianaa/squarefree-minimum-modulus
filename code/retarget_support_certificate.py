"""Reevaluate the same support inequalities on an expanded assignment pool."""
import argparse
import copy
import json
from pathlib import Path
from event_score_majorant import optimized_upper_bound
from extract_weighted_coherent_pool import extract


def retarget(data,pool,optimize_capacity=False,group_bounds=None):
    assert len(data['result']['pricing_history'])==1
    result=extract(data,pool)
    old=data['result'];new=result['result'];c=new['pricing_history'][0]
    for key in ('positive_state_pricing','integer_support_refinement','refined_support_slopes','support_seed'):
        if key in old:new[key]=copy.deepcopy(old[key])
    if optimize_capacity:
        n=len(data['input']['high_primes']);R=data['input']['core_residual']
        for h,row in enumerate(c['group_pricing']):
            shape=[d for d in range(3 if h<n else 1,211) if 210%d==0]
            events=[(e['point'],e['threshold']) for e in row['events']];weights=[e['weight'] for e in row['events']]
            row['capacity_certificate']=optimized_upper_bound(R,shape,events,weights)
    if group_bounds is not None:
        assert group_bounds['input']==data['input'] and group_bounds['complete']
        rows=group_bounds['group_bounds'];assert [r['axis'] for r in rows]==list(range(len(c['group_pricing'])))
        for row,entry in zip(c['group_pricing'],rows):
            if entry['certificate']['bound']<=row['capacity_certificate']['bound']:
                row['capacity_certificate']=copy.deepcopy(entry['certificate'])
    weight_den=c['point_weight_denominator']
    capacity=sum(r['capacity_certificate']['bound'] for r in c['group_pricing'])
    c['capacity_numerator']=sum(c['point_multipliers'])-weight_den*capacity
    coarse=new['positive_state_pricing']
    coarse['pool_numerator']=sum(coarse['point_multipliers'])-weight_den*sum(c['pool_upper_bounds'])
    coarse['unrestricted_capacity_numerator']=sum(coarse['point_multipliers'])-weight_den*capacity
    exact=new.get('integer_support_refinement') or new['refined_support_slopes']
    new['retargeted_bounds']=dict(pool_numerator=sum(r['lower_numerator'] for r in exact['points'])-weight_den*sum(c['pool_upper_bounds']),
        capacity_numerator=sum(r['lower_numerator'] for r in exact['points'])-weight_den*capacity,denominator=c['denominator'],
        configurations=sum(c['pool_sizes']),scope='Same fixed point inequalities; expanded explicit pool and separately certified unrestricted majorants')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--pool',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--optimize-capacity',action='store_true')
    ap.add_argument('--group-bounds')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text());pool=json.loads(Path(args.pool).read_text())
    bounds=json.loads(Path(args.group_bounds).read_text()) if args.group_bounds else None
    result=retarget(data,pool,args.optimize_capacity,bounds);result['provenance']=dict(source=args.source,pool=args.pool,group_bounds=args.group_bounds)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps(result['result']['retargeted_bounds']))
