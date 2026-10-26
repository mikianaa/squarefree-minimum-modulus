"""Integer verification of point reweighting and its explicit scope barrier."""
import argparse
import json
from pathlib import Path
from audit_positive_state_dual import audit as audit_parent
from audit_weighted_coherent_dual import capacity
from audit_refined_support_slopes import audit as audit_refined


def audit(data):
    refined=any(key in data['result'] for key in ('refined_support_slopes','integer_support_refinement'))
    parent=audit_refined(data) if refined else audit_parent(data)
    source,result=data['input'],data['result']
    base=result['pricing_history'][-1];support=result['positive_state_pricing']
    assert support['mode']=='nonpositive_support'
    beta=parent['point_multipliers'] if refined else support['point_multipliers']
    c=result['support_reweighting'];weights=c['point_weights']
    R,n=source['core_residual'],len(source['high_primes'])
    assert len(weights)==len(R) and all(type(w) is int and w>=0 for w in weights)
    scale=c['point_weight_denominator'];assert type(scale) is int and scale>0
    assert c['denominator']==base['denominator']*scale
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(2**n-1)]
    assert len(c['group_pricing'])==len(shapes)
    pool,capacities=[],[]
    for h,(shape,row) in enumerate(zip(shapes,c['group_pricing'])):
        events=[dict(point=e['point'],threshold=e['threshold'],weight=e['weight']*weights[e['point']])
                for e in base['group_pricing'][h]['events'] if weights[e['point']]*e['weight']]
        assert events==row['events']
        capacities.append(capacity(R,shape,events,row['capacity_certificate']))
        values=[]
        for residues in result['group_columns'][h]:
            assert len(residues)==len(shape) and all(type(a) is int and 0<=a<d for a,d in zip(residues,shape))
            counts=[sum(v%d==a for d,a in zip(shape,residues)) for v in R]
            values.append(sum(e['weight'] for e in events if counts[e['point']]>=e['threshold']))
        pool.append(max(values));assert pool[-1]==row['pool_bound']
    demand=sum(w*b for w,b in zip(weights,beta))
    finite=demand-base['point_weight_denominator']*sum(pool)
    unrestricted=demand-base['point_weight_denominator']*sum(capacities)
    assert finite==c['pool_numerator'] and unrestricted==c['capacity_numerator']
    barrier=c.get('barrier_witness');assignments=0;excess=None
    if barrier is not None:
        den=barrier['denominator'];assert type(den) is int and den>0
        assert barrier['score_denominator']==den*base['denominator']
        scores=[0]*len(R);seen=set()
        for group in barrier['groups']:
            h=group['axis'];assert type(h) is int and 0<=h<len(shapes) and h not in seen
            seen.add(h);allowed=set(map(tuple,result['group_columns'][h]))
            assert sum(r['weight'] for r in group['assignments'])==den
            for r in group['assignments']:
                w,residues=r['weight'],r['residues']
                assert type(w) is int and w>0 and tuple(residues) in allowed
                counts=[sum(v%d==a for a,d in zip(residues,shapes[h])) for v in R]
                for e in base['group_pricing'][h]['events']:
                    if counts[e['point']]>=e['threshold']:scores[e['point']]+=w*e['weight']
                assignments+=1
        assert seen==set(range(len(shapes)))
        excess=[v*base['point_weight_denominator']-b*den for v,b in zip(scores,beta)]
        assert excess==barrier['point_score_excess_numerators'] and min(excess)>=0
    return dict(verified=True,pool_numerator=finite,capacity_numerator=unrestricted,denominator=c['denominator'],
                fixed_point_reweighting_barrier_verified=barrier is not None,
                barrier_assignments_checked=assignments,minimum_score_excess_numerator=min(excess) if excess else None,
                score_denominator=barrier['score_denominator'] if barrier else None,
                excludes_stated_pool=finite>0,excludes_unrestricted_core_pattern=unrestricted>0,
                scope='Only combinations of the fixed point event inequalities; no general optimality claim about the method')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args();report=audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
