"""Locate missing local support constraints; this diagnostic is not a proof."""
import argparse
import json
from pathlib import Path
import numpy as np
from optimize_threshold_support import valid_state_pool


def diagnose(candidate,seed,known,share_parity=False):
    source=candidate['input'];assert source==seed['input']==known['input']
    prior=seed['result']['support_seed']+known['result']['point_columns']
    states=valid_state_pool(source,prior,share_parity)
    c=candidate['result']['pricing_history'][-1]
    exact=candidate['result']['integer_support_refinement']['points']
    rewards=[-np.cumsum(np.asarray(axis,dtype=np.int64),axis=1)*c['point_weight_denominator']
             for axis in c['event_multipliers']]
    rows=[]
    for x,point in enumerate(source['core_residual']):
        costs=sum(axis[x,states[x][:,h]] for h,axis in enumerate(rewards))
        index=int(costs.argmin());minimum=int(costs[index])
        priced=next(row for row in exact if row['point_index']==x)
        assert minimum>=priced['lower_numerator']
        rows.append(dict(point_index=x,point=point,prior_state_count=len(states[x]),
            prior_minimum_numerator=minimum,exact_minimum_numerator=priced['lower_numerator'],
            gap_numerator=minimum-priced['lower_numerator'],
            prior_attaining_state=list(map(int,states[x][index])),
            new_attaining_state=priced['attaining_state']))
    rows.sort(key=lambda row:(-row['gap_numerator'],row['point_index']))
    return dict(input=source,denominator=c['denominator'],points=rows,shared_parity_states=share_parity,
        total_gap_numerator=sum(row['gap_numerator'] for row in rows),
        scope='Cost difference for the same saved weights between the supplied prior local pool and the saved exact local minima. Exact minima require their separate audit; no covering or method barrier is claimed.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('candidate','seed','known','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--share-parity',action='store_true')
    args=ap.parse_args()
    result=diagnose(*[json.loads(Path(getattr(args,key)).read_text()) for key in ('candidate','seed','known')],args.share_parity)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(total_gap=result['total_gap_numerator']/result['denominator'],
        largest=[dict(point=row['point'],gap=row['gap_numerator']/result['denominator']) for row in result['points'][:5]])))
