"""Choose a threshold direction that escapes a specific audited mixture.

Escaping this witness is not an obstruction or a proof that the enlarged
curve family can separate every state in the full relaxation.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path


def escape(barrier,report):
    assert report['verified'] and report['curve_cone_barrier_verified']
    source=barrier['input'];n=len(source['high_primes']);G=2**n-1;P=len(source['core_residual'])
    events=[[[0]*(15 if h<n else 17) for x in range(P)] for h in range(G)]
    best={}
    for row in report['threshold_deficits']:
        key=(row['axis'],row['point_index']);gap=Fraction(*row['supply_minus_demand']);assert gap<0
        if key not in best or gap<best[key][0]:best[key]=(gap,row['threshold'])
    assert best;value=Fraction(0)
    for (h,x),(gap,t) in best.items():events[h][x][t]=-1;value-=gap
    assert 0<value<=len(best)
    return dict(input=source,result=dict(pricing_history=[dict(event_multipliers=events)]),
        escape_direction=dict(active_curves=len(best),expected_local_minus_group=[value.numerator,value.denominator],
            scope='Pure threshold curves chosen to defeat this particular audited mixture; not a covering obstruction or complete basis.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','audit','output'):ap.add_argument('--'+key,required=True)
    args=ap.parse_args();result=escape(json.loads(Path(args.source).read_text()),json.loads(Path(args.audit).read_text()))
    result['escape_direction'].update(barrier_source=args.source,audit_source=args.audit)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps(result['escape_direction']))
