"""Certify selected residue cases; these cases need not cover a whole group."""
import argparse
import json
from pathlib import Path
from branched_event_bound import shifted
from event_score_majorant import optimized_upper_bound
from pair_event_majorant import upper_bound


def certify(data,probe,cliques=False):
    source=data['input'];assert source==probe['input']
    events=[row['events'] for row in data['result']['pricing_history'][-1]['group_pricing']]
    result=[];n=len(source['high_primes']);R=source['core_residual']
    for row in probe['comparisons']:
        h=row['axis'];shape=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        remaining,reduced,constant=shifted(R,shape,events[h],row['fixed'])
        pairs=[(e['point'],e['threshold']) for e in reduced];weights=[e['weight'] for e in reduced]
        if cliques:
            from clique_event_majorant import upper_bound as clique_bound
            old=upper_bound(R,remaining,pairs,weights);new=clique_bound(R,remaining,pairs,weights)
        else:
            old=optimized_upper_bound(R,remaining,pairs,weights);new=upper_bound(R,remaining,pairs,weights)
        result.append(dict(axis=h,fixed=row['fixed'],constant=constant,old_certificate=old,new_certificate=new,
                           old_bound=constant+old['bound'],new_bound=constant+new['bound']))
        print(json.dumps(dict(axis=h,case=len(result),old_bound=constant+old['bound'],new_bound=constant+new['bound'])),flush=True)
    return dict(input=source,group_events=events,cases=result,
                scope='Certified upper bounds for the listed partial assignments only; no group-wide coverage claim')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','probe','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--cliques',action='store_true')
    args=ap.parse_args();result=certify(json.loads(Path(args.source).read_text()),json.loads(Path(args.probe).read_text()),args.cliques)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps(dict(cases=len(result['cases']),improved=sum(r['new_bound']<r['old_bound'] for r in result['cases']))))
