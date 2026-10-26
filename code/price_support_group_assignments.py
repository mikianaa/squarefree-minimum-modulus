"""Search omitted assignments against an integer support certificate.

New assignments are concrete counterexamples to finite-pool score bounds.
CP-SAT upper bounds are exploratory and need separate proof verification.
"""
import argparse
import json
from pathlib import Path
from coherent_group_catalogue import Catalogue
from group_event_coherence import maximize


def search(data,seconds=1,checkpoint=None):
    source,result=data['input'],data['result'];base=result['pricing_history'][-1]
    n=len(source['high_primes']);catalogue=Catalogue(source['core_residual'])
    catalogue.add_result(source,result)
    configurations=[[list(r) for r in rows] for rows in result['group_columns']]
    rows=[]
    for h,priced in enumerate(base['group_pricing']):
        shape=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        events=[(e['point'],e['threshold']) for e in priced['events']]
        weights=[e['weight'] for e in priced['events']]
        known=catalogue.price(shape,events,weights)
        row=maximize(source['core_residual'],shape,events,weights,seconds,True,known[0]['residues'])
        candidates=known+row.get('candidate_assignments',[])+([row] if 'residues' in row else [])
        seen={tuple(r) for r in configurations[h]};added=0
        for candidate in candidates:
            residues=candidate['residues'];catalogue.add(shape,residues)
            if tuple(residues) not in seen:seen.add(tuple(residues));configurations[h].append(residues);added+=1
        best=max(candidate['objective'] for candidate in candidates)
        row.update(axis=h,old_pool_bound=base['pool_upper_bounds'][h],best_known_score=best,
                   known_score_gain=best-base['pool_upper_bounds'][h],added_assignments=added)
        rows.append(row)
        output=dict(input=source,result=dict(group_columns=configurations,group_pricing_search=rows,
                    catalogue=catalogue.export(),scope='Explicit new assignments; solver upper bounds are not independently certified'))
        if checkpoint:Path(checkpoint).write_text(json.dumps(output,separators=(',',':'))+'\n')
        print(json.dumps({k:row[k] for k in ('axis','status','known_score_gain','added_assignments')}),flush=True)
    return output


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--seconds',type=float,default=1)
    args=ap.parse_args();data=json.loads(Path(args.source).read_text())
    search(data,args.seconds,args.output)
