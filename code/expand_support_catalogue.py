"""Price new actual core assignments against fixed integer event weights.

Solver upper bounds are not proof inputs. New assignments are explicit and
shared across identical resource families; old column indices are preserved.
"""
import argparse
import copy
import json
from pathlib import Path
from coherent_group_catalogue import Catalogue
from group_event_coherence import maximize
from retarget_support_certificate import retarget


def expand(baseline,certificate,seconds=1):
    assert baseline['input']==certificate['input']
    source=certificate['input'];R=source['core_residual'];n=len(source['high_primes'])
    old=certificate['result'];c=old['pricing_history'][-1]
    assert c['point_weight_denominator']==1 and set(c['point_weights'])=={1}
    catalogue=Catalogue(R);catalogue.add_result(source,old)
    diagnostics=[]
    for h,row in enumerate(c['group_pricing']):
        shape=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        events=[(e['point'],e['threshold']) for e in row['events']];weights=[e['weight'] for e in row['events']]
        known=catalogue.price(shape,events,weights,1)[0]
        suggestion=maximize(R,shape,events,weights,seconds,collect=True,hint=known['residues'])
        added=0
        for candidate in suggestion.get('candidate_assignments',[]):added+=catalogue.add(shape,candidate['residues'])
        current=catalogue.price(shape,events,weights,1)[0]
        diagnostics.append(dict(axis=h,status=suggestion['status'],new_profiles=added,
                                previous_score=known['objective'],found_score=current['objective'],
                                numerical_solver_upper_bound=suggestion.get('upper_bound'),
                                solver_seconds_requested=seconds))
        print(json.dumps(diagnostics[-1]),flush=True)
    columns=[];new_cache=[];attaining=[]
    cache=old.get('threshold_support_optimization',{}).get('group_constraint_rows')
    for h,old_columns in enumerate(old['group_columns']):
        shape=tuple(d for d in range(3 if h<n else 1,211) if 210%d==0)
        rows=copy.deepcopy(old_columns)
        seen={tuple(sum(p%d==a for d,a in zip(shape,column)) for p in R) for column in rows}
        for counts,residues in catalogue.rows[shape]:
            if counts not in seen:rows.append(list(residues));seen.add(counts)
        columns.append(rows)
        if cache:
            new_cache.append([j if j<len(old_columns) else j+len(rows)-len(old_columns) for j in cache[h]])
        row=c['group_pricing'][h]
        attaining.append(catalogue.price(shape,[(e['point'],e['threshold']) for e in row['events']],
                                         [e['weight'] for e in row['events']],1)[0]['residues'])
    pool=dict(input=source,result=dict(group_columns=columns))
    continuation=retarget(certificate,pool)
    if cache:continuation['result']['threshold_support_optimization']=dict(group_constraint_rows=new_cache)
    base=baseline['result']['pricing_history'][-1]
    search_source=dict(input=source,result=dict(group_columns=columns,pricing_history=[dict(
        event_scale=base.get('event_scale',1),event_multipliers=base['event_multipliers'])]),
        scope='Search input preserving original reward directions with an expanded explicit catalogue; not a proof certificate.')
    witness=dict(input=source,denominator=c['denominator'],event_multipliers=c['event_multipliers'],
                 point_states=[dict(point_index=r['point_index'],counts=r['attaining_state'])
                               for r in old['integer_support_refinement']['points']],group_assignments=attaining,
                 claimed_score_difference=continuation['result']['retargeted_bounds']['pool_numerator'])
    report=dict(old_configurations=sum(map(len,old['group_columns'])),new_configurations=sum(map(len,columns)),
                new_distinct_profiles=sum(r['new_profiles'] for r in diagnostics),diagnostics=diagnostics,
                retargeted_bounds=continuation['result']['retargeted_bounds'],
                solver_upper_bounds_used_as_proofs=False,global_problem_settled=False)
    return search_source,continuation,witness,report


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline',required=True);ap.add_argument('--certificate',required=True)
    ap.add_argument('--prefix',required=True);ap.add_argument('--time',type=float,default=1)
    args=ap.parse_args()
    outputs=expand(json.loads(Path(args.baseline).read_text()),json.loads(Path(args.certificate).read_text()),args.time)
    for suffix,obj in zip(('source','continuation','failure-witness','report'),outputs):
        obj['provenance']=dict(baseline=args.baseline,certificate=args.certificate)
        Path(args.prefix+'-'+suffix+'.json').write_text(json.dumps(obj,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in outputs[-1].items() if k!='diagnostics'}))
