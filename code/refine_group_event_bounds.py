"""Repeat proof-checked whole-group separation on a fixed local core model.

Each round saves candidates, proof audits, and the strengthened model. Only
cuts with independently replayed UNSAT proofs enter the model. A positive
local certificate still needs its complete separate dual audit.
"""
import argparse
import json
from pathlib import Path
import time
from group_event_coherence import separate_axis,simplify_candidate
from audit_group_event_cut import verify_bound
from multitail_identity import solve


def save(path,data):
    Path(path).write_text(json.dumps(data,separators=(',',':'))+'\n')


def run(source,prefix,axes,rounds=2,first_round=2,candidates=None,resume=False):
    current = json.loads(Path(source).read_text())
    summary = dict(initial_source=source,rounds=[],global_problem_settled=False)
    start = time.monotonic()
    for number in range(first_round,first_round+rounds):
        if not current['result'].get('point_primal_witness'):
            break
        stem = f'{prefix}-round{number}'
        report = dict(round=number,source=source,attempts=[],certified_cuts=[])
        if resume and Path(stem+'-cuts.json').exists():
            report = json.loads(Path(stem+'-cuts.json').read_text())
            assert json.loads(Path(report['source']).read_text())==current
        initial = json.loads(Path(candidates).read_text()) if candidates and number==first_round else None
        if initial:
            assert json.loads(Path(initial['source']).read_text())==current
        for h in axes:
            saved = next((r for r in report['attempts'] if r['axis']==h),None)
            if saved and 'proof_audit' in saved:
                continue
            separated = next((r for r in initial['results'] if r['axis']==h),None) if initial else None
            if saved:
                separated = saved['separation']
            if separated is None:
                separated = separate_axis(current,h,20)
            item = dict(axis=h,separation=separated)
            if saved:
                report['attempts'].remove(saved)
            report['attempts'].append(item)
            save(stem+'-cuts.json',report)
            candidate = separated.get('cut')
            if candidate is None:
                print(json.dumps(dict(round=number,axis=h,stage='no_candidate')),flush=True)
                continue
            cut = simplify_candidate(current,candidate)
            item['simplified_cut'] = cut
            save(stem+'-cuts.json',report)
            if cut is None:
                print(json.dumps(dict(round=number,axis=h,stage='no_small_cut')),flush=True)
                continue
            shapes = [d for d in range(1 if h>=len(current['input']['high_primes']) else 3,211) if 210%d==0]
            proof_path = Path(stem+f'-axis{h}.rup.gz')
            audit = verify_bound(current['input']['core_residual'],shapes,cut,proof_path,45)
            item['proof_audit'] = audit
            if audit['verified']:
                cut['proof_audit'] = audit
                report['certified_cuts'].append(cut)
            save(stem+'-cuts.json',report)
            print(json.dumps(dict(round=number,axis=h,stage='cut_audit',verified=audit['verified'],
                                  bound=cut['bound'],events=len(cut['events']),seconds=round(time.monotonic()-start,3))),flush=True)
        if not report['certified_cuts']:
            summary['rounds'].append(dict(round=number,status='no_verified_new_cuts',report=stem+'-cuts.json'))
            save(prefix+'-summary.json',summary)
            break
        result = solve(current['input'],seconds=180,seed=current['result']['point_primal_witness'],rank_cuts=True,
                       pair_cuts=True,triple_cuts=True,subset_cuts=True,warm_start=current['result'],
                       group_cuts=report['certified_cuts'])
        current = dict(input=current['input'],result=result)
        source = stem+'-model.json'
        save(source,current)
        row = dict(round=number,status='model_completed',model=source,new_cuts=len(report['certified_cuts']),
                   total_cuts=len(result['group_event_constraints']),
                   negative_primal=bool(result['point_primal_witness']),
                   lower=result['certificate']['integer_lower'] if result['certificate'] else None)
        summary['rounds'].append(row)
        save(prefix+'-summary.json',summary)
        print(json.dumps(row),flush=True)
    return summary


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--prefix',required=True)
    ap.add_argument('--candidates')
    ap.add_argument('--axes',default='0,8,10,16,17')
    ap.add_argument('--rounds',type=int,default=2)
    ap.add_argument('--first-round',type=int,default=2)
    ap.add_argument('--resume',action='store_true')
    args = ap.parse_args()
    run(args.source,args.prefix,list(map(int,args.axes.split(','))),args.rounds,args.first_round,args.candidates,args.resume)
