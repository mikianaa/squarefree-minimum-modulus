"""Apply proved group-event bounds to every group with identical resources.

The SAT formula depends on core points, core moduli, thresholds and weights,
not on higher-prime group labels. Plane clipping only decreases these scores.
"""
import argparse
import copy
import json
from pathlib import Path


def transfer(source,cuts):
    n = len(source['high_primes'])
    result,seen = [],set()
    for cut in cuts:
        old = cut['axis']
        assert type(old) is int and 0<=old<(1<<n)-1
        for h in (range(n) if old<n else range(n,(1<<n)-1)):
            key = (h,cut['bound'],tuple((e['point'],e['threshold'],e['weight']) for e in cut['events']))
            if key in seen:
                continue
            seen.add(key)
            row = copy.deepcopy(cut)
            row['axis'] = h
            row['transferred_from_axis'] = old
            row['proof_audit']['axis'] = h
            # These describe the discovery witness, not the new group.
            for field in ('violation_numerator','denominator','pricing'):
                row.pop(field,None)
            result.append(row)
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    source,result = data['input'],data['result']
    cuts = transfer(source,result['group_event_constraints'])
    witness = result.get('point_primal_witness')
    violations = []
    if witness:
        for i,cut in enumerate(cuts):
            used = sum(r['weight']*e['weight'] for r in witness['support'] for e in cut['events']
                       if r['point_index']==e['point'] and r['counts'][cut['axis']]>=e['threshold'])
            gap = used-cut['bound']*witness['denominator']
            if gap>0:
                violations.append(dict(cut_index=i,axis=cut['axis'],gap=gap))
    report = dict(source=args.source,core_residual=source['core_residual'],cuts=cuts,
                  old_witness_violations=violations,
                  scope='same-resource transfer of previously proved inequalities; no new covering obstruction')
    Path(args.output).write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(dict(cuts=len(cuts),old_witness_violations=len(violations))))
