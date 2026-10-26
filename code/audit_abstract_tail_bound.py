"""Independent integer dual audit with a separately supplied capacity table.

This verifies the local dual and its use of that table. The mathematical
validity of the table itself must be checked by the appropriate separate
capacity audit (exact subsets, compression, or paired slices).
"""
import argparse
import json
from pathlib import Path
from audit_compressed_tail_branch import local_minima


def audit(data,capacity):
    source,result = data['input'],data['result']
    c = result['certificate']
    counts = source['counts']
    H = source['H_upper']
    full = capacity.get('H_upper',capacity.get('two_layer_capacities'))
    assert H==[row[:counts[1]+1] for row in full[:counts[0]+1]]
    assert len(H)==counts[0]+1 and all(len(row)==counts[1]+1 for row in H)
    points,checked = local_minima(c,source['high_pattern'])
    assert points==c['point_multipliers']
    numerator = sum(n*p for n,p in zip(counts,points))
    for row in c['tail_multipliers']:
        t0,t1 = row['thresholds']
        tail = row['axis']>=4
        bound = max(v+(a+b+max(a,b) if tail else 0)-t0*a-t1*b for a,line in enumerate(H) for b,v in enumerate(line))
        assert bound==row['bound']
        numerator += row['multiplier']*bound
    assert numerator==c['numerator']
    lower = -(-numerator//c['denominator'])
    assert lower==c['integer_lower']
    return dict(verified=True,counts=counts,high_pattern=source['high_pattern'],
                exact_numerator=numerator,exact_denominator=c['denominator'],integer_lower=lower,
                positive_obstruction=lower>0,plane_tuples_checked=checked,
                scope='dual inequality conditional on the separately audited capacity table; not a whole-branch proof')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/critical-010-exact-subset-bound.json')
    ap.add_argument('--capacity',default='data/critical-010-exact-subsets.json')
    ap.add_argument('--output',default='data/audit-critical-010-exact-subset-bound.json')
    a = ap.parse_args()
    result = audit(json.loads(Path(a.source).read_text()),json.loads(Path(a.capacity).read_text()))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
