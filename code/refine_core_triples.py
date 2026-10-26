"""Resolve critical optimistic profiles by explicitly placing 105 and 210.

All noncritical pair residuals are already excluded by dominated capacity
profiles. For critical pair residuals, enumerate the actual triple removals.
Redundant triple classes may be moved to gain coverage, so testing all minimal
residuals suffices for cover nonexistence, without a ratio-monotonicity claim
for arbitrary subsets of a residual.
"""
import argparse
import json
from pathlib import Path
import time
from core_residual_geometry import (pair_residuals,optimistic_profiles,dominates_for_cover,
    minimal_triple_residuals,triple_groups,profile,easiest_profiles)
from plane_capacity_bound import bound,refined_bound


def run(pattern):
    start = time.monotonic()
    rows = json.loads(Path('data/core-geometry-profile-audit.json').read_text())
    row = next(x for x in rows if x['pattern']==list(pattern))
    base = next(x for x in json.loads(Path('data/plane-capacity-510510.json').read_text())
                if x['core_pattern']==list(pattern) and x['high_pattern']==[1,1,1])
    refined = [x for x in json.loads(Path('data/refined-plane-capacity-510510.json').read_text())
               if x['core_pattern']==list(pattern) and x['high_pattern']==[1,1,1]]
    bad = {x['profile_index'] for x in refined if not x['excluded']}
    assert bad and bad<=set(base['surviving_indices'])
    assert set(base['surviving_indices'])=={x['profile_index'] for x in refined}
    safe = [tuple(p) for i,p in enumerate(row['easiest_profiles']) if i not in bad]
    points,masks,total = pair_residuals(pattern)
    profiles = optimistic_profiles(points,masks)
    critical = {p for p in profiles if not any(dominates_for_cover(q,p) for q in safe)}
    selected = [m for m in masks if optimistic_profiles(points,[m])&critical]
    groups = triple_groups(points)
    actual_masks = set()
    for mask in selected:
        actual_masks.update(minimal_triple_residuals(mask,groups))
    actual_profiles = {profile(mask,points) for mask in actual_masks}
    unresolved = {p for p in actual_profiles if not any(dominates_for_cover(q,p) for q in safe)}
    remaining = easiest_profiles(unresolved)
    result = dict(core_pattern=pattern,high_pattern=[1,1,1],critical_optimistic_profiles=sorted(critical),
                  critical_pair_masks=sorted(selected),actual_triple_residuals=len(actual_masks),
                  actual_profiles=len(actual_profiles),unresolved_actual_profiles=len(unresolved),
                  remaining_profiles=remaining,remaining_count=len(remaining))
    print(json.dumps({k:v for k,v in result.items() if k not in ('critical_optimistic_profiles','critical_pair_masks','remaining_profiles')}),flush=True)
    tests = []
    for p in remaining:
        simple = bound(p,(11,13,17),(1,1,1))
        test = dict(profile=p,simple=simple,excluded=simple['excluded'])
        if not simple['excluded']:
            test['refined'] = refined_bound(p,(11,13,17),(1,1,1))
            test['excluded'] = test['refined']['excluded']
        tests.append(test)
    result.update(tests=tests,all_excluded=all(t['excluded'] for t in tests),
                  seconds=round(time.monotonic()-start,3),
                  scope='explicit core triple placements, followed by necessary plane/capacity bounds')
    print(json.dumps(dict(pattern=pattern,remaining=len(tests),excluded=sum(t['excluded'] for t in tests),seconds=result['seconds'])),flush=True)
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',required=True,choices=('010','011'))
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    result = run(tuple(map(int,a.pattern)))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
