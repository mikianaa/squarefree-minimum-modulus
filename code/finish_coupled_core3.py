"""Finish the critical 011111 case using shared residue choices modulo 3."""
from collections import Counter
import itertools
import json
from pathlib import Path
import time
from core_residual_geometry import (points_for,triple_groups,minimal_triple_residuals,profile,dominates_for_cover)
from plane_capacity_bound import bound


def run():
    start = time.monotonic()
    refinement = json.loads(Path('data/core-triple-refinement-011.json').read_text())
    rows = json.loads(Path('data/core-geometry-profile-audit.json').read_text())
    row = next(x for x in rows if x['pattern']==[0,1,1])
    previous = json.loads(Path('data/refined-plane-capacity-510510.json').read_text())
    bad = {x['profile_index'] for x in previous if x['core_pattern']==[0,1,1]
           and x['high_pattern']==[1,1,1] and not x['excluded']}
    assert bad=={2}
    safe = [tuple(p) for i,p in enumerate(row['easiest_profiles']) if i not in bad]
    points = points_for((0,1,1))
    groups = triple_groups(points)
    masks = set()
    for m in refinement['critical_pair_masks']:
        masks.update(minimal_triple_residuals(m,groups))
    profiles = {m:profile(m,points) for m in masks}
    unresolved = {p for p in profiles.values() if not any(dominates_for_cover(q,p) for q in safe)}
    enhanced = Counter()
    for m,p in profiles.items():
        if p not in unresolved:
            continue
        # The parity-0 anchor box has only label 2 modulo 3. In parity 1,
        # count the two remaining modulo-3 classes separately, rather than
        # granting their maximum simultaneously with the parity-0 maximum.
        assert all(a[0]==2 for j,(e,a) in enumerate(points) if e==0 and m>>j&1)
        counts = tuple(sum(1 for j,(e,a) in enumerate(points) if m>>j&1 and e==1 and a[0]==v) for v in (1,2))
        assert sum(counts)==p[1]
        enhanced[p,counts] += 1
    tests = []
    for (p,counts),multiplicity in sorted(enhanced.items()):
        cases = []
        for choices in itertools.product((1,2),repeat=3):
            # A modulo-3 choice 0 misses the entire core and is dominated by
            # either remaining choice. Each 3*q may choose 1 or 2 independently.
            capacities = [(0,counts[0]) if c==1 else (p[0],counts[1]) for c in choices]
            result = bound(list(p),(11,13,17),(1,1,1),capacities)
            cases.append(dict(core_mod3_choices=choices,**result))
        tests.append(dict(profile=p,layer1_mod3_counts=counts,actual_residual_sets=multiplicity,
                          cases=cases,all_excluded=all(x['excluded'] for x in cases)))
    result = dict(core_pattern=[0,1,1],high_pattern=[1,1,1],enhanced_profiles=len(tests),
                  actual_residual_sets=sum(enhanced.values()),tests=tests,
                  all_excluded=all(x['all_excluded'] for x in tests),
                  minimum_gap=min(c['gap'] for x in tests for c in x['cases']),
                  seconds=round(time.monotonic()-start,3),
                  scope='exact core-residue coupling plus integer plane-capacity bounds')
    return result


if __name__=='__main__':
    result = run()
    Path('data/coupled-core3-510510.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='tests'}),flush=True)
