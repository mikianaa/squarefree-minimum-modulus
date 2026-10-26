"""Check core-profile symmetry and dominance against unquotiented even lines.

The odd-line quotient is common to both implementations; its separate full
enumeration audit is in test_research_models.py. No solver is used here.
"""
import itertools
import json
from pathlib import Path
import time
from core_residual_geometry import pair_residuals, optimistic_profiles, easiest_profiles, dominates_for_cover


def audit(pattern):
    start = time.monotonic()
    points,masks,total = pair_residuals(pattern)
    profiles = optimistic_profiles(points,masks)
    other_points,other_masks,other_total = pair_residuals(pattern,canonical_even=False)
    other_profiles = optimistic_profiles(other_points,other_masks)
    assert points == other_points and profiles == other_profiles
    easiest = easiest_profiles(profiles)
    assert all(any(dominates_for_cover(e,p) for e in easiest) for p in profiles)
    assert all(not dominates_for_cover(e,p) for e in easiest for p in easiest if e != p)
    expected = {(0,0,0):33,(0,0,1):28,(0,1,0):26,(0,1,1):25}
    assert min(m.bit_count() for m in masks) == expected[pattern]
    assert min(m.bit_count() for m in other_masks) == expected[pattern]
    result = dict(pattern=pattern,canonical_configurations=total,full_even_configurations=other_total,
                  canonical_residuals=len(masks),full_even_residuals=len(other_masks),
                  profiles=len(profiles),easiest_count=len(easiest),easiest_profiles=easiest,
                  minimum_pair_residual=expected[pattern],profile_sets_equal=True,
                  every_profile_dominated=True,seconds=round(time.monotonic()-start,3),
                  scope='optimistic cover-feasibility profiles; not realizable residuals')
    return result


if __name__ == '__main__':
    results = []
    for tail in itertools.product((0,1),repeat=2):
        result = audit((0,)+tail)
        results.append(result)
        print(json.dumps({k:v for k,v in result.items() if k!='easiest_profiles'}),flush=True)
        Path('data/core-geometry-profile-audit.json').write_text(json.dumps(results,indent=2)+'\n')
