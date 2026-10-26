"""Audit the selected 010 profile using ordinary CRT label permutations.

Checks that every ordinary pair residual with the selected histogram profile
is an image of one of the three listed residuals, then audits their paired
capacity tables and positive portable dual inequalities. Uses no optimizer.
"""
import itertools
import json
from pathlib import Path
import time
from audit_joint_tail_envelopes import ordinary_pair_residuals,histogram_profile
from audit_paired_slice_capacity import audit as audit_capacity
from audit_compressed_tail_branch import local_minima


def audit():
    start = time.monotonic()
    read = lambda name: json.loads(Path('data',name).read_text())
    catalogue = read('joint-tail-profiles-010.json')
    source = catalogue['sources'][302]
    target = (*source['counts_before'],tuple(tuple(tuple(pair) for pair in shape) for shape in source['slice_pairs']))
    capacities = [read(f'critical-010-pair{k}-coupled-capacity.json') for k in range(3)]
    crt = {(n%2,n%3,n%5,n%7):n for n in range(210)}
    images = set()
    for cap in capacities:
        assert cap['moduli']==[d for d in range(3,211) if 210%d==0]
        assert histogram_profile(cap['core_residual'])==target
        audit_capacity(cap)
        for a,b in itertools.product(itertools.permutations(range(2,5)),itertools.permutations(range(2,7))):
            p5,p7 = (0,1)+a,(0,1)+b
            images.add(tuple(sorted(crt[n%2,n%3,p5[n%5],p7[n%7]] for n in cap['core_residual'])))
    points,masks = ordinary_pair_residuals([0,1,0])
    even = sum(1<<j for j,n in enumerate(points) if n%2==0)
    matching = 0
    for mask in masks:
        if (mask&even).bit_count()!=target[0] or (mask&~even).bit_count()!=target[1]:
            continue
        R = tuple(n for j,n in enumerate(points) if mask>>j&1)
        if histogram_profile(R)==target:
            assert R in images
            matching += 1
    print(json.dumps(dict(stage='profile_completeness',ordinary_residuals=len(masks),matching=matching,
                          seconds=round(time.monotonic()-start,3))),flush=True)
    screen = read('critical-010-coupled-pair-screen.json')
    pool = screen['pool']
    local_cache = {}
    bounds = []
    for k,cap in enumerate(capacities):
        r0,r1 = cap['layer_sizes']
        for extra in (0,1):
            counts = [max(0,r0-1-(extra==0)),max(0,r1-1-(extra==1))]
            H = [line[:counts[1]+1] for line in cap['H_upper'][:counts[0]+1]]
            K = {}
            for shape in (14,16):
                for t0,t1 in itertools.product(range(17),repeat=2):
                    K[shape,t0,t1] = max(v+(a+b+max(a,b) if shape==16 else 0)-t0*a-t1*b
                                         for a,line in enumerate(H) for b,v in enumerate(line))
            for j,pattern in enumerate(itertools.product((0,1),repeat=4)):
                candidates = []
                for index,item in enumerate(pool):
                    c = item['certificate']
                    points_claim = item['point_multipliers_by_pattern'][j]
                    numerator = sum(n*p for n,p in zip(counts,points_claim))
                    numerator += sum(row['multiplier']*K[16 if row['axis']>=4 else 14,*row['thresholds']]
                                     for row in c['tail_multipliers'])
                    candidates.append(-(-numerator//c['denominator']))
                index = max(range(len(pool)),key=lambda i:candidates[i])
                assert candidates[index]>0,(k,extra,pattern,candidates[index])
                key = index,j
                if key not in local_cache:
                    actual,checked = local_minima(pool[index]['certificate'],pattern)
                    assert actual==pool[index]['point_multipliers_by_pattern'][j]
                    local_cache[key] = checked
                bounds.append(dict(pair=k,extra=extra,high_pattern=pattern,dual_index=index,lower=candidates[index]))
    return dict(verified=True,core_pattern=[0,1,0],histogram_profile_index=302,
                ordinary_pair_residuals=len(masks),matching_profile_residuals=matching,
                anchor_preserving_images=len(images),pair_representatives=3,
                bounds=bounds,uniform_lower=min(r['lower'] for r in bounds),
                local_checks=len(local_cache),plane_tuples_checked=sum(local_cache.values()),
                seconds=round(time.monotonic()-start,3),
                scope='all cores in one specified joint histogram profile; not the entire 010 branch')


if __name__=='__main__':
    result = audit()
    Path('data/audit-critical-010-coupled-profile.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='bounds'}),flush=True)
