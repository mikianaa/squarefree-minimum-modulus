"""Ordinary-residue audit of the seven-prime core catalogue and refinements.

The canonical pair masks come from the previously audited symmetry reduction.
All capacities and final triple placements here use integers modulo 210,
independently of the coordinate-profile and triple-group producer functions.
"""
from collections import Counter
import itertools
import json
from pathlib import Path
import time
from core_residual_geometry import pair_residuals


def easier(q,p):
    for e in (0,1):
        if q[e]==0:
            continue
        if p[e]==0:
            return False
        for j in range(7):
            if q[2+2*j+e]*p[e]<p[2+2*j+e]*q[e]:
                return False
    return True


def statistics(mask,numbers):
    residues = [[n for j,n in enumerate(numbers) if mask>>j&1 and n%2==e] for e in (0,1)]
    counts = tuple(map(len,residues))
    caps = tuple(max(Counter(n%d for n in residues[e]).values(),default=0)
                 for d in (3,5,7,15,21,35,105) for e in (0,1))
    return counts+caps


def optimistic(p):
    result = []
    for extra in (0,1):
        r = [max(0,p[e]-1-int(extra==e)) for e in (0,1)]
        result.append(tuple(r)+tuple(min(r[e],p[2+2*j+e]) for j in range(7) for e in (0,1)))
    return result


def main():
    start = time.monotonic()
    source = json.loads(Path('data/core-geometry-profile-audit.json').read_text())
    planes = json.loads(Path('data/plane-capacity-510510.json').read_text())
    refinements = json.loads(Path('data/refined-plane-capacity-510510.json').read_text())
    coupled = json.loads(Path('data/coupled-core3-510510.json').read_text())
    assert len(source)==4 and len(planes)==32 and len(refinements)==16
    lookup = {(tuple(r['core_pattern']),tuple(r['high_pattern']),r['profile_index']):r for r in refinements}
    bad = {}
    for row in planes:
        key = tuple(row['core_pattern']),tuple(row['high_pattern'])
        bad[key] = {i for i,t in enumerate(row['tests']) if not t['excluded'] and not lookup[key+(i,)]['excluded']}
    assert {k:v for k,v in bad.items() if v}=={((0,1,0),(1,1,1)):{1},((0,1,1),(1,1,1)):{2}}
    by_coordinates = {(n%2,(n%3,n%5,n%7)):n for n in range(210)}
    reports = []
    for row in source:
        pattern = tuple(row['pattern'])
        points,masks,_ = pair_residuals(pattern)
        numbers = [by_coordinates[e,tuple(a)] for e,a in points]
        expected = {n for n in range(210) if all(n%p!=0 and not (n%2==b and n%p==1) for p,b in zip((3,5,7),pattern))}
        assert set(numbers)==expected and len(numbers)==len(expected)
        profiles_by_mask = {m:optimistic(statistics(m,numbers)) for m in masks}
        profiles = set(p for ps in profiles_by_mask.values() for p in ps)
        assert len(profiles)==row['profiles']
        assert all(any(easier(q,p) for q in row['easiest_profiles']) for p in profiles)
        report = dict(pattern=pattern,canonical_pair_residuals=len(masks),profiles=len(profiles),catalogue_complete=True)
        if pattern in ((0,1,0),(0,1,1)):
            data = json.loads(Path('data/core-triple-refinement-'+''.join(map(str,pattern))+'.json').read_text())
            safe = [q for i,q in enumerate(row['easiest_profiles']) if i not in bad[pattern,(1,1,1)]]
            critical = {p for p in profiles if not any(easier(q,p) for q in safe)}
            selected = {m for m,ps in profiles_by_mask.items() if any(p in critical for p in ps)}
            assert selected==set(data['critical_pair_masks'])
            # Actual congruences modulo 105, followed by an arbitrary nonredundant
            # class modulo 210; no coordinate triple-group helper is used.
            final = set()
            for m in selected:
                residues = {numbers[j]%105 for j in range(len(numbers)) if m>>j&1}
                for a in residues:
                    left = m & ~sum(1<<j for j,n in enumerate(numbers) if n%105==a)
                    assert left
                    bits = left
                    while bits:
                        bit = bits & -bits
                        final.add(left ^ bit)
                        bits ^= bit
            assert len(final)==data['actual_triple_residuals']
            actual = {m:statistics(m,numbers) for m in final}
            unresolved = {p for p in actual.values() if not any(easier(q,p) for q in safe)}
            assert len(set(actual.values()))==data['actual_profiles']
            assert len(unresolved)==data['unresolved_actual_profiles']
            assert all(any(easier(q,p) for q in data['remaining_profiles']) for p in unresolved)
            if pattern==(0,1,0):
                assert all(t['excluded'] for t in data['tests'])
            else:
                certified = {(tuple(t['profile']),tuple(t['layer1_mod3_counts'])):t for t in coupled['tests']}
                used = Counter()
                for m,p in actual.items():
                    if p not in unresolved:
                        continue
                    r = [n for j,n in enumerate(numbers) if m>>j&1]
                    assert all(n%3==2 for n in r if n%2==0)
                    counts = tuple(sum(n%2==1 and n%3==a for n in r) for a in (1,2))
                    t = certified[p,counts]
                    assert t['all_excluded'] and len(t['cases'])==8 and all(c['excluded'] for c in t['cases'])
                    used[p,counts] += 1
                assert set(used)==set(certified)
                assert all(used[k]==certified[k]['actual_residual_sets'] for k in used)
                report['coupled_residual_sets_checked'] = sum(used.values())
            report.update(critical_pair_residuals=len(selected),actual_triple_residuals=len(final),refinement_complete=True)
        reports.append(report)
        print(json.dumps(report),flush=True)
    result = dict(verified=True,cases=reports,seconds=round(time.monotonic()-start,3),
                  scope='core-catalogue and refinement completeness; plane arithmetic audited separately')
    Path('data/audit-core-completeness-510510.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)


if __name__=='__main__':
    main()
