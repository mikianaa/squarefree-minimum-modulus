"""Standard-library, ordinary-residue audit of the full 000 eight-prime branch.

Enumerates all nonempty odd/even pair classes without label symmetry, checks
the core count frontier, rebuilds compression capacities from explicit points
and all 128 even-shape allocations, and verifies every reused dual inequality.
"""
import itertools
import json
import math
from pathlib import Path
import time
from audit_conditional_tail import audit as audit_local_certificate


def ordinary_frontier(pattern=(0,0,0)):
    R = [n for n in range(210) if all(n%p!=0 and not(n%2==b and n%p==1) for p,b in zip((3,5,7),pattern))]
    even = sum(1<<j for j,n in enumerate(R) if n%2==0)
    whole = (1<<len(R))-1
    def slices(d):
        return sorted({sum(1<<j for j,n in enumerate(R) if n%d==a) for a in range(d)}-{0})
    odd = [slices(d) for d in (15,21,35)]
    doubled = [slices(d) for d in (30,42,70)]
    bases = {a|b|c for a,b,c in itertools.product(*odd)}
    pairs = set()
    families = 0
    for base in bases:
        first_two = {base|a|b for a,b in itertools.product(*doubled[:2])}
        for union in first_two:
            for last in doubled[2]:
                remaining = whole&~(union|last)
                r0 = (remaining&even).bit_count()
                r1 = remaining.bit_count()-r0
                pairs.add((max(0,r0-2),max(0,r1-1)))
                pairs.add((max(0,r0-1),max(0,r1-2)))
                families += 1
    frontier = sorted([a,b] for a,b in pairs if not any(c<=a and d<=b and (c,d)!=(a,b) for c,d in pairs))
    return frontier,families


def explicit_capacity(pattern=(0,0,0)):
    points = [[n for n in range(210) if n%2==e and all(n%p!=0 and not(e==b and n%p==1) for p,b in zip((3,5,7),pattern))] for e in (0,1)]
    shapes = [3,5,7,15,21,35,105]
    indicators = [[[all(n%p==2 for p in (3,5,7) if d%p==0) for d in shapes] for n in layer] for layer in points]
    H = [[0]*(len(points[1])+1) for _ in range(len(points[0])+1)]
    for assigned in itertools.product((0,1),repeat=7):
        totals = []
        for e in (0,1):
            weights = sorted((sum(int(hit)*(1+int(assigned[j]==e)) for j,hit in enumerate(row)) for row in indicators[e]),reverse=True)
            sums = [0]
            for w in weights:
                sums.append(sums[-1]+w)
            totals.append(sums)
        for a,v in enumerate(totals[0]):
            for b,w in enumerate(totals[1]):
                H[a][b] = max(H[a][b],v+w)
    return H


def main():
    start = time.monotonic()
    read = lambda name: json.loads(Path('data',name).read_text())
    target = read('eight-prime-anchor000-certificate.json')
    assert target['period']==9699690 and target['core_pattern']==[0,0,0]
    frontier,families = ordinary_frontier()
    assert frontier==target['frontier']
    print(json.dumps(dict(stage='ordinary_core_frontier',frontier=frontier,pair_unions_with_last_line_checked=families,seconds=round(time.monotonic()-start,3))),flush=True)
    H = explicit_capacity()
    certificates = read('conditional-tail-core-pattern-certificates.json')['cases']
    concentration = read('core-hit-concentration-7-15.json')
    incumbent = read('core-fibre-eight-plane.json')
    assert len(certificates)==16 and [c['high_pattern'] for c in certificates]==target['high_patterns']
    assert len({tuple(c['high_pattern']) for c in certificates})==16
    # These checks reprice every local state using Python integers. The old
    # fixed-core H is immaterial to the new proof: its RHS is replaced below.
    for c in certificates:
        assert c['aggregate'] and not c['certificate']['resource_multipliers']
        audit_local_certificate(c,concentration,incumbent,check_subsets=False)
    cache = {}
    results = []
    for record in target['cases']:
        a,b = record['core_counts']
        assert [a,b] in frontier
        bounds = []
        for data in certificates:
            c = data['certificate']
            numerator = a*c['point_multipliers'][0]+b*c['point_multipliers'][1]
            for item in c['tail_multipliers']:
                tail = item['axis']>=4
                t0,t1 = item['thresholds']
                key = a,b,tail,t0,t1
                if key not in cache:
                    cache[key] = max(H[m0][m1]+(m0+m1+max(m0,m1) if tail else 0)-t0*m0-t1*m1 for m0 in range(a+1) for m1 in range(b+1))
                assert item['multiplier']<=0
                numerator += item['multiplier']*cache[key]
            bounds.append(-(-numerator//c['denominator']))
        assert bounds==record['bounds'] and min(bounds)>0
        results.append(dict(core_counts=[a,b],bounds=bounds,verified=True))
    assert len(results)==len(frontier)==len({tuple(r['core_counts']) for r in results})
    minimum = min(min(r['bounds']) for r in results)
    assert minimum==target['uniform_lower']
    result = dict(verified=True,period=9699690,core_pattern=[0,0,0],
                  core_frontier=frontier,ordinary_pair_checks=families,
                  high_patterns=target['high_patterns'],cases=results,
                  uniform_uncovered_lower=minimum,seconds=round(time.monotonic()-start,3),
                  scope='all cores in normalized low-anchor branch 000, all 16 higher patterns; other low-anchor branches remain open')
    Path('data/audit-eight-prime-anchor000.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','high_patterns','core_frontier')}),flush=True)


if __name__=='__main__':
    main()
