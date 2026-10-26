"""Independent integer audit of the quantitative core-capacity frontiers.

Uses ordinary residues for statistics, explicit subset sums for even resources,
and the independent histogram audit for the second prime. The first core's
symmetry reduction was independently checked by audit_core_geometry.py.
"""
from collections import Counter
import itertools
import json
from pathlib import Path
import time
from core_residual_geometry import pair_residuals
from audit_core_completeness import statistics,optimistic
from audit_compressed_parity import boundary,audit_extension


def profile_rows(pattern):
    points,masks,_ = pair_residuals(tuple(pattern))
    by_crt = {(n%2,(n%3,n%5,n%7)):n for n in range(210)}
    numbers = [by_crt[e,tuple(a)] for e,a in points]
    rows = set()
    for m in masks:
        p = statistics(m,numbers)
        counts = [sum(m>>j&1 for j,n in enumerate(numbers) if n%2==1 and n%3==v) for v in (1,2)]
        for op in optimistic(p):
            rows.add(op+tuple(min(op[1],v) for v in counts))
    return sorted(rows)


def direct_frontier(rows,q,b):
    best = {}
    for row in rows:
        r = row[:2]
        caps = list(zip(row[2:16:2],row[3:16:2]))
        even = {(0,0)}
        for c0,c1 in caps:
            even = {(x+c0,y) for x,y in even}|{(x,y+c1) for x,y in even}
        other_odd = [sum(c[e] for c in caps[1:]) for e in (0,1)]
        for odd3 in ((0,row[16]),(r[0],row[17])):
            needed = [(q-1-(b==e))*r[e]-other_odd[e]-odd3[e] for e in (0,1)]
            for x,y in even:
                a,c = max(0,needed[0]-x),max(0,needed[1]-y)
                best[a] = min(best.get(a,c),c)
    return boundary(best.items())


def main():
    start = time.monotonic()
    data = json.loads(Path('data/core-capacity-frontiers.json').read_text())
    previous = json.loads(Path('data/compressed-parity-frontiers.json').read_text())
    optima = {tuple(r['pattern']):r['uncovered'] for r in
              json.loads(Path('data/min-uncovered-2310-five-patterns-optimal.json').read_text())}
    reports = []
    for tail in itertools.product((0,1),repeat=2):
        pattern = [0]+list(tail)
        rows = profile_rows(pattern)
        for b in (0,1):
            p = pattern+[b]
            target = next(r for r in data['at2310'] if r['pattern']==p)
            assert len(rows)==target['enhanced_profile_count']
            raw = direct_frontier(rows,11,b)
            assert raw==target['raw_frontier']
            old = next(r for r in previous['at2310'] if r['pattern']==p)['frontier']
            front = boundary((max(a,c),max(b,d)) for a,b in raw for c,d in old)
            assert front==target['frontier'] and min(map(sum,front))==target['uncovered_lower']
            report = dict(period=2310,pattern=p,verified=True,lower=target['uncovered_lower'])
            reports.append(report)
            print(json.dumps(report),flush=True)
    for target in data['at30030']:
        p = target['pattern']
        base = next(r for r in data['at2310'] if r['pattern']==p[:4])['frontier']
        if target['status']=='including_solver':
            lower = optima.get(tuple(p[:4]),0)
            base = boundary((a+t,max(b,lower-a-t)) for a,b in base for t in range(max(0,lower-a-b)+1))
        actual = audit_extension((3,5,7,11),p[:4],base,13,p[4])
        assert actual==target['frontier'] and min(map(sum,actual))==target['uncovered_lower']
        report = dict(period=30030,pattern=p,status=target['status'],verified=True,lower=target['uncovered_lower'])
        reports.append(report)
        print(json.dumps(report),flush=True)
    result = dict(verified=True,cases=reports,seconds=round(time.monotonic()-start,3),
                  scope='quantitative capacity arithmetic and propagation; solver-dependent inputs retain their dependence')
    Path('data/audit-core-capacity-frontiers.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
