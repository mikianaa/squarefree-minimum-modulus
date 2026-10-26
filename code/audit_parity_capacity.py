"""Independent audit of the 2310 frontier and 30030 obstruction certificates.

Uses brute integer residual ranges instead of convex minimization, and exhaustive
even-resource assignments instead of the knapsack feasibility routine. The 210
input is separately checked by the full unquotiented axis-line enumeration in
test_research_models.py. Standard library only.
"""
import itertools
import json
import math
from pathlib import Path


def pareto(points):
    # Independent quadratic dominance check; sets here are small after grouping.
    by_first = {}
    for a,b in points:
        by_first[a] = min(b,by_first.get(a,b))
    return sorted((a,b) for a,b in by_first.items()
                  if not any(x < a and y <= b for x,y in by_first.items()))


def cap_table(ps,bs):
    modulus = math.prod(ps)
    ds = [d for d in range(2,modulus+1) if modulus % d == 0]
    dims = [[p-1-int(b == e) for p,b in zip(ps,bs)] for e in (0,1)]
    sizes = [math.prod(row) for row in dims]
    caps = [[sizes[e]//math.prod(dims[e][i] for i,p in enumerate(ps) if d % p == 0)
             for d in ds] for e in (0,1)]
    return sizes,caps


def brute_extension(ps,bs,frontier,q,bq):
    sizes,caps = cap_table(ps,bs)
    points = []
    for allocation in itertools.product((0,1),repeat=len(caps[0])):
        functions = []
        for e in (0,1):
            coeff = q-1-int(bq == e)
            resources = caps[e]+[c for c,a in zip(caps[e],allocation) if a == e]
            values = [max(0,coeff*r-sum(min(r,c) for c in resources)) for r in range(sizes[e]+1)]
            # Suffix minima give the best value for every r >= the lower bound.
            for r in range(len(values)-2,-1,-1):
                values[r] = min(values[r],values[r+1])
            functions.append(values)
        points.extend((functions[0][r],functions[1][s]) for r,s in frontier)
    return pareto(points)


def brute_cover_margin(ps,bs,r,q,bq):
    _,caps = cap_table(ps,bs)
    c0 = [min(r[0],c) for c in caps[0]]
    c1 = [min(r[1],c) for c in caps[1]]
    # All odd shapes contribute in both layers. Initially allocate every even
    # shape to layer 1, then enumerate every subset transferred to layer 0.
    margin0 = sum(c0)-(q-1-int(bq == 0))*r[0]
    margin1 = 2*sum(c1)-(q-1-int(bq == 1))*r[1]
    totals0,totals1 = [0],[0]
    for a,b in zip(c0,c1):
        totals0 += [x+a for x in totals0]
        totals1 += [x+b for x in totals1]
    return max(min(margin0+a,margin1-b) for a,b in zip(totals0,totals1))


def audit_simplified():
    root = Path(__file__).resolve().parents[1]
    certificate = json.loads((root/'data/six-prime-obstruction.json').read_text())
    results = []
    for row in certificate:
        pattern = row['pattern']
        sizes,caps = cap_table((3,5,7),pattern[:3])
        lower = math.inf
        for assignment in itertools.product((0,1),repeat=len(caps[0])):
            deficits = []
            for e in (0,1):
                resources = caps[e]+[c for c,a in zip(caps[e],assignment) if a == e]
                deficits.append([max(0,(10-int(pattern[3] == e))*r-sum(min(r,c) for c in resources))
                                 for r in range(sizes[e]+1)])
            suffix = deficits[1][:]
            for i in range(len(suffix)-2,-1,-1):
                suffix[i] = min(suffix[i],suffix[i+1])
            for r0,g0 in enumerate(deficits[0]):
                r1 = max(0,row['uncovered_210_lower']-r0)
                if r1 <= sizes[1]:
                    lower = min(lower,g0+suffix[r1])
        assert lower >= row['uncovered_2310_lower']
        # Stronger independent check: inspect ALL feasible count pairs in the
        # anchor boxes, not just pairs on the lower-sum boundary.
        sizes,caps = cap_table((3,5,7,11),pattern)
        truncated = [[[min(r,c) for c in caps[e]] for r in range(sizes[e]+1)] for e in (0,1)]
        margins = [math.inf,math.inf]
        for r0,c0 in enumerate(truncated[0]):
            for r1 in range(max(0,row['uncovered_2310_lower']-r0),sizes[1]+1):
                c1 = truncated[1][r1]
                deficit = 12*(r0+r1)-sum(c0)-sum(c1)-sum(max(a,b) for a,b in zip(c0,c1))
                margins[0] = min(margins[0],deficit-r0)
                margins[1] = min(margins[1],deficit-r1)
        assert min(margins) > 0
        result = dict(pattern=pattern,brute_uncovered_2310_lower=lower,
                      all_count_pairs_deficits=margins,verified=True)
        results.append(result)
        print(json.dumps(result),flush=True)
    (root/'data/audit-six-prime-obstruction.json').write_text(json.dumps(results,indent=2)+'\n')


def main():
    root = Path(__file__).resolve().parents[1]
    base = json.loads((root/'data/three-axis-210-certificate.json').read_text())
    target = json.loads((root/'data/parity-capacity-30030-certificate.json').read_text())
    results = []
    for row in target:
        pattern = row['pattern']
        source = next(r for r in base['cases'] if r['pattern'] == pattern[:3])
        actual = brute_extension((3,5,7),pattern[:3],source['residual_pair_frontier'],11,pattern[3])
        assert actual == [tuple(v) for v in row['frontier']]
        margins = [max(brute_cover_margin((3,5,7,11),pattern,r,13,b) for r in actual) for b in (0,1)]
        assert all(v < 0 for v in margins), (pattern,margins)
        result = dict(pattern=pattern,frontier_size=len(actual),minimum_uncovered=min(sum(r) for r in actual),
                      maximal_cover_margins=margins,excluded_both_extensions=True)
        results.append(result)
        print(json.dumps(result),flush=True)
    (root/'data/audit-parity-capacity-30030.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    main()
    audit_simplified()
