"""Integer-only audit of event bounds from weighted residue capacities.

This verifier neither imports an optimizer nor reads a SAT proof. Every hit
vector rounds down to a tuple of event thresholds without changing its score.
Nonnegative capacity weights therefore suffice to exclude each excessive tuple.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def audit(certificate):
    points, moduli = certificate['points'], certificate['moduli']
    events, inequalities = certificate['events'], certificate['inequalities']
    assert len(points)==len(set(points)) and all(type(x) is int for x in points)
    assert len(moduli)==len(set(moduli)) and all(type(d) is int and d>0 for d in moduli)
    bound = certificate['bound']
    assert type(bound) is int and bound>=0
    levels = [{0} for x in points]
    for e in events:
        x,t,w = e['point'],e['threshold'],e['weight']
        assert type(x) is int and 0<=x<len(points)
        assert type(t) is int and 1<=t<=len(moduli)
        assert type(w) is int and w>0
        levels[x].add(t)
    for row in inequalities:
        weights = row['weights']
        assert len(weights)==len(points) and all(type(w) is int and w>=0 for w in weights)
        capacity = 0
        for d in moduli:
            # Enumerating all residues is intentionally independent of the
            # generator's grouping over only occupied residue classes.
            capacity += max(sum(w for x,w in zip(points,weights) if x%d==a) for a in range(d))
        assert type(row['capacity']) is int and row['capacity']==capacity
    checked = excessive = feasible_max = 0
    for thresholds in itertools.product(*(sorted(line) for line in levels)):
        checked += 1
        score = sum(e['weight'] for e in events if thresholds[e['point']]>=e['threshold'])
        excluded = any(sum(w*t for w,t in zip(row['weights'],thresholds))>row['capacity']
                       for row in inequalities)
        if score>bound:
            excessive += 1
            assert excluded, (thresholds,score)
        if not excluded:
            feasible_max = max(feasible_max,score)
    residues = certificate.get('attaining_residues')
    attained = None
    if residues is not None:
        assert len(residues)==len(moduli)
        assert all(type(a) is int and 0<=a<d for a,d in zip(residues,moduli))
        hits = [sum(x%d==a for a,d in zip(residues,moduli)) for x in points]
        attained = sum(e['weight'] for e in events if hits[e['point']]>=e['threshold'])
        assert attained==bound
    return dict(verified=True,bound=bound,inequalities_checked=len(inequalities),
                threshold_profiles_checked=checked,excessive_profiles_excluded=excessive,
                relaxation_maximum=feasible_max,attained_score=attained,
                scope='one fixed finite core event inequality; not a global covering obstruction')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    raw = Path(args.source).read_bytes()
    report = audit(json.loads(raw))
    report.update(source=args.source,source_sha256=hashlib.sha256(raw).hexdigest())
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
