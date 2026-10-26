"""Independent integer audit of local tail models retaining core identities."""
import argparse
import itertools
import json
import math
from collections import Counter
import gzip
import hashlib
from pathlib import Path
from audit_multitail import local_minima,audit as audit_aggregate


def audit_group_cut(source,cut,proof_cache=None):
    h = cut['axis']
    assert type(h) is int and 0<=h<(1<<len(source['high_primes']))-1
    assert type(cut['bound']) is int and cut['bound']>=0
    shapes = [d for d in range(1 if h>=len(source['high_primes']) else 3,211) if 210%d==0]
    if cut['proof_audit'].get('kind')=='weighted_capacity':
        # Reconstruct every resource maximum using only integer arithmetic;
        # neither the generator nor its saved verified flag is trusted.
        events = cut['events']
        weights = cut['capacity_certificate']['point_weights']
        assert len(events)>=2 and len(events)==len(weights)
        assert cut['bound']==len(events)-1
        assert all(type(w) is int and w>0 for w in weights)
        points = [e['point'] for e in events]
        assert len(points)==len(set(points))
        assert all(type(x) is int and 0<=x<len(source['core_residual']) for x in points)
        assert all(type(e['weight']) is int and e['weight']==1 for e in events)
        assert all(type(e['threshold']) is int and 1<=e['threshold']<=len(shapes) for e in events)
        maxima = []
        for d in shapes:
            buckets = [0]*d
            for x,w in zip(points,weights):
                buckets[source['core_residual'][x]%d] += w
            maxima.append(max(buckets))
        assert maxima==cut['capacity_certificate']['resource_capacities']
        assert sum(e['threshold']*w for e,w in zip(events,weights))>sum(maxima)
        return 0
    # Other whole-group cuts require a separately replayed propositional proof.
    report = cut['proof_audit']
    assert report['axis']==h and report['bound']==cut['bound']
    # A call-local cache reuses a successfully replayed proof only when its
    # entire mathematical input and proof identity agree. Group labels may
    # differ, but the actual core resource set must be identical.
    key = json.dumps([source['core_residual'],shapes,cut['events'],cut['bound'],
                      report['cnf_sha256'],report['proof_sha256'],report['proof_path']],
                     sort_keys=True,separators=(',',':'))
    if proof_cache is not None and key in proof_cache:
        return 0
    from audit_group_event_cut import build
    from sat_uncovered_audit import check_rup
    clauses,_ = build(source['core_residual'],shapes,cut['events'],cut['bound'])
    assert hashlib.sha256(json.dumps(clauses,separators=(',',':')).encode()).hexdigest()==report['cnf_sha256']
    with gzip.open(report['proof_path'],'rb') as f:
        raw = f.read()
    assert hashlib.sha256(raw).hexdigest()==report['proof_sha256']
    audit = check_rup(clauses,raw.decode().splitlines())
    assert audit['verified']
    if proof_cache is not None:
        proof_cache.add(key)
    return audit['checked_additions']


def audit(data):
    source,result = data['input'],data['result']
    proof_cache = set()
    R,high,pattern,H = [source[k] for k in ('core_residual','high_primes','high_pattern','H_upper')]
    n = len(high)
    assert len(R)==len(set(R)) and all(type(v) is int and 0<=v<210 for v in R)
    counts = [sum(v%2==e for v in R) for e in (0,1)]
    assert counts==source['counts']
    assert len(H)==counts[0]+1 and all(len(line)==counts[1]+1 for line in H)
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    core_moduli = [d for d in range(3,211) if 210%d==0]
    shapes = [core_moduli if h<n else [1,2]+core_moduli for h in range(len(groups))]
    witness = result.get('point_primal_witness')
    if witness:
        scale = witness['denominator']
        assert type(scale) is int and scale>0
        rhs = [[0]*len(groups) for x in R]
        expected = {(h,d) for h,moduli in enumerate(shapes) for d in moduli}
        seen = set()
        for row in witness['projections']:
            h,d = row['axis'],row['modulus']
            assert (h,d) in expected and (h,d) not in seen
            seen.add((h,d))
            assert sum(a['weight'] for a in row['assignments'])==scale
            residues = set()
            for assignment in row['assignments']:
                a,w = assignment['residue'],assignment['weight']
                assert type(a) is int and 0<=a<d and a not in residues
                assert type(w) is int and w>0
                residues.add(a)
                for x,v in enumerate(R):
                    if v%d==a:
                        rhs[x][h] += w
        assert seen==expected
        masses = [0]*len(R)
        lhs = [[0]*len(groups) for x in R]
        hist = [[[0]*(len(shape)+1) for shape in shapes] for e in (0,1)]
        point_hist = [[[0]*(len(shape)+1) for shape in shapes] for x in R]
        objective = 0
        for row in witness['support']:
            x,state,w = row['point_index'],row['counts'],row['weight']
            assert type(x) is int and 0<=x<len(R) and type(w) is int and w>0
            assert len(state)==len(groups)
            e = R[x]%2
            D = [q-1-(b==e) for q,b in zip(high,pattern)]
            masses[x] += w
            for h,v in enumerate(state):
                assert type(v) is int and 0<=v<=(min(14,D[h]) if h<n else 16)
                lhs[x][h] += w*v
                hist[e][h][v] += w
                point_hist[x][h][v] += w
            free = [d-v for d,v in zip(D,state[:n])]
            objective += w*(math.prod(free)-sum(state[h]*math.prod(free[i] for i in range(n) if i not in S)
                                               for h,S in enumerate(groups[n:],n)))
        assert masses==[scale]*len(R)
        assert all(a<=b for la,lb in zip(lhs,rhs) for a,b in zip(la,lb))
        checked = 0
        for h,shape in enumerate(shapes):
            for t0,t1 in itertools.product(range(len(shape)+1),repeat=2):
                used = sum(w*max(0,v-t) for e,t in enumerate((t0,t1)) for v,w in enumerate(hist[e][h]))
                bound = max(v+(a+b+max(a,b) if h>=n else 0)-t0*a-t1*b for a,line in enumerate(H) for b,v in enumerate(line))
                assert used<=scale*bound
                checked += 1
        assert objective==witness['objective_numerator'] and objective<0
        rank_checked = 0
        if result.get('rank_cuts_enabled'):
            aggregated = dict(denominator=scale,objective_numerator=objective,
                              support=[dict(parity=R[r['point_index']]%2,counts=r['counts'],weight=r['weight']) for r in witness['support']])
            extra = audit_aggregate(dict(input=source,result=dict(primal_relaxation_witness=aggregated,rank_cuts_enabled=True)))
            rank_checked = extra['rank_inequalities_checked']
        pair_checked = 0
        if result.get('pair_cuts_enabled'):
            for h,shape in enumerate(shapes):
                limit = len(shape)
                for x,y in itertools.combinations(range(len(R)),2):
                    cap = limit+sum(R[x]%d==R[y]%d for d in shape)
                    for t in range(1,limit+1):
                        u = cap+1-t
                        if 1<=u<=limit:
                            assert sum(point_hist[x][h][t:])+sum(point_hist[y][h][u:])<=scale
                            pair_checked += 1
        triple_checked = 0
        if result.get('triple_cuts_enabled'):
            tails = [[[sum(line[t:]) for t in range(len(line))] for line in point] for point in point_hist]
            for h,shape in enumerate(shapes):
                limit = len(shape)
                for selected in itertools.combinations(range(len(R)),3):
                    cap = sum(max(Counter(R[x]%d for x in selected).values()) for d in shape)
                    x,y,z = selected
                    for t,u in itertools.product(range(1,limit+1),repeat=2):
                        v = max(1,cap+1-t-u)
                        if v<=limit:
                            assert tails[x][h][t]+tails[y][h][u]+tails[z][h][v]<=2*scale
                            triple_checked += 1
        subset_checked = 0
        for h,selected,thresholds in result.get('subset_constraints',[]):
            assert type(h) is int and 0<=h<len(groups)
            assert len(selected)>=3 and len(selected)==len(set(selected))==len(thresholds)
            assert all(type(x) is int and 0<=x<len(R) for x in selected)
            assert all(type(t) is int and 1<=t<=len(shapes[h]) for t in thresholds)
            cap = sum(max(Counter(R[x]%d for x in selected).values()) for d in shapes[h])
            assert sum(thresholds)>cap
            assert sum(sum(point_hist[x][h][t:]) for x,t in zip(selected,thresholds))<=(len(selected)-1)*scale
            subset_checked += 1
        group_checked = group_proof_additions = 0
        for cut in result.get('group_event_constraints',[]):
            used = sum(event['weight']*sum(point_hist[event['point']][cut['axis']][event['threshold']:]) for event in cut['events'])
            assert used<=scale*cut['bound']
            group_proof_additions += audit_group_cut(source,cut,proof_cache)
            group_checked += 1
        return dict(verified=True,counts=counts,high_primes=high,high_pattern=pattern,
                    objective_numerator=objective,denominator=scale,projection_resources_checked=len(expected),
                    point_link_inequalities_checked=len(R)*len(groups),concentration_inequalities_checked=checked,
                    rank_inequalities_checked=rank_checked,
                    pair_inequalities_checked=pair_checked,
                    triple_inequalities_checked=triple_checked,
                    stored_subset_inequalities_checked=subset_checked,
                    group_event_inequalities_checked=group_checked,group_event_rup_additions_checked=group_proof_additions,
                    support_size=len(witness['support']),
                    scope='feasible negative-cost point-identity relaxation distribution; not an actual covering system')
    certificate = result['certificate']
    links = certificate['link_multipliers']
    assert len(links)==len(R) and all(len(line)==len(groups) and all(type(v) is int and v<=0 for v in line) for line in links)
    points,checked = [],0
    pair_penalties = [[] for x in R]
    for row in certificate.get('pair_multipliers',[]):
        h = row['axis']
        x,y = row['points']
        t,u = row['thresholds']
        coefficient = row['multiplier']
        assert type(h) is int and 0<=h<len(groups)
        assert type(x) is int and type(y) is int and 0<=x<y<len(R)
        assert type(t) is int and type(u) is int and 1<=t<=len(shapes[h]) and 1<=u<=len(shapes[h])
        assert type(coefficient) is int and coefficient<=0
        cap = len(shapes[h])+sum(R[x]%d==R[y]%d for d in shapes[h])
        assert t+u>cap
        pair_penalties[x].append(dict(axis=h,thresholds=[t,t],weights=[1,1],multiplier=coefficient))
        pair_penalties[y].append(dict(axis=h,thresholds=[u,u],weights=[1,1],multiplier=coefficient))
    for row in certificate.get('subset_multipliers',[]):
        h,selected,thresholds,coefficient = [row[k] for k in ('axis','points','thresholds','multiplier')]
        assert type(h) is int and 0<=h<len(groups)
        assert len(selected)>=3 and len(selected)==len(set(selected))==len(thresholds)
        assert all(type(x) is int and 0<=x<len(R) for x in selected)
        assert all(type(t) is int and 1<=t<=len(shapes[h]) for t in thresholds)
        assert type(coefficient) is int and coefficient<=0
        cap = sum(max(Counter(R[x]%d for x in selected).values()) for d in shapes[h])
        assert sum(thresholds)>cap
        for x,t in zip(selected,thresholds):
            pair_penalties[x].append(dict(axis=h,thresholds=[t,t],weights=[1,1],multiplier=coefficient))
    for cut in certificate.get('group_event_multipliers',[]):
        audit_group_cut(source,cut,proof_cache)
        coefficient = cut['multiplier']
        assert type(coefficient) is int and coefficient<=0
        for event in cut['events']:
            x,t,w = [event[k] for k in ('point','threshold','weight')]
            pair_penalties[x].append(dict(axis=cut['axis'],thresholds=[t,t],weights=[w,w],multiplier=coefficient))
    for x,v in enumerate(R):
        augmented = dict(denominator=certificate['denominator'],rank_multipliers=certificate.get('rank_multipliers',[])+pair_penalties[x],tail_multipliers=certificate['tail_multipliers']+
                         [dict(axis=h,thresholds=[0,0],multiplier=y) for h,y in enumerate(links[x]) if y])
        values,work = local_minima(source,augmented,parity_only=v%2)
        assert values[0]==certificate['point_multipliers'][x]
        points.append(values[0])
        checked += work
    terms = [min(sum(links[x][h] for x,v in enumerate(R) if v%d==a) for a in {v%d for v in R})
             for h,moduli in enumerate(shapes) for d in moduli]
    assert terms==certificate['resource_multipliers']
    numerator = sum(points)+sum(terms)
    numerator += sum(row['multiplier'] for row in certificate.get('pair_multipliers',[]))
    numerator += sum(row['multiplier']*(len(row['points'])-1) for row in certificate.get('subset_multipliers',[]))
    numerator += sum(row['multiplier']*row['bound'] for row in certificate.get('group_event_multipliers',[]))
    for row in certificate['tail_multipliers']:
        t0,t1 = row['thresholds']
        bound = max(v+(a+b+max(a,b) if row['axis']>=n else 0)-t0*a-t1*b for a,line in enumerate(H) for b,v in enumerate(line))
        assert bound==row['bound']
        numerator += row['multiplier']*bound
    for row in certificate.get('rank_multipliers',[]):
        t0,t1 = row['thresholds']
        w0,w1 = row['weights']
        feasible = [(a,b) for a,line in enumerate(H) for b,v in enumerate(line)
                    if t0*a+t1*b<=v+(a+b+max(a,b) if row['axis']>=n else 0)]
        bound = max(w0*a+w1*b for a,b in feasible)
        assert bound==row['bound']
        numerator += row['multiplier']*bound
    assert numerator==certificate['numerator']
    lower = -(-numerator//certificate['denominator'])
    assert lower==certificate['integer_lower']
    return dict(verified=True,counts=counts,high_primes=high,high_pattern=pattern,integer_lower=lower,
                positive_obstruction=lower>0,plane_tuples_checked=checked,
                scope='one specified core and pattern, conditional on separately audited capacity H')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    result = audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
