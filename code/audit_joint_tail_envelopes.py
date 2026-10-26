"""Independent ordinary-residue audit of count-indexed joint H envelopes.

Enumerates all even pair classes, using only the explicitly checked two-label
normalization of the odd pair classes. Histogram counts and explicit-point
compression replace the producer's coordinate masks and weight tables.
NumPy accelerates exact integer arrays; no optimization solver is used.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from audit_eight_prime_anchor_branch import explicit_capacity


def normalized_pairs(p):
    representatives = [(1,1),(1,2),(2,1),(2,2)]+([(2,3)] if p>3 else [])
    for pair in itertools.product(range(1,p),repeat=2):
        found = False
        for perm in itertools.permutations(range(2,p)):
            image = {1:1,**dict(zip(range(2,p),perm))}
            if tuple(image[x] for x in pair) in representatives:
                found = True
                break
        assert found
    return representatives


def ordinary_pair_residuals(pattern):
    points = [n for n in range(210) if all(n%p!=0 and not(n%2==b and n%p==1) for p,b in zip((3,5,7),pattern))]
    whole = (1<<len(points))-1
    masks = {d:[sum(1<<j for j,n in enumerate(points) if n%d==a) for a in range(d)] for d in (15,21,35,30,42,70)}
    def crt(p,a,q,b):
        return next(n for n in range(p*q) if n%p==a and n%q==b)
    bases = set()
    for a,b,c in itertools.product(*(normalized_pairs(p) for p in (3,5,7))):
        residues = [crt(3,a[0],5,b[0]),crt(3,a[1],7,c[0]),crt(5,b[1],7,c[1])]
        bases.add(masks[15][residues[0]]|masks[21][residues[1]]|masks[35][residues[2]])
    even = [set(masks[d])-{0} for d in (30,42,70)]
    residuals = set()
    for base in bases:
        first = {base|a|b for a,b in itertools.product(*even[:2])}
        residuals.update(whole&~(u|v) for u,v in itertools.product(first,even[2]))
    return points,residuals


def histogram_profile(numbers):
    layers = [np.array([n for n in numbers if n%2==e],dtype=np.int64) for e in (0,1)]
    profiles = []
    for d in (3,5,7,15,21,35,105):
        histograms = [np.bincount(layer%d,minlength=d) for layer in layers]
        pairs = sorted(set(zip(*[h.tolist() for h in histograms])),reverse=True)
        frontier = []
        maximum_second = -1
        for a,b in pairs:
            if b>maximum_second:
                frontier.append((a,b))
                maximum_second = b
        profiles.append(tuple(frontier))
    return (len(layers[0]),len(layers[1]),tuple(profiles))


def from_profile(profile,compression):
    r0,r1,shapes = profile
    a = np.arange(r0+1)[:,None]
    b = np.arange(r1+1)[None,:]
    H = np.zeros((r0+1,r1+1),dtype=np.int64)
    for pairs in shapes:
        odd = np.zeros_like(H)
        for x,y in pairs:
            odd = np.maximum(odd,np.minimum(a,x)+np.minimum(b,y))
        even = np.maximum(np.minimum(a,max(x for x,y in pairs)),np.minimum(b,max(y for x,y in pairs)))
        H += odd+even
    return np.minimum(H,compression[:r0+1,:r1+1])


def audit(data,profile_data=None):
    start = time.monotonic()
    pattern = data['core_pattern']
    assert pattern in ([0,0,0],[0,0,1],[0,1,0],[0,1,1])
    points,masks = ordinary_pair_residuals(pattern)
    print(json.dumps(dict(stage='ordinary_pair_enumeration',pattern=pattern,masks=len(masks),seconds=round(time.monotonic()-start,3))),flush=True)
    compression = np.array(explicit_capacity(pattern),dtype=np.int64)
    cache = {}
    envelopes = {}
    for index,mask in enumerate(masks):
        numbers = [n for j,n in enumerate(points) if mask>>j&1]
        profile = histogram_profile(numbers)
        if profile not in cache:
            cache[profile] = from_profile(profile,compression)
        H = cache[profile]
        r0,r1 = profile[:2]
        for extra in (0,1):
            a,b = max(0,r0-1-(extra==0)),max(0,r1-1-(extra==1))
            key = a,b
            if key not in envelopes:
                envelopes[key] = H[:a+1,:b+1].copy()
            else:
                np.maximum(envelopes[key],H[:a+1,:b+1],out=envelopes[key])
        if index and index%25000==0:
            print(json.dumps(dict(stage='histograms',processed=index,total=len(masks),profiles=len(cache),seconds=round(time.monotonic()-start,3))),flush=True)
    expected = {tuple(r['counts']):np.array(r['H_upper'],dtype=np.int64) for r in data['cases']}
    assert set(expected)==set(envelopes)
    for key,H in envelopes.items():
        assert np.array_equal(H,expected[key]),key
    profile_report = {}
    if profile_data is not None:
        assert profile_data['core_pattern']==pattern
        source_keys = {(r['counts_before'][0],r['counts_before'][1],tuple(tuple(tuple(p) for p in pairs) for pairs in r['slice_pairs'])) for r in profile_data['sources']}
        assert source_keys==set(cache)
        tables = profile_data['cases']
        used = set()
        for row in profile_data['sources']:
            key = (row['counts_before'][0],row['counts_before'][1],tuple(tuple(tuple(p) for p in pairs) for pairs in row['slice_pairs']))
            assert histogram_profile(row['representative_residual'])==key
            H = cache[key]
            for extra,index in enumerate(row['target_cases']):
                a,b = max(0,key[0]-1-(extra==0)),max(0,key[1]-1-(extra==1))
                assert tables[index]['counts']==[a,b]
                assert np.array_equal(np.array(tables[index]['H_upper']),H[:a+1,:b+1])
                used.add(index)
        assert used==set(range(len(tables)))
        basis = set(profile_data['maximal_case_indices'])
        assert len(profile_data['dominator'])==len(tables)
        for i,j in enumerate(profile_data['dominator']):
            assert j in basis and tables[i]['counts']==tables[j]['counts']
            assert np.all(np.array(tables[i]['H_upper'])<=np.array(tables[j]['H_upper']))
        profile_report = dict(profile_catalogue_verified=True,final_count_tables=len(tables),maximal_tables=len(basis))
    return dict(verified=True,core_pattern=pattern,ordinary_pair_residuals=len(masks),
                distinct_histogram_profiles=len(cache),count_envelopes=len(envelopes),
                all_entries_match=True,seconds=round(time.monotonic()-start,3),
                **profile_report,scope='complete ordinary-residue audit of joint concentration upper bounds; no covering exclusion by itself')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='010')
    ap.add_argument('--profiles',action='store_true')
    a = ap.parse_args()
    data = json.loads(Path(f'data/joint-tail-envelopes-{a.pattern}.json').read_text())
    profiles = json.loads(Path(f'data/joint-tail-profiles-{a.pattern}.json').read_text()) if a.profiles else None
    result = audit(data,profiles)
    kind = 'profiles' if a.profiles else 'envelopes'
    Path(f'data/audit-joint-tail-{kind}-{a.pattern}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
