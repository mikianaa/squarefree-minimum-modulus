"""Independent arithmetic/completeness audit of the whole 010 branch.

Ordinary residues and explicit CRT permutations certify coverage of the
refined profiles. Weighted-point pair tables certify capacity upper bounds.
Exact local minimization checks only the duals actually used in the proof.
"""
import itertools
import argparse
import json
from pathlib import Path
import time
import numpy as np
from audit_joint_tail_envelopes import audit as audit_catalogue
from audit_joint_tail_envelopes import ordinary_pair_residuals,histogram_profile,from_profile
from audit_eight_prime_anchor_branch import explicit_capacity
from audit_paired_slice_capacity import audit as audit_pair_capacity
from audit_compressed_tail_branch import local_minima


def profile_key(source):
    return (*source['counts_before'],tuple(tuple(tuple(p) for p in shape) for shape in source['slice_pairs']))


def direct_concentration(H):
    a,b = np.indices(H.shape,dtype=np.int64)
    result = np.empty((2,17,17),dtype=np.int64)
    for kind in (0,1):
        table = H+(a+b+np.maximum(a,b) if kind else 0)
        for t0,t1 in itertools.product(range(17),repeat=2):
            result[kind,t0,t1] = np.max(table-t0*a-t1*b)
    return result


def audit(pattern=(0,1,0)):
    start = time.monotonic()
    pattern = list(pattern)
    assert pattern in ([0,1,0],[0,1,1])
    tag = ''.join(map(str,pattern))
    read = lambda name: json.loads(Path('data',name).read_text())
    cat = read(f'joint-tail-profiles-{tag}.json')
    report = audit_catalogue(read(f'joint-tail-envelopes-{tag}.json'),cat)
    data = read(f'paired-tail-refinement-{tag}.json')
    assert data['core_pattern']==pattern
    pool = data['pool']
    compression = np.array(explicit_capacity(pattern),dtype=np.int64)
    keys = {profile_key(source):i for i,source in enumerate(cat['sources'])}
    selected = set(data['source_profiles'])
    assert len(selected)==len(data['source_profiles'])
    crt = {(n%2,n%3,n%5,n%7):n for n in range(210)}
    images = {i:set() for i in selected}
    cap_tables = []
    for index,cap in enumerate(data['capacities']):
        assert cap['moduli']==[d for d in range(3,211) if 210%d==0]
        audit_pair_capacity(cap)
        key = histogram_profile(cap['core_residual'])
        source_index = keys[key]
        assert source_index in selected and source_index==cap['profile_index']
        joint = from_profile(key,compression)
        assert np.array_equal(joint,np.array(cap['joint_upper']))
        cap_tables.append(np.minimum(joint,np.array(cap['H_upper'],dtype=np.int64)))
        for a,b in itertools.product(itertools.permutations(range(2,5)),itertools.permutations(range(2,7))):
            p5,p7 = (0,1)+a,(0,1)+b
            image = tuple(sorted(crt[n%2,n%3,p5[n%5],p7[n%7]] for n in cap['core_residual']))
            images[source_index].add(image)
        if (index+1)%20==0:
            print(json.dumps(dict(stage='paired_capacities',done=index+1,total=len(data['capacities']),seconds=round(time.monotonic()-start,2))),flush=True)
    points,masks = ordinary_pair_residuals(pattern)
    matching = 0
    selected_counts = {tuple(cat['sources'][i]['counts_before']) for i in selected}
    even = sum(1<<j for j,n in enumerate(points) if n%2==0)
    for mask in masks:
        a = (mask&even).bit_count()
        if (a,mask.bit_count()-a) not in selected_counts:
            continue
        R = tuple(n for j,n in enumerate(points) if mask>>j&1)
        source_index = keys[histogram_profile(R)]
        if source_index in selected:
            assert R in images[source_index],(source_index,R)
            matching += 1
    print(json.dumps(dict(stage='refined_completeness',ordinary_residuals=len(masks),matching=matching,seconds=round(time.monotonic()-start,2))),flush=True)

    # Check every original count table, allowing its independently checked
    # pointwise dominator to supply a stronger computationally reusable proof.
    tables = [np.array(row['H_upper'],dtype=np.int64) for row in cat['cases']]
    counts = [row['counts'] for row in cat['cases']]
    pair_indices = {}
    for k,H in enumerate(cap_tables):
        r0,r1 = [v-1 for v in H.shape]
        for extra in (0,1):
            a,b = max(0,r0-1-(extra==0)),max(0,r1-1-(extra==1))
            pair_indices[k,extra] = len(tables)
            tables.append(H[:a+1,:b+1])
            counts.append([a,b])
    # This screening computation does not establish local dual validity;
    # that is checked below for the finally selected certificates.
    weights = np.zeros((len(pool),2,17,17),dtype=np.int64)
    for k,item in enumerate(pool):
        c = item['certificate']
        assert type(c['denominator']) is int and c['denominator']>0
        for row in c['tail_multipliers']:
            assert type(row['multiplier']) is int and row['multiplier']<=0
            assert 0<=row['axis']<15
            t0,t1 = row['thresholds']
            assert 0<=t0<=16 and 0<=t1<=16
            weights[k,int(row['axis']>=4),t0,t1] += row['multiplier']
    # Cache transforms by literal table to avoid repeated calculations for
    # duplicate profiles. Each transform uses the defining direct maximum.
    cache = {}
    constants = []
    for table_index,H in enumerate(tables):
        key = H.shape,H.tobytes()
        if key not in cache:
            cache[key] = direct_concentration(H)
        constants.append(cache[key])
        if (table_index+1)%2000==0:
            print(json.dumps(dict(stage='concentration_tables',done=table_index+1,total=len(tables),seconds=round(time.monotonic()-start,2))),flush=True)
    constants = np.array(constants)
    values = []
    for k,item in enumerate(pool):
        c = item['certificate']
        common = np.sum(constants*weights[k][None,:,:,:],axis=(1,2,3))
        numerators = np.array(counts,dtype=np.int64)@np.array(item['point_multipliers_by_pattern'],dtype=np.int64).T+common[:,None]
        values.append(-((-numerators)//c['denominator']))
    values = np.array(values)
    best = values.max(axis=0)
    required = set()
    for source_index,source in enumerate(cat['sources']):
        if source_index in selected:
            continue
        for i in source['target_cases']:
            d = cat['dominator'][i]
            for j in range(16):
                target = d if best[d,j]>0 else i
                assert best[target,j]>0,(source_index,i,j,best[target,j])
                required.add((target,j))
    for i in pair_indices.values():
        for j in range(16):
            assert best[i,j]>0,(i,j,best[i,j])
            required.add((i,j))
    print(json.dumps(dict(stage='positive_screen',required_inequalities=len(required),unique_capacity_tables=len(cache),seconds=round(time.monotonic()-start,2))),flush=True)
    # Cover the required inequalities greedily to reduce independent local
    # state enumeration. Every selected local inequality is recomputed.
    proof = []
    local_checks = 0
    tuples = 0
    for j,high_pattern in enumerate(itertools.product((0,1),repeat=4)):
        pending = {i for i,p in required if p==j}
        while pending:
            indices = sorted(pending)
            gain = np.sum(values[:,indices,j]>0,axis=1)
            k = int(gain.argmax())
            assert gain[k]>0
            actual,checked = local_minima(pool[k]['certificate'],high_pattern)
            assert actual==pool[k]['point_multipliers_by_pattern'][j]
            local_checks += 1
            tuples += checked
            covered = [i for i in indices if values[k,i,j]>0]
            # Repeat the final arithmetic with Python integers, independently
            # of the array matrix calculation (and without overflow).
            c = pool[k]['certificate']
            for i in covered:
                numerator = sum(n*p for n,p in zip(counts[i],actual))
                numerator += sum(row['multiplier']*int(constants[i,int(row['axis']>=4),*row['thresholds']]) for row in c['tail_multipliers'])
                bound = -(-numerator//c['denominator'])
                assert bound==int(values[k,i,j]) and bound>0
                proof.append(dict(table_index=i,high_pattern=high_pattern,dual_index=k,lower=bound))
            pending.difference_update(covered)
    return dict(verified=True,core_pattern=pattern,period=9699690,
                catalogue_audit=report,refined_profiles=len(selected),refined_pair_orbits=len(cap_tables),
                ordinary_refined_residuals=matching,paired_capacity_tables_checked=len(cap_tables),
                inequalities_checked=len(proof),local_checks=local_checks,plane_tuples_checked=tuples,
                uniform_lower=min(r['lower'] for r in proof),proof_selection=proof,
                seconds=round(time.monotonic()-start,3),
                scope=f'entire normalized low-anchor branch {tag} for squarefree moduli dividing 9699690; no assertion about other branches or larger supports')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='010')
    a = ap.parse_args()
    result = audit(list(map(int,a.pattern)))
    Path(f'data/audit-paired-tail-branch-{a.pattern}.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='proof_selection'}),flush=True)
