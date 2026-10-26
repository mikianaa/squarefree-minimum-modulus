"""Keep actual joint-slice profiles instead of mixing them in one envelope.

At each fixed final count pair, retain pointwise maximal H tables. This is a
quantitative dominance rule: every concentration constant for a smaller table
is no larger, and every dual multiplier on that constant is nonpositive.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from core_residual_geometry import pair_residuals
from portable_tail_certificates import CoreCapacityUpper,concentration_bounds
from compressed_tail_branch import evaluate


def joint_profile(model,mask):
    counts = [(mask&model.even).bit_count(),(mask&~model.even).bit_count()]
    profiles = []
    moduli = [d for d in range(3,211) if 210%d==0]
    for d,slices in zip(moduli,model.slices):
        if d%2==0:
            continue
        pairs = sorted({((s&mask&model.even).bit_count(),(s&mask&~model.even).bit_count()) for s in slices},reverse=True)
        kept = []
        second = -1
        for a,b in pairs:
            if b>second:
                kept.append((a,b))
                second = b
        profiles.append(tuple(kept))
    return tuple(counts)+(tuple(profiles),)


def maximal_tables(cases):
    by_counts = {}
    for i,row in enumerate(cases):
        by_counts.setdefault(tuple(row['counts']),[]).append(i)
    kept = []
    dominator = [None]*len(cases)
    for counts,indices in sorted(by_counts.items()):
        arrays = {i:np.array(cases[i]['H_upper'],dtype=np.int64).ravel() for i in indices}
        selected = []
        for i in sorted(indices,key=lambda i:(int(arrays[i].sum()),i),reverse=True):
            source = next((j for j in selected if np.all(arrays[j]>=arrays[i])),None)
            if source is None:
                selected.append(i)
                source = i
            dominator[i] = source
        kept.extend(selected)
    return kept,dominator


def run(pattern):
    start = time.monotonic()
    model = CoreCapacityUpper(pattern)
    _,masks,_ = pair_residuals(pattern)
    examples = {}
    for mask in masks:
        key = joint_profile(model,mask)
        if key not in examples or mask<examples[key]:
            examples[key] = mask
    sources = []
    cases = []
    lookup = {}
    for key,mask in sorted(examples.items()):
        H = model.upper(mask)
        indices = []
        for extra in (0,1):
            a,b = max(0,key[0]-1-(extra==0)),max(0,key[1]-1-(extra==1))
            target = H[:a+1,:b+1].copy()
            table_key = (a,b,target.tobytes())
            if table_key not in lookup:
                lookup[table_key] = len(cases)
                cases.append(dict(counts=[a,b],H_upper=target.tolist()))
            indices.append(lookup[table_key])
        sources.append(dict(counts_before=list(key[:2]),slice_pairs=key[2],target_cases=indices,
                            representative_residual=[n for j,n in enumerate(model.numbers) if mask>>j&1]))
    kept,dominator = maximal_tables(cases)
    pool = json.loads(Path('data/joint-tail-branch-010.json').read_text())['pool']
    counts = np.array([cases[i]['counts'] for i in kept])
    Ks = np.array([np.stack([concentration_bounds(np.array(cases[i]['H_upper']),s) for s in (14,16)]) for i in kept])
    best = np.array([evaluate(item,counts,Ks) for item in pool]).max(axis=0)
    result = dict(core_pattern=pattern,raw_pair_residuals=len(masks),sources=sources,cases=cases,
                  maximal_case_indices=kept,dominator=dominator,best_initial_bounds=best.tolist(),
                  remaining=int(np.sum(best<=0)),independent_audit=False,seconds=round(time.monotonic()-start,3),
                  scope='complete generated joint-slice profile catalogue with pointwise H dominance at fixed counts; requires independent audit')
    print(json.dumps(dict(pattern=pattern,raw=len(masks),profiles=len(sources),tables=len(cases),maximal=len(kept),remaining=result['remaining'],minimum=int(best.min()),seconds=result['seconds'])),flush=True)
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='010')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    Path(a.output).write_text(json.dumps(run(list(map(int,a.pattern))),separators=(',',':'))+'\n')
