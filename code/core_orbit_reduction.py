"""Exact CRT-label orbits and inclusion reduction for normalized 210 cores.

The 3-coordinate has only label 2 movable-free after anchoring. Labels
2,3,4 of the 5-coordinate and 2,...,6 of the 7-coordinate may be permuted.
Represent a residual by four rows (5-label) and six columns (7-label),
with three channel bits for (parity,3-label)=(0,2),(1,1),(1,2).
Enumerate the 120 column permutations and sort the three free rows. This
is the exact minimum over all 720 anchor-preserving label permutations.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from core_residual_geometry import pair_residuals,points_for


class CoreOrbits:
    def __init__(self,pattern):
        self.pattern = tuple(pattern)
        self.points = points_for(pattern)
        self.column_permutations = [(0,)+p for p in itertools.permutations(range(1,6))]
        self.row_permutations = [(0,)+p for p in itertools.permutations(range(1,4))]
        self.locations = [(a[1]-1,3*(a[2]-1)+(0 if e==0 else a[0])) for e,a in self.points]
        self.row_cache = {}
        self.placed_row_cache = [{} for _ in range(4)]
        self.placed_row_bits = [[(bit,1<<j) for j,(r,bit) in enumerate(self.locations) if r==row] for row in range(4)]

    def rows(self,mask):
        result = [0]*4
        for j,(row,bit) in enumerate(self.locations):
            if mask>>j&1:
                result[row] |= 1<<bit
        return tuple(result)

    def mask(self,rows):
        result = 0
        for row,value in enumerate(rows):
            cache = self.placed_row_cache[row]
            if value not in cache:
                cache[value] = sum(target for bit,target in self.placed_row_bits[row] if value>>bit&1)
            result |= cache[value]
        return result

    def permuted_row(self,row):
        if row not in self.row_cache:
            cells = [(row>>(3*j))&7 for j in range(6)]
            self.row_cache[row] = np.array([sum(cells[p[j]]<<(3*j) for j in range(6)) for p in self.column_permutations],dtype=np.uint32)
        return self.row_cache[row]

    def canonical(self,masks,batch=2000):
        masks = list(masks)
        result = []
        for start in range(0,len(masks),batch):
            rows = [self.rows(int(m)) for m in masks[start:start+batch]]
            candidates = np.stack([np.stack([self.permuted_row(r) for r in row],axis=1) for row in rows])
            candidates[:,:,1:].sort(axis=2)
            allowed = np.ones(candidates.shape[:2],dtype=bool)
            for column in range(4):
                minimum = np.where(allowed,candidates[:,:,column],2**32-1).min(axis=1)
                allowed &= candidates[:,:,column]==minimum[:,None]
            chosen = candidates[np.arange(len(rows)),allowed.argmax(axis=1)]
            result.extend(self.mask(tuple(map(int,row))) for row in chosen)
        return result

    def orbit(self,mask):
        rows = self.rows(mask)
        columns = [self.permuted_row(r) for r in rows]
        images = set()
        for j in range(len(self.column_permutations)):
            row = [int(c[j]) for c in columns]
            for p in self.row_permutations:
                images.add(self.mask([row[p[i]] for i in range(4)]))
        return images


class SubsetIndex:
    """Bitset index: a stored mask is contained in R iff it omits R's zeros."""
    def __init__(self,points):
        self.absent = [0]*points
        self.count = 0
        self.order = list(range(points))

    def extend(self,masks):
        masks = np.asarray(masks,dtype=np.uint64)
        for j in range(len(self.absent)):
            packed = np.packbits(((masks>>j)&1)==0,bitorder='little').tobytes()
            self.absent[j] |= int.from_bytes(packed,'little')<<self.count
        self.count += len(masks)
        self.order = sorted(range(len(self.absent)),key=lambda j:self.absent[j].bit_count())

    def contains_subset(self,mask):
        candidates = None
        for j in self.order:
            if mask>>j&1:
                continue
            candidates = self.absent[j] if candidates is None else candidates&self.absent[j]
            if not candidates:
                return False
        return self.count>0 if candidates is None else bool(candidates)


def minimal_orbits(orbits,canonical_masks):
    """Retain representatives minimal under inclusion after relabelling.

    Every discarded R contains an image of a retained E. If a later
    congruence family covers R, it also covers that E, and applying the
    inverse relabelling covers the representative. Thus it is enough to
    exclude extensions of every retained representative.
    """
    masks = sorted(set(canonical_masks),key=lambda m:(m.bit_count(),m))
    kept = []
    index = SubsetIndex(len(orbits.points))
    current_size = -1
    pending = []
    for mask in masks:
        size = mask.bit_count()
        if size!=current_size:
            if pending:
                index.extend(np.concatenate(pending))
                pending = []
            current_size = size
        if index.contains_subset(mask):
            continue
        kept.append(mask)
        pending.append(np.fromiter(orbits.orbit(mask),dtype=np.uint64))
    return kept


def run(full_even=False):
    start = time.monotonic()
    results = []
    for suffix in itertools.product((0,1),repeat=2):
        pattern = (0,)+suffix
        points,masks,total = pair_residuals(pattern,canonical_even=not full_even)
        orbits = CoreOrbits(pattern)
        canonical = sorted(set(orbits.canonical(masks)))
        print(json.dumps(dict(stage='canonical',pattern=pattern,raw=len(masks),orbits=len(canonical),seconds=round(time.monotonic()-start,3))),flush=True)
        minimal = minimal_orbits(orbits,canonical)
        row = dict(pattern=list(pattern),raw_masks=len(masks),canonical_orbits=len(canonical),
                   minimal_orbits=len(minimal),canonical_masks=canonical,minimal_masks=minimal,
                   size_histogram={str(s):sum(m.bit_count()==s for m in minimal) for s in sorted({m.bit_count() for m in minimal})})
        results.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ('canonical_masks','minimal_masks')}|{'seconds':round(time.monotonic()-start,3)}),flush=True)
    return dict(full_even=full_even,cases=results,seconds=round(time.monotonic()-start,3),
                scope='pair-residual core orbits before classes 105 and 210; no new covering obstruction')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full-even',action='store_true')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    Path(a.output).write_text(json.dumps(run(a.full_even),indent=2)+'\n')
