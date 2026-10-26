"""Exact subset capacities controlling plane-projection hit concentration.

For a fixed core residual R and the 14 core shapes d>=3 dividing 210, set
H(m)=max_{T subset R, |T|=m} sum_d max_a |T intersect (a mod d)|.
Every choice of one projection per shape, with point hit counts n_x, satisfies
sum of the m largest n_x <= H(m), and hence
sum_x (n_x-t)_+ <= max_m(H(m)-t*m).
The finite enumeration below computes H, not an actual higher-prime cover.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np


def enumerate_capacities(residual,moduli,chunk_bits=18):
    assert len(residual)<=24, 'explicit exponential-work guard'
    placements = [[sum(1<<j for j,x in enumerate(residual) if x%d==a) for a in sorted({x%d for x in residual})]
                  for d in moduli]
    best = np.full(len(residual)+1,-1,dtype=np.int64)
    witnesses = [None]*(len(residual)+1)
    layer_mask = sum(1<<j for j,x in enumerate(residual) if x%2==0)
    layer_sizes = [sum(x%2==e for x in residual) for e in (0,1)]
    two_layer = np.full((layer_sizes[0]+1,layer_sizes[1]+1),-1,dtype=np.int64)
    for start in range(0,1<<len(residual),1<<chunk_bits):
        masks = np.arange(start,min(1<<len(residual),start+(1<<chunk_bits)),dtype=np.uint32)
        sizes = np.bitwise_count(masks)
        even_sizes = np.bitwise_count(masks&layer_mask)
        scores = np.zeros(len(masks),dtype=np.int16)
        for choices in placements:
            maximum = np.zeros(len(masks),dtype=np.int16)
            for placement in choices:
                maximum = np.maximum(maximum,np.bitwise_count(masks&placement))
            scores += maximum
        np.maximum.at(two_layer,(even_sizes,sizes-even_sizes),scores)
        for size in range(len(residual)+1):
            indices = np.flatnonzero(sizes==size)
            if not len(indices):
                continue
            i = indices[int(scores[indices].argmax())]
            if int(scores[i])>best[size]:
                best[size] = int(scores[i])
                witnesses[size] = int(masks[i])
    H = best.tolist()
    tails = [max(h-t*m for m,h in enumerate(H)) for t in range(len(moduli)+1)]
    return dict(capacities=H,layer_sizes=layer_sizes,two_layer_capacities=two_layer.tolist(),
                excess_bounds=tails,attaining_subset_masks=witnesses,subsets_checked=1<<len(residual))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/min-uncovered-9699690-witness.json')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    start = time.monotonic()
    data = json.loads(Path(a.source).read_text())
    core = {r['modulus']:r['residue'] for r in data['system'] if 210%r['modulus']==0}
    residual = [x for x in range(210) if not any(x%d==v for d,v in core.items())]
    result = dict(core_residual=residual,moduli=sorted(core),**enumerate_capacities(residual,sorted(core)))
    result['seconds'] = round(time.monotonic()-start,3)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='attaining_subset_masks'}))


if __name__=='__main__':
    main()
