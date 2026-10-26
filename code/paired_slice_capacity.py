"""Couple core slice capacities two moduli at a time.

For two projection resources A,B and a subset T of prescribed parity sizes,
the maximum of |T intersect A|+|T intersect B| is the sum of the largest point
weights 1_A+1_B in the two layers. Enumerate their residue choices exactly.
Every perfect matching of the 14 shapes supplies an upper bound on total
capacity. Minimize over matchings, recording a witnessing matching per cell.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import numpy as np


def paired_capacity(R,moduli):
    assert len(moduli)%2==0 and len(moduli)<=16
    parity = [sum(1<<j for j,n in enumerate(R) if n%2==e) for e in (0,1)]
    sizes = [p.bit_count() for p in parity]
    m0 = np.arange(sizes[0]+1)[:,None]
    m1 = np.arange(sizes[1]+1)[None,:]
    slices = [[sum(1<<j for j,n in enumerate(R) if n%d==a) for a in sorted({n%d for n in R})] for d in moduli]
    pair_tables = {}
    for i in range(len(moduli)):
        for j in range(i+1,len(moduli)):
            statistics = set()
            for A in slices[i]:
                for B in slices[j]:
                    U,I = A|B,A&B
                    statistics.add(tuple(v for p in parity for v in ((U&p).bit_count(),(I&p).bit_count())))
            table = np.zeros((sizes[0]+1,sizes[1]+1),dtype=np.int64)
            for u0,i0,u1,i1 in statistics:
                table = np.maximum(table,np.minimum(m0,u0)+np.minimum(m0,i0)+np.minimum(m1,u1)+np.minimum(m1,i1))
            pair_tables[i,j] = table
    decisions = {}
    @lru_cache(None)
    def visit(mask):
        if not mask:
            return np.zeros((sizes[0]+1,sizes[1]+1),dtype=np.int64)
        i = (mask&-mask).bit_length()-1
        rest = mask^(1<<i)
        best = np.full((sizes[0]+1,sizes[1]+1),10**9,dtype=np.int64)
        choice = np.full(best.shape,-1,dtype=np.int8)
        for j in range(i+1,len(moduli)):
            if rest>>j&1:
                candidate = pair_tables[i,j]+visit(rest^(1<<j))
                better = candidate<best
                best = np.minimum(best,candidate)
                choice[better] = j
        decisions[mask] = choice
        return best
    full = (1<<len(moduli))-1
    H = visit(full)
    matchings = []
    for a in range(sizes[0]+1):
        line = []
        for b in range(sizes[1]+1):
            mask = full
            pairs = []
            while mask:
                i = (mask&-mask).bit_length()-1
                j = int(decisions[mask][a,b])
                pairs.append([i,j])
                mask ^= (1<<i)|(1<<j)
            line.append(pairs)
        matchings.append(line)
    return dict(core_residual=R,moduli=moduli,layer_sizes=sizes,H_upper=H.tolist(),
                matching_witnesses=matchings,matching_dp_states=visit.cache_info().currsize,
                scope='upper subset-capacity bound from two-resource couplings; not a covering construction')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/critical-010-exact-subsets.json')
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    source = json.loads(Path(a.source).read_text())
    result = paired_capacity(source['core_residual'],source['moduli'])
    Path(a.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('core_residual','H_upper','matching_witnesses')}))
