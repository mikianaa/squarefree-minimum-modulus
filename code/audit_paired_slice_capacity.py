"""Independent weighted-point audit of paired slice capacity certificates."""
import argparse
import itertools
import json
from pathlib import Path
import time


def audit(data):
    start = time.monotonic()
    R,moduli = data['core_residual'],data['moduli']
    assert len(set(R))==len(R) and len(moduli)%2==0
    layers = [[n for n in R if n%2==e] for e in (0,1)]
    sizes = list(map(len,layers))
    assert sizes==data['layer_sizes']
    pair_tables = {}
    for i,d in enumerate(moduli):
        for j in range(i+1,len(moduli)):
            q = moduli[j]
            table = [[0]*(sizes[1]+1) for _ in range(sizes[0]+1)]
            for a,b in itertools.product(sorted({n%d for n in R}),sorted({n%q for n in R})):
                partial = []
                for layer in layers:
                    weights = sorted((int(n%d==a)+int(n%q==b) for n in layer),reverse=True)
                    sums = [0]
                    for w in weights:
                        sums.append(sums[-1]+w)
                    partial.append(sums)
                for m0,x in enumerate(partial[0]):
                    for m1,y in enumerate(partial[1]):
                        table[m0][m1] = max(table[m0][m1],x+y)
            pair_tables[i,j] = table
    assert len(data['H_upper'])==len(data['matching_witnesses'])==sizes[0]+1
    for a in range(sizes[0]+1):
        assert len(data['H_upper'][a])==len(data['matching_witnesses'][a])==sizes[1]+1
        for b in range(sizes[1]+1):
            matching = data['matching_witnesses'][a][b]
            assert sorted(v for pair in matching for v in pair)==list(range(len(moduli)))
            value = sum(pair_tables[tuple(sorted(pair))][a][b] for pair in matching)
            assert value==data['H_upper'][a][b]
    return dict(verified=True,core_size=len(R),layer_sizes=sizes,pair_tables_checked=len(pair_tables),
                matching_witnesses_checked=(sizes[0]+1)*(sizes[1]+1),seconds=round(time.monotonic()-start,3),
                scope='valid upper bound on subset projection capacities; matching optimality not required or asserted by this audit')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',default='data/critical-010-paired-capacity.json')
    ap.add_argument('--output',default='data/audit-critical-010-paired-capacity.json')
    a = ap.parse_args()
    result = audit(json.loads(Path(a.source).read_text()))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
