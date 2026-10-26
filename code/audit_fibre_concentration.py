"""Standard-library audits of all core subsets and exact LP dual inequalities."""
import itertools
import json
import math
from pathlib import Path
import time


def gray_capacities(residual,moduli):
    sizes = [sum(x%2==e for x in residual) for e in (0,1)]
    H = [[-1]*(sizes[1]+1) for _ in range(sizes[0]+1)]
    H[0][0] = 0
    residues = [[x%d for d in moduli] for x in residual]
    counts = [[0]*d for d in moduli]
    histograms = [[d]+[0]*len(residual) for d in moduli]
    maxima = [0]*len(moduli)
    parity_counts = [0,0]
    total = 0
    for index in range(1,1<<len(residual)):
        j = (index&-index).bit_length()-1
        adding = bool((index^(index>>1))>>j&1)
        parity_counts[residual[j]%2] += 1 if adding else -1
        for k,a in enumerate(residues[j]):
            old = counts[k][a]
            new = old+(1 if adding else -1)
            counts[k][a] = new
            histograms[k][old] -= 1
            histograms[k][new] += 1
            if new>maxima[k]:
                maxima[k] = new
                total += 1
            elif old==maxima[k] and histograms[k][old]==0:
                maxima[k] -= 1
                total -= 1
        a,b = parity_counts
        H[a][b] = max(H[a][b],total)
    return H


def audit_dual(data,concentration,incumbent):
    residual = concentration['core_residual']
    assert data['core_residual']==residual
    high,pattern = data['high_primes'],data['high_pattern']
    certificate = data['certificate']
    scale = certificate['denominator']
    assert scale>0
    H = concentration['two_layer_capacities']
    maximum = len(concentration['moduli'])
    tables = [[[0]*(maximum+1) for _ in high] for _ in (0,1)]
    tail_term = 0
    for row in certificate['tail_multipliers']:
        axis = row['axis']
        t0,t1 = row['thresholds']
        bound = max(v-t0*a-t1*b for a,line in enumerate(H) for b,v in enumerate(line))
        assert bound==row['bound'] and row['multiplier']<=0
        tail_term += bound*row['multiplier']
        for e,t in enumerate((t0,t1)):
            for n in range(maximum+1):
                tables[e][axis][n] += row['multiplier']*max(0,n-t)
    links = certificate['link_multipliers']
    assert len(links)==len(residual) and all(len(row)==len(high) and all(v<=0 for v in row) for row in links)
    checked = 0
    for x,value in enumerate(residual):
        e = value%2
        D = [q-1-(b==e) for q,b in zip(high,pattern)]
        potentials = [[tables[e][i][n]+links[x][i]*n for n in range(min(maximum,D[i])+1)] for i in range(len(high))]
        actual = min(scale*math.prod(d-n for d,n in zip(D,ns))-sum(potentials[i][n] for i,n in enumerate(ns))
                     for ns in itertools.product(*(range(len(t)) for t in potentials)))
        assert actual==certificate['point_multipliers'][x]
        checked += math.prod(map(len,potentials))
    moduli = [r['modulus'] for r in incumbent['core_system']]
    resources = []
    for i,q in enumerate(high):
        for d in sorted(moduli):
            resources.append(min(sum(links[x][i] for x,n in enumerate(residual) if n%d==a) for a in {n%d for n in residual}))
    assert resources==certificate['resource_multipliers']
    numerator = sum(certificate['point_multipliers'])+sum(resources)+tail_term
    assert numerator==certificate['numerator']
    lower = -(-numerator//scale)
    assert lower==certificate['integer_lower']
    return dict(verified=True,integer_lower=lower,hit_count_columns_checked=checked,
                exact_numerator=numerator,exact_denominator=scale)


def main():
    start = time.monotonic()
    concentration = json.loads(Path('data/core-hit-concentration-7-15.json').read_text())
    incumbent = json.loads(Path('data/core-fibre-eight-plane.json').read_text())
    assert gray_capacities(concentration['core_residual'],concentration['moduli'])==concentration['two_layer_capacities']
    print(json.dumps(dict(stage='all_subsets',verified=True,count=1<<len(concentration['core_residual']),seconds=round(time.monotonic()-start,3))),flush=True)
    dual = audit_dual(json.loads(Path('data/identified-fibre-eight-lp.json').read_text()),concentration,incumbent)
    # Check the tail capacity separately from ordinary residue histograms.
    R = concentration['core_residual']
    high,pattern = incumbent['high_primes'],incumbent['high_pattern']
    tail = 0
    for mask in range(1,1<<len(high)):
        if mask.bit_count()<2:
            continue
        volumes = [math.prod(q-1-(b==e) for i,(q,b) in enumerate(zip(high,pattern)) if not mask>>i&1) for e in (0,1)]
        for d in range(1,106):
            if 105%d:
                continue
            odd = max(sum(volumes[n%2] for n in R if n%d==a) for a in range(d))
            even = max(sum(volumes[n%2] for n in R if n%(2*d)==a) for a in range(2*d))
            tail += odd+even
    assert tail==incumbent['tail_capacity']==132663
    out = dict(verified=True,all_core_subsets=1<<len(R),dual=dual,tail_capacity=tail,
               excludes_this_core=dual['integer_lower']>tail,seconds=round(time.monotonic()-start,3),
               scope='bounds for one specified core residual only; no global covering obstruction obtained')
    Path('data/audit-fibre-concentration-eight.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':
    main()
