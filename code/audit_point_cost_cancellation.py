"""Exact point-cost diagnostics for previously audited rational witnesses.

This checks local arithmetic and masses, not the resource/event constraints.
It is not an audit of covering feasibility or a lower-bound certificate.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def diagnose(data):
    source = data['input']
    witness = data['result'].get('point_primal_witness') or data['result']['coherent_primal_witness']
    R,high,pattern = [source[k] for k in ('core_residual','high_primes','high_pattern')]
    n,den = len(high),witness['denominator']
    groups = [S for k in range(1,n+1) for S in itertools.combinations(range(n),k)]
    masses,costs,positive = [0]*len(R),[0]*len(R),[0]*len(R)
    signs = dict(positive=0,negative=0,zero=0)
    for row in witness['support']:
        x,state,w = [row[k] for k in ('point_index','counts','weight')]
        assert type(w) is int and w>0 and len(state)==len(groups)
        D = [q-1-(b==R[x]%2) for q,b in zip(high,pattern)]
        assert all(type(v) is int and 0<=v<=(min(14,D[h]) if h<n else 16) for h,v in enumerate(state))
        f = [d-v for d,v in zip(D,state)]
        value = math.prod(f)-sum(state[h]*math.prod(f[i] for i in range(n) if i not in S)
                                for h,S in enumerate(groups[n:],n))
        masses[x] += w;costs[x] += w*value;positive[x] += w*max(0,value)
        signs['positive' if value>0 else 'negative' if value<0 else 'zero'] += 1
    assert masses==[den]*len(R)
    assert sum(costs)==witness['objective_numerator']
    return dict(verified_local_arithmetic=True,denominator=den,raw_numerator=sum(costs),
                positive_point_means_numerator=sum(max(0,v) for v in costs),
                positive_state_costs_numerator=sum(positive),state_signs=signs,
                positive_mean_points=sum(v>0 for v in costs),
                point_rows=[dict(point_index=x,residue=v,cost_numerator=costs[x],
                                 positive_state_costs_numerator=positive[x]) for x,v in enumerate(R)],
                scope='evaluation of a fixed witness, not optimization of either strengthened model')


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args = ap.parse_args()
    report = diagnose(json.loads(Path(args.source).read_text()))
    report['source'] = args.source
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='point_rows'}))
